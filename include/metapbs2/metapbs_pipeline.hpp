#pragma once
// =============================================================
// metapbs_pipeline.hpp — Meta-PBS Algorithm 1 (Exact Pipeline)
//
// Mirrors tfhe-go/metapbs/metapbs_d1.go:
//   BuildMetaPBSTV[T](eval, f) -> LookUpTable
// and tfhe-go/metapbs/blind_rotate_glwe.go:
//   BlindRotateGLWEAssign[T](eval, ctLWE, ctAcc, ctOut)
// and tfhe-go/metapbs/br_quotient_fix.go:
//   BlindRotateGLWEFromQuotient[T](eval, cquoK, ctAcc, ctOut, currentMod)
//   RunAlgorithm1QuotientFix(ct, f, brEval, trEvals, cfg, params)
// and tfhe-go/metapbs/tv_nonredundant.go:
//   runExplicitWithTV(ct, f, brEval, trEvals, cfg, tv, params)
//
// Algorithm 1 (paper §5):
//   Input: LWE ciphertext c encrypting m ∈ Z_t
//   Step 0: Build TV = paper-scale LUT of f
//   Step 1: C0 = BlindRotate(c, TV)      -- GLWE ciphertext encoding f(m)
//   Step 1': (cquo0, crem0) = HomDivRem(c, Q/(2N))
//   For k = 1..K:
//     Step 2k: C'_k = HomTruncRepeatShifted(C_{k-1}, [a_k, b_k], β_k, δ_k)
//     Step 2k': (cquo_k, crem_k) = HomDivRemAtScale(crem_{k-1}, currentMod, β_k)
//     Step 3k: C_k = BlindRotateGLWEFromQuotient(cquo_k, C'_k)
//   Output: Cout = SampleExtract(C_K, 0)
// =============================================================

#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "gatebootstrapping.hpp"
#include "detwfa.hpp"
#include "params.hpp"
#include "tlwe.hpp"
#include "trlwe.hpp"
#include "trgsw.hpp"

#include "metapbs2/sym_range.hpp"
#include "metapbs2/homdivrem.hpp"
#include "metapbs2/trunc_repeat.hpp"

namespace MetaPBS2 {

template <typename T>
inline std::make_signed_t<T> SignedTorus(T x) {
    return static_cast<std::make_signed_t<T>>(x);
}

template <class targetP, typename Int>
inline int Mod2N(Int x) {
    const Int twoN = static_cast<Int>(2 * targetP::n);
    Int r = x % twoN;
    if (r < 0) r += twoN;
    return static_cast<int>(r);
}

// =============================================================
// RoundConfig: parameters for one Algorithm 1 round
// Mirrors: type RoundParams struct { Beta, T int }
// =============================================================
struct RoundConfig {
    int beta;       // Expansion factor (β_k)
    int T;          // Number of redundant slots after expansion
    int delta_bound; // Pre-computed delta bound (optional, for validation)
};

// =============================================================
// Algorithm1Config: full configuration for Algorithm 1
// Mirrors: Algorithm1ExplicitConfig in tfhe-go
// =============================================================
struct Algorithm1Config {
    int K;          // Number of expansion rounds
    int t;          // Message modulus (= MessageModulus)
    std::vector<RoundConfig> rounds;   // length K
    // TR (TruncRepeat) gadget parameters baked into key generation
};

struct BlindRotatePruneStats {
    std::uint64_t total = 0;    // CMUXes that would run without periodic pruning.
    std::uint64_t skipped = 0;  // CMUXes skipped by period invariance.

    std::uint64_t executed() const { return total - skipped; }
    double prune_rate() const {
        return total == 0 ? 0.0 : static_cast<double>(skipped) / static_cast<double>(total);
    }
};

// =============================================================
// BuildMetaPBSTV: build the test vector (paper-scale LUT)
//
// Paper encoding uses scale = Q/t (NOT Q/(2t) as in standard TFHE).
// For r0 = 2N/t, places Q/t * f(i) at position i*r0 in [r0]_sym slots.
// The second half Z_t is handled negacyclically via X^N = -1.
//
// Mirrors: func BuildMetaPBSTV[T TorusInt](eval, f) LookUpTable[T]
// =============================================================
template <class targetP>
TFHEpp::Polynomial<targetP> BuildMetaPBSTV(
    std::function<int(int)> f, int t) {
    constexpr int N = targetP::n;
    if (t <= 0 || (2 * N) % t != 0 || t % 2 != 0)
        throw std::invalid_argument("Meta-PBS TV requires even t dividing 2N");

    // paper scale = Q / t = 2^64 / t
    const typename targetP::T scale =
        ((~typename targetP::T(0)) / static_cast<typename targetP::T>(t)) + 1;

    TFHEpp::Polynomial<targetP> tv = {};
    const int r0 = 2 * N / t;
    auto [lo, hi] = SymRange(r0);

    for (int i = 0; i < t / 2; i++) {
        int value = ((f(i) % t) + t) % t;
        typename targetP::T encoded = static_cast<typename targetP::T>(value) * scale;
        for (int j = lo; j <= hi; j++) {
            auto [sign, idx] = ReduceLaurentExponent(i * r0 + j, N);
            if (sign == +1)
                tv[idx] += encoded;
            else
                tv[idx] -= encoded;
        }
    }
    return tv;
}

// =============================================================
// BlindRotateGLWE: BR with GLWE ciphertext as initial accumulator
//
// Mirrors: func BlindRotateGLWEAssign[T](eval, ctLWE, ctAcc, ctOut)
//
// Standard BlindRotate starts from a polynomial TV.
// Here, ctAcc is a GLWE ciphertext used as the initial accumulator.
// The rotation is computed using the GLWE key's bootstrapping key.
//
// For each coordinate a_i of the LWE ciphertext:
//   ctOut = CMUX(BK[i], X^{a_i} ctOut, ctOut)
// with initial state: ctOut = X^{-b} * ctAcc
// =============================================================
template <class brP>
void BlindRotateGLWE(
    TFHEpp::TRLWE<typename brP::targetP>& ctOut,
    const TFHEpp::TLWE<typename brP::domainP>& ctLWE,
    const TFHEpp::TRLWE<typename brP::targetP>& ctAcc,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    constexpr int n = domP::k * domP::n;

    // Use TFHEpp's exact ModSwitch to accumulate rounding errors correctly
    TFHEpp::ModswitchTLWE<domP> moded;
    TFHEpp::BRModSwitch<brP, tgtP::nbit>(moded, ctLWE);

    uint32_t b_bar = moded[n];

    // Initial state: ctOut = X^{b_bar} * ctAcc
    for (int k = 0; k <= (int)tgtP::k; k++)
        TFHEpp::PolynomialMulByXai<tgtP>(ctOut[k], ctAcc[k], b_bar);

    // CMUX for each coefficient
    for (int i = 0; i < n; i++) {
        uint32_t a_bar = moded[i];
        if (a_bar == 0) continue;
        TFHEpp::CMUXwithPolynomialMulByXaiMinusOne<brP>(ctOut, bkfft[i], a_bar);
    }
}

// =============================================================
// BlindRotateGLWEFromQuotient: BR on a quotient LWE ciphertext
//
// Quotient coordinates are small integers from HomDivRem.
// BlindRotate's ModSwitch(x) ≈ 0 for small x, so we rescale:
//   rescaled[i] = cquo[i] * (Q / (2N))
// Then ModSwitch(rescaled[i]) = cquo[i], standard BR works.
//
// NOTE: rescale = Q/(2N), NOT Q/currentMod.
//
// Mirrors: func BlindRotateGLWEFromQuotient[T](eval, cquoK, ctAcc, ctOut, currentMod)
// =============================================================
template <class brP>
void BlindRotateGLWEFromQuotient(
    TFHEpp::TRLWE<typename brP::targetP>& ctOut,
    const TFHEpp::TLWE<typename brP::domainP>& cquo,
    const TFHEpp::TRLWE<typename brP::targetP>& ctAcc,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    int N,
    int period = 0,
    BlindRotatePruneStats* stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    constexpr int n = domP::k * domP::n;
    const int twoN = 2 * tgtP::n;
    (void)N;

    if (period > 1 && twoN % period != 0)
        throw std::invalid_argument("BlindRotateGLWEFromQuotient: period must divide 2N");

    int bbar = Mod2N<tgtP>(-SignedTorus(cquo[n]));
    for (int k = 0; k <= static_cast<int>(tgtP::k); k++)
        TFHEpp::PolynomialMulByXai<tgtP>(
            ctOut[k], ctAcc[k], static_cast<typename tgtP::T>(bbar));

    for (int i = 0; i < n; i++) {
        int abar = Mod2N<tgtP>(SignedTorus(cquo[i]));
        if (abar == 0) continue;
        if (stats) stats->total++;
        if (period > 1 && abar % period == 0) {
            if (stats) stats->skipped++;
            continue;
        }
        TFHEpp::CMUXwithPolynomialMulByXaiMinusOne<brP>(ctOut, bkfft[i], abar);
    }
}

// =============================================================
// SampleExtractIndex0: extract coefficient 0 from TRLWE
// (identical to TFHEpp::SampleExtractIndex with index=0)
// =============================================================
template <class P>
void SampleExtractIndex0(TFHEpp::TLWE<P>& tlwe, const TFHEpp::TRLWE<P>& trlwe) {
    TFHEpp::SampleExtractIndex<P>(tlwe, trlwe, 0);
}

template <class P>
TFHEpp::TRLWE<P> MakeAccumulator(const TFHEpp::Polynomial<P>& tv) {
    TFHEpp::TRLWE<P> acc{};
    acc[P::k] = tv;
    return acc;
}

// =============================================================
// RunAlgorithm1: complete Meta-PBS Algorithm 1 pipeline
//
// brP:   BootstrappingKey parameter (domainP=LWE domain, targetP=GLWE target)
// tgtP:  TRLWE target parameter (= brP::targetP)
//
// Mirrors: func RunAlgorithm1QuotientFix(ct, f, brEval, trEvals, cfg, params)
//
// Parameters:
//   ct      - input LWE ciphertext encrypting m ∈ Z_t
//   f       - function to evaluate: m → f(m) ∈ Z_t
//   bkfft   - bootstrapping key FFT
//   trkeys  - TruncRepeat keys for each round (length K)
//   cfg     - Algorithm 1 config (rounds, betas, t)
//   N       - ring degree of TRLWE (= brP::targetP::n)
// =============================================================
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
RunAlgorithm1(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    std::function<int(int)> f,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    int first_blind_rotate_period = 0,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    constexpr int N = tgtP::n;

    // Compute round parameters: a_k, b_k, delta_k
    // r_0 = 2N / t,  r_k = r_{k-1} * beta_k
    // [r]_sym = [lo, hi],  a_k = lo(r_{k-1}),  b_k = hi(r_{k-1})
    // delta_k = DeltaOffset(r_{k-1}, beta_k)
    std::vector<int> redundancies(cfg.K + 1);
    std::vector<int> deltas(cfg.K);
    redundancies[0] = 2 * N / cfg.t;
    for (int k = 0; k < cfg.K; k++) {
        deltas[k] = DeltaOffset(redundancies[k], cfg.rounds[k].beta);
        redundancies[k + 1] = redundancies[k] * cfg.rounds[k].beta;
    }

    // ---- Step 0: build test vector ----
    auto tv = BuildMetaPBSTV<tgtP>(f, cfg.t);

    // ---- Step 1': initial HomDivRem on ct ----
    TFHEpp::TLWE<domP> cquo0, crem0;
    HomDivRemLWE<domP>(cquo0, crem0, ct, static_cast<typename domP::T>(2 * N));

    // ---- Step 1: initial BlindRotate from the paper quotient ----
    TFHEpp::TRLWE<tgtP> C0;
    BlindRotateGLWEFromQuotient<brP>(
        C0, cquo0, MakeAccumulator<tgtP>(tv), bkfft, N,
        first_blind_rotate_period, prune_stats);

    // ---- Expansion rounds ----
    int current_mod = 2 * N;
    TFHEpp::TRLWE<tgtP> prev_GLWE = C0;
    TFHEpp::TLWE<domP> prev_crem = crem0;

    for (int k = 0; k < cfg.K; k++) {
        const auto& rnd = cfg.rounds[k];
        int a_k = -rnd.T;
        int b_k = rnd.T;
        int delta_k = deltas[k];

        // Step 2k: TruncRepeat + delta shift
        TFHEpp::TRLWE<tgtP> Ck_prime;
        HomTruncRepeatShifted<tgtP>(Ck_prime, prev_GLWE,
                                     a_k, b_k, rnd.beta, delta_k,
                                     trkeys[k]);

        // Step 2k': HomDivRem on remainder
        TFHEpp::TLWE<domP> cquo_k, crem_k;
        HomDivRemAtScale<domP>(cquo_k, crem_k, prev_crem, current_mod, rnd.beta);

        // Step 3k: BlindRotate from quotient
        TFHEpp::TRLWE<tgtP> Ck;
        BlindRotateGLWEFromQuotient<brP>(Ck, cquo_k, Ck_prime, bkfft, N);

        current_mod *= rnd.beta;
        prev_GLWE = Ck;
        prev_crem = crem_k;
    }

    // ---- Extract output: SampleExtract coefficient 0 ----
    TFHEpp::TLWE<tgtP> cout;
    SampleExtractIndex0<tgtP>(cout, prev_GLWE);
    return cout;
}

}  // namespace MetaPBS2
