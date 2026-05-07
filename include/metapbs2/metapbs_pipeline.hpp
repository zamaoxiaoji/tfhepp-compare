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
#include <cassert>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
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
    std::uint64_t pbs_calls = 0;  // PBS / blind-rotation calls observed.
    std::uint64_t total = 0;    // CMUXes that would run without periodic pruning.
    std::uint64_t cmux_calls = 0;  // CMUXes actually executed.
    std::uint64_t skipped = 0;  // CMUXes skipped by period invariance.
    std::uint64_t key_switch_count = 0;
    std::uint64_t pbs_count_reducer = 0;
    std::uint64_t pbs_count_gapmsb = 0;
    std::uint64_t pbs_count_bit_extract = 0;
    std::uint64_t pbs_count_bool_to_weight = 0;
    std::uint64_t pbs_count_final_msb = 0;
    std::uint64_t pbs_count_recursive_bit_extract = 0;
    std::uint64_t pbs_count_recursive_bool_to_weight = 0;
    std::vector<int> periods;   // Exact slot-domain periods used by each LUT.
    std::vector<std::uint64_t> total_by_pbs;
    std::vector<std::uint64_t> cmux_by_pbs;
    std::vector<std::uint64_t> skipped_by_pbs;

    std::uint64_t executed() const { return total - skipped; }
    double prune_rate() const {
        return total == 0 ? 0.0 : static_cast<double>(skipped) / static_cast<double>(total);
    }
};

inline bool IsPowerOfTwo(int x) {
    return x > 0 && (x & (x - 1)) == 0;
}

inline bool IsMultipleOfPeriod(std::uint32_t x, int period) {
    if (period <= 1) return true;
    const auto p = static_cast<std::uint32_t>(period);
    return IsPowerOfTwo(period) ? ((x & (p - 1)) == 0) : ((x % p) == 0);
}

inline bool IsMultipleOfPeriod(int x, int period) {
    if (period <= 1) return true;
    return IsPowerOfTwo(period) ? ((x & (period - 1)) == 0) : ((x % period) == 0);
}

template <class P>
inline typename P::T NegacyclicExtendedCoeff(const TFHEpp::Polynomial<P>& tv, int j) {
    constexpr int N = P::n;
    const int twoN = 2 * N;
    int idx = j % twoN;
    if (idx < 0) idx += twoN;
    if (idx < N) return tv[idx];
    return typename P::T(0) - tv[idx - N];
}

template <class P>
bool NegacyclicRotationInvariant(const TFHEpp::Polynomial<P>& tv, int shift) {
    constexpr int N = P::n;
    const int twoN = 2 * N;
    int s = shift % twoN;
    if (s < 0) s += twoN;
    for (int j = 0; j < twoN; j++)
        if (NegacyclicExtendedCoeff<P>(tv, j) !=
            NegacyclicExtendedCoeff<P>(tv, j + s))
            return false;
    return true;
}

template <class P>
std::size_t HashPolynomial(const TFHEpp::Polynomial<P>& tv) {
    std::size_t h = 1469598103934665603ull;
    for (auto x : tv) {
        h ^= static_cast<std::size_t>(x);
        h *= 1099511628211ull;
    }
    h ^= static_cast<std::size_t>(P::n);
    h *= 1099511628211ull;
    return h;
}

template <class P>
int ExactNegacyclicPeriodUncached(const TFHEpp::Polynomial<P>& tv) {
    constexpr int N = P::n;
    const int twoN = 2 * N;
    for (int p = 1; p <= twoN; p++)
        if (twoN % p == 0 && NegacyclicRotationInvariant<P>(tv, p))
            return p;
    return twoN;
}

template <class P>
int ExactNegacyclicPeriod(const TFHEpp::Polynomial<P>& tv) {
    struct CacheEntry {
        TFHEpp::Polynomial<P> tv;
        int period;
    };
    static std::unordered_map<std::size_t, std::vector<CacheEntry>> cache;
    const auto key = HashPolynomial<P>(tv);
    auto it = cache.find(key);
    if (it != cache.end()) {
        for (const auto& entry : it->second)
            if (entry.tv == tv) return entry.period;
    }
    const int period = ExactNegacyclicPeriodUncached<P>(tv);
    cache[key].push_back(CacheEntry{tv, period});
    return period;
}

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
    TFHEpp::BRModSwitch<brP, 1>(moded, ctLWE);

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
// BlindRotateTLWEWithPeriod: standard first-round BR on the original TLWE.
//
// This matches tfhe-go's runExplicitWithTV first round: C0 = BlindRotate(ct, TV).
// Optional period pruning is only sound for the caller-provided periodic LUT
// and only skips rotations whose mod-switched exponent is a multiple of the
// LUT period.
// =============================================================
template <class brP>
void BlindRotateTLWEWithPeriod(
    TFHEpp::TRLWE<typename brP::targetP>& ctOut,
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::Polynomial<typename brP::targetP>& tv,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    int period = 0,
    BlindRotatePruneStats* stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    constexpr int n = domP::k * domP::n;
    const int twoN = 2 * tgtP::n;

    if (period > 1 && twoN % period != 0)
        throw std::invalid_argument("BlindRotateTLWEWithPeriod: period must divide 2N");
    if (stats) {
        stats->pbs_calls++;
        stats->key_switch_count++;
        stats->periods.push_back(period > 0 ? period : twoN);
    }

#ifdef USE_KEY_BUNDLE
    if (stats) {
        TFHEpp::ModswitchTLWE<domP> moded;
        TFHEpp::BRModSwitch<brP, 1>(moded, ct);
        std::uint64_t local_total = 0;
        for (int i = 0; i < n; i++)
            if (moded[i] != 0) local_total++;
        stats->total += local_total;
        stats->cmux_calls += local_total;
        stats->total_by_pbs.push_back(local_total);
        stats->cmux_by_pbs.push_back(local_total);
        stats->skipped_by_pbs.push_back(0);
    }
    TFHEpp::BlindRotate<brP>(ctOut, ct, bkfft, tv);
#else
    TFHEpp::ModswitchTLWE<domP> moded;
    TFHEpp::BRModSwitch<brP, 1>(moded, ct);

    ctOut = {};
    TFHEpp::PolynomialMulByXai<tgtP>(ctOut[tgtP::k], tv, moded[n]);

    std::uint64_t local_total = 0;
    std::uint64_t local_skipped = 0;
    std::uint64_t local_cmux = 0;
    for (int i = 0; i < n; i++) {
        uint32_t abar = moded[i];
        if (abar == 0) continue;
        local_total++;
        if (period > 1 && IsMultipleOfPeriod(abar, period)) {
#ifndef NDEBUG
            assert(NegacyclicRotationInvariant<tgtP>(
                tv, static_cast<int>(abar % static_cast<std::uint32_t>(twoN))));
#endif
            local_skipped++;
            continue;
        }
        local_cmux++;
        TFHEpp::CMUXwithPolynomialMulByXaiMinusOne<brP>(ctOut, bkfft[i], abar);
    }
    if (stats) {
        stats->total += local_total;
        stats->skipped += local_skipped;
        stats->cmux_calls += local_cmux;
        stats->total_by_pbs.push_back(local_total);
        stats->cmux_by_pbs.push_back(local_cmux);
        stats->skipped_by_pbs.push_back(local_skipped);
    }
#endif
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
    if (stats) {
        stats->pbs_calls++;
        stats->key_switch_count++;
        stats->periods.push_back(period > 0 ? period : twoN);
    }

    int bbar = Mod2N<tgtP>(-SignedTorus(cquo[n]));
    for (int k = 0; k <= static_cast<int>(tgtP::k); k++)
        TFHEpp::PolynomialMulByXai<tgtP>(
            ctOut[k], ctAcc[k], static_cast<typename tgtP::T>(bbar));

    std::uint64_t local_total = 0;
    std::uint64_t local_skipped = 0;
    std::uint64_t local_cmux = 0;
    for (int i = 0; i < n; i++) {
        int abar = Mod2N<tgtP>(SignedTorus(cquo[i]));
        if (abar == 0) continue;
        local_total++;
        if (period > 1 && IsMultipleOfPeriod(abar, period)) {
            local_skipped++;
            continue;
        }
        local_cmux++;
        TFHEpp::CMUXwithPolynomialMulByXaiMinusOne<brP>(ctOut, bkfft[i], abar);
    }
    if (stats) {
        stats->total += local_total;
        stats->skipped += local_skipped;
        stats->cmux_calls += local_cmux;
        stats->total_by_pbs.push_back(local_total);
        stats->cmux_by_pbs.push_back(local_cmux);
        stats->skipped_by_pbs.push_back(local_skipped);
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

template <class brP>
void FirstRoundRemainderFromBRModSwitch(
    TFHEpp::TLWE<typename brP::domainP>& crem,
    const TFHEpp::TLWE<typename brP::domainP>& ct) {
    using domP = typename brP::domainP;
    using tgtP = typename brP::targetP;
    constexpr int n = domP::k * domP::n;
    constexpr int twoN = 2 * tgtP::n;
    constexpr int shift =
        std::numeric_limits<typename domP::T>::digits - 1 - tgtP::nbit;
    const auto scale = static_cast<typename domP::T>(typename domP::T(1) << shift);

    TFHEpp::ModswitchTLWE<domP> moded;
    TFHEpp::BRModSwitch<brP, 1>(moded, ct);

    for (int i = 0; i < n; i++)
        crem[i] = ct[i] - static_cast<typename domP::T>(moded[i]) * scale;

    const uint32_t exponent = moded[n] % static_cast<uint32_t>(twoN);
    const uint32_t q_body =
        (static_cast<uint32_t>(twoN) - exponent) % static_cast<uint32_t>(twoN);
    crem[n] = ct[n] - static_cast<typename domP::T>(q_body) * scale;
}

template <class brP>
TFHEpp::TLWE<typename brP::targetP>
RunAlgorithm1WithTV(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::Polynomial<typename brP::targetP>& tv,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    int first_blind_rotate_period = 0,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    constexpr int N = tgtP::n;

    if (cfg.K < 0 || static_cast<int>(cfg.rounds.size()) != cfg.K)
        throw std::invalid_argument("RunAlgorithm1WithTV requires cfg.rounds.size() == cfg.K");
    if (static_cast<int>(trkeys.size()) < cfg.K)
        throw std::invalid_argument("RunAlgorithm1WithTV requires one TruncRepeat key per round");
    if (cfg.t <= 0 || (2 * N) % cfg.t != 0)
        throw std::invalid_argument("RunAlgorithm1WithTV requires cfg.t dividing 2N");

    std::vector<int> redundancies(cfg.K + 1);
    std::vector<int> deltas(cfg.K);
    redundancies[0] = 2 * N / cfg.t;
    for (int k = 0; k < cfg.K; k++) {
        deltas[k] = DeltaOffset(redundancies[k], cfg.rounds[k].beta);
        redundancies[k + 1] = redundancies[k] * cfg.rounds[k].beta;
    }

    TFHEpp::TRLWE<tgtP> prev_GLWE;
    BlindRotateTLWEWithPeriod<brP>(
        prev_GLWE, ct, tv, bkfft, first_blind_rotate_period, prune_stats);

    TFHEpp::TLWE<domP> prev_crem;
    FirstRoundRemainderFromBRModSwitch<brP>(prev_crem, ct);

    int current_mod = 2 * N;
    for (int k = 0; k < cfg.K; k++) {
        const auto& rnd = cfg.rounds[k];

        TFHEpp::TRLWE<tgtP> Ck_prime;
        HomTruncRepeatShifted<tgtP>(Ck_prime, prev_GLWE,
                                     -rnd.T, rnd.T, rnd.beta, deltas[k],
                                     trkeys[k]);

        TFHEpp::TLWE<domP> cquo_k, crem_k;
        HomDivRemAtScale<domP>(cquo_k, crem_k, prev_crem, current_mod, rnd.beta);

        TFHEpp::TRLWE<tgtP> Ck;
        BlindRotateGLWEFromQuotient<brP>(
            Ck, cquo_k, Ck_prime, bkfft, N, 0, prune_stats);

        current_mod *= rnd.beta;
        prev_GLWE = Ck;
        prev_crem = crem_k;
    }

    TFHEpp::TLWE<tgtP> cout;
    SampleExtractIndex0<tgtP>(cout, prev_GLWE);
    return cout;
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
    auto tv = BuildMetaPBSTV<tgtP>(f, cfg.t);
    return RunAlgorithm1WithTV<brP>(
        ct, tv, bkfft, trkeys, cfg, first_blind_rotate_period, prune_stats);
}

}  // namespace MetaPBS2
