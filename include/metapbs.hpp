#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "params.hpp"
#include "gatebootstrapping.hpp"

namespace TFHEpp::metapbs {

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  Prune statistics (public — useful for benchmarking)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

struct BlindRotatePruneStats {
    std::uint64_t total = 0;    // CMUXes that would have run (a != 0).
    std::uint64_t skipped = 0;  // CMUXes skipped because the rotation is a no-op.
};

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  Internal implementation details
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

namespace detail {

// ── Math helpers ────────────────────────────────────────────────────────────

constexpr int sym_min(const int r) { return (r == 2) ? 0 : -(r / 2); }
constexpr int sym_max(const int r) { return sym_min(r) + r - 1; }

template <class Int>
constexpr Int centered_mod(Int x, const Int m)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    Int r = x % m;
    const Int half_floor = m / 2;
    const Int half_ceil_minus1 = (m - 1) / 2;
    if (r < -half_floor) r += m;
    if (r > half_ceil_minus1) r -= m;
    return r;
}

template <class TargetP, class Int>
constexpr int mod2N(const Int x)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    constexpr Int twoN = static_cast<Int>(2) * static_cast<Int>(TargetP::n);
    Int r = x % twoN;
    if (r < 0) r += twoN;
    return static_cast<int>(r);
}

constexpr int DeltaRB(const int r, const int B)
{
    return sym_min(r * B) - B * sym_min(r) - sym_min(B);
}

// ── LWE integer representation ─────────────────────────────────────────────

template <class DomainP, class Int = std::int64_t>
struct LWE {
    std::array<Int, DomainP::k * DomainP::n + 1> c{};
    Int modulus = 0;
};

template <class DomainP, class Int = std::int64_t>
LWE<DomainP, Int> LiftTLWEToInt(const TLWE<DomainP> &in)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    constexpr int digits = std::numeric_limits<typename DomainP::T>::digits;
    static_assert(digits < 63,
                  "LiftTLWEToInt assumes modulus fits into signed 64-bit");
    LWE<DomainP, Int> out;
    out.modulus = static_cast<Int>(1) << digits;
    for (size_t i = 0; i < in.size(); i++) {
        using SignedT = std::make_signed_t<typename DomainP::T>;
        out.c[i] = static_cast<Int>(static_cast<SignedT>(in[i]));
    }
    return out;
}

// ── HomDivRem (Definition 1) ───────────────────────────────────────────────

template <class DomainP, class Int = std::int64_t>
std::pair<LWE<DomainP, Int>, LWE<DomainP, Int>> HomDivRem(
    const LWE<DomainP, Int> &in, const Int q_quo)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    const Int q = in.modulus;
    if (q_quo <= 0 || q <= 0 || (q % q_quo) != 0)
        throw std::invalid_argument("HomDivRem: require q_quo | q and q_quo>0");
    const Int q_rem = q / q_quo;

    LWE<DomainP, Int> quo;
    quo.modulus = q_quo;
    LWE<DomainP, Int> rem;
    rem.modulus = q_rem;

    for (size_t i = 0; i < in.c.size(); i++) {
        const Int r = centered_mod<Int>(in.c[i], q_rem);
        rem.c[i] = r;
        quo.c[i] = (in.c[i] - r) / q_rem;
    }
    return {quo, rem};
}

// ── TLWE / TRLWE helpers ───────────────────────────────────────────────────

template <class P>
inline void TLWENegate(TLWE<P> &c)
{
    for (auto &x : c) x = -x;
}

template <class P>
inline void TRLWEMulByXai(TRLWE<P> &out, const TRLWE<P> &in, const int a)
{
    const int twoN = 2 * static_cast<int>(P::n);
    int aa = a % twoN;
    if (aa < 0) aa += twoN;
    for (int k = 0; k < P::k + 1; k++)
        PolynomialMulByXai<P>(out[k], in[k], static_cast<typename P::T>(aa));
}

// ── TruncRepeat (Definition 5) ─────────────────────────────────────────────

template <class P>
void TruncRepeat(TRLWE<P> &out, const TRLWE<P> &in, const int T, const int B,
                 const AnnihilateKey<P> &ahk)
{
    static_assert(P::k == 1, "TruncRepeat currently assumes k=1");
    if (T < 0 || B <= 0) throw std::invalid_argument("TruncRepeat: bad T/B");
    if ((2 * T + 1) * B > static_cast<int>(P::n))
        throw std::invalid_argument("TruncRepeat: (2T+1)*B must be <= N");

    const int kmin = sym_min(B);
    const int kmax = sym_max(B);
    const int out_min = -T * B + kmin;
    const int out_max = T * B + kmax;
    const int pos_len = out_max >= 0 ? out_max + 1 : 0;
    const int neg_len = out_min < 0 ? -out_min : 0;

    std::vector<TLWE<P>> coeffs(2 * T + 1);
    for (int d = -T; d <= T; d++) {
        const int idx = (d >= 0) ? d : (static_cast<int>(P::n) + d);
        SampleExtractIndex<P>(coeffs[d + T], in, idx);
        if (d < 0) TLWENegate<P>(coeffs[d + T]);
    }

    std::vector<TLWE<P>> pos(pos_len);
    std::vector<TLWE<P>> neg(neg_len);
    for (auto &c : pos) c = {};
    for (auto &c : neg) c = {};

    for (int d = -T; d <= T; d++) {
        const TLWE<P> &cd = coeffs[d + T];
        for (int k = kmin; k <= kmax; k++) {
            const int out_deg = d * B + k;
            if (out_deg >= 0) {
                pos[out_deg] = cd;
            }
            else {
                const int p = out_deg + neg_len;
                TLWE<P> tmp = cd;
                TLWENegate<P>(tmp);
                neg[p] = tmp;
            }
        }
    }

    TRLWE<P> pos_trlwe;
    TLWE2TRLWEPacking<P>(pos_trlwe, pos, ahk);
    if (neg_len == 0) {
        out = pos_trlwe;
        return;
    }

    TRLWE<P> neg_trlwe;
    TLWE2TRLWEPacking<P>(neg_trlwe, neg, ahk);

    TRLWE<P> neg_rot;
    TRLWEMulByXai<P>(neg_rot, neg_trlwe, static_cast<int>(P::n) - neg_len);
    TRLWEAdd<P>(out, pos_trlwe, neg_rot);
}

// ── Torus encoding ─────────────────────────────────────────────────────────

template <class TargetP, class SignedInt = std::int64_t>
typename TargetP::T EncodeTorus(const SignedInt x, const std::int64_t t)
{
    static_assert(std::is_signed_v<SignedInt>, "x must be signed");
    constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
    static_assert(digits <= 64,
                  "EncodeTorus currently supports up to 64-bit torus types");
    const __int128_t q = static_cast<__int128_t>(1) << digits;
    const __int128_t val = (static_cast<__int128_t>(x) * q) / t;
    return static_cast<typename TargetP::T>(val);
}

// ── Test vector construction ───────────────────────────────────────────────

template <class TargetP, class SignedInt = std::int64_t>
Polynomial<TargetP> GenerateTestVectorNegacyclic(
    const std::vector<SignedInt> &f_half, const std::int64_t t)
{
    constexpr int N = static_cast<int>(TargetP::n);
    const int twoN = 2 * N;
    if ((t & 1) != 0 || t <= 0 || t > twoN || (twoN % t) != 0)
        throw std::invalid_argument("GenerateTestVectorNegacyclic: bad t");
    if (static_cast<std::int64_t>(f_half.size()) != (t / 2))
        throw std::invalid_argument(
            "GenerateTestVectorNegacyclic: f_half size must be t/2");

    const int block = twoN / static_cast<int>(t);
    const int shift = block / 2;

    Polynomial<TargetP> tv{};
    for (auto &x : tv) x = 0;

    for (int i = 0; i < static_cast<int>(t / 2); i++) {
        const auto y = EncodeTorus<TargetP, SignedInt>(f_half[i], t);
        for (int j = 0; j < block; j++) {
            int e = i * block + j - shift;
            e %= twoN;
            if (e < 0) e += twoN;
            if (e >= N)
                tv[e - N] -= y;
            else
                tv[e] += y;
        }
    }
    return tv;
}

template <class TargetP>
Polynomial<TargetP> GenerateTestVectorNegacyclicTorus(
    const std::vector<typename TargetP::T> &f_half_torus, const std::int64_t t)
{
    constexpr int N = static_cast<int>(TargetP::n);
    const int twoN = 2 * N;
    if ((t & 1) != 0 || t <= 0 || t > twoN || (twoN % t) != 0)
        throw std::invalid_argument("GenerateTestVectorNegacyclicTorus: bad t");
    if (static_cast<std::int64_t>(f_half_torus.size()) != (t / 2))
        throw std::invalid_argument(
            "GenerateTestVectorNegacyclicTorus: f_half_torus size must be t/2");

    const int block = twoN / static_cast<int>(t);
    const int shift = block / 2;

    Polynomial<TargetP> tv{};
    tv.fill(0);

    for (int i = 0; i < static_cast<int>(t / 2); i++) {
        const auto y = f_half_torus[i];
        for (int j = 0; j < block; j++) {
            int e = i * block + j - shift;
            e %= twoN;
            if (e < 0) e += twoN;
            if (e >= N)
                tv[e - N] -= y;
            else
                tv[e] += y;
        }
    }
    return tv;
}

// ── Bit extraction test vectors ────────────────────────────────────────────

template <class TargetP>
const Polynomial<TargetP> &BitExtractTestVector2N(const int bit)
{
    constexpr int N = static_cast<int>(TargetP::n);
    if (bit < 0 || bit >= static_cast<int>(TargetP::nbit))
        throw std::invalid_argument("BitExtractTestVector2N: bad bit index");

    struct Cache {
        std::array<Polynomial<TargetP>, TargetP::nbit> tv{};
        std::array<bool, TargetP::nbit> inited{};
    };
    static Cache cache;

    if (!cache.inited[bit]) {
        Polynomial<TargetP> tv{};
        tv.fill(0);
        constexpr auto half_q =
            static_cast<typename TargetP::T>(1)
            << (std::numeric_limits<typename TargetP::T>::digits - 1);
        for (int i = 0; i < N; i++)
            tv[i] = ((i >> bit) & 1) ? half_q : 0;
        cache.tv[bit] = tv;
        cache.inited[bit] = true;
    }
    return cache.tv[bit];
}

// ── BlindRotate (Meta-PBS variant, consumes integer quotient) ───────────────

template <class brP, class Int = std::int64_t, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename brP::targetP> &res,
                 const LWE<typename brP::domainP, Int> &tlwe_quo,
                 const BootstrappingKeyFFT<brP> &bkfft,
                 const Polynomial<typename brP::targetP> &testvector)
{
    static_assert(brP::targetP::k == 1,
                  "Meta-PBS currently assumes TLWE/TRLWE k=1");

    const int bbar = mod2N<typename brP::targetP, Int>(-tlwe_quo.c[brP::domainP::k * brP::domainP::n]);

    res = {};
    PolynomialMulByXai<typename brP::targetP>(
        res[brP::targetP::k], testvector,
        static_cast<typename brP::targetP::T>(bbar));

#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < brP::domainP::k * brP::domainP::n / brP::Addends; i++) {
        std::array<typename brP::domainP::T, brP::Addends> bara{};
        for (int j = 0; j < brP::Addends; j++)
            bara[j] = static_cast<typename brP::domainP::T>(mod2N<typename brP::targetP, Int>(
                tlwe_quo.c[brP::Addends * i + j]));
        alignas(64) TRGSWFFT<typename brP::targetP> BKadded;
        KeyBundleFFT<brP>(BKadded, bkfft[i], bara);
        ExternalProduct<typename brP::targetP>(res, res, BKadded);
    }
#else
    for (int i = 0; i < brP::domainP::k * brP::domainP::n; i++) {
        const int a = mod2N<typename brP::targetP, Int>(tlwe_quo.c[i]);
        if (a == 0) continue;
        CMUXwithPolynomialMulByXaiMinusOne<brP>(res, bkfft[i], a);
    }
#endif
}

template <class brP, class Int = std::int64_t, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename brP::targetP> &res,
                 const LWE<typename brP::domainP, Int> &tlwe_quo,
                 const BootstrappingKeyFFT<brP> &bkfft,
                 const TRLWE<typename brP::targetP> &testvector)
{
    static_assert(brP::targetP::k == 1,
                  "Meta-PBS currently assumes TLWE/TRLWE k=1");

    const int bbar = mod2N<typename brP::targetP, Int>(-tlwe_quo.c[brP::domainP::k * brP::domainP::n]);

    for (int k = 0; k < brP::targetP::k + 1; k++)
        PolynomialMulByXai<typename brP::targetP>(
            res[k], testvector[k],
            static_cast<typename brP::targetP::T>(bbar));

#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < brP::domainP::k * brP::domainP::n / brP::Addends; i++) {
        std::array<typename brP::domainP::T, brP::Addends> bara{};
        for (int j = 0; j < brP::Addends; j++)
            bara[j] = static_cast<typename brP::domainP::T>(mod2N<typename brP::targetP, Int>(
                tlwe_quo.c[brP::Addends * i + j]));
        alignas(64) TRGSWFFT<typename brP::targetP> BKadded;
        KeyBundleFFT<brP>(BKadded, bkfft[i], bara);
        ExternalProduct<typename brP::targetP>(res, res, BKadded);
    }
#else
    for (int i = 0; i < brP::domainP::k * brP::domainP::n; i++) {
        const int a = mod2N<typename brP::targetP, Int>(tlwe_quo.c[i]);
        if (a == 0) continue;
        CMUXwithPolynomialMulByXaiMinusOne<brP>(res, bkfft[i], a);
    }
#endif
}

// ── BlindRotatePeriodic (with CMUX pruning) ────────────────────────────────

template <class brP, class Int>
void BlindRotatePeriodic(TRLWE<typename brP::targetP> &res,
                         const LWE<typename brP::domainP, Int> &tlwe_quo,
                         const BootstrappingKeyFFT<brP> &bkfft,
                         const Polynomial<typename brP::targetP> &testvector,
                         const int period,
                         BlindRotatePruneStats *stats)
{
    static_assert(brP::targetP::k == 1,
                  "Meta-PBS currently assumes TLWE/TRLWE k=1");

    if (period <= 1) {
        BlindRotate<brP, Int>(res, tlwe_quo, bkfft, testvector);
        return;
    }
    constexpr int N = static_cast<int>(brP::targetP::n);
    const int twoN = 2 * N;
    if ((twoN % period) != 0)
        throw std::invalid_argument("BlindRotatePeriodic: period must divide 2N");

    const int bbar = mod2N<typename brP::targetP, Int>(
        -tlwe_quo.c[brP::domainP::k * brP::domainP::n]);

    res = {};
    PolynomialMulByXai<typename brP::targetP>(
        res[brP::targetP::k], testvector,
        static_cast<typename brP::targetP::T>(bbar));

#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < brP::domainP::k * brP::domainP::n / brP::Addends; i++) {
        std::array<typename brP::domainP::T, brP::Addends> bara{};
        bool all_zero = true;
        for (int j = 0; j < brP::Addends; j++) {
            const int a = mod2N<typename brP::targetP, Int>(
                tlwe_quo.c[brP::Addends * i + j]);
            int apruned = a;
            if (a != 0 && (a % period) == 0) apruned = 0;

            bara[j] = static_cast<typename brP::domainP::T>(apruned);
            if (stats && a != 0) {
                stats->total++;
                if (apruned == 0) stats->skipped++;
            }
            if (apruned != 0) all_zero = false;
        }
        if (all_zero) continue;
        alignas(64) TRGSWFFT<typename brP::targetP> BKadded;
        KeyBundleFFT<brP>(BKadded, bkfft[i], bara);
        ExternalProduct<typename brP::targetP>(res, res, BKadded);
    }
#else
    for (int i = 0; i < brP::domainP::k * brP::domainP::n; i++) {
        const int a = mod2N<typename brP::targetP, Int>(tlwe_quo.c[i]);
        if (a == 0) continue;
        if (stats) stats->total++;
        if ((a % period) == 0) {
            if (stats) stats->skipped++;
            continue;
        }
        CMUXwithPolynomialMulByXaiMinusOne<brP>(res, bkfft[i], a);
    }
#endif
}

// ── MetaPBSExtractBit2N (internal bit extraction helper) ───────────────────

template <class brP, class Int = std::int64_t>
void MetaPBSExtractBit2N(TLWE<typename brP::targetP> &out,
                         const TLWE<typename brP::domainP> &in,
                         const BootstrappingKeyFFT<brP> &bkfft,
                         const AnnihilateKey<typename brP::targetP> &ahk,
                         const int bit, const std::vector<int> &betas,
                         const std::vector<int> &Ts,
                         BlindRotatePruneStats *prune_stats = nullptr)
{
    if (betas.size() != Ts.size())
        throw std::invalid_argument("MetaPBSExtractBit2N: betas/Ts size mismatch");

    constexpr int N = static_cast<int>(brP::targetP::n);
    const int twoN = 2 * N;

    const auto &tv = BitExtractTestVector2N<typename brP::targetP>(bit);
    const int period = 1 << (bit + 1);

    int r = 1;

    const auto lifted = LiftTLWEToInt<typename brP::domainP, Int>(in);
    auto [cquo, crem] =
        HomDivRem<typename brP::domainP, Int>(lifted, static_cast<Int>(twoN));

    TRLWE<typename brP::targetP> C;
    BlindRotatePeriodic<brP, Int>(C, cquo, bkfft, tv, period, prune_stats);

    for (size_t i = 0; i < betas.size(); i++) {
        const int beta = betas[i];
        const int T = Ts[i];
        if (beta < 2)
            throw std::invalid_argument("MetaPBSExtractBit2N: beta must be >=2");

        TRLWE<typename brP::targetP> Ctr;
        TruncRepeat<typename brP::targetP>(Ctr, C, T, beta, ahk);
        TRLWE<typename brP::targetP> Cprime;
        TRLWEMulByXai<typename brP::targetP>(Cprime, Ctr, DeltaRB(r, beta));

        if ((crem.modulus % static_cast<Int>(beta)) != 0)
            throw std::invalid_argument(
                "MetaPBSExtractBit2N: beta must divide remainder modulus");
        auto [cquo_next, crem_next] =
            HomDivRem<typename brP::domainP, Int>(crem, static_cast<Int>(beta));

        TRLWE<typename brP::targetP> Cnext;
        BlindRotate<brP, Int>(Cnext, cquo_next, bkfft, Cprime);

        C = Cnext;
        crem = crem_next;
        r *= beta;
    }

    SampleExtractIndex<typename brP::targetP>(out, C, 0);
}

// Default iteration schedule overload.
template <class brP, class Int = std::int64_t>
void MetaPBSExtractBit2N(TLWE<typename brP::targetP> &out,
                         const TLWE<typename brP::domainP> &in,
                         const BootstrappingKeyFFT<brP> &bkfft,
                         const AnnihilateKey<typename brP::targetP> &ahk,
                         const int bit,
                         BlindRotatePruneStats *prune_stats = nullptr)
{
    constexpr int N = static_cast<int>(brP::targetP::n);
    static const std::vector<int> betas = {8, 8};
    static const std::vector<int> Ts = {
        (N / betas[0] - 1) / 2,
        (N / betas[1] - 1) / 2,
    };
    MetaPBSExtractBit2N<brP, Int>(out, in, bkfft, ahk, bit, betas, Ts,
                                 prune_stats);
}

}  // namespace detail

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  Public helper: test vector construction (needed by advanced users)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

using detail::GenerateTestVectorNegacyclic;
using detail::GenerateTestVectorNegacyclicTorus;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  Public Interface 1: MetaPBS
//
//  Algorithm 1 (Meta-PBS) for a negacyclic LUT (single-output), with optional
//  periodic-invariance CMUX pruning in the *first* BlindRotate.
//
//  Parameters:
//    out       — output TLWE ciphertext
//    in        — input TLWE ciphertext
//    bkfft     — bootstrapping key (FFT domain)
//    ahk       — annihilate key (for TruncRepeat / LWE-to-RLWE packing)
//    testvector— the negacyclic LUT to evaluate
//    betas     — β parameters for each TruncRepeat iteration
//    Ts        — T parameters for each TruncRepeat iteration
//    t         — plaintext space size (must be even, divide 2N, and > 0)
//    period    — TV periodicity for CMUX pruning (default 0 = no pruning).
//                When period > 1, the first BlindRotate skips CMUXes whose
//                rotation amount is a multiple of period.
//    stats     — optional pointer to collect pruning statistics
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

template <class brP, class Int = std::int64_t>
void MetaPBS(TLWE<typename brP::targetP> &out,
             const TLWE<typename brP::domainP> &in,
             const BootstrappingKeyFFT<brP> &bkfft,
             const AnnihilateKey<typename brP::targetP> &ahk,
             const Polynomial<typename brP::targetP> &testvector,
             const std::vector<int> &betas, const std::vector<int> &Ts,
             const std::int64_t t,
             const int period = 0,
             BlindRotatePruneStats *stats = nullptr)
{
    if (betas.size() != Ts.size())
        throw std::invalid_argument("MetaPBS: betas/Ts size mismatch");

    constexpr int N = static_cast<int>(brP::targetP::n);
    const int twoN = 2 * N;
    if ((t & 1) != 0 || t <= 0 || t > twoN || (twoN % t) != 0)
        throw std::invalid_argument("MetaPBS: bad t");

    int r = twoN / static_cast<int>(t);

    const auto lifted = detail::LiftTLWEToInt<typename brP::domainP, Int>(in);
    auto [cquo, crem] =
        detail::HomDivRem<typename brP::domainP, Int>(lifted, static_cast<Int>(twoN));

    // First BlindRotate: use periodic pruning if period > 1.
    TRLWE<typename brP::targetP> C;
    if (period > 1) {
        detail::BlindRotatePeriodic<brP, Int>(C, cquo, bkfft, testvector,
                                             period, stats);
    } else {
        detail::BlindRotate<brP, Int>(C, cquo, bkfft, testvector);
    }

    for (size_t i = 0; i < betas.size(); i++) {
        const int beta = betas[i];
        const int T = Ts[i];
        if (beta < 2) throw std::invalid_argument("MetaPBS: beta must be >=2");

        TRLWE<typename brP::targetP> Ctr;
        detail::TruncRepeat<typename brP::targetP>(Ctr, C, T, beta, ahk);
        TRLWE<typename brP::targetP> Cprime;
        detail::TRLWEMulByXai<typename brP::targetP>(Cprime, Ctr,
                                                     detail::DeltaRB(r, beta));

        if ((crem.modulus % static_cast<Int>(beta)) != 0)
            throw std::invalid_argument("MetaPBS: beta must divide remainder modulus");
        auto [cquo_next, crem_next] =
            detail::HomDivRem<typename brP::domainP, Int>(crem, static_cast<Int>(beta));

        TRLWE<typename brP::targetP> Cnext;
        detail::BlindRotate<brP, Int>(Cnext, cquo_next, bkfft, Cprime);

        C = Cnext;
        crem = crem_next;
        r *= beta;
    }

    SampleExtractIndex<typename brP::targetP>(out, C, 0);
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  Public Interface 2: ExtractBit
//
//  Extract a single bit from an input ciphertext in Z_{2N} and put it back
//  into its original position (in-place semantics).
//
//  Output: LWE ciphertext encrypting (bit_value * 2^bit_lsb) in Z_{2N},
//  where bit_lsb = nbit - 1 - bit_msb.
//
//  Implementation (high-precision via Lvl2):
//    1) MetaPBSExtractBit2N: periodic non-redundant LUT → {0, q/2}
//    2) IdentityKeySwitch: lvl1 → brP_weight::domainP
//    3) GateBootstrapping to lvl2 with sign LUT → {-w/2, +w/2} + shift → {0, w}
//    4) IdentityKeySwitch: lvl2 → lvl1
//
//  Parameters:
//    out          — output TLWE ciphertext (at brP_extract::targetP)
//    in           — input TLWE ciphertext (at brP_extract::domainP)
//    bkfft_extract— bootstrapping key for the extraction stage
//    ahk          — annihilate key
//    ksk_to_dom   — key switching key: targetP → brP_weight::domainP
//    bkfft_weight — bootstrapping key for the weighting stage (→ lvl2)
//    ksk_down     — key switching key: lvl2 → targetP
//    bit_msb      — bit index (MSB=0, i.e., 0 = most significant bit)
//    prune_stats  — optional pointer to collect pruning statistics
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

template <class brP_extract, class brP_weight, class ksP_to_dom, class ksP_down,
          class Int = std::int64_t>
void ExtractBit(TLWE<typename brP_extract::targetP> &out,
                const TLWE<typename brP_extract::domainP> &in,
                const BootstrappingKeyFFT<brP_extract> &bkfft_extract,
                const AnnihilateKey<typename brP_extract::targetP> &ahk,
                const KeySwitchingKey<ksP_to_dom> &ksk_to_dom,
                const BootstrappingKeyFFT<brP_weight> &bkfft_weight,
                const KeySwitchingKey<ksP_down> &ksk_down,
                const int bit_msb,
                BlindRotatePruneStats *prune_stats = nullptr)
{
    using OutP = typename brP_extract::targetP;
    using DomP = typename brP_weight::domainP;
    using WeightP = typename brP_weight::targetP;

    static_assert(std::is_same_v<typename ksP_to_dom::domainP, OutP>,
                  "ExtractBit: ksP_to_dom::domainP must be targetP");
    static_assert(std::is_same_v<typename ksP_to_dom::targetP, DomP>,
                  "ExtractBit: ksP_to_dom::targetP must match brP_weight::domainP");
    static_assert(std::is_same_v<typename ksP_down::domainP, WeightP>,
                  "ExtractBit: ksP_down::domainP must match brP_weight::targetP");
    static_assert(std::is_same_v<typename ksP_down::targetP, OutP>,
                  "ExtractBit: ksP_down::targetP must be targetP");

    constexpr int nbit = static_cast<int>(OutP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument("ExtractBit: bad bit index (MSB=0)");

    const int bit_lsb = (nbit - 1) - bit_msb;
    const std::uint32_t weight = 1u << bit_lsb;

    // Stage 1: periodic bit extraction => {0, q/2} in targetP.
    TLWE<OutP> bit_ct{};
    detail::MetaPBSExtractBit2N<brP_extract, Int>(bit_ct, in, bkfft_extract, ahk,
                                                  bit_lsb, prune_stats);

    // Stage 2: key switch to the weight-bootstrap domain.
    TLWE<DomP> bit_dom{};
    IdentityKeySwitch<ksP_to_dom>(bit_dom, bit_ct, ksk_to_dom);

    // Stage 2.5: center away from the torus wrap boundary at 0.
    //   0   -> -q/4
    //   q/2 -> +q/4
    constexpr int dom_digits = std::numeric_limits<typename DomP::T>::digits;
    constexpr auto q_over_4 =
        static_cast<typename DomP::T>(static_cast<typename DomP::T>(1)
                                      << (dom_digits - 2));
    bit_dom[DomP::k * DomP::n] -= q_over_4;

    // Stage 3: bootstrap to lvl2 with a sign LUT producing {-w/2,+w/2},
    // then shift by +w/2 to get {0,w}.
    constexpr int out_plain_bits = nbit + 1;
    constexpr int w_digits = std::numeric_limits<typename WeightP::T>::digits;
    constexpr int delta_shift = w_digits - out_plain_bits;
    static_assert(delta_shift >= 1,
                  "ExtractBit: invalid torus width / nbit");

    const auto w_half = static_cast<typename WeightP::T>(
        static_cast<typename WeightP::T>(weight) << (delta_shift - 1));

    Polynomial<WeightP> tv{};
    tv.fill(w_half);

    TLWE<WeightP> w_ct{};
    GateBootstrappingTLWE2TLWE<brP_weight>(w_ct, bit_dom, bkfft_weight, tv);
    w_ct[WeightP::k * WeightP::n] += w_half;

    // Stage 4: key switch back down to targetP.
    IdentityKeySwitch<ksP_down>(out, w_ct, ksk_down);
}

}  // namespace TFHEpp::metapbs
