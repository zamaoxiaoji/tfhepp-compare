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

// Symmetric interval [r]_sym used in Meta-PBS paper:
//   [r]_sym = Z ∩ [-floor(r/2), ceil(r/2)-1], with the special case [2]_sym={0,1}.
// Examples:
//   r=8 -> {-4,-3,-2,-1,0,1,2,3}, r=1 -> {0}, r=2 -> {0,1}.
constexpr int sym_min(const int r) { return (r == 2) ? 0 : -(r / 2); }
constexpr int sym_max(const int r) { return sym_min(r) + r - 1; }

// Centered remainder: returns x mod m in [-floor(m/2), floor((m-1)/2)].
template <class Int>
constexpr Int centered_mod(Int x, const Int m)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    // C++ % keeps the sign of x; adjust into the centered interval.
    Int r = x % m;
    const Int half_floor = m / 2;
    const Int half_ceil_minus1 = (m - 1) / 2;
    if (r < -half_floor) r += m;
    if (r > half_ceil_minus1) r -= m;
    return r;
}

template <class DomainP, class Int = std::int64_t>
struct LWE {
    std::array<Int, DomainP::k * DomainP::n + 1> c{};
    Int modulus = 0;  // Current modulus q (positive).
};

// Lift a torus TLWE ciphertext (mod 2^digits) into a signed integer vector in
// [-q/2, q/2).
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

// Homomorphic approximate Euclidean division (Definition 1 in Meta-PBS paper):
// for input modulus q and quotient modulus q_quo (q_quo | q),
// returns (c_quo, c_rem) s.t. c = (q/q_quo)*c_quo + c_rem and
// c_rem coefficients are centered mod (q/q_quo).
template <class DomainP, class Int = std::int64_t>
std::pair<LWE<DomainP, Int>, LWE<DomainP, Int>> HomDivRem(
    const LWE<DomainP, Int> &in, const Int q_quo)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    // Require q_quo | q.
    // (Paper uses q' | q; here q_quo plays the role of q'.)
    const Int q = in.modulus;
    if (q_quo <= 0 || q <= 0 || (q % q_quo) != 0) {
        throw std::invalid_argument("HomDivRem: require q_quo | q and q_quo>0");
    }
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

template <class TargetP, class Int>
constexpr int mod2N(const Int x)
{
    static_assert(std::is_signed_v<Int>, "Int must be signed");
    constexpr Int twoN = static_cast<Int>(2) * static_cast<Int>(TargetP::n);
    Int r = x % twoN;
    if (r < 0) r += twoN;
    return static_cast<int>(r);
}

// Δ_{r,B} from Theorem 1: used to re-center the duplicated coefficients.
constexpr int DeltaRB(const int r, const int B)
{
    return sym_min(r * B) - B * sym_min(r) - sym_min(B);
}

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

struct BlindRotatePruneStats;

template <class brP, class Int>
void BlindRotatePeriodic(TRLWE<typename brP::targetP> &res,
                         const LWE<typename brP::domainP, Int> &tlwe_quo,
                         const BootstrappingKeyFFT<brP> &bkfft,
                         const Polynomial<typename brP::targetP> &testvector,
                         const int period, BlindRotatePruneStats *stats);

// Homomorphic truncRepeat(·, [-T,T], B) (Definition 5) via
// SampleExtract + LWE-to-RLWE packing.
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

    // 1) SampleExtract the Laurent coefficients in [-T, T].
    std::vector<TLWE<P>> coeffs(2 * T + 1);
    for (int d = -T; d <= T; d++) {
        const int idx = (d >= 0) ? d : (static_cast<int>(P::n) + d);
        SampleExtractIndex<P>(coeffs[d + T], in, idx);
        if (d < 0) TLWENegate<P>(coeffs[d + T]);  // Lift: coeff(-d) = -M[N-d]
    }

    // 2) Build two contiguous coefficient segments: degrees >=0 (pos) and <0
    // (neg) after Embed.
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
                const int p = out_deg + neg_len;  // out_deg in [-neg_len, -1]
                TLWE<P> tmp = cd;
                TLWENegate<P>(tmp);  // Embed: X^{out_deg} -> -X^{N+out_deg}
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

// Generate a TFHE-style test vector for a negacyclic function f over Z_t.
// The caller provides f(i) for i in [0, t/2), represented as signed integers
// in [-t/2, t/2) (centered mod t).
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

    const int block = twoN / static_cast<int>(t);  // 2N/t
    const int shift = block / 2;                   // ≈ N/t, integer for even block

    Polynomial<TargetP> tv{};
    for (auto &x : tv) x = 0;

    for (int i = 0; i < static_cast<int>(t / 2); i++) {
        const auto y = EncodeTorus<TargetP, SignedInt>(f_half[i], t);
        for (int j = 0; j < block; j++) {
            int e = i * block + j - shift;
            e %= twoN;
            if (e < 0) e += twoN;
            if (e >= N)
                tv[e - N] -= y;  // X^{N+k} = -X^k
            else
                tv[e] += y;
        }
    }
    return tv;
}

// Same as GenerateTestVectorNegacyclic, but the caller provides the first-half
// LUT values directly as torus elements (no EncodeTorus scaling).
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

    const int block = twoN / static_cast<int>(t);  // 2N/t
    const int shift = block / 2;                   // ≈ N/t, integer for even block

    Polynomial<TargetP> tv{};
    tv.fill(0);

    for (int i = 0; i < static_cast<int>(t / 2); i++) {
        const auto y = f_half_torus[i];
        for (int j = 0; j < block; j++) {
            int e = i * block + j - shift;
            e %= twoN;
            if (e < 0) e += twoN;
            if (e >= N)
                tv[e - N] -= y;  // X^{N+k} = -X^k
            else
                tv[e] += y;
        }
    }
    return tv;
}

// Build a non-redundant (t=2N) negacyclic test vector for extracting bit "bit"
// (LSB=0) from m in Z_{2N}. Output encoding is 0 or q/2 (torus half-turn).
//
// NOTE: This only supports bit < log2(N) (i.e., bit < TargetP::nbit). The MSB
// that distinguishes [0,N) vs [N,2N) is not negacyclic in general.
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

// Build a non-redundant (t=2N) negacyclic test vector for extracting bit
// "bit_lsb" (LSB=0) and putting it back in-place, i.e., output plaintext is
// either 0 or 2^{bit_lsb}.
//
// IMPORTANT: Because this is a negacyclic LUT, this matches "in-place bit"
// semantics for messages in [0, N). Outside that range the second half follows
// the forced relation f(x+N) = -f(x).
template <class TargetP>
const Polynomial<TargetP> &BitExtractInPlaceTestVector2N(const int bit_lsb)
{
    constexpr int nbit = static_cast<int>(TargetP::nbit);
    if (bit_lsb < 0 || bit_lsb >= nbit)
        throw std::invalid_argument(
            "BitExtractInPlaceTestVector2N: bad bit index");

    struct Cache {
        std::array<Polynomial<TargetP>, TargetP::nbit> tv{};
        std::array<bool, TargetP::nbit> inited{};
    };
    static Cache cache;

    if (!cache.inited[bit_lsb]) {
        constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
        constexpr int delta_shift = digits - (nbit + 1);  // q / (2N)
        static_assert(delta_shift >= 1,
                      "BitExtractInPlaceTestVector2N: invalid torus width / nbit");

        Polynomial<TargetP> tv{};
        tv.fill(0);

        const auto w_torus = static_cast<typename TargetP::T>(
            static_cast<typename TargetP::T>(1)
            << (delta_shift + bit_lsb));  // (2^bit_lsb) * q / (2N)

        constexpr int N = static_cast<int>(TargetP::n);
        for (int i = 0; i < N; i++)
            tv[i] = ((i >> bit_lsb) & 1) ? w_torus : 0;

        cache.tv[bit_lsb] = tv;
        cache.inited[bit_lsb] = true;
    }
    return cache.tv[bit_lsb];
}

// Algorithm 1 (Meta-PBS) for a negacyclic LUT (single-output).
template <class brP, class Int = std::int64_t>
void MetaPBS(TLWE<typename brP::targetP> &out,
             const TLWE<typename brP::domainP> &in,
             const BootstrappingKeyFFT<brP> &bkfft,
             const AnnihilateKey<typename brP::targetP> &ahk,
             const Polynomial<typename brP::targetP> &testvector,
             const std::vector<int> &betas, const std::vector<int> &Ts,
             const std::int64_t t)
{
    if (betas.size() != Ts.size())
        throw std::invalid_argument("MetaPBS: betas/Ts size mismatch");

    constexpr int N = static_cast<int>(brP::targetP::n);
    const int twoN = 2 * N;
    if ((t & 1) != 0 || t <= 0 || t > twoN || (twoN % t) != 0)
        throw std::invalid_argument("MetaPBS: bad t");

    // r0 = 2N/t (redundancy of the initial test vector around the target entry).
    int r = twoN / static_cast<int>(t);

    // (c_quo,0, c_rem,0) <- HomDivRem(c, 2N)
    const auto lifted = LiftTLWEToInt<typename brP::domainP, Int>(in);
    auto [cquo, crem] =
        HomDivRem<typename brP::domainP, Int>(lifted, static_cast<Int>(twoN));

    // C0 <- BlindRotate(c_quo,0, TV)
    TRLWE<typename brP::targetP> C;
    BlindRotate<brP, Int>(C, cquo, bkfft, testvector);

    for (size_t i = 0; i < betas.size(); i++) {
        const int beta = betas[i];
        const int T = Ts[i];
        if (beta < 2) throw std::invalid_argument("MetaPBS: beta must be >=2");

        // C'_i <- TruncRepeat(C_i, [-T_i,T_i], beta_i) * X^{Δ_{r_i,beta_i}}
        TRLWE<typename brP::targetP> Ctr;
        TruncRepeat<typename brP::targetP>(Ctr, C, T, beta, ahk);
        TRLWE<typename brP::targetP> Cprime;
        TRLWEMulByXai<typename brP::targetP>(Cprime, Ctr, DeltaRB(r, beta));

        // c_quo,i+1, c_rem,i+1 <- HomDivRem(c_rem,i, beta_i)
        if ((crem.modulus % static_cast<Int>(beta)) != 0)
            throw std::invalid_argument("MetaPBS: beta must divide remainder modulus");
        auto [cquo_next, crem_next] =
            HomDivRem<typename brP::domainP, Int>(crem, static_cast<Int>(beta));

        // C_{i+1} <- BlindRotate(c_quo,i+1, C'_i)
        TRLWE<typename brP::targetP> Cnext;
        BlindRotate<brP, Int>(Cnext, cquo_next, bkfft, Cprime);

        C = Cnext;
        crem = crem_next;
        r *= beta;
    }

    SampleExtractIndex<typename brP::targetP>(out, C, 0);
}

// Meta-PBS specialized for bit extraction from m in Z_{2N} using a periodic
// LUT (period 2^{bit+1}) with outputs in {0, q/2}.
//
// This wrapper enforces t=2N and constructs the corresponding non-redundant
// test vector internally. It also applies the periodic-CMUX pruning in the
// *first* BlindRotate (since the plaintext test vector is periodic), optionally
// collecting stats.
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

    // r0 = 2N/t = 1 for t=2N.
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

// Convenience overload with a sane default iteration schedule for t=2N.
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

// One-PBS, no-KS bit extraction "put back in-place" (MSB=0 indexing).
//
// Implementation:
// 1) Run the periodic t=2N extraction once: output in {0, q/2}.
// 2) Apply a centered power-of-two downscale to map q/2 -> weight_torus.
//
// For bit_lsb = nbit-1-bit_msb, we need:
//   weight_torus = 2^(digits-(nbit+1)+bit_lsb)
// while q/2 is 2^(digits-1), so shift = nbit-bit_lsb.
//
// This keeps the path "one PBS + linear ops", i.e., no additional PBS and no
// key-switching.
template <class brP, class Int = std::int64_t>
void ExtractBitInPlaceOnePBS(
    TLWE<typename brP::targetP> &out, const TLWE<typename brP::domainP> &in,
    const BootstrappingKeyFFT<brP> &bkfft,
    const AnnihilateKey<typename brP::targetP> &ahk, const int bit_msb)
{
    using TargetP = typename brP::targetP;
    using TorusT = typename TargetP::T;
    using SignedT = std::make_signed_t<TorusT>;

    constexpr int nbit = static_cast<int>(TargetP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument(
            "ExtractBitInPlaceOnePBS: bad bit index (MSB=0)");

    const int bit_lsb = (nbit - 1) - bit_msb;
    TLWE<TargetP> half_turn{};
    MetaPBSExtractBit2N<brP, Int>(half_turn, in, bkfft, ahk, bit_lsb);

    const int shift = nbit - bit_lsb;
    for (size_t i = 0; i < out.size(); i++) {
        using WideSignedT = std::conditional_t<
            (std::numeric_limits<TorusT>::digits <= 32), std::int64_t,
            __int128_t>;

        const WideSignedT sx =
            static_cast<WideSignedT>(static_cast<SignedT>(half_turn[i]));
        const WideSignedT add = static_cast<WideSignedT>(1) << (shift - 1);
        const WideSignedT sy =
            (sx >= 0) ? ((sx + add) >> shift) : -(((-sx) + add) >> shift);

        out[i] = static_cast<TorusT>(static_cast<SignedT>(-sy));
    }
}

// Rounds x / 2^shift to the nearest integer (ties round up).
// Implemented without relying on wider intermediates to avoid overflow for
// large torus types.
template <class TorusT>
constexpr TorusT RoundDivPow2(const TorusT x, const int shift)
{
    static_assert(std::is_unsigned_v<TorusT>, "TorusT must be unsigned");
    if (shift <= 0) return x;
    if (shift >= std::numeric_limits<TorusT>::digits) return 0;
    const TorusT hi = x >> shift;
    const TorusT round = (x >> (shift - 1)) & static_cast<TorusT>(1);
    return hi + round;
}

// Rounds x / 2^shift in the *centered* torus representation (two's complement),
// i.e., interpret x as a signed integer in [-q/2, q/2), then divide with
// rounding and map back to the torus.
template <class TorusT>
constexpr TorusT RoundDivPow2Centered(const TorusT x, const int shift)
{
    static_assert(std::is_unsigned_v<TorusT>, "TorusT must be unsigned");
    if (shift <= 0) return x;

    using SignedT = std::make_signed_t<TorusT>;
    // Wide enough for int32/int64 torus types.
    using WideSignedT = std::conditional_t<
        (std::numeric_limits<TorusT>::digits <= 32), std::int64_t, __int128_t>;
    const WideSignedT sx =
        static_cast<WideSignedT>(static_cast<SignedT>(x));
    const WideSignedT add =
        static_cast<WideSignedT>(1) << (shift - 1);

    WideSignedT sy = 0;
    if (sx >= 0)
        sy = (sx + add) >> shift;
    else
        sy = -(((-sx) + add) >> shift);

    return static_cast<TorusT>(static_cast<SignedT>(sy));
}

// Extract a plaintext bit (indexed from MSB=0) from an input ciphertext in
// Z_{2N} using the non-redundant periodic LUT. The output encoding is always
// {0, q/2} (torus half-turn).
//
// How "weight" applies:
// - Do NOT scale the periodic LUT coefficients by weight, since that breaks the
//   periodic no-op rotations across the X^N = -1 sign flip.
// - Instead, interpret the returned TLWE as an integer ciphertext modulo
//   t_out = 2*weight. Under TFHEpp's centered integer decode,
//   q/2 corresponds to -t_out/2 = -weight, whose (k+1)-bit two's-complement bit
//   pattern is exactly `1 << k` (e.g., 10000 for weight=1<<4).
//
// In other words, this interface returns the extracted bit as the MSB of a
// (k+1)-bit signed integer ring, without changing the periodic LUT.
template <class brP, class Int = std::int64_t>
void ExtractBitHalfTurn(TLWE<typename brP::targetP> &out,
                        const TLWE<typename brP::domainP> &in,
                        const BootstrappingKeyFFT<brP> &bkfft,
                        const AnnihilateKey<typename brP::targetP> &ahk,
                        const int bit_msb, const std::uint32_t weight,
                        BlindRotatePruneStats *prune_stats = nullptr)
{
    constexpr int nbit = static_cast<int>(brP::targetP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument("ExtractBitHalfTurn: bad bit index (MSB=0)");
    if (weight == 0 || (weight & (weight - 1)) != 0)
        throw std::invalid_argument(
            "ExtractBitHalfTurn: weight must be a non-zero power of two");
    if (weight > (1u << (nbit - 1)))
        throw std::invalid_argument("ExtractBitHalfTurn: weight too large");

    const int bit_lsb = (nbit - 1) - bit_msb;
    MetaPBSExtractBit2N<brP, Int>(out, in, bkfft, ahk, bit_lsb, prune_stats);
}

// Weighted ExtractBit interface (hi-precision post-step).
//
// For the periodic non-redundant extraction, the output is {0,q/2}. Converting
// this to an *exact* 11-bit integer (mod 2N) at lvl1 can be tight under typical
// TFHE bootstrapping noise. This helper does the weighting step by
// bootstrapping to a higher-precision target (typically lvl2param / 64-bit
// torus), then key-switches back to lvl1.
//
// Requirements on parameters:
// - brP_extract: domainP -> lvl1 (extract stage, Meta-PBS periodic LUT)
// - ksP_to_dom: lvl1 -> brP_weight::domainP (to feed the weight bootstrap)
// - brP_weight: brP_weight::domainP -> lvl2 (weight bootstrap stage)
// - ksP_down: lvl2 -> lvl1 (final key switch)
template <class brP_extract, class brP_weight, class ksP_to_dom, class ksP_down,
          class Int = std::int64_t>
void ExtractBitViaLvl2(TLWE<typename brP_extract::targetP> &out,
                       const TLWE<typename brP_extract::domainP> &in,
                       const BootstrappingKeyFFT<brP_extract> &bkfft_extract,
                       const AnnihilateKey<typename brP_extract::targetP> &ahk,
                       const KeySwitchingKey<ksP_to_dom> &ksk_to_dom,
                       const BootstrappingKeyFFT<brP_weight> &bkfft_weight,
                       const KeySwitchingKey<ksP_down> &ksk_down,
                       const int bit_msb, const std::uint32_t weight,
                       BlindRotatePruneStats *prune_stats = nullptr)
{
    using OutP = typename brP_extract::targetP;
    using DomP = typename brP_weight::domainP;
    using WeightP = typename brP_weight::targetP;

    static_assert(std::is_same_v<typename ksP_to_dom::domainP, OutP>,
                  "ExtractBitViaLvl2: ksP_to_dom::domainP must be lvl1");
    static_assert(std::is_same_v<typename ksP_to_dom::targetP, DomP>,
                  "ExtractBitViaLvl2: ksP_to_dom::targetP must match brP_weight::domainP");
    static_assert(std::is_same_v<typename ksP_down::domainP, WeightP>,
                  "ExtractBitViaLvl2: ksP_down::domainP must match brP_weight::targetP");
    static_assert(std::is_same_v<typename ksP_down::targetP, OutP>,
                  "ExtractBitViaLvl2: ksP_down::targetP must be lvl1");

    constexpr int nbit = static_cast<int>(OutP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument("ExtractBitViaLvl2: bad bit index (MSB=0)");
    if (weight == 0 || (weight & (weight - 1)) != 0)
        throw std::invalid_argument(
            "ExtractBitViaLvl2: weight must be a non-zero power of two");
    if (weight > (1u << (nbit - 1)))
        throw std::invalid_argument("ExtractBitViaLvl2: weight too large");

    // Stage 1: periodic bit extraction => {0, q/2} in lvl1.
    TLWE<OutP> bit_ct{};
    ExtractBitHalfTurn<brP_extract, Int>(bit_ct, in, bkfft_extract, ahk, bit_msb,
                                         weight, prune_stats);

    // Stage 2: key switch to the weight-bootstrap domain (typically lvlhalf).
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

    // Stage 3: bootstrap to lvl2 with a sign LUT producing {-w/2,+w/2}, then
    // shift by +w/2 to get {0,w}.
    constexpr int out_plain_bits = nbit + 1;  // plain_modulus = 2N
    constexpr int w_digits = std::numeric_limits<typename WeightP::T>::digits;
    constexpr int delta_shift = w_digits - out_plain_bits;
    static_assert(delta_shift >= 1,
                  "ExtractBitViaLvl2: invalid torus width / nbit");

    const auto w_half = static_cast<typename WeightP::T>(
        static_cast<typename WeightP::T>(weight) << (delta_shift - 1));

    Polynomial<WeightP> tv{};
    tv.fill(w_half);

    TLWE<WeightP> w_ct{};
    GateBootstrappingTLWE2TLWE<brP_weight>(w_ct, bit_dom, bkfft_weight, tv);
    w_ct[WeightP::k * WeightP::n] += w_half;

    // Stage 4: key switch back down to lvl1.
    IdentityKeySwitch<ksP_down>(out, w_ct, ksk_down);
}

// Convenience wrapper: put extracted bit back into its original position.
template <class brP_extract, class brP_weight, class ksP_to_dom, class ksP_down,
          class Int = std::int64_t>
void ExtractBitInPlaceViaLvl2(TLWE<typename brP_extract::targetP> &out,
                              const TLWE<typename brP_extract::domainP> &in,
                              const BootstrappingKeyFFT<brP_extract> &bkfft_extract,
                              const AnnihilateKey<typename brP_extract::targetP> &ahk,
                              const KeySwitchingKey<ksP_to_dom> &ksk_to_dom,
                              const BootstrappingKeyFFT<brP_weight> &bkfft_weight,
                              const KeySwitchingKey<ksP_down> &ksk_down,
                              const int bit_msb,
                              BlindRotatePruneStats *prune_stats = nullptr)
{
    constexpr int nbit = static_cast<int>(brP_extract::targetP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument(
            "ExtractBitInPlaceViaLvl2: bad bit index (MSB=0)");

    const int bit_lsb = (nbit - 1) - bit_msb;
    const std::uint32_t weight = 1u << bit_lsb;

    ExtractBitViaLvl2<brP_extract, brP_weight, ksP_to_dom, ksP_down, Int>(
        out, in, bkfft_extract, ahk, ksk_to_dom, bkfft_weight, ksk_down,
        bit_msb, weight, prune_stats);
}

// Weighted ExtractBit interface.
//
// This returns an LWE ciphertext encrypting (bit * weight) in Z_{2N}, i.e., the
// same integer encoding used by tlweSymIntEncrypt/tlweSymIntDecrypt with
// plain_modulus = 2N (TargetP::nbit+1 bits).
//
// Internals:
// 1) ExtractBitHalfTurn: periodic non-redundant LUT => {0, q/2}.
// 2) IdentityKeySwitch: bring the extracted {0,q/2} bit to the chosen bootstrap
//    domain for the weighting stage.
// 3) Meta-PBS with t=2N that maps {0,q/2} -> {-w/2,+w/2}, then add +w/2
//    (and a small bias) to shift to {0, w}, where w = weight * (q / 2N).
//
// Notes:
// - Do NOT scale the periodic LUT coefficients by weight; the post-PBS weighting
//   avoids breaking the periodic pruning across the X^N = -1 sign flip.
// - "ksP" must satisfy: ksP::domainP == brP_extract::targetP and
//   ksP::targetP == brP_weight::domainP.
template <class brP_extract, class brP_weight, class ksP,
          class Int = std::int64_t>
void ExtractBit(TLWE<typename brP_extract::targetP> &out,
                const TLWE<typename brP_extract::domainP> &in,
                const BootstrappingKeyFFT<brP_extract> &bkfft_extract,
                const AnnihilateKey<typename brP_extract::targetP> &ahk,
                const KeySwitchingKey<ksP> &ksk,
                const BootstrappingKeyFFT<brP_weight> &bkfft_weight,
                const int bit_msb, const std::uint32_t weight,
                BlindRotatePruneStats *prune_stats = nullptr)
{
    static_assert(std::is_same_v<typename brP_extract::targetP,
                                 typename brP_weight::targetP>,
                  "ExtractBit: targetP mismatch between extract and weight stages");
    static_assert(std::is_same_v<typename ksP::domainP,
                                 typename brP_extract::targetP>,
                  "ExtractBit: ksP::domainP must match brP_extract::targetP");
    static_assert(std::is_same_v<typename ksP::targetP,
                                 typename brP_weight::domainP>,
                  "ExtractBit: ksP::targetP must match brP_weight::domainP");

    constexpr int nbit = static_cast<int>(brP_extract::targetP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument("ExtractBit: bad bit index (MSB=0)");
    if (weight == 0 || (weight & (weight - 1)) != 0)
        throw std::invalid_argument(
            "ExtractBit: weight must be a non-zero power of two");
    if (weight > (1u << (nbit - 1)))
        throw std::invalid_argument("ExtractBit: weight too large");

    // Stage 1: periodic bit extraction => {0, q/2}.
    TLWE<typename brP_extract::targetP> bit_ct{};
    ExtractBitHalfTurn<brP_extract, Int>(bit_ct, in, bkfft_extract, ahk,
                                         bit_msb, weight, prune_stats);

    // Stage 2: key switch bit_ct to the bootstrap domain for the t=2N stage.
    TLWE<typename brP_weight::domainP> bit_dom{};
    IdentityKeySwitch<ksP>(bit_dom, bit_ct, ksk);

    // Stage 2.5: center {0,q/2} away from the torus wrap boundary at 0.
    //   0      -> -q/4
    //   q/2    -> +q/4
    // This makes the next LUT robust: small noise around 0 would otherwise
    // wrap to q-ε and look like an index near 2N.
    {
        using DomP = typename brP_weight::domainP;
        constexpr int dom_digits =
            std::numeric_limits<typename DomP::T>::digits;
        constexpr auto q_over_4 =
            static_cast<typename DomP::T>(static_cast<typename DomP::T>(1)
                                          << (dom_digits - 2));
        bit_dom[DomP::k * DomP::n] -= q_over_4;
    }

    // Stage 3: bootstrap (full 2N LUT) to map {-q/4,+q/4} -> {-w/2,+w/2},
    // then shift by +w/2.
    //
    // Rationale:
    // - Stage-1 returns a torus half-turn: {0, q/2}.
    //   Under integer decoding modulo 2N, this is {0, -N} == {0, N}.
    // - Using a t=2 PBS would only see parity, but both 0 and N are even.
    // - With t=2N, the two possible inputs differ by N, so a negacyclic LUT
    //   can naturally output opposite signs (-w/2 vs +w/2), which we then shift
    //   into {0, w}.
    using TargetP = typename brP_extract::targetP;
    constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
    constexpr int delta_shift = digits - (static_cast<int>(TargetP::nbit) + 1);
    static_assert(delta_shift >= 1,
                  "ExtractBit: invalid target torus width / nbit");

    const auto w_half = static_cast<typename TargetP::T>(
        static_cast<typename TargetP::T>(weight) << (delta_shift - 1));

    // g(x) = +w_half for x in [0,N), and g(x) = -w_half for x in [N,2N).
    // With the -q/4 centering above, bit=1 lands in the first half (q/4) and
    // bit=0 lands in the second half (3q/4), so this produces {-w_half,+w_half}
    // for {0,1} respectively.
    constexpr std::int64_t t_nr = 2 * static_cast<std::int64_t>(TargetP::n);

    Polynomial<TargetP> tv{};
    tv.fill(w_half);

    // Iteration schedule:
    // We need betas that divide the current remainder modulus (q / 2N) at each
    // step. For a 16-bit domain torus (lvl0), q/2N is only 32, so {8,4} works
    // but {8,8} would not. For 32-bit domains, {8,8} matches the defaults used
    // elsewhere in this file.
    std::vector<int> betas;
    std::vector<int> Ts;
    {
        using DomP = typename brP_weight::domainP;
        constexpr int dom_digits =
            std::numeric_limits<typename DomP::T>::digits;
        constexpr int rem_pow0 =
            dom_digits - (static_cast<int>(TargetP::nbit) + 1);
        static_assert(rem_pow0 >= 0, "ExtractBit: bad domain/target bit widths");

        if constexpr (rem_pow0 >= 6) {
            betas = {8, 8};
        }
        else if constexpr (rem_pow0 == 5) {
            betas = {8, 4};
        }
        else if constexpr (rem_pow0 == 4) {
            betas = {4, 4};
        }
        else if constexpr (rem_pow0 == 3) {
            betas = {8};
        }
        else if constexpr (rem_pow0 == 2) {
            betas = {4};
        }
        else if constexpr (rem_pow0 == 1) {
            betas = {2};
        }
        else {
            betas = {};
        }

        Ts.reserve(betas.size());
        for (const int beta : betas)
            Ts.push_back((static_cast<int>(TargetP::n) / beta - 1) / 2);
    }

    MetaPBS<brP_weight, Int>(out, bit_dom, bkfft_weight, ahk, tv, betas, Ts,
                             /*t=*/t_nr);
    // Shift {-w_half,+w_half} -> {0,2*w_half} which equals {0, weight} under
    // integer encoding modulo 2N.
    out[TargetP::k * TargetP::n] += w_half;
}

// Convenience wrapper for "put the extracted bit back to its original place":
// weight is chosen as 1<<bit_lsb (with bit indexed from MSB=0).
template <class brP_extract, class brP_weight, class ksP,
          class Int = std::int64_t>
void ExtractBitInPlace(TLWE<typename brP_extract::targetP> &out,
                       const TLWE<typename brP_extract::domainP> &in,
                       const BootstrappingKeyFFT<brP_extract> &bkfft_extract,
                       const AnnihilateKey<typename brP_extract::targetP> &ahk,
                       const KeySwitchingKey<ksP> &ksk,
                       const BootstrappingKeyFFT<brP_weight> &bkfft_weight,
                       const int bit_msb,
                       BlindRotatePruneStats *prune_stats = nullptr)
{
    constexpr int nbit = static_cast<int>(brP_extract::targetP::nbit);
    if (bit_msb < 0 || bit_msb >= nbit)
        throw std::invalid_argument(
            "ExtractBitInPlace: bad bit index (MSB=0)");

    const int bit_lsb = (nbit - 1) - bit_msb;
    const std::uint32_t weight = 1u << bit_lsb;

    ExtractBit<brP_extract, brP_weight, ksP, Int>(
        out, in, bkfft_extract, ahk, ksk, bkfft_weight, bit_msb, weight,
        prune_stats);
}

// Meta-PBS BlindRotate that consumes the quotient ciphertext from HomDivRem,
// i.e., coefficients are already in Z (no torus modulus switching needed).
// The exponent is interpreted modulo 2N (N = targetP::n).
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

struct BlindRotatePruneStats {
    std::uint64_t total = 0;    // CMUXes that would have run (a != 0).
    std::uint64_t skipped = 0;  // CMUXes skipped because the rotation is a no-op.
};

// BlindRotate variant that can prune CMUXes when the (signed) 2N coefficient
// representation of the test vector is periodic with period "period".
// In that case, Rot_a(TV) == TV for any a ≡ 0 (mod period), so CMUX becomes
// an unconditional identity map and can be skipped.
//
// NOTE: This is only correct if the caller guarantees the periodicity condition
// for the provided plaintext testvector in R_q = Z_q[X]/(X^N+1).
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

}  // namespace TFHEpp::metapbs
