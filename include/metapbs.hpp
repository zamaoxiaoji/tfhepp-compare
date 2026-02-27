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

}  // namespace TFHEpp::metapbs
