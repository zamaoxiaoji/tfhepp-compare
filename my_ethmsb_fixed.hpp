#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

#ifdef USE_KEY_BUNDLE
#include "keybundle.hpp"
#endif
#include "evalkeygens.hpp"
#include "gatebootstrapping.hpp"
#include "key.hpp"
#include "params.hpp"
#include "tlwe.hpp"

namespace my_ethmsb {

using P = TFHEpp::lvl2param;
using Torus = typename P::T;
using Wide = unsigned __int128;
using SignedWide = __int128;

static_assert(std::numeric_limits<Torus>::digits == 64);

struct my_lvl22pbsparam {
    using domainP = P;
    using targetP = P;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

inline constexpr Torus BOOL_ONE = Torus{1} << 61;

inline void require_k(const int k)
{
    if (k < 1 || k > 63)
        throw std::invalid_argument("ETHMSB supports 1 <= k <= 63");
}

inline void require_kappa(const int kappa)
{
    if (kappa < 1 || kappa > 62)
        throw std::invalid_argument("ETHMSB supports 1 <= kappa <= 62");
}

inline Torus delta(const int k)
{
    require_k(k);
    return Torus{1} << (64 - k);
}

inline Torus mask_bits(const int k)
{
    require_k(k);
    return (Torus{1} << k) - Torus{1};
}

inline Torus encode_unsigned(const Torus m, const int k)
{
    require_k(k);
    return static_cast<Torus>(Wide{m} * Wide{delta(k)});
}

inline Torus decode_unsigned_phase(const Torus phase, const int k)
{
    require_k(k);
    const Torus rounded = phase + delta(k) / 2;
    return (rounded >> (64 - k)) & mask_bits(k);
}

inline Torus torus_abs_centered(const Torus x)
{
    constexpr Torus half = Torus{1} << 63;
    return x <= half ? x : Torus{0} - x;
}

template <class Param>
inline bool decrypt_arith_bit(const TFHEpp::TLWE<Param>& ct,
                              const TFHEpp::SecretKey& sk,
                              const typename Param::T out_value)
{
    static_assert(std::is_same_v<Param, P>,
                  "my_ethmsb arithmetic-bit decryption is lvl2-only");
    const Torus phase = TFHEpp::tlweSymPhase<Param>(ct, sk.key.get<Param>());
    const Torus dist_one = torus_abs_centered(phase - out_value);
    const Torus dist_zero = torus_abs_centered(phase);
    return dist_one < dist_zero;
}

inline bool decrypt_arith_bit(const TFHEpp::TLWE<P>& ct,
                              const TFHEpp::SecretKey& sk,
                              const Torus out_value = BOOL_ONE)
{
    return decrypt_arith_bit<P>(ct, sk, out_value);
}

template <class Param>
inline void add_const_inplace(TFHEpp::TLWE<Param>& ct,
                              const typename Param::T c)
{
    ct[Param::k * Param::n] += c;
}

template <class Param>
inline TFHEpp::TLWE<Param> add_const(TFHEpp::TLWE<Param> ct,
                                     const typename Param::T c)
{
    add_const_inplace<Param>(ct, c);
    return ct;
}

template <class Param>
inline void sub(TFHEpp::TLWE<Param>& out, const TFHEpp::TLWE<Param>& a,
                const TFHEpp::TLWE<Param>& b)
{
    for (int i = 0; i <= Param::k * Param::n; ++i) out[i] = a[i] - b[i];
}

template <class Param>
inline TFHEpp::TLWE<Param> sub(const TFHEpp::TLWE<Param>& a,
                               const TFHEpp::TLWE<Param>& b)
{
    TFHEpp::TLWE<Param> out;
    sub<Param>(out, a, b);
    return out;
}

template <class Param>
inline void add(TFHEpp::TLWE<Param>& out, const TFHEpp::TLWE<Param>& a,
                const TFHEpp::TLWE<Param>& b)
{
    for (int i = 0; i <= Param::k * Param::n; ++i) out[i] = a[i] + b[i];
}

template <class Param>
inline TFHEpp::TLWE<Param> add(const TFHEpp::TLWE<Param>& a,
                               const TFHEpp::TLWE<Param>& b)
{
    TFHEpp::TLWE<Param> out;
    add<Param>(out, a, b);
    return out;
}

template <class Param>
inline void scalar_mul_pow2(TFHEpp::TLWE<Param>& out,
                            const TFHEpp::TLWE<Param>& a, const int s)
{
    if (s < 0 || s >= std::numeric_limits<typename Param::T>::digits)
        throw std::invalid_argument("invalid torus power-of-two scale");
    for (int i = 0; i <= Param::k * Param::n; ++i) out[i] = a[i] << s;
}

template <class Param>
inline TFHEpp::TLWE<Param> scalar_mul_pow2(const TFHEpp::TLWE<Param>& a,
                                           const int s)
{
    TFHEpp::TLWE<Param> out;
    scalar_mul_pow2<Param>(out, a, s);
    return out;
}

template <class Param>
inline TFHEpp::TLWE<Param> trivial_constant(const typename Param::T value)
{
    TFHEpp::TLWE<Param> out = {};
    out[Param::k * Param::n] = value;
    return out;
}

template <class BRP, uint32_t num_out = 1>
void br_mod_switch_safe(TFHEpp::ModswitchTLWE<typename BRP::domainP>& moded,
                        const TFHEpp::TLWE<typename BRP::domainP>& tlwe)
{
    static_assert(std::is_same_v<typename BRP::domainP, P>);
    static_assert(std::is_same_v<typename BRP::targetP, P>);
    constexpr uint32_t bitwidth = TFHEpp::bits_needed<num_out - 1>();
    constexpr int phase_shift =
        std::numeric_limits<typename BRP::domainP::T>::digits - 1 -
        BRP::targetP::nbit + bitwidth;
    constexpr typename BRP::domainP::T roundoffset =
        typename BRP::domainP::T{1}
        << (std::numeric_limits<typename BRP::domainP::T>::digits - 2 -
            BRP::targetP::nbit + bitwidth);
    constexpr int recon_shift =
        std::numeric_limits<typename BRP::domainP::T>::digits - 1 -
        BRP::targetP::nbit;

    SignedWide correction = 0;
    for (int i = 0; i < BRP::domainP::k * BRP::domainP::n; ++i) {
        moded[i] = ((tlwe[i] + roundoffset) >> phase_shift) << bitwidth;
        const typename BRP::domainP::T recon = moded[i] << recon_shift;
        correction += static_cast<int64_t>(tlwe[i] - recon);
    }

    typename BRP::domainP::T corrected_body =
        tlwe[BRP::domainP::k * BRP::domainP::n];
    if constexpr (BRP::domainP::key_value_min == 0 &&
                  BRP::domainP::key_value_max == 1) {
        corrected_body -=
            static_cast<typename BRP::domainP::T>(correction / 2);
    }
    moded[BRP::domainP::k * BRP::domainP::n] =
        2 * BRP::targetP::n -
        (((corrected_body + roundoffset) >> phase_shift) << bitwidth);
}

template <class BRP, uint32_t num_out = 1>
void blind_rotate_safe(TFHEpp::TRLWE<typename BRP::targetP>& res,
                       const TFHEpp::TLWE<typename BRP::domainP>& tlwe,
                       const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
                       const TFHEpp::Polynomial<typename BRP::targetP>&
                           testvector)
{
    TFHEpp::ModswitchTLWE<typename BRP::domainP> moded;
    br_mod_switch_safe<BRP, num_out>(moded, tlwe);
    res = {};
    TFHEpp::PolynomialMulByXai<typename BRP::targetP>(
        res[BRP::targetP::k], testvector,
        moded[BRP::domainP::k * BRP::domainP::n]);
#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < BRP::domainP::k * BRP::domainP::n / BRP::Addends;
         ++i) {
        alignas(64) TFHEpp::TRGSWFFT<typename BRP::targetP> BKadded;
        TFHEpp::KeyBundleFFT<BRP>(
            BKadded, bkfft[i],
            std::span(moded)
                .subspan(BRP::Addends * i, BRP::Addends)
                .template first<BRP::Addends>());
        TFHEpp::trgswfftExternalProduct<typename BRP::targetP>(
            res, res, BKadded);
    }
#else
    for (int i = 0; i < BRP::domainP::k * BRP::domainP::n; ++i) {
        if (moded[i] == 0) continue;
        TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<BRP>(
            res, bkfft[i], static_cast<int>(moded[i]));
    }
#endif
}

template <class BRP>
void gate_bootstrapping_tlwe2tlwefft_safe(
    TFHEpp::TLWE<typename BRP::targetP>& res,
    const TFHEpp::TLWE<typename BRP::domainP>& tlwe,
    const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
    const TFHEpp::Polynomial<typename BRP::targetP>& testvector)
{
    alignas(64) TFHEpp::TRLWE<typename BRP::targetP> acc;
    blind_rotate_safe<BRP>(acc, tlwe, bkfft, testvector);
    TFHEpp::SampleExtractIndex<typename BRP::targetP>(res, acc, 0);
}

template <class BRP>
void pbs_msb_value(TFHEpp::TLWE<typename BRP::targetP>& out,
                   TFHEpp::TLWE<typename BRP::domainP> in,
                   const int input_bits, const Torus offset,
                   const Torus out_value,
                   const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft)
{
    static_assert(std::is_same_v<typename BRP::domainP, P>);
    static_assert(std::is_same_v<typename BRP::targetP, P>);
    require_k(input_bits);
    if ((out_value & Torus{1}) != 0)
        throw std::invalid_argument("pbs_msb_value requires even out_value");

    add_const_inplace<typename BRP::domainP>(in, offset);
    const Torus half = out_value / 2;
    TFHEpp::Polynomial<typename BRP::targetP> testvector;
    testvector.fill(Torus{0} - half);

    alignas(64) TFHEpp::TLWE<typename BRP::targetP> tmp;
#ifdef MY_ETHMSB_USE_CUSTOM_L22_BR
    gate_bootstrapping_tlwe2tlwefft_safe<BRP>(tmp, in, bkfft, testvector);
#else
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BRP>(tmp, in, bkfft, testvector);
#endif
    out = tmp;
    add_const_inplace<typename BRP::targetP>(out, half);
}

template <class BRP>
void ethmsb_value(TFHEpp::TLWE<typename BRP::targetP>& out,
                  const TFHEpp::TLWE<typename BRP::domainP>& ct, const int k,
                  const int kappa, const Torus out_value,
                  const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft)
{
    static_assert(std::is_same_v<typename BRP::domainP, P>);
    static_assert(std::is_same_v<typename BRP::targetP, P>);
    require_k(k);
    require_kappa(kappa);
    if (out_value == 0) throw std::invalid_argument("out_value must be nonzero");

    if (k <= kappa) {
        pbs_msb_value<BRP>(out, ct, k, delta(k) / 2, out_value, bkfft);
        return;
    }

    alignas(64) TFHEpp::TLWE<typename BRP::domainP> shifted;
    scalar_mul_pow2<typename BRP::domainP>(shifted, ct, kappa);

    const int suffix_bits = k - kappa;
    const int w = k - kappa - 1;
    const Torus guard_weight = Torus{1} << w;
    const Torus guard_value =
        static_cast<Torus>(Wide{delta(k)} * Wide{guard_weight});

    alignas(64) TFHEpp::TLWE<typename BRP::targetP> guard;
    ethmsb_value<BRP>(guard, shifted, suffix_bits, kappa, guard_value, bkfft);

    alignas(64) TFHEpp::TLWE<typename BRP::domainP> guarded;
    sub<typename BRP::domainP>(guarded, ct, guard);

    const Torus final_offset =
        static_cast<Torus>((Wide{(Torus{1} << w) + Torus{1}} *
                            Wide{delta(k)}) /
                           Wide{2});
    pbs_msb_value<BRP>(out, guarded, k, final_offset, out_value, bkfft);
}

inline int ethmsb_pbs_count(const int k, const int kappa)
{
    require_k(k);
    require_kappa(kappa);
    return k <= kappa ? 1 : 1 + ethmsb_pbs_count(k - kappa, kappa);
}

inline TFHEpp::TLWE<P> encrypt_unsigned_for_compare(
    const Torus x, const int t, const TFHEpp::SecretKey& sk)
{
    if (t < 1 || t > 62)
        throw std::invalid_argument("comparison supports 1 <= t <= 62");
    if (x >= (Torus{1} << t))
        throw std::invalid_argument("unsigned input is outside t-bit range");
    const int k = t + 1;
    return TFHEpp::tlweSymEncrypt<P>(encode_unsigned(x, k), sk.key.get<P>());
}

template <class BRP = my_lvl22pbsparam>
void lt(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
        const TFHEpp::TLWE<P>& b, const int t, const int kappa,
        const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
        const Torus out_value = BOOL_ONE)
{
    const int k = t + 1;
    alignas(64) TFHEpp::TLWE<P> diff;
    sub<P>(diff, a, b);
    ethmsb_value<BRP>(out, diff, k, kappa, out_value, bkfft);
}

template <class BRP = my_lvl22pbsparam>
void gt(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
        const TFHEpp::TLWE<P>& b, const int t, const int kappa,
        const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
        const Torus out_value = BOOL_ONE)
{
    const int k = t + 1;
    alignas(64) TFHEpp::TLWE<P> diff;
    sub<P>(diff, b, a);
    ethmsb_value<BRP>(out, diff, k, kappa, out_value, bkfft);
}

template <class BRP = my_lvl22pbsparam>
void ge(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
        const TFHEpp::TLWE<P>& b, const int t, const int kappa,
        const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
        const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P> is_lt;
    lt<BRP>(is_lt, a, b, t, kappa, bkfft, out_value);
    const TFHEpp::TLWE<P> one = trivial_constant<P>(out_value);
    sub<P>(out, one, is_lt);
}

template <class BRP = my_lvl22pbsparam>
void le(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
        const TFHEpp::TLWE<P>& b, const int t, const int kappa,
        const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
        const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P> is_gt;
    gt<BRP>(is_gt, a, b, t, kappa, bkfft, out_value);
    const TFHEpp::TLWE<P> one = trivial_constant<P>(out_value);
    sub<P>(out, one, is_gt);
}

template <class BRP = my_lvl22pbsparam>
void neq(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
         const TFHEpp::TLWE<P>& b, const int t, const int kappa,
         const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
         const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P> is_lt;
    alignas(64) TFHEpp::TLWE<P> is_gt;
    lt<BRP>(is_lt, a, b, t, kappa, bkfft, out_value);
    gt<BRP>(is_gt, a, b, t, kappa, bkfft, out_value);
    add<P>(out, is_lt, is_gt);
}

template <class BRP = my_lvl22pbsparam>
void eq(TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
        const TFHEpp::TLWE<P>& b, const int t, const int kappa,
        const TFHEpp::BootstrappingKeyFFT<BRP>& bkfft,
        const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P> is_neq;
    neq<BRP>(is_neq, a, b, t, kappa, bkfft, out_value);
    const TFHEpp::TLWE<P> one = trivial_constant<P>(out_value);
    sub<P>(out, one, is_neq);
}

}  // namespace my_ethmsb
