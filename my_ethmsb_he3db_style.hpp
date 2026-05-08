#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "evalkeygens.hpp"
#include "gatebootstrapping.hpp"
#include "key.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "params.hpp"
#include "tlwe.hpp"

namespace my_ethmsb_he3db_style {

using P2 = TFHEpp::lvl2param;
using P0 = TFHEpp::lvl0param;
using KS20 = TFHEpp::lvl20param;
using BR02 = TFHEpp::lvl02param;
using Torus = typename P2::T;
using Wide = unsigned __int128;

static_assert(std::numeric_limits<Torus>::digits == 64);

inline constexpr Torus BOOL_ONE_L2 = Torus{1} << 61;

inline Torus delta(const int k) { return my_ethmsb::delta(k); }

inline Torus encode_unsigned(const Torus m, const int k)
{
    return my_ethmsb::encode_unsigned(m, k);
}

inline bool decrypt_arith_bit_l2(const TFHEpp::TLWE<P2>& ct,
                                 const TFHEpp::SecretKey& sk,
                                 const Torus out_value = BOOL_ONE_L2)
{
    return my_ethmsb::decrypt_arith_bit<P2>(ct, sk, out_value);
}

inline TFHEpp::TLWE<P2> trivial_constant_l2(const Torus value)
{
    return my_ethmsb::trivial_constant<P2>(value);
}

inline TFHEpp::TLWE<P2> encrypt_unsigned_for_compare_l2(
    const Torus x, const int t, const TFHEpp::SecretKey& sk)
{
    if (t < 1 || t > 62)
        throw std::invalid_argument("comparison supports 1 <= t <= 62");
    if (x >= (Torus{1} << t))
        throw std::invalid_argument("unsigned input is outside t-bit range");
    const int k = t + 1;
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned(x, k),
                                      sk.key.get<P2>());
}

inline void pbs_msb_value_l2_to_l2(
    TFHEpp::TLWE<P2>& out, TFHEpp::TLWE<P2> in, const int input_bits,
    const Torus offset_l2, const Torus out_value_l2,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
{
    my_ethmsb::require_k(input_bits);
    if ((out_value_l2 & Torus{1}) != 0)
        throw std::invalid_argument(
            "pbs_msb_value_l2_to_l2 requires even out_value");

    my_ethmsb::add_const_inplace<P2>(in, offset_l2);

    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20);

    const Torus half = out_value_l2 / 2;
    TFHEpp::Polynomial<P2> testvector;
    testvector.fill(Torus{0} - half);

    alignas(64) TFHEpp::TLWE<P2> tmp;
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(tmp, in0, bkfft02,
                                                testvector);
    out = tmp;
    my_ethmsb::add_const_inplace<P2>(out, half);
}

inline void ethmsb_value_l2_to_l2(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& ct, const int k,
    const int kappa, const Torus out_value_l2,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
{
    my_ethmsb::require_k(k);
    my_ethmsb::require_kappa(kappa);
    if (out_value_l2 == 0)
        throw std::invalid_argument("out_value_l2 must be nonzero");

    if (k <= kappa) {
        pbs_msb_value_l2_to_l2(out, ct, k,
                               my_ethmsb::base_offset_for_current_layer(k),
                               out_value_l2, iksk20, bkfft02);
        return;
    }

    alignas(64) TFHEpp::TLWE<P2> shifted;
    my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);

    const int suffix_bits = k - kappa;
    const Torus guard_value =
        my_ethmsb::guard_value_for_parent_scale(k, kappa);

    alignas(64) TFHEpp::TLWE<P2> guard;
    ethmsb_value_l2_to_l2(guard, shifted, suffix_bits, kappa, guard_value,
                          iksk20, bkfft02);

    alignas(64) TFHEpp::TLWE<P2> guarded;
    my_ethmsb::sub<P2>(guarded, ct, guard);

    pbs_msb_value_l2_to_l2(out, guarded, k,
                           my_ethmsb::gap_offset_for_current_layer(k, kappa),
                           out_value_l2, iksk20, bkfft02);
}

inline int ethmsb_pbs_count(const int k, const int kappa)
{
    return my_ethmsb::ethmsb_pbs_count(k, kappa);
}

inline void lt(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int t, const int kappa,
               const TFHEpp::KeySwitchingKey<KS20>& iksk20,
               const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
               const Torus out_value_l2 = BOOL_ONE_L2)
{
    alignas(64) TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, a, b);
    ethmsb_value_l2_to_l2(out, diff, t + 1, kappa, out_value_l2, iksk20,
                          bkfft02);
}

inline void gt(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int t, const int kappa,
               const TFHEpp::KeySwitchingKey<KS20>& iksk20,
               const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
               const Torus out_value_l2 = BOOL_ONE_L2)
{
    alignas(64) TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, b, a);
    ethmsb_value_l2_to_l2(out, diff, t + 1, kappa, out_value_l2, iksk20,
                          bkfft02);
}

inline void le(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int t, const int kappa,
               const TFHEpp::KeySwitchingKey<KS20>& iksk20,
               const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
               const Torus out_value_l2 = BOOL_ONE_L2)
{
    alignas(64) TFHEpp::TLWE<P2> is_gt;
    gt(is_gt, a, b, t, kappa, iksk20, bkfft02, out_value_l2);
    my_ethmsb::sub<P2>(out, trivial_constant_l2(out_value_l2), is_gt);
}

inline void ge(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int t, const int kappa,
               const TFHEpp::KeySwitchingKey<KS20>& iksk20,
               const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
               const Torus out_value_l2 = BOOL_ONE_L2)
{
    alignas(64) TFHEpp::TLWE<P2> is_lt;
    lt(is_lt, a, b, t, kappa, iksk20, bkfft02, out_value_l2);
    my_ethmsb::sub<P2>(out, trivial_constant_l2(out_value_l2), is_lt);
}

inline void neq(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
                const TFHEpp::TLWE<P2>& b, const int t, const int kappa,
                const TFHEpp::KeySwitchingKey<KS20>& iksk20,
                const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
                const Torus out_value_l2 = BOOL_ONE_L2)
{
    alignas(64) TFHEpp::TLWE<P2> is_lt;
    alignas(64) TFHEpp::TLWE<P2> is_gt;
    lt(is_lt, a, b, t, kappa, iksk20, bkfft02, out_value_l2);
    gt(is_gt, a, b, t, kappa, iksk20, bkfft02, out_value_l2);
    my_ethmsb::add<P2>(out, is_lt, is_gt);
}

inline void eq(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int t, const int kappa,
               const TFHEpp::KeySwitchingKey<KS20>& iksk20,
               const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
               const Torus out_value_l2 = BOOL_ONE_L2)
{
    alignas(64) TFHEpp::TLWE<P2> is_neq;
    neq(is_neq, a, b, t, kappa, iksk20, bkfft02, out_value_l2);
    my_ethmsb::sub<P2>(out, trivial_constant_l2(out_value_l2), is_neq);
}

}  // namespace my_ethmsb_he3db_style
