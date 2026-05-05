#include "ethmsb.h"
#include "extract_msb.h"
#include <algorithm>
#include <stdexcept>
using namespace TFHEpp;

namespace HEDB
{

namespace
{

static constexpr uint32_t KAPPA = 5;

uint32_t default_guard_bit(uint32_t plain_bits)
{
    if (plain_bits <= KAPPA)
        return 0;
    return plain_bits - KAPPA; // MSB-first index; weight index is KAPPA - 1.
}

uint32_t window_guard_bit(uint32_t plain_bits)
{
    const uint32_t weight_index =
        std::min<uint32_t>(KAPPA - 1, plain_bits - 1);
    return plain_bits - 1 - weight_index;
}

template <class P>
void validate_plain_bits(uint32_t plain_bits)
{
    constexpr uint32_t q = std::numeric_limits<typename P::T>::digits;
    if (plain_bits == 0 || plain_bits >= q)
        throw std::invalid_argument("ETHMSB plain bits must leave torus scaling margin.");
}

void final_gap_msb(TLWELvl1 &res, const TLWELvl1 &tlwe,
                   uint32_t plain_bits,
                   const TFHEEvalKey &ek, bool result_type)
{
    const uint32_t guard_bit = window_guard_bit(plain_bits);
    TFHEpp::GapMSBGateBootstrapping(res, tlwe, ek, result_type,
                                    guard_bit, plain_bits);
}

void eth_extract_msb5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                      uint32_t plain_bits,
                      const TFHEEvalKey &ek, bool result_type)
{
    final_gap_msb(res, tlwe, plain_bits, ek, result_type);
}

void eth_extract_msb9(TLWELvl1 &res, const TLWELvl1 &tlwe,
                      uint32_t plain_bits, const TFHEEvalKey &ek,
                      bool result_type)
{
    TLWELvl1 shift_tlwe, sign_tlwe5;
    const uint32_t scale_bits =
        std::numeric_limits<Lvl1::T>::digits - plain_bits;

    for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
        shift_tlwe[i] = tlwe[i] << (plain_bits - KAPPA);

    MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC);

    for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
        shift_tlwe[i] = shift_tlwe[i] - sign_tlwe5[i];

    IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);

    for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
        res[i] = tlwe[i] - shift_tlwe[i];

    eth_extract_msb5(res, res, plain_bits, ek, result_type);
}

void eth_im_extract_msb9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                         uint32_t plain_bits, const TFHEEvalKey &ek,
                         bool result_type)
{
    IdentityKeySwitch<lvl21param>(res, tlwe, *ek.iksklvl21);
    eth_extract_msb9(res, res, plain_bits, ek, result_type);
}

void recurse_lvl2_step(TLWELvl2 &out, const TLWELvl2 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek)
{
    TLWELvl2 shift_tlwe, sign_tlwe6;
    const uint32_t scale_bits =
        std::numeric_limits<Lvl2::T>::digits - plain_bits;

    for (size_t i = 0; i <= Lvl2::k * Lvl2::n; i++)
        shift_tlwe[i] = tlwe[i] << (plain_bits - 6);

    MSBGateBootstrapping(sign_tlwe6, shift_tlwe, ek, ARITHMETIC);

    for (size_t i = 0; i <= Lvl2::k * Lvl2::n; i++)
        shift_tlwe[i] = shift_tlwe[i] - sign_tlwe6[i];

    IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);

    for (size_t i = 0; i <= Lvl2::k * Lvl2::n; i++)
        out[i] = tlwe[i] - shift_tlwe[i];
}

void eth_im_extract_msb14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type)
{
    TLWELvl2 reduced;
    recurse_lvl2_step(reduced, tlwe, plain_bits, ek);
    eth_im_extract_msb9(res, reduced, plain_bits - KAPPA, ek, result_type);
}

void eth_im_extract_msb19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type)
{
    TLWELvl2 reduced;
    recurse_lvl2_step(reduced, tlwe, plain_bits, ek);
    eth_im_extract_msb14(res, reduced, plain_bits - KAPPA, ek, result_type);
}

void eth_im_extract_msb24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type)
{
    TLWELvl2 reduced;
    recurse_lvl2_step(reduced, tlwe, plain_bits, ek);
    eth_im_extract_msb19(res, reduced, plain_bits - KAPPA, ek, result_type);
}

void eth_im_extract_msb29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type)
{
    TLWELvl2 reduced;
    recurse_lvl2_step(reduced, tlwe, plain_bits, ek);
    eth_im_extract_msb24(res, reduced, plain_bits - KAPPA, ek, result_type);
}

void eth_im_extract_msb33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type)
{
    TLWELvl2 reduced;
    recurse_lvl2_step(reduced, tlwe, plain_bits, ek);
    eth_im_extract_msb29(res, reduced, plain_bits - KAPPA, ek, result_type);
}

} // namespace

uint32_t DefaultETHMSBGuardBit(uint32_t plain_bits)
{
    return default_guard_bit(plain_bits);
}

void ETHMSB_lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                 uint32_t plain_bits, const TFHEEvalKey &ek,
                 bool result_type)
{
    validate_plain_bits<Lvl1>(plain_bits);

    if (plain_bits <= KAPPA) {
        HomMSB(res, tlwe, plain_bits, ek, result_type);
        return;
    }
    if (plain_bits <= 9) {
        eth_extract_msb9(res, tlwe, plain_bits, ek, result_type);
        return;
    }

    throw std::invalid_argument("ETHMSB lvl1 supports up to 9 plaintext bits.");
}

void ETHMSB_lvl2(TLWELvl1 &res, const TLWELvl2 &tlwe,
                 uint32_t plain_bits, const TFHEEvalKey &ek,
                 bool result_type)
{
    validate_plain_bits<Lvl2>(plain_bits);

    if (plain_bits <= KAPPA) {
        HomMSB(res, tlwe, plain_bits, ek, result_type);
        return;
    }

    if (plain_bits <= 9)
        eth_im_extract_msb9(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 14)
        eth_im_extract_msb14(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 19)
        eth_im_extract_msb19(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 24)
        eth_im_extract_msb24(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 29)
        eth_im_extract_msb29(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 33)
        eth_im_extract_msb33(res, tlwe, plain_bits, ek, result_type);
    else
        throw std::invalid_argument("ETHMSB lvl2 supports up to 33 plaintext bits.");
}

void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
               uint32_t plain_bits, const TFHEEvalKey &ek,
               bool result_type)
{
    ETHMSB_lvl1(res, tlwe, plain_bits, ek, result_type);
}

void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
               uint32_t plain_bits, const TFHEEvalKey &ek,
               bool result_type)
{
    ETHMSB_lvl2(res, tlwe, plain_bits, ek, result_type);
}

} // namespace HEDB
