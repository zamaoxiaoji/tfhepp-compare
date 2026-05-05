#pragma once
#include "HEDB/utils/types.h"
#include "tfhepp_utils.h"

namespace HEDB
{
    // ETHMSB offset-only path. This mirrors the HE3DB recursive MSB
    // structure and replaces only the final small-window MSB offset by
    // the guard-gap midpoint offset from samplepaper.tex.

    uint32_t DefaultETHMSBGuardBit(uint32_t plain_bits);

    void ETHMSB_lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     uint32_t plain_bits, const TFHEEvalKey &ek,
                     bool result_type);

    void ETHMSB_lvl2(TLWELvl1 &res, const TLWELvl2 &tlwe,
                     uint32_t plain_bits, const TFHEEvalKey &ek,
                     bool result_type);

    void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type);

    void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type);

} // namespace HEDB
