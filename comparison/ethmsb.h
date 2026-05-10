#pragma once
/**
 * @file ethmsb.h
 * @brief Algorithm 1 — ETHMSB with gap offset (samplepaper.tex).
 *
 * Recursive HE3DB-style HomMSB pipeline that extracts the sign bit of a
 * (plain_bits)-bit ciphertext using one guard-bit-extraction PBS per level
 * plus a final MSB-BS, with a refined gap offset on every internal level.
 *
 * Public API mirrors HE3DB's extract_msb.h: a family of fixed-width
 * ExtractMSB{N} entry points and a HomMSB dispatcher.
 */
#include "tfhepp_utils.h"

namespace tfhepp_compare::ethmsb
{
    using namespace tfhepp_compare;

    // ── Lvl1 path: plain_bits in (0, 10] ──
    void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     const TFHEEvalKey &ek, bool result_type);

    void ExtractMSB10(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                      const TFHEEvalKey &ek, bool result_type);

    // ── Lvl2 → Lvl1 path: plain_bits in (0, 33] ──
    void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type);

    void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type);

    void ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type);

    void ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type);

    void ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type);

    void ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type);

    void ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type);

    // Public dispatch: select the correct ExtractMSB{N} based on plain_bits.
    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type);

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type);

} // namespace tfhepp_compare::ethmsb
