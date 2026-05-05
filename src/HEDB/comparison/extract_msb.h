#pragma once
#include "tfhepp_utils.h"
#include "HEDB/utils/types.h"

namespace HEDB
{
    // ============================================================
    //  Standard MSB extraction  (HE3DB recursive approach)
    // ============================================================

    // 1–5 bits: single gate bootstrap
    void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     const TFHEEvalKey &ek, bool result_type);

    // 6–9 bits: shift + MSB + IDE + recurse
    void ExtractMSB9(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     uint32_t plain_bits, const TFHEEvalKey &ek,
                     bool result_type);

    // 11 bits (lvl2)
    void ExtractMSB11(TLWELvl2 &res, const TLWELvl2 &tlwe,
                      uint32_t plain_bits, const TFHEEvalKey &ek,
                      bool result_type);

    // Key-switch + recursive  (lvl2 → lvl1 output)
    void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type);

    void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type);

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

    // ============================================================
    //  Unified dispatcher
    // ============================================================

    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                uint32_t plain_bits, const TFHEEvalKey &ek,
                bool result_type);

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                uint32_t plain_bits, const TFHEEvalKey &ek,
                bool result_type);

    // ============================================================
    //  Pruned MSB extraction  (with periodic CMUX skip)
    // ============================================================
    //
    // These variants add periodic pruning to the first-round
    // MSBGateBootstrapping and use gap offset on the final MSB.
    //
    // guard_k: the bit index to clear (configurable, 1 ≤ guard_k ≤ p-1)
    //          determines both the LUT period for the bit-extraction
    //          step and the gap width for the final MSB.

    void PrunedExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                           uint32_t plain_bits, const TFHEEvalKey &ek,
                           bool result_type, uint32_t guard_k);

    void PrunedExtractMSB9(TLWELvl1 &res, const TLWELvl1 &tlwe,
                           uint32_t plain_bits, const TFHEEvalKey &ek,
                           bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                             uint32_t plain_bits, const TFHEEvalKey &ek,
                             bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                             uint32_t plain_bits, const TFHEEvalKey &ek,
                             bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type, uint32_t guard_k);

    void PrunedImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type, uint32_t guard_k);

    // Unified pruned dispatcher
    void PrunedHomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                      uint32_t plain_bits, const TFHEEvalKey &ek,
                      bool result_type, uint32_t guard_k);

    void PrunedHomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                      uint32_t plain_bits, const TFHEEvalKey &ek,
                      bool result_type, uint32_t guard_k);

} // namespace HEDB
