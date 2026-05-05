#pragma once
#include "tfhepp_utils.h"
#include "HEDB/utils/types.h"

namespace HEDB
{
    // ============================================================
    //  ETHMSB: Error-Truncating Homomorphic MSB Extraction
    //  (Algorithm 2 from the paper)
    //
    //  k ≤ κ:  direct PBS-MSB with offset = Δ/2
    //  k > κ:  shift → extract b_κ with scaled LUT → subtract → recurse(k-κ)
    //
    //  κ = 5 for lvl0→lvl1 (N=1024), κ = 6 for lvl0→lvl2 (N=2048)
    // ============================================================

    // lvl1 version (κ=5): input and output are lvl1
    void ETHMSB_lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     uint32_t k, const TFHEEvalKey &ek,
                     bool result_type);

    // lvl2→lvl1 version: input lvl2, output lvl1
    void ETHMSB_lvl2(TLWELvl1 &res, const TLWELvl2 &tlwe,
                     uint32_t k, const TFHEEvalKey &ek,
                     bool result_type);

    // Unified dispatcher
    void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                   uint32_t k, const TFHEEvalKey &ek,
                   bool result_type);

    void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                   uint32_t k, const TFHEEvalKey &ek,
                   bool result_type);

} // namespace HEDB
