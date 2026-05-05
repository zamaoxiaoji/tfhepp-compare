#include "extract_msb.h"
#include "HEDB/utils/utils.h"
#include "tfhepp_utils.h"
using namespace TFHEpp;

namespace HEDB
{

// ================================================================
//  Standard MSB extraction  (identical to HE3DB)
// ================================================================

void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                 const TFHEEvalKey &ek, bool result_type)
{
    MSBGateBootstrapping(res, tlwe, ek, result_type);
}

void ExtractMSB9(TLWELvl1 &res, const TLWELvl1 &tlwe,
                 uint32_t plain_bits, const TFHEEvalKey &ek,
                 bool result_type)
{
    TLWELvl1 shift_tlwe, sign_tlwe5;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits;
    for (size_t i = 0; i <= Lvl1::n; i++)
        shift_tlwe[i] = tlwe[i] << (plain_bits - 5);
    MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC);
    for (size_t i = 0; i <= Lvl1::n; i++)
        shift_tlwe[i] = shift_tlwe[i] - sign_tlwe5[i];
    IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);
    for (size_t i = 0; i <= Lvl1::n; i++)
        res[i] = tlwe[i] - shift_tlwe[i];
    ExtractMSB5(res, res, ek, result_type);
}

void ExtractMSB11(TLWELvl2 &res, const TLWELvl2 &tlwe,
                  uint32_t plain_bits, const TFHEEvalKey &ek,
                  bool result_type)
{
    TLWELvl2 shift_tlwe, sign_tlwe6;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - plain_bits;
    for (size_t i = 0; i <= Lvl2::n; i++)
        shift_tlwe[i] = tlwe[i] << (plain_bits - 6);
    MSBGateBootstrapping(sign_tlwe6, shift_tlwe, ek, ARITHMETIC);
    for (size_t i = 0; i <= Lvl2::n; i++)
        shift_tlwe[i] = shift_tlwe[i] - sign_tlwe6[i];
    IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);
    for (size_t i = 0; i <= Lvl2::n; i++)
        res[i] = tlwe[i] - shift_tlwe[i];
    MSBGateBootstrapping(res, res, ek, result_type);
}

void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type)
{
    IdentityKeySwitch<lvl21param>(res, tlwe, *ek.iksklvl21);
    ExtractMSB5(res, res, ek, result_type);
}

void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type)
{
    IdentityKeySwitch<lvl21param>(res, tlwe, *ek.iksklvl21);
    ExtractMSB9(res, res, plain_bits, ek, result_type);
}

// Helper: one recursion level for lvl2 (shift→MSB→IDE→subtract)
static void RecurseLvl2Step(TLWELvl2 &out, const TLWELvl2 &tlwe,
                            uint32_t plain_bits, const TFHEEvalKey &ek)
{
    TLWELvl2 shift_tlwe, sign_tlwe6;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - plain_bits;
    for (size_t i = 0; i <= Lvl2::n; i++)
        shift_tlwe[i] = tlwe[i] << (plain_bits - 6);
    MSBGateBootstrapping(sign_tlwe6, shift_tlwe, ek, ARITHMETIC);
    for (size_t i = 0; i <= Lvl2::n; i++)
        shift_tlwe[i] = shift_tlwe[i] - sign_tlwe6[i];
    IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);
    for (size_t i = 0; i <= Lvl2::n; i++)
        out[i] = tlwe[i] - shift_tlwe[i];
}

void ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                    uint32_t plain_bits, const TFHEEvalKey &ek,
                    bool result_type)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    IdentityKeySwitch<lvl21param>(res, reduced, *ek.iksklvl21);
    ExtractMSB9(res, res, plain_bits - 5, ek, result_type);
}

void ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                    uint32_t plain_bits, const TFHEEvalKey &ek,
                    bool result_type)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    ImExtractMSB14(res, reduced, plain_bits - 5, ek, result_type);
}

void ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                    uint32_t plain_bits, const TFHEEvalKey &ek,
                    bool result_type)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    ImExtractMSB19(res, reduced, plain_bits - 5, ek, result_type);
}

void ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                    uint32_t plain_bits, const TFHEEvalKey &ek,
                    bool result_type)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    ImExtractMSB24(res, reduced, plain_bits - 5, ek, result_type);
}

void ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                    uint32_t plain_bits, const TFHEEvalKey &ek,
                    bool result_type)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    ImExtractMSB29(res, reduced, plain_bits - 5, ek, result_type);
}

// ================================================================
//  Unified dispatcher  (standard)
// ================================================================

void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
            uint32_t plain_bits, const TFHEEvalKey &ek, bool result_type)
{
    if (plain_bits <= 5)       ExtractMSB5(res, tlwe, ek, result_type);
    else if (plain_bits <= 9)  ExtractMSB9(res, tlwe, plain_bits, ek, result_type);
    else throw std::invalid_argument("Plain bits out of range for lvl1.");
}

void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
            uint32_t plain_bits, const TFHEEvalKey &ek, bool result_type)
{
    if (plain_bits <= 5)       ImExtractMSB5(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 9)  ImExtractMSB9(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 14) ImExtractMSB14(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 19) ImExtractMSB19(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 24) ImExtractMSB24(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 29) ImExtractMSB29(res, tlwe, plain_bits, ek, result_type);
    else if (plain_bits <= 33) ImExtractMSB33(res, tlwe, plain_bits, ek, result_type);
    else throw std::invalid_argument("Plain bits out of range for lvl2.");
}

// ================================================================
//  Pruned MSB extraction  (GapMSB approach)
// ================================================================
//
//  The algorithm from chapter3revised:
//
//  Step 1: Extract bit guard_k from ct using a periodic LUT
//          with period M_k = 2^(p-guard_k) in mod-2N space.
//          The pruned blind rotation skips CMUX when ā ≡ 0 (mod M_k).
//          Output: ct_k ∈ {0, Q/2} encoding bit_k.
//
//  Step 2: Clear bit guard_k:
//          ct_gap = ct - ct_k >> (guard_k)
//          where >> means right-shift all TLWE components to rescale
//          from Q/2 encoding to w_k·Δ encoding.
//          BUT: right-shifting introduces truncation noise.
//
//  Alternative for Step 2 (what Tang Li does):
//          Since MSB weight = Q/2 in any encoding, we can use
//          the trick: shift ct to bring bit guard_k into MSB
//          position, extract MSB (prunable since LUT is periodic),
//          then subtract. The MSB output at ARITHMETIC scale = Q/2
//          exactly equals the weight of the now-MSB-position bit.
//
//  Step 3: Do standard HE3DB recursive MSB extraction on ct_gap,
//          with gap offset in the final MSB round.
//
//  Concrete implementation: We adopt the same recursive approach
//  as HE3DB/Tang Li. The key insight is:
//  - The FIRST round's MSBGateBootstrapping extracts the MSB of
//    the shifted ciphertext. After shifting, the effective LUT
//    in the mod-2N domain IS periodic (with period related to the
//    shift amount), enabling pruning.
//  - ALL subsequent rounds use the standard MSB LUT (no pruning).
//  - The FINAL MSB round uses gap offset for increased tolerance.

// --- Pruned versions: use standard recursive structure,
//     but with GapMSB on the final round ---

void PrunedExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type, uint32_t guard_k)
{
    // Final round: use gap MSB with offset
    GapMSBGateBootstrapping(res, tlwe, ek, result_type,
                            guard_k, plain_bits);
}

void PrunedExtractMSB9(TLWELvl1 &res, const TLWELvl1 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type, uint32_t guard_k)
{
    TLWELvl1 shift_tlwe, sign_tlwe5;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits;
    for (size_t i = 0; i <= Lvl1::n; i++)
        shift_tlwe[i] = tlwe[i] << (plain_bits - 5);
    // Standard MSB for intermediate round
    MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC);
    for (size_t i = 0; i <= Lvl1::n; i++)
        shift_tlwe[i] = shift_tlwe[i] - sign_tlwe5[i];
    IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);
    for (size_t i = 0; i <= Lvl1::n; i++)
        res[i] = tlwe[i] - shift_tlwe[i];
    // Final round: gap MSB
    PrunedExtractMSB5(res, res, plain_bits - 4, ek, result_type, guard_k);
}

void PrunedImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                         uint32_t plain_bits, const TFHEEvalKey &ek,
                         bool result_type, uint32_t guard_k)
{
    IdentityKeySwitch<lvl21param>(res, tlwe, *ek.iksklvl21);
    PrunedExtractMSB5(res, res, plain_bits, ek, result_type, guard_k);
}

void PrunedImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                         uint32_t plain_bits, const TFHEEvalKey &ek,
                         bool result_type, uint32_t guard_k)
{
    IdentityKeySwitch<lvl21param>(res, tlwe, *ek.iksklvl21);
    PrunedExtractMSB9(res, res, plain_bits, ek, result_type, guard_k);
}

void PrunedImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type, uint32_t guard_k)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    IdentityKeySwitch<lvl21param>(res, reduced, *ek.iksklvl21);
    PrunedExtractMSB9(res, res, plain_bits - 5, ek, result_type, guard_k);
}

void PrunedImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type, uint32_t guard_k)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    PrunedImExtractMSB14(res, reduced, plain_bits - 5, ek,
                         result_type, guard_k);
}

void PrunedImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type, uint32_t guard_k)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    PrunedImExtractMSB19(res, reduced, plain_bits - 5, ek,
                         result_type, guard_k);
}

void PrunedImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type, uint32_t guard_k)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    PrunedImExtractMSB24(res, reduced, plain_bits - 5, ek,
                         result_type, guard_k);
}

void PrunedImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                          uint32_t plain_bits, const TFHEEvalKey &ek,
                          bool result_type, uint32_t guard_k)
{
    TLWELvl2 reduced;
    RecurseLvl2Step(reduced, tlwe, plain_bits, ek);
    PrunedImExtractMSB29(res, reduced, plain_bits - 5, ek,
                         result_type, guard_k);
}

// ================================================================
//  Unified pruned dispatcher
// ================================================================

void PrunedHomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                  uint32_t plain_bits, const TFHEEvalKey &ek,
                  bool result_type, uint32_t guard_k)
{
    if (plain_bits <= 5)
        PrunedExtractMSB5(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 9)
        PrunedExtractMSB9(res, tlwe, plain_bits, ek, result_type, guard_k);
    else
        throw std::invalid_argument("Plain bits out of range for lvl1.");
}

void PrunedHomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                  uint32_t plain_bits, const TFHEEvalKey &ek,
                  bool result_type, uint32_t guard_k)
{
    if (plain_bits <= 5)
        PrunedImExtractMSB5(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 9)
        PrunedImExtractMSB9(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 14)
        PrunedImExtractMSB14(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 19)
        PrunedImExtractMSB19(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 24)
        PrunedImExtractMSB24(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 29)
        PrunedImExtractMSB29(res, tlwe, plain_bits, ek, result_type, guard_k);
    else if (plain_bits <= 33)
        PrunedImExtractMSB33(res, tlwe, plain_bits, ek, result_type, guard_k);
    else
        throw std::invalid_argument("Plain bits out of range for lvl2.");
}

} // namespace HEDB
