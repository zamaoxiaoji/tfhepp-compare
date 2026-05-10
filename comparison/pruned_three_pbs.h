#pragma once
/**
 * @file pruned_three_pbs.h
 * @brief Algorithm 2 — Periodic-pruned BitExtract + Boolean-to-Arithmetic +
 *        ETHMSB (Chapter 3 of Tang Li's thesis, Algorithm
 *        alg:full_periodic_pruned_iterative_pbs and alg:et-hmsb-expanded).
 *
 * Each recursion level performs:
 *   1. BitExtract (Q/2 LUT, periodic-pruned blind rotation) → boolean Enc(b_κ)
 *   2. B2A weight PBS                                       → arithmetic Enc(b_κ · 2^{k-κ-1})
 *   3. Subtract from the original ciphertext to clear b_κ
 * The base case is one direct PBS-MSB (κ=5).
 *
 * For a (5+m·κ)-bit input this uses 2m PBS calls plus the final base-case PBS.
 */
#include "tfhepp_utils.h"

namespace tfhepp_compare::three_pbs
{
    using namespace tfhepp_compare;

    // ── Lvl1 path ──
    void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     const TFHEEvalKey &ek, bool result_type);

    void ExtractMSB10(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                      const TFHEEvalKey &ek, bool result_type);

    // ── Lvl2 → Lvl1 path ──
    void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type);

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

    // Public dispatch.
    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type);

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type);

} // namespace tfhepp_compare::three_pbs
