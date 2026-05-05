#pragma once
#include "extract_msb.h"
#include "operators.h"
#include "tfhepp_utils.h"

namespace HEDB
{
    // ============================================================
    //  Standard comparison  (HE3DB original)
    // ============================================================

    template <typename P>
    void greater_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                      TLWELvl1 &res, uint32_t plain_bits,
                      TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void less_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                   TLWELvl1 &res, uint32_t plain_bits,
                   TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void greater_than_equal(TFHEpp::TLWE<P> &cipher1,
                            TFHEpp::TLWE<P> &cipher2,
                            TLWELvl1 &res, uint32_t plain_bits,
                            TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<TFHEpp::lvl1param>(res, res);
        if (IS_ARITHMETIC(result_type))
            TFHEpp::LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void less_than_equal(TFHEpp::TLWE<P> &cipher1,
                         TFHEpp::TLWE<P> &cipher2,
                         TLWELvl1 &res, uint32_t plain_bits,
                         TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<TFHEpp::lvl1param>(res, res);
        if (IS_ARITHMETIC(result_type))
            TFHEpp::LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
               TLWELvl1 &res, uint32_t plain_bits,
               TFHEEvalKey &ek, bool result_type)
    {
        TLWELvl1 gt_tlwe, lt_tlwe;
        greater_than_equal<P>(cipher1, cipher2, gt_tlwe,
                              plain_bits, ek, LOGIC);
        less_than_equal<P>(cipher1, cipher2, lt_tlwe,
                           plain_bits, ek, LOGIC);
        HomAND(res, gt_tlwe, lt_tlwe, ek, result_type);
    }

    // ============================================================
    //  GapMSB comparison  (offset-shifted decision boundary)
    // ============================================================
    //
    //  Key idea from Chapter 3:
    //  The constraint bit_k(m) = 0 creates an unreachable interval
    //  of width g(k) = (w_k + 1) * Δ near the MSB threshold.
    //
    //  We add offset' = (w_k + 1)/2 * Δ to the difference ciphertext
    //  BEFORE feeding it to the standard MSB extraction pipeline.
    //  This shifts the MSB decision boundary to the midpoint of the
    //  gap, maximizing the noise tolerance margin.
    //
    //  The standard HE3DB MSBGateBootstrapping already has its own
    //  internal offset (Q/64 for lvl1, Q/128 for lvl2) for rounding
    //  correction. Our gap offset is **added to the input ciphertext**
    //  at the HomCompare level, so both offsets work together:
    //    total_shift = gap_offset + internal_rounding_offset
    //
    //  guard_k: configurable bit index (1 ≤ guard_k ≤ p-2).
    //    p_eff = plain_bits + 1 (bit width of the difference)
    //    w_k = 2^(p_eff - 1 - guard_k)
    //    Δ = Q / 2^p_eff
    //    gap_offset = (w_k + 1) / 2 * Δ

    template <typename P>
    static void add_gap_offset(TFHEpp::TLWE<P> &tlwe,
                               uint32_t plain_bits, uint32_t guard_k)
    {
        const uint32_t p_eff = plain_bits + 1;
        constexpr uint32_t q = std::numeric_limits<typename P::T>::digits;
        const uint64_t w_k = 1ULL << (p_eff - 1 - guard_k);
        const uint64_t delta = 1ULL << (q - p_eff);
        const typename P::T gap_offset =
            static_cast<typename P::T>(((w_k + 1) * delta) / 2);
        tlwe[P::k * P::n] += gap_offset;
    }

    template <typename P>
    void gap_greater_than(TFHEpp::TLWE<P> &cipher1,
                          TFHEpp::TLWE<P> &cipher2,
                          TLWELvl1 &res, uint32_t plain_bits,
                          TFHEEvalKey &ek, bool result_type,
                          uint32_t guard_k)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        add_gap_offset<P>(sub_tlwe, plain_bits, guard_k);
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void gap_less_than(TFHEpp::TLWE<P> &cipher1,
                       TFHEpp::TLWE<P> &cipher2,
                       TLWELvl1 &res, uint32_t plain_bits,
                       TFHEEvalKey &ek, bool result_type,
                       uint32_t guard_k)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        add_gap_offset<P>(sub_tlwe, plain_bits, guard_k);
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void gap_greater_than_equal(TFHEpp::TLWE<P> &cipher1,
                                TFHEpp::TLWE<P> &cipher2,
                                TLWELvl1 &res, uint32_t plain_bits,
                                TFHEEvalKey &ek, bool result_type,
                                uint32_t guard_k)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        add_gap_offset<P>(sub_tlwe, plain_bits, guard_k);
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<TFHEpp::lvl1param>(res, res);
        if (IS_ARITHMETIC(result_type))
            TFHEpp::LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void gap_less_than_equal(TFHEpp::TLWE<P> &cipher1,
                             TFHEpp::TLWE<P> &cipher2,
                             TLWELvl1 &res, uint32_t plain_bits,
                             TFHEEvalKey &ek, bool result_type,
                             uint32_t guard_k)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        add_gap_offset<P>(sub_tlwe, plain_bits, guard_k);
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<TFHEpp::lvl1param>(res, res);
        if (IS_ARITHMETIC(result_type))
            TFHEpp::LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void gap_equal(TFHEpp::TLWE<P> &cipher1,
                   TFHEpp::TLWE<P> &cipher2,
                   TLWELvl1 &res, uint32_t plain_bits,
                   TFHEEvalKey &ek, bool result_type,
                   uint32_t guard_k)
    {
        TLWELvl1 gt_tlwe, lt_tlwe;
        gap_greater_than_equal<P>(cipher1, cipher2, gt_tlwe,
                                  plain_bits, ek, LOGIC, guard_k);
        gap_less_than_equal<P>(cipher1, cipher2, lt_tlwe,
                               plain_bits, ek, LOGIC, guard_k);
        HomAND(res, gt_tlwe, lt_tlwe, ek, result_type);
    }

} // namespace HEDB
