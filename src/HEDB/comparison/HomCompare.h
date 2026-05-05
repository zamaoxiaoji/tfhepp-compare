#pragma once
#include "extract_msb.h"
#include "ethmsb.h"
#include "operators.h"
#include "tfhepp_utils.h"

namespace HEDB
{
    // ============================================================
    //  Standard comparison  (HE3DB original via HomMSB)
    // ============================================================

    template <typename P>
    void greater_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                      TLWELvl1 &res, uint32_t plain_bits,
                      TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = cipher2[i] - cipher1[i];
        HomMSB(res, sub, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void less_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                   TLWELvl1 &res, uint32_t plain_bits,
                   TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = cipher1[i] - cipher2[i];
        HomMSB(res, sub, plain_bits + 1, ek, result_type);
    }

    // ============================================================
    //  ETHMSB comparison  (Algorithm 2 from the paper)
    //  k ≤ κ: identical to standard PBS-MSB (no gap offset)
    //  k > κ: guard-bit injection + gap offset
    // ============================================================

    template <typename P>
    void ethmsb_greater_than(TFHEpp::TLWE<P> &cipher1,
                             TFHEpp::TLWE<P> &cipher2,
                             TLWELvl1 &res, uint32_t plain_bits,
                             TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = cipher2[i] - cipher1[i];
        HomETHMSB(res, sub, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void ethmsb_less_than(TFHEpp::TLWE<P> &cipher1,
                          TFHEpp::TLWE<P> &cipher2,
                          TLWELvl1 &res, uint32_t plain_bits,
                          TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = cipher1[i] - cipher2[i];
        HomETHMSB(res, sub, plain_bits + 1, ek, result_type);
    }

} // namespace HEDB
