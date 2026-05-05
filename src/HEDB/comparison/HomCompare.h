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

    // a >= b  ⟺  NOT(MSB(a - b))
    template <typename P>
    void ethmsb_greater_than_equal(TFHEpp::TLWE<P> &cipher1,
                                   TFHEpp::TLWE<P> &cipher2,
                                   TLWELvl1 &res, uint32_t plain_bits,
                                   TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = cipher1[i] - cipher2[i];
        HomETHMSB(res, sub, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek);
    }

    // a <= b  ⟺  NOT(MSB(b - a))
    template <typename P>
    void ethmsb_less_than_equal(TFHEpp::TLWE<P> &cipher1,
                                TFHEpp::TLWE<P> &cipher2,
                                TLWELvl1 &res, uint32_t plain_bits,
                                TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = cipher2[i] - cipher1[i];
        HomETHMSB(res, sub, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek);
    }

    // a == b  ⟺  (a >= b) AND (a <= b)
    template <typename P>
    void ethmsb_equal(TFHEpp::TLWE<P> &cipher1,
                      TFHEpp::TLWE<P> &cipher2,
                      TLWELvl1 &res, uint32_t plain_bits,
                      TFHEEvalKey &ek, bool result_type)
    {
        TLWELvl1 ge, le;
        ethmsb_greater_than_equal<P>(cipher1, cipher2, ge, plain_bits, ek, LOGIC);
        ethmsb_less_than_equal<P>(cipher1, cipher2, le, plain_bits, ek, LOGIC);
        HomAND(res, ge, le, ek, result_type);
    }

} // namespace HEDB
