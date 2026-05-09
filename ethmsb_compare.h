#pragma once
/**
 * @file ethmsb_compare.h
 * @brief Homomorphic comparison operators using ETHMSB
 *        (mirrors HE3DB's HomCompare.h)
 */
#include "ethmsb.h"

namespace ETHMSB_NS
{
    // cipher1 > cipher2 ⟺ msb(cipher2 - cipher1) == 1
    template <typename P>
    void greater_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                      TLWELvl1 &res, uint32_t plain_bits,
                      TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomETHMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void greater_than_pruned_fast(TFHEpp::TLWE<P> &cipher1,
                                  TFHEpp::TLWE<P> &cipher2,
                                  TLWELvl1 &res, uint32_t plain_bits,
                                  TFHEEvalKey &ek, bool result_type,
                                  const PrunedETHMSBOptions &options,
                                  PruneStats *stats = nullptr)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomETHMSBPrunedFast(res, sub_tlwe, plain_bits + 1, ek, result_type,
                            options, stats);
    }

    // cipher1 >= cipher2 ⟺ NOT(msb(cipher1 - cipher2))
    template <typename P>
    void greater_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                            TLWELvl1 &res, uint32_t plain_bits,
                            TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        HomETHMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) ETHMSB_NS::LOG_to_ARI(res, res, ek);
    }

    // cipher1 < cipher2 ⟺ msb(cipher1 - cipher2) == 1
    template <typename P>
    void less_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                   TLWELvl1 &res, uint32_t plain_bits,
                   TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        HomETHMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    // cipher1 <= cipher2 ⟺ NOT(msb(cipher2 - cipher1))
    template <typename P>
    void less_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                         TLWELvl1 &res, uint32_t plain_bits,
                         TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomETHMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) ETHMSB_NS::LOG_to_ARI(res, res, ek);
    }

    // cipher1 == cipher2 ⟺ (cipher1 >= cipher2) AND (cipher1 <= cipher2)
    template <typename P>
    void equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
               TLWELvl1 &res, uint32_t plain_bits,
               TFHEEvalKey &ek, bool result_type)
    {
        TLWELvl1 greater_tlwe, less_tlwe;
        greater_than_equal<P>(cipher1, cipher2, greater_tlwe, plain_bits, ek, LOGIC);
        less_than_equal<P>(cipher1, cipher2, less_tlwe, plain_bits, ek, LOGIC);
        HomAND(res, greater_tlwe, less_tlwe, ek, result_type);
    }

} // namespace ETHMSB_NS
