#pragma once
/**
 * @file HomCompare.h
 * @brief Comparison primitives (>, ≥, <, ≤, =) layered on top of a
 *        HomMSB-like sign-extraction routine.
 *
 * Mirrors HE3DB's HEDB/comparison/HomCompare.h. To avoid duplicating the
 * five operators per algorithm, each algorithm namespace just imports the
 * generic templates and binds them to its own HomMSB.
 *
 * Usage:
 *     using namespace tfhepp_compare::ethmsb;       // 1-PBS guard-bit ETHMSB
 *     greater_than<Lvl1>(c0, c1, res, plain_bits, ek, ARITHMETIC);
 *
 *     using namespace tfhepp_compare::three_pbs;    // BitExtract + B2A + ETHMSB
 *     greater_than<Lvl1>(c0, c1, res, plain_bits, ek, ARITHMETIC);
 */
#include "ethmsb.h"
#include "pruned_three_pbs.h"

namespace tfhepp_compare::detail
{
    using namespace tfhepp_compare;

    // a > b  ⇔  msb(b - a) == 1
    template <auto HomMSBFn, typename P>
    inline void greater_than_impl(TFHEpp::TLWE<P> &cipher1,
                                  TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                                  uint32_t plain_bits, TFHEEvalKey &ek,
                                  bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomMSBFn(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    // a ≥ b  ⇔  NOT(msb(a - b))
    template <auto HomMSBFn, typename P>
    inline void greater_than_equal_impl(TFHEpp::TLWE<P> &cipher1,
                                        TFHEpp::TLWE<P> &cipher2,
                                        TLWELvl1 &res, uint32_t plain_bits,
                                        TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        HomMSBFn(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) LOG_to_ARI(res, res, ek);
    }

    // a < b  ⇔  msb(a - b) == 1
    template <auto HomMSBFn, typename P>
    inline void less_than_impl(TFHEpp::TLWE<P> &cipher1,
                               TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                               uint32_t plain_bits, TFHEEvalKey &ek,
                               bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        HomMSBFn(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    // a ≤ b  ⇔  NOT(msb(b - a))
    template <auto HomMSBFn, typename P>
    inline void less_than_equal_impl(TFHEpp::TLWE<P> &cipher1,
                                     TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                                     uint32_t plain_bits, TFHEEvalKey &ek,
                                     bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        HomMSBFn(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) LOG_to_ARI(res, res, ek);
    }

    template <auto HomMSBFn, typename P>
    inline void equal_impl(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2,
                           TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek,
                           bool result_type)
    {
        TLWELvl1 ge_tlwe, le_tlwe;
        greater_than_equal_impl<HomMSBFn, P>(cipher1, cipher2, ge_tlwe,
                                             plain_bits, ek, LOGIC);
        less_than_equal_impl<HomMSBFn, P>(cipher1, cipher2, le_tlwe, plain_bits,
                                          ek, LOGIC);
        HomAND(res, ge_tlwe, le_tlwe, ek, result_type);
    }
} // namespace tfhepp_compare::detail

// ── Bind the templates to each algorithm's HomMSB. ──
//
// HomMSB has two overloads (Lvl1 input, Lvl2 input). To pick the right one
// inside a generic template we wrap each as a lambda that infers from the P
// template parameter via the input ciphertext type.

#define TFHEPP_COMPARE_DEFINE_OPS(NS)                                        \
    namespace tfhepp_compare::NS                                              \
    {                                                                         \
        template <typename P>                                                 \
        inline void greater_than(TFHEpp::TLWE<P> &c1, TFHEpp::TLWE<P> &c2,    \
                                 TLWELvl1 &res, uint32_t plain_bits,         \
                                 TFHEEvalKey &ek, bool result_type)          \
        {                                                                     \
            tfhepp_compare::detail::greater_than_impl<                        \
                static_cast<void (*)(TLWELvl1 &, const TFHEpp::TLWE<P> &,    \
                                     uint32_t, const TFHEEvalKey &, bool)>(  \
                    &NS::HomMSB),                                             \
                P>(c1, c2, res, plain_bits, ek, result_type);                \
        }                                                                     \
        template <typename P>                                                 \
        inline void greater_than_equal(TFHEpp::TLWE<P> &c1,                  \
                                       TFHEpp::TLWE<P> &c2, TLWELvl1 &res,   \
                                       uint32_t plain_bits, TFHEEvalKey &ek, \
                                       bool result_type)                      \
        {                                                                     \
            tfhepp_compare::detail::greater_than_equal_impl<                  \
                static_cast<void (*)(TLWELvl1 &, const TFHEpp::TLWE<P> &,    \
                                     uint32_t, const TFHEEvalKey &, bool)>(  \
                    &NS::HomMSB),                                             \
                P>(c1, c2, res, plain_bits, ek, result_type);                \
        }                                                                     \
        template <typename P>                                                 \
        inline void less_than(TFHEpp::TLWE<P> &c1, TFHEpp::TLWE<P> &c2,       \
                              TLWELvl1 &res, uint32_t plain_bits,            \
                              TFHEEvalKey &ek, bool result_type)              \
        {                                                                     \
            tfhepp_compare::detail::less_than_impl<                           \
                static_cast<void (*)(TLWELvl1 &, const TFHEpp::TLWE<P> &,    \
                                     uint32_t, const TFHEEvalKey &, bool)>(  \
                    &NS::HomMSB),                                             \
                P>(c1, c2, res, plain_bits, ek, result_type);                \
        }                                                                     \
        template <typename P>                                                 \
        inline void less_than_equal(TFHEpp::TLWE<P> &c1, TFHEpp::TLWE<P> &c2, \
                                    TLWELvl1 &res, uint32_t plain_bits,      \
                                    TFHEEvalKey &ek, bool result_type)        \
        {                                                                     \
            tfhepp_compare::detail::less_than_equal_impl<                     \
                static_cast<void (*)(TLWELvl1 &, const TFHEpp::TLWE<P> &,    \
                                     uint32_t, const TFHEEvalKey &, bool)>(  \
                    &NS::HomMSB),                                             \
                P>(c1, c2, res, plain_bits, ek, result_type);                \
        }                                                                     \
        template <typename P>                                                 \
        inline void equal(TFHEpp::TLWE<P> &c1, TFHEpp::TLWE<P> &c2,           \
                          TLWELvl1 &res, uint32_t plain_bits,                \
                          TFHEEvalKey &ek, bool result_type)                  \
        {                                                                     \
            tfhepp_compare::detail::equal_impl<                               \
                static_cast<void (*)(TLWELvl1 &, const TFHEpp::TLWE<P> &,    \
                                     uint32_t, const TFHEEvalKey &, bool)>(  \
                    &NS::HomMSB),                                             \
                P>(c1, c2, res, plain_bits, ek, result_type);                \
        }                                                                     \
    }

TFHEPP_COMPARE_DEFINE_OPS(ethmsb)
TFHEPP_COMPARE_DEFINE_OPS(three_pbs)

#undef TFHEPP_COMPARE_DEFINE_OPS
