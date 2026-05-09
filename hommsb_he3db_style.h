#pragma once
/**
 * @file hommsb_he3db_style.h
 * @brief Re-implementation of HE3DB's HomMSB algorithm using our TFHEpp params.
 *        Follows HE3DB's exact algorithm: shift → MSB-BS → subtract → Identity-BS → subtract → recurse
 */

#include "ethmsb.h"  // reuse base types, μ_polygen, encrypt/decrypt

namespace HomMSB_NS
{
    using namespace ETHMSB_NS;  // import base types

    // ── gpolygen: identity LUT polynomial (same as HE3DB) ──
    template <class P>
    TFHEpp::Polynomial<P> gpolygen(uint32_t plain_bits, uint32_t scale_bits)
    {
        TFHEpp::Polynomial<P> poly;
        uint32_t padding_bits = P::nbit - plain_bits;
        for (int i = 0; i < (int)P::n; i++)
            poly[i] = (1ULL << scale_bits) * (i >> padding_bits);
        return poly;
    }

    // ── IdeGateBootstrapping (same as HE3DB) ──
    // Identity bootstrapping: extracts the message and re-encodes it
    inline void IdeGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                     uint32_t scale_bits, const TFHEEvalKey &ek)
    {
        constexpr uint64_t offset =
            1ULL << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;
        constexpr uint32_t plain_bits = 4;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01,
            gpolygen<Lvl1>(plain_bits, scale_bits));
    }

    inline void IdeGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                     uint32_t scale_bits, const TFHEEvalKey &ek)
    {
        constexpr uint64_t offset =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - 7);
        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;
        constexpr uint32_t plain_bits = 5;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02,
            gpolygen<Lvl2>(plain_bits, scale_bits));
    }

    // ══════════════════════════════════════════════════════════════
    //  HE3DB-style HomMSB: Lvl1 path
    // ══════════════════════════════════════════════════════════════

    // Base case: ≤5 bits → direct MSB-BS
    inline void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                            const TFHEEvalKey &ek, bool result_type)
    {
        ETHMSB_NS::MSBGateBootstrapping(res, tlwe, ek, result_type);
    }

    // 5 < plain_bits ≤ 9: one HE3DB recursion level on Lvl1
    inline void ExtractMSB9(TLWELvl1 &res, const TLWELvl1 &tlwe,
                            uint32_t plain_bits, const TFHEEvalKey &ek,
                            bool result_type)
    {
        TLWELvl1 shift_tlwe, sign_tlwe5;
        uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits;
        // Step 1: left-shift to align 5-bit window
        for (size_t i = 0; i <= Lvl1::n; i++)
            shift_tlwe[i] = tlwe[i] << (plain_bits - 5);
        // Step 2: extract MSB of 5-bit window
        ETHMSB_NS::MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC);
        // Step 3: subtract MSB from shifted
        for (size_t i = 0; i <= Lvl1::n; i++)
            shift_tlwe[i] = shift_tlwe[i] - sign_tlwe5[i];
        // Step 4: Identity BS to round
        IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);
        // Step 5: subtract from original
        for (size_t i = 0; i <= Lvl1::n; i++)
            res[i] = tlwe[i] - shift_tlwe[i];
        // Step 6: final MSB extraction
        ExtractMSB5(res, res, ek, result_type);
    }

    // ══════════════════════════════════════════════════════════════
    //  HE3DB-style HomMSB: Lvl2 path (output to Lvl1)
    // ══════════════════════════════════════════════════════════════

    // Helper: one HE3DB level on Lvl2 (shift→MSB-BS→sub→Ide-BS→sub)
    static inline void he3db_lvl2_one_level(TLWELvl2 &res_lvl2,
                                            const TLWELvl2 &tlwe,
                                            uint32_t plain_bits,
                                            const TFHEEvalKey &ek)
    {
        TLWELvl2 shift_tlwe, sign_tlwe;
        uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - plain_bits;
        for (size_t i = 0; i <= Lvl2::n; i++)
            shift_tlwe[i] = tlwe[i] << (plain_bits - 6);
        ETHMSB_NS::MSBGateBootstrapping(sign_tlwe, shift_tlwe, ek, ARITHMETIC);
        for (size_t i = 0; i <= Lvl2::n; i++)
            shift_tlwe[i] = shift_tlwe[i] - sign_tlwe[i];
        IdeGateBootstrapping(shift_tlwe, shift_tlwe, scale_bits, ek);
        for (size_t i = 0; i <= Lvl2::n; i++)
            res_lvl2[i] = tlwe[i] - shift_tlwe[i];
    }

    // ImExtractMSB5: Lvl2 → IKS21 → Lvl1 → MSB5
    inline void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              const TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        ExtractMSB5(res, res, ek, result_type);
    }

    // ImExtractMSB9: Lvl2 → IKS21 → Lvl1 → ExtractMSB9
    inline void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type)
    {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        ExtractMSB9(res, res, plain_bits, ek, result_type);
    }

    // Recursive Lvl2 levels (each peels 5 bits via HE3DB's 2-PBS approach)
    inline void ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 residual;
        he3db_lvl2_one_level(residual, tlwe, plain_bits, ek);
        ImExtractMSB9(res, residual, plain_bits - 5, ek, result_type);
    }

    inline void ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 residual;
        he3db_lvl2_one_level(residual, tlwe, plain_bits, ek);
        ImExtractMSB14(res, residual, plain_bits - 5, ek, result_type);
    }

    inline void ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 residual;
        he3db_lvl2_one_level(residual, tlwe, plain_bits, ek);
        ImExtractMSB19(res, residual, plain_bits - 5, ek, result_type);
    }

    inline void ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 residual;
        he3db_lvl2_one_level(residual, tlwe, plain_bits, ek);
        ImExtractMSB24(res, residual, plain_bits - 5, ek, result_type);
    }

    inline void ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 residual;
        he3db_lvl2_one_level(residual, tlwe, plain_bits, ek);
        ImExtractMSB29(res, residual, plain_bits - 5, ek, result_type);
    }

    // ── Dispatch ──
    inline void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type)
    {
        if (plain_bits <= 5) ExtractMSB5(res, tlwe, ek, result_type);
        else if (plain_bits <= 9)
            ExtractMSB9(res, tlwe, plain_bits, ek, result_type);
        else throw std::invalid_argument("Lvl1 plain_bits out of range.");
    }

    inline void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type)
    {
        if (plain_bits <= 5)
            ImExtractMSB5(res, tlwe, ek, result_type);
        else if (plain_bits <= 9)
            ImExtractMSB9(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 14)
            ImExtractMSB14(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 19)
            ImExtractMSB19(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 24)
            ImExtractMSB24(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 29)
            ImExtractMSB29(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 33)
            ImExtractMSB33(res, tlwe, plain_bits, ek, result_type);
        else throw std::invalid_argument("Plain bits out of range.");
    }

    // ── Comparison operators ──
    template <typename P>
    void greater_than(TFHEpp::TLWE<P> &c1, TFHEpp::TLWE<P> &c2,
                      TLWELvl1 &res, uint32_t plain_bits,
                      TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = c2[i] - c1[i];
        HomMSB(res, sub, plain_bits + 1, ek, result_type);
    }

    template <typename P>
    void less_than(TFHEpp::TLWE<P> &c1, TFHEpp::TLWE<P> &c2,
                   TLWELvl1 &res, uint32_t plain_bits,
                   TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub[i] = c1[i] - c2[i];
        HomMSB(res, sub, plain_bits + 1, ek, result_type);
    }

} // namespace HomMSB_NS
