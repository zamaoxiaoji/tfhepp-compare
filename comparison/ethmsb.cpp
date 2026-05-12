#include "ethmsb.h"

namespace tfhepp_compare::ethmsb
{
    // ══════════════════════════════════════════════════════════════
    //  Internal helpers
    // ══════════════════════════════════════════════════════════════

    // Refined gap offset (Eq. (offset₂) of samplepaper). For the standard
    // κ=5 reliable window, this becomes (2^{κ-1}+1)/2 · Δ when the parent
    // level provides a gap; otherwise we fall back to the κ=5 default.
    static Lvl1::T GapOffsetLvl1(uint32_t guarded_plain_bits)
    {
        constexpr uint32_t kappa  = 5;
        const uint32_t     digits = std::numeric_limits<Lvl1::T>::digits;
        if (guarded_plain_bits <= kappa + 1 || guarded_plain_bits >= digits)
            return Lvl1::T(1) << (digits - 6);
        return (Lvl1::T(1) << (digits - kappa - 2)) +
               (Lvl1::T(1) << (digits - guarded_plain_bits - 1));
    }

    // PBS-MSB on a guarded ciphertext, using the gap-aware offset above.
    static void GapMSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                        uint32_t guarded_plain_bits,
                                        const TFHEEvalKey &ek,
                                        bool result_type)
    {
        Lvl1::T μ = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += GapOffsetLvl1(guarded_plain_bits);

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(μ));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += μ;
    }

    // Guard-bit extraction PBS at level 1: a single PBS that emits an
    // arithmetic-encoded ±weight encrypting the boundary bit at its original
    // weight (Algorithm 2 of samplepaper).
    static void GuardBitExtractBS_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                       Lvl1::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl1::T offset =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(weight));
        res[Lvl1::k * Lvl1::n] += weight;
    }

    static void GuardBitExtractBS_Lvl2(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                       Lvl2::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl2::T offset =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - 7);
        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(weight));
        res[Lvl2::k * Lvl2::n] += weight;
    }

    // One ETHMSB recursion level on Lvl1: shift → guard-bit PBS → subtract.
    // gap_parent_bits=0 indicates the very first call (no gap yet).
    static void ethmsb_lvl1_recursive(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                      uint32_t plain_bits,
                                      const TFHEEvalKey &ek, bool result_type,
                                      uint32_t gap_parent_bits)
    {
        constexpr uint32_t kappa = 5;
        if (plain_bits <= kappa) {
            if (gap_parent_bits == 0)
                MSBGateBootstrapping(res, tlwe, ek, result_type);
            else
                GapMSBGateBootstrapping(res, tlwe, gap_parent_bits, ek,
                                        result_type);
            return;
        }

        // Step 1: left-shift to align b_κ to the MSB of the κ-bit window.
        TLWELvl1 shift_tlwe;
        for (size_t i = 0; i <= Lvl1::n; i++) shift_tlwe[i] = tlwe[i] << kappa;

        // Step 2: single guard-bit extraction PBS, weight 2^{digits-κ-2}.
        constexpr Lvl1::T weight =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - kappa - 2);
        TLWELvl1 guard_bit;
        GuardBitExtractBS_Lvl1(guard_bit, shift_tlwe, weight, ek);

        // Step 3: subtract from the original ciphertext to clear b_κ.
        TLWELvl1 guarded;
        for (size_t i = 0; i <= Lvl1::n; i++)
            guarded[i] = tlwe[i] - guard_bit[i];

        // Step 4: recurse, telling the next level that this guard came from a
        // (plain_bits)-bit window.
        ethmsb_lvl1_recursive(res, guarded, plain_bits - kappa, ek, result_type,
                              plain_bits);
    }

    // One ETHMSB recursion level on Lvl2 (output stays in Lvl2 for further
    // recursion; the final descent step switches down to Lvl1).
    static void ethmsb_lvl2_one_level(TLWELvl2 &out, const TLWELvl2 &tlwe,
                                      const TFHEEvalKey &ek)
    {
        constexpr uint32_t kappa = 5;
        TLWELvl2 shift_tlwe;
        for (size_t i = 0; i <= Lvl2::n; i++) shift_tlwe[i] = tlwe[i] << kappa;

        constexpr Lvl2::T weight =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - kappa - 2);
        TLWELvl2 guard_bit;
        GuardBitExtractBS_Lvl2(guard_bit, shift_tlwe, weight, ek);

        for (size_t i = 0; i <= Lvl2::n; i++)
            out[i] = tlwe[i] - guard_bit[i];
    }

    // Top-of-Lvl2 dispatcher: peel one Lvl2 guard level until the residual
    // fits into Lvl1's reliable range, then switch down to Lvl1 and continue.
    static void ethmsb_lvl2_dispatch(TLWELvl1 &res, const TLWELvl2 &tlwe,
                                     uint32_t plain_bits,
                                     const TFHEEvalKey &ek, bool result_type,
                                     uint32_t gap_parent_bits)
    {
        if (plain_bits <= 6) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            if (gap_parent_bits == 0)
                MSBGateBootstrapping(res, res, plain_bits, ek, result_type);
            else
                GapMSBGateBootstrapping(res, res, gap_parent_bits, ek,
                                        result_type);
            return;
        }
        if (plain_bits <= 9) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            ethmsb_lvl1_recursive(res, res, plain_bits, ek, result_type,
                                  gap_parent_bits);
            return;
        }
        TLWELvl2 guarded;
        ethmsb_lvl2_one_level(guarded, tlwe, ek);
        ethmsb_lvl2_dispatch(res, guarded, plain_bits - 5, ek, result_type,
                             plain_bits);
    }

    // ══════════════════════════════════════════════════════════════
    //  Public API
    // ══════════════════════════════════════════════════════════════

    void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     const TFHEEvalKey &ek, bool result_type)
    {
        MSBGateBootstrapping(res, tlwe, 5, ek, result_type);
    }

    void ExtractMSB10(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                      const TFHEEvalKey &ek, bool result_type)
    {
        ethmsb_lvl1_recursive(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                       uint32_t plain_bits, const TFHEEvalKey &ek,
                       bool result_type)
    {
        (void) plain_bits;
        ethmsb_lvl2_dispatch(res, tlwe, 5, ek, result_type, 0);
    }

    void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type)
    {
        ethmsb_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        ethmsb_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        ethmsb_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        ethmsb_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        ethmsb_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        ethmsb_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type, 0);
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 6)
            MSBGateBootstrapping(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 10)
            ExtractMSB10(res, tlwe, plain_bits, ek, result_type);
        else
            throw std::invalid_argument(
                "ETHMSB: Lvl1 plain_bits out of range (max 10).");
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 6)
            ImExtractMSB5(res, tlwe, plain_bits, ek, result_type);
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
        else
            throw std::invalid_argument(
                "ETHMSB: Lvl2 plain_bits out of range (max 33).");
    }

} // namespace tfhepp_compare::ethmsb
