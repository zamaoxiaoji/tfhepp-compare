#include "tfhepp_utils.h"

namespace tfhepp_compare
{
    // ── MSBGateBootstrapping ──
    void MSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                              const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T μ = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        constexpr uint64_t offset =
            1ULL << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(μ));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += μ;
    }

    void MSBGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                              const TFHEEvalKey &ek, bool result_type)
    {
        Lvl2::T μ = Lvl2::μ;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        constexpr uint64_t offset =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - 7);
        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(μ));
        if (IS_ARITHMETIC(result_type)) res[Lvl2::k * Lvl2::n] += μ;
    }

    // ── IdeGateBootstrapping (identity rounding) ──
    void IdeGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
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
            res, tlwelvl0, *ek.bkfftlvl01, g_polygen<Lvl1>(plain_bits, scale_bits));
    }

    void IdeGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
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
            res, tlwelvl0, *ek.bkfftlvl02, g_polygen<Lvl2>(plain_bits, scale_bits));
    }

    // ── ARI ↔ LOG conversion at level 1 ──
    void ARI_to_LOG(TLWELvl1 &res, const TLWELvl1 &tlwe, const TFHEEvalKey &ek)
    {
        Lvl1::T μ = Lvl1::μ;
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += μ;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(μ));
    }

    void LOG_to_ARI(TLWELvl1 &res, const TLWELvl1 &tlwe, const TFHEEvalKey &ek)
    {
        Lvl1::T μ = Lvl1::μ << 1;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwe, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-μ));
        res[Lvl1::k * Lvl1::n] += μ;
    }

    // ── Gate operators ──
    void HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
                const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T offset = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++) res[i] = ca[i] + cb[i];
        res[Lvl1::k * Lvl1::n] -= Lvl1::μ >> 1;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-offset));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
    }

    void HomOR(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
               const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T offset = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++) res[i] = ca[i] + cb[i];
        res[Lvl1::k * Lvl1::n] += (Lvl1::μ >> 1);
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-offset));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
    }

} // namespace tfhepp_compare
