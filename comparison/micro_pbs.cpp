#include "micro_pbs.h"

namespace tfhepp_compare::micro_pbs
{
    const char *ConversionModeName(ConversionMode mode)
    {
        switch (mode) {
            case ConversionMode::BoolArithExisting:
                return "boolarith_existing";
            case ConversionMode::MicroLutLvl01Compat:
                return "micro_lut_lvl01_compat";
            case ConversionMode::MicroUnsafeDirect:
                return "micro_unsafe_direct";
            case ConversionMode::MicroUnsafeKsPre:
                return "micro_unsafe_ks_pre";
            case ConversionMode::MicroUnsafeBest:
                return "micro_unsafe_best";
        }
        return "unknown";
    }

    MicroUnsafeBestPlan SelectMicroUnsafeBestPlan(uint32_t window_local_k)
    {
        MicroUnsafeBestPlan plan;
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        if (window_local_k < digits) {
            plan.guard_value = Lvl1::T(1) << (digits - 1 - window_local_k);
            plan.A = plan.guard_value >> 1;
            plan.margin_log2 =
                plan.A == 0 ? 0.0 : std::log2(static_cast<double>(plan.A));
        }
        if (window_local_k >= 1 && window_local_k <= 6) {
            plan.supported = true;
            plan.status = "selected";
            plan.classification = "OK";
            plan.offset_name = "Q/32";
            plan.offset_value = Lvl1::T(1) << (digits - 5);
            plan.expected_noise_floor_log2 = 22.4;
        }
        else if (window_local_k <= 8) {
            plan.supported = true;
            plan.status = "selected";
            plan.classification = "OK";
            plan.offset_name = "Q/16";
            plan.offset_value = Lvl1::T(1) << (digits - 4);
            plan.expected_noise_floor_log2 = 22.4;
        }
        return plan;
    }

    // Build a negacyclic LUT polynomial that returns -A on coefficients
    // [0, N) and (implicitly) +A after the X^N wrap. Concretely all 2N slots
    // get coefficient -A; the negacyclic extraction at index ≥ N flips the
    // sign for free, exactly satisfying F(x + Q/2) = -F(x) in the boolean
    // Q/2 encoding because Q/2 in the phase domain maps to N in the 2N
    // domain.
    template <class P>
    static TFHEpp::Polynomial<P> MinusA_Polygen(typename P::T A)
    {
        TFHEpp::Polynomial<P> poly;
        for (typename P::T &p : poly) p = -A;
        return poly;
    }

    void MicroPBSQHalfToGuardValue_Lvl1Compat(TLWELvl1 &ct_guard,
                                              const TLWELvl1 &bit_qhalf,
                                              Lvl1::T guard_value,
                                              const TFHEEvalKey &ek)
    {
        // A = round_torus(guard_value / 2). guard_value is a power of two
        // multiple of 2 by construction (= 2^{digits-1-k}, k ≥ 0), so the
        // division is exact and rounding_error = 0.
        const Lvl1::T A = guard_value >> 1;

        // The boolean ciphertext sits at phase ∈ {0, Q/2} + noise. Add the
        // standard Q/64 LUT-segment-midpoint offset used by the rest of the
        // pipeline so the mod-switched representative falls in the LUT
        // segment centre.
        constexpr Lvl1::T offset =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 6);

        TLWELvl1 tlweoffset = bit_qhalf;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);

        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            ct_guard, tlwelvl0, *ek.bkfftlvl01, MinusA_Polygen<Lvl1>(A));

        // Add A homomorphically to map {-A, +A} → {0, 2A = guard_value}.
        ct_guard[Lvl1::k * Lvl1::n] += A;
    }

    void MicroPBSQHalfToGuardValue_Lvl1(TLWELvl1 &ct_guard,
                                        const TLWELvl1 &bit_qhalf,
                                        Lvl1::T guard_value,
                                        const TFHEEvalKey &ek)
    {
        MicroPBSQHalfToGuardValue_Lvl1Compat(ct_guard, bit_qhalf, guard_value,
                                             ek);
    }

    void MicroPBSQHalfToGuardValue_Lvl2(TLWELvl2 &ct_guard,
                                        const TLWELvl2 &bit_qhalf,
                                        Lvl2::T guard_value,
                                        const TFHEEvalKey &ek)
    {
        const Lvl2::T A = guard_value >> 1;

        constexpr Lvl2::T offset =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - 7);

        TLWELvl2 tlweoffset = bit_qhalf;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);

        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            ct_guard, tlwelvl0, *ek.bkfftlvl02, MinusA_Polygen<Lvl2>(A));

        ct_guard[Lvl2::k * Lvl2::n] += A;
    }

    // Legacy Boolean2Weight (kept for A/B). Mathematically equivalent: same
    // LUT shape (μ_polygen returns the same all-(-A) polynomial when called
    // with A), same offset, same chain. The function exists so the isolated
    // test can call something with the historical name without depending on
    // the static helper inside pruned_three_pbs.cpp.
    void LegacyBool2WeightBaseline_Lvl1(TLWELvl1 &ct_guard,
                                        const TLWELvl1 &bit_qhalf,
                                        Lvl1::T guard_value,
                                        const TFHEEvalKey &ek)
    {
        const Lvl1::T A = guard_value >> 1;

        constexpr Lvl1::T offset =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 6);

        TLWELvl1 tlweoffset = bit_qhalf;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);

        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            ct_guard, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(A));

        ct_guard[Lvl1::k * Lvl1::n] += A;
    }

} // namespace tfhepp_compare::micro_pbs
