#pragma once
/**
 * @file micro_pbs.h
 * @brief Boolean Q/2 -> guard-value conversion entry points.
 *
 *  Mathematics (negacyclic LUT):
 *      A = round_torus(guard_value / 2)
 *      F(0)   = -A
 *      F(Q/2) = +A          // satisfies F(x + Q/2) = -F(x) automatically
 *      ct_pm   = MicroPBS(bit_qhalf, F)
 *      ct_guard = ct_pm + A             // {0 → 0,  Q/2 → guard_value}
 *
 *  Lvl1Compat parameter source:
 *      Reuses TFHEpp::lvl01param (Lvl0 -> Lvl1 PBS). This is not a true
 *      micro-parameter path; it is preserved as the compatibility baseline
 *      and named explicitly as `micro_lut_lvl01_compat`.
 *
 *  This is an UNSAFE_PERFORMANCE_ONLY skeleton: correctness of the LUT
 *  semantics is the only thing checked against by the isolated test.
 */
#include "tfhepp_utils.h"

namespace tfhepp_compare::micro_pbs
{
    using namespace tfhepp_compare;

    // Conversion-mode selector for the Boolean -> guard-value step in the
    // 3-PBS pipeline and isolated sweep.
    enum class ConversionMode {
        BoolArithExisting,  // legacy Boolean2Weight (μ_polygen + +w)
        MicroLutLvl01Compat,
        MicroUnsafeDirect,
        MicroUnsafeKsPre,
        MicroUnsafeBest
    };

    const char *ConversionModeName(ConversionMode mode);

    struct MicroUnsafeBestPlan {
        bool supported = false;
        const char *status = "MICRO_CONVERSION_UNSUPPORTED_K";
        const char *classification = "MICRO_CONVERSION_UNSUPPORTED_K";
        const char *conversion_mode = "micro_unsafe_best";
        const char *candidate = "micro_n64_N1024_l3_b6";
        const char *offset_name = "unsupported";
        Lvl1::T offset_value = 0;
        Lvl1::T guard_value = 0;
        Lvl1::T A = 0;
        double margin_log2 = 0.0;
        double expected_noise_floor_log2 = 22.4;
    };

    MicroUnsafeBestPlan SelectMicroUnsafeBestPlan(uint32_t window_local_k);

    // Output: ct_guard ∈ Enc(0) when bit_qhalf decrypts to 0,
    //         ct_guard ∈ Enc(guard_value) when bit_qhalf decrypts to Q/2.
    //
    // guard_value must be a torus integer of the form 2^(digits-1-k) where k
    // is the bit index being cleared at this recursion level.
    //
    // Compatibility implementation. It reuses the Lvl10 IKS + Lvl01 PBS chain
    // already present in TFHEEvalKey. Latency equality with
    // LegacyBool2WeightBaseline_Lvl1 is expected because both are the same
    // conversion chain with equivalent LUTs.
    void MicroPBSQHalfToGuardValue_Lvl1Compat(TLWELvl1 &ct_guard,
                                              const TLWELvl1 &bit_qhalf,
                                              Lvl1::T guard_value,
                                              const TFHEEvalKey &ek);

    // Backward-compatible old API name. Keep callers building, but do not use
    // this name in new reports because it can be mistaken for a real micro
    // parameter path.
    void MicroPBSQHalfToGuardValue_Lvl1(TLWELvl1 &ct_guard,
                                        const TLWELvl1 &bit_qhalf,
                                        Lvl1::T guard_value,
                                        const TFHEEvalKey &ek);

    void MicroPBSQHalfToGuardValue_Lvl2(TLWELvl2 &ct_guard,
                                        const TLWELvl2 &bit_qhalf,
                                        Lvl2::T guard_value,
                                        const TFHEEvalKey &ek);

    // Legacy Boolean2Weight implementation (μ_polygen + +w on the same Lvl01
    // PBS) exposed here for A/B comparison in the isolated test. This is the
    // BoolArithExisting baseline; the math is mathematically equivalent to
    // MicroPBSQHalfToGuardValue_Lvl1, so the isolated test mainly verifies
    // they agree end-to-end and measures latency overhead (which should be
    // close to zero since both run one PBS on the same Lvl01 chain).
    void LegacyBool2WeightBaseline_Lvl1(TLWELvl1 &ct_guard,
                                        const TLWELvl1 &bit_qhalf,
                                        Lvl1::T guard_value,
                                        const TFHEEvalKey &ek);

} // namespace tfhepp_compare::micro_pbs
