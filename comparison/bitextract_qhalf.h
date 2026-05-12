#pragma once
/**
 * @file bitextract_qhalf.h
 * @brief Window-local BitExtract to Boolean Q/2 encoding with trace counters.
 */

#include <cstdint>
#include <limits>

#include "tfhepp_utils.h"

namespace tfhepp_compare::bitextract_qhalf
{
    using namespace tfhepp_compare;

    enum class FailureClass {
        None,
        PlaintextLutBug,
        ModswitchRepresentativeBug,
        CenterOffsetBug,
        PeriodicSkipBug,
        IKSRoundingBug,
        EncryptedNoiseBug,
        BitextractNoiseExceedsCellMargin,
        ScaleOrIndexingBug,
        ExactBitextractNotCertifiable,
        Unknown,
    };

    const char *FailureClassName(FailureClass c);

    struct Trace {
        uint32_t br_index = 0;
        uint32_t modswitch_index = 0;
        uint32_t lut_index = 0;
        uint32_t skip_count = 0;
        uint32_t cmux_count = 0;
        uint32_t period = 0;
        Lvl1::T offset_used = 0;
        bool bitextract_centering_enabled = false;
        Lvl1::T center_offset = 0;
        Lvl1::T original_phase = 0;
        Lvl1::T centered_phase_for_bitextract = 0;
        uint64_t expected_nearest_message = 0;
        uint64_t distance_to_nearest_bit_boundary_in_torus = 0;
        double distance_to_nearest_bit_boundary_in_delta_units = 0.0;
        Lvl1::T phase_before_iks = 0;
        Lvl0::T phase_after_iks_before_offset = 0;
        Lvl1::T phase_after_centering = 0;
        Lvl0::T phase_after_centering_and_iks = 0;
        int64_t signed_error_to_expected_center = 0;
        double signed_error_in_delta_units = 0.0;
        const char *representative_mode = "floor_via_negative_roundoffset";
    };

    constexpr uint32_t WindowLocalPeriodLvl1(uint32_t window_local_k)
    {
        return uint32_t(1) << (Lvl1::nbit + 1 - window_local_k);
    }

    constexpr Lvl1::T Lvl01BRRoundOffset()
    {
        return Lvl1::T(1) << (std::numeric_limits<Lvl0::T>::digits - 2 -
                              Lvl1::nbit);
    }

    constexpr Lvl1::T BitExtractFloorOffsetLvl1()
    {
#ifdef USE_HE3DB_COMPAT
        return Lvl1::T(0);
#else
        return Lvl1::T(0) - Lvl01BRRoundOffset();
#endif
    }

    Lvl1::T ExpectedQHalfBit(uint64_t message, uint32_t p,
                             uint32_t window_local_k);

    uint64_t WindowWeight(uint32_t p, uint32_t window_local_k);

    Lvl1::T PhaseDelta(uint32_t p);

    Lvl1::T CenterOffset(uint32_t p);

    double DistanceToNearestBitBoundaryDeltaUnits(uint64_t message, uint32_t p,
                                                  uint32_t window_local_k,
                                                  bool centered_cell);

    uint64_t DistanceToNearestBitBoundaryTorus(uint64_t message, uint32_t p,
                                               uint32_t window_local_k,
                                               bool centered_cell);

    bool IsValidWindowLocalK(uint32_t window_local_k);

    void BitExtractQHalf_NoPrune_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                      uint32_t window_local_k,
                                      const TFHEEvalKey &ek,
                                      Trace *trace = nullptr);

    void BitExtractQHalf_Pruned_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                     uint32_t window_local_k,
                                     const TFHEEvalKey &ek,
                                     Trace *trace = nullptr);

    void BitExtractQHalf_NoPrune_Lvl1_CenteredCell(
        TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t p,
        uint32_t window_local_k, const TFHEEvalKey &ek,
        Trace *trace = nullptr);

    void BitExtractQHalf_Pruned_Lvl1_CenteredCell(
        TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t p,
        uint32_t window_local_k, const TFHEEvalKey &ek,
        Trace *trace = nullptr);

    // Pure representative oracle used by the diagnostic target. Returns the
    // LUT output phase for an exact plaintext message with no LWE key terms.
    Lvl1::T QExactOraclePhase_Lvl1(uint64_t message, uint32_t p,
                                   uint32_t window_local_k,
                                   Trace *trace = nullptr);

    Lvl1::T QExactOraclePhase_Lvl1_CenteredCell(uint64_t message, uint32_t p,
                                                uint32_t window_local_k,
                                                Trace *trace = nullptr);

} // namespace tfhepp_compare::bitextract_qhalf
