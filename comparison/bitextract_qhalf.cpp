#include "bitextract_qhalf.h"

#include <algorithm>
#include <cmath>

#include "gatebootstrapping.hpp"

namespace tfhepp_compare::bitextract_qhalf
{
    const char *FailureClassName(FailureClass c)
    {
        switch (c) {
            case FailureClass::None:
                return "NONE";
            case FailureClass::PlaintextLutBug:
                return "PLAINTEXT_LUT_BUG";
            case FailureClass::ModswitchRepresentativeBug:
                return "MODSWITCH_REPRESENTATIVE_BUG";
            case FailureClass::CenterOffsetBug:
                return "CENTER_OFFSET_BUG";
            case FailureClass::PeriodicSkipBug:
                return "PERIODIC_SKIP_BUG";
            case FailureClass::IKSRoundingBug:
                return "IKS_ROUNDING_BUG";
            case FailureClass::EncryptedNoiseBug:
                return "ENCRYPTED_NOISE_BUG";
            case FailureClass::BitextractNoiseExceedsCellMargin:
                return "BITEXTRACT_NOISE_EXCEEDS_CELL_MARGIN";
            case FailureClass::ScaleOrIndexingBug:
                return "SCALE_OR_INDEXING_BUG";
            case FailureClass::ExactBitextractNotCertifiable:
                return "EXACT_BITEXTRACT_NOT_CERTIFIABLE";
            case FailureClass::Unknown:
                return "UNKNOWN";
        }
        return "UNKNOWN";
    }

    uint64_t WindowWeight(uint32_t p, uint32_t window_local_k)
    {
        if (window_local_k >= p) return 0;
        return uint64_t(1) << (p - 1 - window_local_k);
    }

    Lvl1::T ExpectedQHalfBit(uint64_t message, uint32_t p,
                             uint32_t window_local_k)
    {
        const uint64_t W = WindowWeight(p, window_local_k);
        const uint64_t bit = W == 0 ? 0 : ((message / W) & 1ULL);
        return bit ? (Lvl1::T(1)
                      << (std::numeric_limits<Lvl1::T>::digits - 1))
                   : Lvl1::T(0);
    }

    Lvl1::T PhaseDelta(uint32_t p)
    {
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        return p >= digits ? Lvl1::T(1) : (Lvl1::T(1) << (digits - p));
    }

    Lvl1::T CenterOffset(uint32_t p)
    {
        return PhaseDelta(p) >> 1;
    }

    double DistanceToNearestBitBoundaryDeltaUnits(uint64_t message, uint32_t p,
                                                  uint32_t window_local_k,
                                                  bool centered_cell)
    {
        (void) p;
        const uint64_t W = WindowWeight(p, window_local_k);
        if (W == 0) return 0.0;
        const double x = static_cast<double>(message % W) +
                         (centered_cell ? 0.5 : 0.0);
        return std::min(x, static_cast<double>(W) - x);
    }

    uint64_t DistanceToNearestBitBoundaryTorus(uint64_t message, uint32_t p,
                                               uint32_t window_local_k,
                                               bool centered_cell)
    {
        const double units = DistanceToNearestBitBoundaryDeltaUnits(
            message, p, window_local_k, centered_cell);
        return static_cast<uint64_t>(
            std::llround(units * static_cast<double>(PhaseDelta(p))));
    }

    bool IsValidWindowLocalK(uint32_t window_local_k)
    {
        return window_local_k >= 1 && window_local_k <= 8 &&
               window_local_k <= Lvl1::nbit;
    }

    template <class P>
    static TFHEpp::Polynomial<P> QHalfBitPolygen(uint32_t period)
    {
        TFHEpp::Polynomial<P> poly = {};
        const uint32_t half_period = period / 2;
        const typename P::T q_half =
            typename P::T(1) << (std::numeric_limits<typename P::T>::digits - 1);
        for (uint32_t i = 0; i < P::n; i++)
            poly[i] = ((i % period) >= half_period) ? q_half : 0;
        return poly;
    }

    template <class P>
    static typename P::targetP::T LutValueFromBRIndex(
        uint32_t br_index,
        const TFHEpp::Polynomial<typename P::targetP> &testvector)
    {
        alignas(64) TFHEpp::TRLWE<typename P::targetP> acc = {};
        TFHEpp::PolynomialMulByXai<typename P::targetP>(
            acc[P::targetP::k], testvector, br_index);
        return acc[P::targetP::k][0];
    }

    template <class P>
    static void BlindRotateTrace(
        TFHEpp::TRLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t period, bool pruned, Trace *trace)
    {
        res = {};
        constexpr uint32_t bitwidth = TFHEpp::bits_needed<0>();
        const uint32_t br_index =
            2 * P::targetP::n -
            ((tlwe[P::domainP::k * P::domainP::n] >>
              (std::numeric_limits<typename P::domainP::T>::digits - 1 -
               P::targetP::nbit + bitwidth))
             << bitwidth);
        TFHEpp::PolynomialMulByXai<typename P::targetP>(
            res[P::targetP::k], testvector, br_index);
        uint32_t skip_count = 0;
        uint32_t cmux_count = 0;
        for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
            constexpr typename P::domainP::T roundoffset =
                typename P::domainP::T(1)
                << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                    P::targetP::nbit + bitwidth);
            const uint32_t a_bar =
                (tlwe[i] + roundoffset) >>
                (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                 P::targetP::nbit + bitwidth)
                    << bitwidth;
            if (a_bar == 0) continue;
            if (pruned && period > 1 && (a_bar % period == 0)) {
                skip_count++;
                continue;
            }
            cmux_count++;
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(
                res, bkfft[i], a_bar);
        }
        if (trace != nullptr) {
            trace->br_index = br_index;
            trace->modswitch_index = (2 * P::targetP::n - br_index) %
                                     (2 * P::targetP::n);
            trace->lut_index = trace->modswitch_index % period;
            trace->skip_count = skip_count;
            trace->cmux_count = cmux_count;
            trace->period = period;
            trace->offset_used = BitExtractFloorOffsetLvl1();
            trace->representative_mode = "he3db_compat_b_trunc_a_round";
        }
    }

    template <class P>
    static void GateBootstrappingTrace(
        TFHEpp::TLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t period, bool pruned, Trace *trace)
    {
        alignas(64) TFHEpp::TRLWE<typename P::targetP> acc;
        BlindRotateTrace<P>(acc, tlwe, bkfft, testvector, period, pruned,
                            trace);
        TFHEpp::SampleExtractIndex<typename P::targetP>(res, acc, 0);
    }

    static void BitExtractQHalfImpl_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                         uint32_t p,
                                         uint32_t window_local_k,
                                         const TFHEEvalKey &ek, bool pruned,
                                         bool centered_cell, Trace *trace)
    {
        const uint32_t period = WindowLocalPeriodLvl1(window_local_k);
        const Lvl1::T center = centered_cell ? CenterOffset(p) : Lvl1::T(0);
        TLWELvl1 tlwecentered = tlwe;
        tlwecentered[Lvl1::k * Lvl1::n] += center;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwecentered,
                                         *ek.iksklvl10);
        tlwelvl0[Lvl0::k * Lvl0::n] +=
            static_cast<Lvl0::T>(BitExtractFloorOffsetLvl1());
        GateBootstrappingTrace<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, QHalfBitPolygen<Lvl1>(period),
            period, pruned, trace);
        if (trace != nullptr) {
            trace->bitextract_centering_enabled = centered_cell;
            trace->center_offset = center;
            trace->period = period;
            trace->offset_used = BitExtractFloorOffsetLvl1();
        }
    }

    void BitExtractQHalf_NoPrune_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                      uint32_t window_local_k,
                                      const TFHEEvalKey &ek, Trace *trace)
    {
        BitExtractQHalfImpl_Lvl1(res, tlwe, 0, window_local_k, ek, false,
                                 false, trace);
    }

    void BitExtractQHalf_Pruned_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                     uint32_t window_local_k,
                                     const TFHEEvalKey &ek, Trace *trace)
    {
        BitExtractQHalfImpl_Lvl1(res, tlwe, 0, window_local_k, ek, true,
                                 false, trace);
    }

    void BitExtractQHalf_NoPrune_Lvl1_CenteredCell(
        TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t p,
        uint32_t window_local_k, const TFHEEvalKey &ek, Trace *trace)
    {
        BitExtractQHalfImpl_Lvl1(res, tlwe, p, window_local_k, ek, false,
                                 true, trace);
    }

    void BitExtractQHalf_Pruned_Lvl1_CenteredCell(
        TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t p,
        uint32_t window_local_k, const TFHEEvalKey &ek, Trace *trace)
    {
        BitExtractQHalfImpl_Lvl1(res, tlwe, p, window_local_k, ek, true,
                                 true, trace);
    }

    static Lvl1::T QExactOraclePhaseImpl_Lvl1(uint64_t message, uint32_t p,
                                              uint32_t window_local_k,
                                              bool centered_cell, Trace *trace)
    {
        const uint32_t period = WindowLocalPeriodLvl1(window_local_k);
        const Lvl1::T delta = PhaseDelta(p);
        const Lvl1::T center = centered_cell ? CenterOffset(p) : Lvl1::T(0);
        TFHEpp::TLWE<Lvl0> tlwe = {};
        tlwe[Lvl0::k * Lvl0::n] =
            static_cast<Lvl0::T>(static_cast<Lvl1::T>(message) * delta +
                                 center +
                                 BitExtractFloorOffsetLvl1());
        constexpr uint32_t bitwidth = TFHEpp::bits_needed<0>();
        const uint32_t br_index =
            2 * Lvl1::n -
            ((tlwe[Lvl0::k * Lvl0::n] >>
              (std::numeric_limits<Lvl0::T>::digits - 1 - Lvl1::nbit +
               bitwidth))
             << bitwidth);
        if (trace != nullptr) {
            trace->br_index = br_index;
            trace->modswitch_index = (2 * Lvl1::n - br_index) % (2 * Lvl1::n);
            trace->lut_index = trace->modswitch_index % period;
            trace->period = period;
            trace->offset_used = BitExtractFloorOffsetLvl1();
            trace->representative_mode = "he3db_compat_b_trunc_a_round";
            trace->bitextract_centering_enabled = centered_cell;
            trace->center_offset = center;
            trace->expected_nearest_message = message;
            trace->distance_to_nearest_bit_boundary_in_delta_units =
                DistanceToNearestBitBoundaryDeltaUnits(message, p,
                                                       window_local_k,
                                                       centered_cell);
            trace->distance_to_nearest_bit_boundary_in_torus =
                DistanceToNearestBitBoundaryTorus(message, p, window_local_k,
                                                  centered_cell);
            trace->skip_count = 0;
            trace->cmux_count = 0;
        }
        return LutValueFromBRIndex<Lvl01>(br_index,
                                          QHalfBitPolygen<Lvl1>(period));
    }

    Lvl1::T QExactOraclePhase_Lvl1(uint64_t message, uint32_t p,
                                   uint32_t window_local_k, Trace *trace)
    {
        return QExactOraclePhaseImpl_Lvl1(message, p, window_local_k, false,
                                          trace);
    }

    Lvl1::T QExactOraclePhase_Lvl1_CenteredCell(uint64_t message, uint32_t p,
                                                uint32_t window_local_k,
                                                Trace *trace)
    {
        return QExactOraclePhaseImpl_Lvl1(message, p, window_local_k, true,
                                          trace);
    }

} // namespace tfhepp_compare::bitextract_qhalf
