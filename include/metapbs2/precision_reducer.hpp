#pragma once

#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "metapbs2/bit_extraction.hpp"
#include "metapbs2/metapbs_pipeline.hpp"
#include "tlwe.hpp"

namespace MetaPBS2 {

template <typename T>
inline T ReducerTorusScaleForPrecision(int p) {
    static_assert(std::is_unsigned_v<T>, "ReducerTorusScaleForPrecision expects unsigned torus");
    const int digits = std::numeric_limits<T>::digits;
    if (p <= 0 || p >= digits)
        throw std::invalid_argument("message precision out of torus range");
    return T(1) << (digits - p);
}

struct PrecisionReduceStep {
    int p_in = 0;
    int p_out = 0;
    int removed_bits = 0;
    std::string kind;
};

struct PrecisionReduceSchedule {
    int p_original = 0;
    int p_final = 0;
    std::vector<PrecisionReduceStep> steps;
};

template <class P>
struct PrecisionReduceRoundOutput {
    int round = 0;
    int p_in = 0;
    int p_out = 0;
    TFHEpp::TLWE<P> ct_out{};
};

template <class P>
struct PrecisionReduceResult {
    TFHEpp::TLWE<P> ct_out{};
    int p_out = 0;
    typename P::T delta_out = 0;
    std::string semantic;
    PrecisionReduceSchedule schedule;
};

inline PrecisionReduceSchedule MakeHE3DBStylePrecisionSchedule(
    int p_original, int p_native = 12) {
    if (p_original <= 0)
        throw std::invalid_argument("precision reducer requires positive p_original");
    PrecisionReduceSchedule schedule;
    schedule.p_original = p_original;
    int p = p_original;
    while (p > p_native) {
        if (p <= 5)
            throw std::invalid_argument("precision reducer cannot remove five bits");
        schedule.steps.push_back(PrecisionReduceStep{
            .p_in = p,
            .p_out = p - 5,
            .removed_bits = 5,
            .kind = "HE3DBStyleLowFiveResidueNormalization",
        });
        p -= 5;
    }
    schedule.p_final = p;
    return schedule;
}

inline std::string ScheduleString(const PrecisionReduceSchedule& schedule,
                                  int p_work = 0) {
    std::ostringstream os;
    os << schedule.p_original;
    for (const auto& step : schedule.steps) os << "->" << step.p_out;
    if (p_work > 0 && schedule.p_final < p_work) os << "->zeroextend" << p_work;
    return os.str();
}

inline std::uint64_t SimulateRemoveLowFiveBitsRound(std::uint64_t m) {
    return m >> 5;
}

inline std::uint64_t SimulateHE3DBStylePrecisionReduce(
    std::uint64_t m, int p_original, int& p_out, int p_native = 12) {
    const auto schedule = MakeHE3DBStylePrecisionSchedule(p_original, p_native);
    for (const auto& step : schedule.steps) {
        (void)step;
        m = SimulateRemoveLowFiveBitsRound(m);
    }
    p_out = schedule.p_final;
    return m;
}

template <typename T>
inline bool CheckNegacyclicCompatibility(const std::vector<T>& ext, int N) {
    if (N <= 0 || static_cast<int>(ext.size()) != 2 * N) return false;
    for (int j = 0; j < N; j++)
        if (ext[j + N] != T(0) - ext[j]) return false;
    return true;
}

template <class P>
inline std::vector<typename P::T> NegacyclicExtension(
    const TFHEpp::Polynomial<P>& tv) {
    std::vector<typename P::T> ext(2 * P::n);
    for (int j = 0; j < P::n; j++) {
        ext[j] = tv[j];
        ext[j + P::n] = typename P::T(0) - tv[j];
    }
    return ext;
}

template <class P>
TFHEpp::Polynomial<P> HE3DBMuPoly(typename P::T mu) {
    TFHEpp::Polynomial<P> poly{};
    poly.fill(typename P::T(0) - mu);
    return poly;
}

template <class P>
TFHEpp::Polynomial<P> HE3DBGPoly(int plain_bits, int scale_bits) {
    if (plain_bits <= 0 || plain_bits > P::nbit)
        throw std::invalid_argument("gpolygen plain_bits out of range");
    if (scale_bits < 0 ||
        scale_bits >= std::numeric_limits<typename P::T>::digits)
        throw std::invalid_argument("gpolygen scale_bits out of range");
    TFHEpp::Polynomial<P> poly{};
    const int padding_bits = P::nbit - plain_bits;
    const auto scale =
        static_cast<typename P::T>(typename P::T(1) << scale_bits);
    for (int i = 0; i < P::n; i++)
        poly[i] = scale * static_cast<typename P::T>(i >> padding_bits);
    return poly;
}

template <class inP, class midP, class outP>
inline void RecordReducerPBSStats(const TFHEpp::Polynomial<outP>& tv,
                                  BlindRotatePruneStats* stats) {
    if (!stats) return;
    stats->pbs_calls++;
    stats->pbs_count_reducer++;
    stats->key_switch_count++;
    stats->total += midP::k * midP::n;
    stats->cmux_calls += midP::k * midP::n;
    stats->periods.push_back(ExactNegacyclicPeriod<outP>(tv));
    stats->total_by_pbs.push_back(midP::k * midP::n);
    stats->cmux_by_pbs.push_back(midP::k * midP::n);
    stats->skipped_by_pbs.push_back(0);
    (void)sizeof(inP);
}

// HE3DB-style arithmetic sign PBS for the low-six-bit window. The returned
// ciphertext is arithmetic Boolean encoded as 0 or Q/2 at outP.
template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
PBS_MSB_Weighted_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    typename brP_logari::targetP::T weight_torus,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match reducer PBS domain");
    static_assert(std::is_same_v<inP, outP>,
                  "precision reducer expects lvl2 -> lvl0 -> lvl2");

    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    TFHEpp::TLWE<inP> shifted = ct;
    shifted[inP::k * inP::n] +=
        static_cast<typename inP::T>(typename inP::T(1) << (digits - 7));

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    const auto mu = static_cast<typename outP::T>(weight_torus >> 1);
    if (mu == 0)
        throw std::invalid_argument("weighted MSB PBS requires nonzero weight");
    auto tv = HE3DBMuPoly<outP>(mu);
    RecordReducerPBSStats<inP, midP, outP>(tv, stats);

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += mu;
    return out;
}

template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
PBS_MSB_Arithmetic_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using outP = typename brP_logari::targetP;
    return PBS_MSB_Weighted_Lvl02<iksP, brP_logari>(
        ct, static_cast<typename outP::T>(outP::μ << 2),
        iksk, bkfft_logari, stats);
}

// Identity-style PBS mapping a low-five-bit residue r to r * Q/2^p_in.
template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
PBS_IdentityLowFiveResidue_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& residue_phase,
    int p_in,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match reducer PBS domain");
    static_assert(std::is_same_v<inP, outP>,
                  "precision reducer expects lvl2 -> lvl0 -> lvl2");

    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    if (p_in <= 5 || p_in >= digits)
        throw std::invalid_argument("residue identity PBS requires 6 <= p_in < q");
    TFHEpp::TLWE<inP> shifted = residue_phase;
    shifted[inP::k * inP::n] +=
        static_cast<typename inP::T>(typename inP::T(1) << (digits - 7));

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    auto tv = HE3DBGPoly<outP>(5, digits - p_in);
    const auto ext = NegacyclicExtension<outP>(tv);
    if (!CheckNegacyclicCompatibility(ext, outP::n))
        throw std::logic_error("residue identity LUT failed negacyclic check");
    RecordReducerPBSStats<inP, midP, outP>(tv, stats);

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    return out;
}

template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
RemoveLowFiveBitsRound(
    const TFHEpp::TLWE<typename iksP::domainP>& ct_in,
    int p_in,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<inP, outP>,
                  "precision reducer round must return the input key type");
    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    if (p_in <= 5 || p_in >= digits)
        throw std::invalid_argument("RemoveLowFiveBitsRound requires 6 <= p_in < q");

    // HE3DB's signed low-five-residue PBS is centered: without this public
    // half-bin shift it computes round(m / 32), which wraps Q-epsilon inputs
    // to zero.  Subtracting 16 * Delta_in first turns the centered operation
    // into exact floor(m / 32) modulo the p_out domain.
    TFHEpp::TLWE<inP> adjusted = ct_in;
    const auto half_bin =
        static_cast<typename inP::T>(typename inP::T(1)
                                     << (digits - p_in + 4));
    adjusted[inP::k * inP::n] -= half_bin;

    TFHEpp::TLWE<inP> shift_tlwe{};
    const int shift_bits = p_in - 6;
    for (std::size_t i = 0; i < ct_in.size(); i++)
        shift_tlwe[i] = static_cast<typename inP::T>(adjusted[i] << shift_bits);

    auto sign6 = PBS_MSB_Arithmetic_Lvl02<iksP, brP_logari>(
        shift_tlwe, iksk, bkfft_logari, stats);
    for (std::size_t i = 0; i < shift_tlwe.size(); i++)
        shift_tlwe[i] -= sign6[i];

    auto cleaned = PBS_IdentityLowFiveResidue_Lvl02<iksP, brP_logari>(
        shift_tlwe, p_in, iksk, bkfft_logari, stats);

    TFHEpp::TLWE<outP> out{};
    for (std::size_t i = 0; i < out.size(); i++)
        out[i] = adjusted[i] - cleaned[i];
    return out;
}

template <class iksP, class brP_logari>
PrecisionReduceResult<typename brP_logari::targetP>
HE3DBStylePrecisionReduceToNative(
    const TFHEpp::TLWE<typename iksP::domainP>& ct_in,
    int p_original,
    int p_native,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr,
    std::vector<PrecisionReduceRoundOutput<typename brP_logari::targetP>>*
        round_outputs = nullptr) {
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, outP>,
                  "precision reducer output type must match input type");
    const auto schedule = MakeHE3DBStylePrecisionSchedule(p_original, p_native);

    TFHEpp::TLWE<outP> current = ct_in;
    int p = p_original;
    int round = 0;
    for (const auto& step : schedule.steps) {
        current = RemoveLowFiveBitsRound<iksP, brP_logari>(
            current, step.p_in, iksk, bkfft_logari, stats);
        p = step.p_out;
        if (round_outputs) {
            round_outputs->push_back(PrecisionReduceRoundOutput<outP>{
                .round = round,
                .p_in = step.p_in,
                .p_out = step.p_out,
                .ct_out = current,
            });
        }
        round++;
    }

    std::ostringstream semantic;
    semantic << "Enc_{Q/2^" << p << "}(floor(m/2^"
             << (p_original - p) << "))";
    return PrecisionReduceResult<outP>{
        .ct_out = current,
        .p_out = p,
        .delta_out = ReducerTorusScaleForPrecision<typename outP::T>(p),
        .semantic = semantic.str(),
        .schedule = schedule,
    };
}

}  // namespace MetaPBS2
