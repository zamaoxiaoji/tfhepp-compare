#pragma once
// =============================================================
// gap_msb.hpp - GapMSB built on the MetaPBS2 BitExtract wrapper.
//
// This follows Chapter_3_revised (1).tex, alg:et-hmsb-expanded:
//   1. ct_bool <- BitExtractBoolPruned(ct, BK, p, k), with bit output
//      encoded as 0/Q/2 under the target/original key.
//   2. Convert that BoolHalf ciphertext to arithmetic weight
//      bit_k(m) * w_k * Delta, where w_k = 2^{p-1-k}, Delta = Q/2^p.
//   3. Clear the bit by ct_gap <- ct - weighted_bit.
//   4. Add offset' = (w_k + 1) Delta / 2 to the body.
//   5. Run a sign/MSB PBS.
//
// The logical-to-arithmetic conversion plays the same role as HE3DB's
// LOG_to_ARI, but this implementation uses a direct LUT to produce the
// exact weighted arithmetic contribution needed by GapMSB. It does not use
// HE3DB's identity-bootstrap right-shift reduction path.
// =============================================================

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "metapbs2/bit_extraction.hpp"
#include "metapbs2/precision_reducer.hpp"
#include "tlwe.hpp"

namespace MetaPBS2 {

enum class WeightedBitMode {
    CheckedDirectThenLogical,
    DirectOnly,
    LogicalOnly,
};

struct GapMSBOptions {
    int p;       // Message precision, cfg.t must be 2^p.
    int k;       // Chapter-3/MSB-side bit to clear. Current binary encoding requires 1 <= k < p.
    int kappa = 5;  // Current one-shot GapMSB supports the final kappa-bit window.
    bool enable_periodic_pruning = true;
    int period = 0;  // Optional override for the first BitExtract blind-rotate period.
    WeightedBitMode weighted_mode = WeightedBitMode::CheckedDirectThenLogical;
};

struct HomMSBOptions {
    int kappa = 5;
    bool enable_periodic_pruning = true;
    int period = 0;
    WeightedBitMode weighted_mode = WeightedBitMode::CheckedDirectThenLogical;
};

template <typename T>
struct PrecisionState {
    int p = 0;      // Current effective message precision.
    T delta = 0;    // Current server-side torus scale, Q / 2^p.
    int round = 0;  // Recursive/MetaPBS round that produced this state.
};

template <typename T>
struct GapMSBRoundParams {
    int p_original = 0;
    PrecisionState<T> current{};
    PrecisionState<T> next{};
    int k_current = 0;
    std::uint64_t w_k = 0;
    T A_k = 0;
    T offset = 0;
};

template <typename T>
inline T TorusScaleForPrecision(int p) {
    static_assert(std::is_unsigned_v<T>, "TorusScaleForPrecision expects unsigned torus");
    const int digits = std::numeric_limits<T>::digits;
    if (p <= 0 || p >= digits)
        throw std::invalid_argument("message precision out of torus range");
    return T(1) << (digits - p);
}

inline int DefaultGapBitForPrecision(int p) {
    if (p <= 1)
        throw std::invalid_argument("GapMSB requires p >= 2");
    const int k = p >= 6 ? p - 5 : std::max(1, p - 2);
    if (k <= 0 || k >= p)
        throw std::invalid_argument("default GapMSB bit is outside current precision");
    return k;
}

template <typename T>
inline T StandardMSBOffsetForPrecision(int p) {
    static_assert(std::is_unsigned_v<T>, "StandardMSBOffsetForPrecision expects unsigned torus");
    return TorusScaleForPrecision<T>(p) >> 1;
}

template <typename T>
inline T GapExtraOffsetForChapterBit(int p, int k) {
    static_assert(std::is_unsigned_v<T>, "GapExtraOffsetForChapterBit expects unsigned torus");
    if (k <= 0)
        throw std::invalid_argument("GapMSB cannot clear the Chapter/MSB top bit");
    const int lsb_k = LSBIndexFromChapterBit(p, k);
    const unsigned __int128 delta = TorusScaleForPrecision<T>(p);
    const unsigned __int128 w = static_cast<unsigned __int128>(1) << lsb_k;
    return static_cast<T>((w * delta) >> 1);
}

template <typename T>
inline T HalfGapOffsetForChapterBit(int p, int k) {
    static_assert(std::is_unsigned_v<T>, "HalfGapOffsetForChapterBit expects unsigned torus");
    return StandardMSBOffsetForPrecision<T>(p) +
           GapExtraOffsetForChapterBit<T>(p, k);
}

template <typename T>
inline T ArithmeticWeightForChapterBit(int p, int k) {
    static_assert(std::is_unsigned_v<T>, "ArithmeticWeightForChapterBit expects unsigned torus");
    if (k <= 0)
        throw std::invalid_argument("GapMSB cannot clear the Chapter/MSB top bit");
    const int lsb_k = LSBIndexFromChapterBit(p, k);
    const unsigned __int128 delta = TorusScaleForPrecision<T>(p);
    const unsigned __int128 w = static_cast<unsigned __int128>(1) << lsb_k;
    return static_cast<T>(w * delta);
}

template <typename T>
inline GapMSBRoundParams<T> MakeGapMSBRoundParams(
    int p_original, int p_current, int round, int k_current = 0,
    int p_next = 0) {
    static_assert(std::is_unsigned_v<T>, "MakeGapMSBRoundParams expects unsigned torus");
    if (p_current <= 1)
        throw std::invalid_argument("current precision must be at least two bits");
    if (k_current == 0)
        k_current = DefaultGapBitForPrecision(p_current);
    if (k_current <= 0 || k_current >= p_current)
        throw std::invalid_argument("GapMSB requires 1 <= k_current < p_current");
    if (p_next == 0)
        p_next = p_current;
    const int lsb_k = LSBIndexFromChapterBit(p_current, k_current);
    const auto delta_current = TorusScaleForPrecision<T>(p_current);
    const auto delta_next = TorusScaleForPrecision<T>(p_next);
    const std::uint64_t w_k = std::uint64_t{1} << lsb_k;
    return GapMSBRoundParams<T>{
        .p_original = p_original,
        .current = PrecisionState<T>{.p = p_current, .delta = delta_current, .round = round},
        .next = PrecisionState<T>{.p = p_next, .delta = delta_next, .round = round + 1},
        .k_current = k_current,
        .w_k = w_k,
        .A_k = ArithmeticWeightForChapterBit<T>(p_current, k_current),
        .offset = HalfGapOffsetForChapterBit<T>(p_current, k_current),
    };
}

template <class P>
inline int DecodeSignPhase(typename P::T phase) {
    return static_cast<std::make_signed_t<typename P::T>>(phase) < 0 ? 1 : 0;
}

template <class P>
int DecodeSignCout(
    const TFHEpp::TLWE<P>& cout,
    const TFHEpp::Key<P>& key) {
    return DecodeSignPhase<P>(TFHEpp::tlweSymPhase<P>(cout, key));
}

// Convert a BoolHalf ciphertext encoded as 0/Q/2 into a WeightedBit ciphertext
// encoded as 0 or weight_torus. The public +Q/4 shift places bit 0 at Q/4 and
// bit 1 at 3Q/4, both Q/4 away from the sign boundaries. The sign LUT is
// negacyclic-valid and maps [0,Q/2) to -A/2 and [Q/2,Q) to +A/2; adding A/2
// after key-compatible PBS yields 0/A.
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
BoolToWeightPBS(
    const TFHEpp::TLWE<typename brP::domainP>& bit_ct,
    typename brP::targetP::T weight_torus,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    BlindRotatePruneStats* stats = nullptr) {
    using domP = typename brP::domainP;
    using tgtP = typename brP::targetP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "BoolToWeightPBS currently expects same-domain bootstrap parameters");

    TFHEpp::TLWE<domP> shifted = bit_ct;
    constexpr int digits = std::numeric_limits<typename domP::T>::digits;
    const auto q_over_4 =
        static_cast<typename domP::T>(typename domP::T(1) << (digits - 2));
    shifted[domP::k * domP::n] += q_over_4;

    auto half_weight = static_cast<typename tgtP::T>(weight_torus >> 1);
    TFHEpp::Polynomial<tgtP> tv{};
    tv.fill(typename tgtP::T(0) - half_weight);
    if (stats) {
        stats->pbs_calls++;
        stats->pbs_count_gapmsb++;
        stats->pbs_count_bool_to_weight++;
        stats->key_switch_count++;
        stats->total += domP::k * domP::n;
        stats->cmux_calls += domP::k * domP::n;
        stats->periods.push_back(ExactNegacyclicPeriod<tgtP>(tv));
        stats->total_by_pbs.push_back(domP::k * domP::n);
        stats->cmux_by_pbs.push_back(domP::k * domP::n);
        stats->skipped_by_pbs.push_back(0);
    }

    TFHEpp::TLWE<tgtP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP>(out, shifted, bkfft, tv);
    out[tgtP::k * tgtP::n] += half_weight;
    return out;
}

template <class brP>
[[deprecated("Use BoolToWeightPBS for the explicit BoolHalf -> WeightedBit conversion")]]
TFHEpp::TLWE<typename brP::targetP>
LogicalBitToArithmeticWeight(
    const TFHEpp::TLWE<typename brP::domainP>& bit_ct,
    typename brP::targetP::T weight_torus,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft) {
    return BoolToWeightPBS<brP>(bit_ct, weight_torus, bkfft);
}

template <class brP>
TFHEpp::TLWE<typename brP::targetP>
ExtractWeightedChapterBit(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const GapMSBOptions& options,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "ExtractWeightedChapterBit currently expects same-domain bootstrap parameters");

    const int lsb_k = LSBIndexFromChapterBit(options.p, options.k);
    const auto weight =
        ArithmeticWeightForChapterBit<typename tgtP::T>(options.p, options.k);
    const BitExtractOptions bit_options{
        .p = options.p,
        .k = options.k,
        .enable_periodic_pruning = options.enable_periodic_pruning,
        .period = options.period,
    };

    const bool direct_ok =
        WeightedBitNegacyclicCompatible<typename tgtP::T>(lsb_k, tgtP::n, weight);
    if (options.weighted_mode != WeightedBitMode::LogicalOnly) {
        if (direct_ok)
            return BitExtractWeighted<brP>(
                ct, bkfft, trkeys, cfg, bit_options, weight, prune_stats);
        if (options.weighted_mode == WeightedBitMode::DirectOnly)
            throw std::invalid_argument("direct weighted bit LUT violates negacyclic encoding");
    }

    const auto before_bitextract = prune_stats ? prune_stats->pbs_calls : 0;
    const auto bit_ct = BitExtractBoolPruned<brP>(
        ct, bkfft, trkeys, cfg, bit_options, prune_stats);
    if (prune_stats) {
        const auto delta = prune_stats->pbs_calls - before_bitextract;
        prune_stats->pbs_count_gapmsb += delta;
        prune_stats->pbs_count_bit_extract += delta;
    }

    return BoolToWeightPBS<brP>(bit_ct, weight, bkfft, prune_stats);
}

template <class brP>
void ClearChapterBitAssign(
    TFHEpp::TLWE<typename brP::targetP>& ct_gap,
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::TLWE<typename brP::targetP>& weighted_bit) {
    using domP = typename brP::domainP;
    using tgtP = typename brP::targetP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "ClearChapterBitAssign currently expects same-domain bootstrap parameters");
    for (size_t i = 0; i < ct_gap.size(); i++)
        ct_gap[i] = ct[i] - weighted_bit[i];
}

template <class brP>
TFHEpp::TLWE<typename brP::targetP>
SignMetaPBSWithOffset(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    typename brP::domainP::T offset,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using domP = typename brP::domainP;
    using tgtP = typename brP::targetP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "SignMetaPBSWithOffset currently expects same-domain bootstrap parameters");
    constexpr int N = tgtP::n;
    if (cfg.t != 2 * N)
        throw std::invalid_argument("SignMetaPBSWithOffset currently supports only cfg.t=2N");
    if (static_cast<int>(trkeys.size()) < cfg.K)
        throw std::invalid_argument("SignMetaPBSWithOffset requires one TruncRepeat key per round");

    TFHEpp::TLWE<domP> shifted = ct;
    shifted[domP::k * domP::n] += offset;

    TFHEpp::Polynomial<tgtP> tv{};
    tv.fill(static_cast<typename tgtP::T>(tgtP::μ));

    const auto before = prune_stats ? prune_stats->pbs_calls : 0;
    auto out = RunAlgorithm1WithTV<brP>(shifted, tv, bkfft, trkeys, cfg, 0, prune_stats);
    if (prune_stats) {
        const auto delta = prune_stats->pbs_calls - before;
        prune_stats->pbs_count_gapmsb += delta;
        prune_stats->pbs_count_final_msb += delta;
    }
    return out;
}

template <class brP>
TFHEpp::TLWE<typename brP::targetP>
GapMSB(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const GapMSBOptions& options,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "GapMSB currently expects same-domain bootstrap parameters");

    if (cfg.t != CheckedPowerOfTwo(options.p))
        throw std::invalid_argument("GapMSB requires cfg.t == 2^p");
    if (options.k <= 0 || options.k >= options.p)
        throw std::invalid_argument("GapMSB requires 1 <= k < p");
    if (options.kappa <= 0 || options.k < options.p - options.kappa)
        throw std::invalid_argument("GapMSB requires k inside the final kappa-bit window");

    auto weighted_bit =
        ExtractWeightedChapterBit<brP>(ct, bkfft, trkeys, cfg, options, prune_stats);

    TFHEpp::TLWE<tgtP> ct_gap{};
    ClearChapterBitAssign<brP>(ct_gap, ct, weighted_bit);
    return SignMetaPBSWithOffset<brP>(
        ct_gap,
        HalfGapOffsetForChapterBit<typename tgtP::T>(options.p, options.k),
        bkfft, trkeys, cfg, prune_stats);
}

// Compatibility wrapper for older tests. New recursive comparisons should call
// HomMSB, which returns the HE3DB-style lvl1 sign ciphertext via a base PBS.
template <class brP>
[[deprecated("Use HomMSB for recursive comparison")]]
TFHEpp::TLWE<typename brP::targetP>
RecursiveGapMSB(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const GapMSBOptions& options,
    int rounds,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    using domP = typename brP::domainP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "RecursiveGapMSB currently expects same-domain bootstrap parameters");
    if (rounds <= 0)
        return GapMSB<brP>(ct, bkfft, trkeys, cfg, options, prune_stats);
    (void)bkfft;
    (void)trkeys;
    (void)cfg;
    (void)options;
    (void)prune_stats;
    throw std::invalid_argument(
        "RecursiveGapMSB is disabled because it cleared bits while keeping the "
        "ciphertext at the original torus scale. Delta-prime must be produced "
        "by a real PBS/MetaPBS LUT re-encoding step before it can be used.");
}

// =============================================================
// NaiveSignPBS_Lvl01: base-case MSB extraction at lvl0→lvl1.
//
// Matches HE3DB's MSBGateBootstrapping(TLWELvl1, TLWELvl2, ek):
//   1. Add offset to shift decision boundary
//   2. IKS lvl2→lvl0  (dimension reduction: 2048 → 636)
//   3. GateBootstrapping lvl0→lvl1  (N=1024, fast!)
//
// Output: TLWE<lvl1param> sign ciphertext.
// This avoids the expensive lvl2→lvl2 PBS (N=2048) that was 12x slower.
// =============================================================
template <class iksP, class brP_base>
TFHEpp::TLWE<typename brP_base::targetP>
NaiveSignPBS_Lvl01(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    typename iksP::domainP::T offset,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base) {
    using inP = typename iksP::domainP;      // lvl2param
    using midP = typename iksP::targetP;     // lvl0param
    using outP = typename brP_base::targetP; // lvl1param
    static_assert(std::is_same_v<midP, typename brP_base::domainP>,
                  "IKS target must match base BK domain");

    TFHEpp::TLWE<inP> shifted = ct;
    shifted[inP::k * inP::n] += offset;

    // IKS: lvl2 → lvl0
    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    TFHEpp::Polynomial<outP> tv{};
    tv.fill(static_cast<typename outP::T>(outP::μ));

    // GateBootstrapping: lvl0 → lvl1 (N=1024)
    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_base>(out, tlwe_mid, bkfft_base, tv);
    return out;
}

template <class iksP, class brP_base>
TFHEpp::TLWE<typename brP_base::targetP>
NaiveSignPBS_Lvl01(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    int p,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base) {
    using inP = typename iksP::domainP;
    constexpr int digits_in = std::numeric_limits<typename inP::T>::digits;
    if (p <= 0 || p >= digits_in)
        throw std::invalid_argument("NaiveSignPBS_Lvl01 precision out of range");
    const auto offset = static_cast<typename inP::T>(
        typename inP::T(1) << (digits_in - p - 1));
    return NaiveSignPBS_Lvl01<iksP, brP_base>(
        ct, offset, iksk, bkfft_base);
}

// =============================================================
// BoolToWeightPBS_Lvl02: lightweight BoolHalf -> WeightedBit conversion.
//
// IKS lvl2→lvl0 + GateBS lvl0→lvl2 (636 CMUXes, not 2048!)
// Output stays at lvl2 for direct subtraction.
// =============================================================
template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
BoolToWeightPBS_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& bit_ct,
    typename brP_logari::targetP::T weight_torus,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match LOG_to_ARI BK domain");

    // Shift: 0/Q/2 -> Q/4/3Q/4, both Q/4 away from sign boundaries.
    TFHEpp::TLWE<inP> shifted = bit_ct;
    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    shifted[inP::k * inP::n] +=
        static_cast<typename inP::T>(typename inP::T(1) << (digits - 2));

    // IKS: lvl2 → lvl0
    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    // Sign TV: [0,Q/2) -> -A/2, [Q/2,Q) -> +A/2.
    auto half_w = static_cast<typename outP::T>(weight_torus >> 1);
    TFHEpp::Polynomial<outP> tv{};
    tv.fill(typename outP::T(0) - half_w);
    if (stats) {
        stats->pbs_calls++;
        stats->pbs_count_gapmsb++;
        stats->pbs_count_bool_to_weight++;
        stats->key_switch_count++;
        stats->total += midP::k * midP::n;
        stats->cmux_calls += midP::k * midP::n;
        stats->periods.push_back(ExactNegacyclicPeriod<outP>(tv));
        stats->total_by_pbs.push_back(midP::k * midP::n);
        stats->cmux_by_pbs.push_back(midP::k * midP::n);
        stats->skipped_by_pbs.push_back(0);
    }

    // GateBS: lvl0 → lvl2 (636 CMUXes)
    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += half_w;
    return out;
}

template <class iksP, class brP_logari>
[[deprecated("Use BoolToWeightPBS_Lvl02 for the explicit BoolHalf -> WeightedBit conversion")]]
TFHEpp::TLWE<typename brP_logari::targetP>
LogicalBitToArithmeticWeight_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& bit_ct,
    typename brP_logari::targetP::T weight_torus,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari) {
    return BoolToWeightPBS_Lvl02<iksP, brP_logari>(
        bit_ct, weight_torus, iksk, bkfft_logari);
}

template <class brP_metapbs, class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_metapbs::targetP>
ExtractWeightedChapterBitFast(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& ct,
    int encoding_p,
    int chapter_bit,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const HomMSBOptions& options,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using domP = typename brP_metapbs::domainP;
    using tgtP = typename brP_metapbs::targetP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "ExtractWeightedChapterBitFast expects same-domain MetaPBS parameters");
    static_assert(std::is_same_v<typename iksP::domainP, tgtP>,
                  "LOG_to_ARI IKS domain must match MetaPBS output");
    static_assert(std::is_same_v<typename iksP::targetP, typename brP_logari::domainP>,
                  "LOG_to_ARI IKS target must match LOG_to_ARI BK domain");
    static_assert(std::is_same_v<typename brP_logari::targetP, tgtP>,
                  "LOG_to_ARI BK target must match MetaPBS output");

    if (options.kappa <= 0)
        throw std::invalid_argument("HomMSBOptions.kappa must be positive");

    const int p_metapbs = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    if (chapter_bit <= 0 || chapter_bit >= p_metapbs)
        throw std::invalid_argument("HomMSB chapter bit is outside the MetaPBS window");
    const int lsb_k = LSBIndexFromChapterBit(p_metapbs, chapter_bit);
    const auto weight =
        ArithmeticWeightForChapterBit<typename tgtP::T>(encoding_p, chapter_bit);
    const BitExtractOptions bit_options{
        .p = p_metapbs,
        .k = chapter_bit,
        .enable_periodic_pruning = options.enable_periodic_pruning,
        .period = options.period,
    };

    const bool direct_ok =
        WeightedBitNegacyclicCompatible<typename tgtP::T>(lsb_k, tgtP::n, weight);
    if (options.weighted_mode != WeightedBitMode::LogicalOnly) {
        if (direct_ok)
            return BitExtractWeighted<brP_metapbs>(
                ct, bkfft, trkeys, cfg, bit_options, weight, prune_stats);
        if (options.weighted_mode == WeightedBitMode::DirectOnly)
            throw std::invalid_argument("direct weighted bit LUT violates negacyclic encoding");
    }

    auto bit_ct = BitExtractBoolPruned<brP_metapbs>(
        ct, bkfft, trkeys, cfg, bit_options, prune_stats);
    return BoolToWeightPBS_Lvl02<iksP, brP_logari>(
        bit_ct, weight, iksk, bkfft_logari, prune_stats);
}

// =============================================================
// HomGapMSB: recursive MSB with lightweight LOG_to_ARI + base case.
//
// brP_metapbs: lvl2→lvl2 (MetaPBS BitExtract, 2048 CMUXes/round)
// brP_logari:  lvl0→lvl2 (LOG_to_ARI, 636 CMUXes — 3.2x faster)
// brP_base:    lvl0→lvl1 (base case, N=1024)
// iksP:        lvl2→lvl0 (dimension reduction)
// =============================================================
template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomGapMSB(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& ct,
    int p,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP_metapbs::targetP;
    using domP = typename brP_metapbs::domainP;
    static_assert(std::is_same_v<domP, tgtP>);
    static_assert(std::is_same_v<typename brP_logari::targetP, tgtP>);

    if (options.kappa <= 0)
        throw std::invalid_argument("HomMSBOptions.kappa must be positive");
    if (p <= 0)
        throw std::invalid_argument("HomMSB requires positive plaintext precision");
    if (p <= options.kappa) {
        const auto offset = StandardMSBOffsetForPrecision<typename domP::T>(p);
        return NaiveSignPBS_Lvl01<iksP, brP_base>(
            ct, offset, iksk, bkfft_base);
    }
    (void)bkfft;
    (void)trkeys;
    (void)cfg;
    (void)bkfft_logari;
    (void)prune_stats;
    throw std::invalid_argument(
        "HomGapMSB recursive Delta-prime path is disabled: the previous "
        "implementation reduced the logical precision p without a PBS LUT "
        "that re-encoded the ciphertext to Q/2^p_next. Use the one-shot "
        "HomMSB/GapMSB path with cfg.t == 2^p, or add an explicit "
        "PrecisionState re-encoding round.");
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomGapMSBAtOriginalScale(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& ct,
    int encoding_p,
    int remaining_p,
    int chapter_bit,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP_metapbs::targetP;
    using domP = typename brP_metapbs::domainP;
    static_assert(std::is_same_v<domP, tgtP>);
    static_assert(std::is_same_v<typename brP_logari::targetP, tgtP>);

    (void)ct;
    (void)encoding_p;
    (void)remaining_p;
    (void)chapter_bit;
    (void)bkfft;
    (void)trkeys;
    (void)cfg;
    (void)bkfft_logari;
    (void)bkfft_base;
    (void)iksk;
    (void)options;
    (void)prune_stats;
    throw std::invalid_argument(
        "HomGapMSBAtOriginalScale is disabled: it intentionally kept the "
        "original Delta while reducing the logical remaining precision, which "
        "is not a valid Delta-prime re-encoding semantics.");
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomMSB(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& ct,
    int plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options = HomMSBOptions{},
    BlindRotatePruneStats* prune_stats = nullptr) {
    using inP = typename brP_metapbs::domainP;
    static_assert(std::is_same_v<inP, typename brP_metapbs::targetP>,
                  "HomMSB safe path expects same-domain MetaPBS parameters");
    if (options.kappa <= 0)
        throw std::invalid_argument("HomMSBOptions.kappa must be positive");
    if (plain_bits <= 0)
        throw std::invalid_argument("HomMSB requires positive plaintext precision");

    const int p_native = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    TFHEpp::TLWE<inP> ct_work = ct;
    bool invert_after_gap = false;
    if (plain_bits > p_native) {
        // Move the signed MSB boundary away from the torus wrap before the
        // HE3DB-style precision reducer.  The final GapMSB computes the MSB of
        // this shifted value at p_native, then the lvl1 Boolean is inverted to
        // recover the original signed MSB.
        ct_work[inP::k * inP::n] +=
            static_cast<typename inP::T>(
                typename inP::T(1)
                << (std::numeric_limits<typename inP::T>::digits - 1));
        auto reduced = HE3DBStylePrecisionReduceToNative<iksP, brP_logari>(
            ct_work, plain_bits, p_native, iksk, bkfft_logari, prune_stats);
        ct_work = reduced.ct_out;
        invert_after_gap = true;
    }

    // For plain_bits < p_native, and for reducer outputs with p_final <
    // p_native, the physical torus phase is unchanged:
    // m * Q/2^p == (m << (p_native - p)) * Q/2^p_native.
    // This is legal zero-extension. The final Chapter-3 GapMSB is always
    // defined at the native/current precision, not at the original precision.
    const int p_work = p_native;
    const int chapter_bit = p_work - options.kappa;
    if (chapter_bit <= 0 || chapter_bit >= p_work)
        throw std::invalid_argument("HomMSBOptions.kappa selects an invalid native GapMSB bit");
    GapMSBOptions gap_options{
        .p = p_work,
        .k = chapter_bit,
        .kappa = options.kappa,
        .enable_periodic_pruning = options.enable_periodic_pruning,
        .period = options.period,
        .weighted_mode = options.weighted_mode,
    };
    auto lvl2_sign = GapMSB<brP_metapbs>(
        ct_work, bkfft, trkeys, cfg, gap_options, prune_stats);

    // GapMSB already performed the aligned sign decision and returns a lvl2 sign
    // ciphertext. This final base PBS is only a key/level-compatible refresh to
    // the public HomMSB lvl1 output type; no extra arithmetic offset is applied.
    auto out = NaiveSignPBS_Lvl01<iksP, brP_base>(
        lvl2_sign, typename inP::T(0), iksk, bkfft_base);
    if (invert_after_gap) {
        using outP = typename brP_base::targetP;
        out[outP::k * outP::n] +=
            static_cast<typename outP::T>(
                typename outP::T(1)
                << (std::numeric_limits<typename outP::T>::digits - 1));
    }
    return out;
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomMSB(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& ct,
    int plain_bits, int kappa,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    BlindRotatePruneStats* prune_stats = nullptr) {
    HomMSBOptions options;
    options.kappa = kappa;
    return HomMSB<brP_metapbs, brP_logari, brP_base, iksP>(
        ct, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

}  // namespace MetaPBS2
