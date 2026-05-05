#pragma once
// =============================================================
// gap_msb.hpp - GapMSB built on the MetaPBS2 BitExtract wrapper.
//
// This follows Chapter_3_revised (1).tex, alg:et-hmsb-expanded:
//   1. ct_k <- BitExtract(ct, BK, p, k), with bit output encoded as 0/Q/2.
//   2. Convert that logical bit ciphertext to arithmetic weight
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
inline T TorusScaleForPrecision(int p) {
    static_assert(std::is_unsigned_v<T>, "TorusScaleForPrecision expects unsigned torus");
    const int digits = std::numeric_limits<T>::digits;
    if (p <= 0 || p >= digits)
        throw std::invalid_argument("message precision out of torus range");
    return T(1) << (digits - p);
}

template <typename T>
inline T HalfGapOffsetForChapterBit(int p, int k) {
    static_assert(std::is_unsigned_v<T>, "HalfGapOffsetForChapterBit expects unsigned torus");
    if (k <= 0)
        throw std::invalid_argument("GapMSB cannot clear the Chapter/MSB top bit");
    const int lsb_k = LSBIndexFromChapterBit(p, k);
    const T delta = TorusScaleForPrecision<T>(p);
    const T w = T(1) << lsb_k;
    return (w + T(1)) * (delta >> 1);
}

template <typename T>
inline T ArithmeticWeightForChapterBit(int p, int k) {
    static_assert(std::is_unsigned_v<T>, "ArithmeticWeightForChapterBit expects unsigned torus");
    if (k <= 0)
        throw std::invalid_argument("GapMSB cannot clear the Chapter/MSB top bit");
    const int lsb_k = LSBIndexFromChapterBit(p, k);
    const T delta = TorusScaleForPrecision<T>(p);
    return (T(1) << lsb_k) * delta;
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

// Convert a logical bit ciphertext encoded as 0/Q/2 into an arithmetic
// contribution encoded as 0 or weight_torus. The input is first centered by
// subtracting Q/4, then a sign PBS maps the two stable regions to 0/weight.
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
LogicalBitToArithmeticWeight(
    const TFHEpp::TLWE<typename brP::domainP>& bit_ct,
    typename brP::targetP::T weight_torus,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft) {
    using domP = typename brP::domainP;
    using tgtP = typename brP::targetP;
    static_assert(std::is_same_v<domP, tgtP>,
                  "LogicalBitToArithmeticWeight currently expects same-domain bootstrap parameters");

    TFHEpp::TLWE<domP> centered = bit_ct;
    constexpr int digits = std::numeric_limits<typename domP::T>::digits;
    const auto q_over_4 =
        static_cast<typename domP::T>(typename domP::T(1) << (digits - 2));
    centered[domP::k * domP::n] -= q_over_4;

    auto half_weight = static_cast<typename tgtP::T>(weight_torus >> 1);
    TFHEpp::Polynomial<tgtP> tv{};
    tv.fill(half_weight);

    TFHEpp::TLWE<tgtP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP>(out, centered, bkfft, tv);
    out[tgtP::k * tgtP::n] += half_weight;
    return out;
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

    const auto bit_ct = BitExtract<brP>(
        ct, bkfft, trkeys, cfg, bit_options, prune_stats);

    return LogicalBitToArithmeticWeight<brP>(bit_ct, weight, bkfft);
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
    const Algorithm1Config& cfg) {
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

    return RunAlgorithm1WithTV<brP>(shifted, tv, bkfft, trkeys, cfg);
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

    auto weighted_bit =
        ExtractWeightedChapterBit<brP>(ct, bkfft, trkeys, cfg, options, prune_stats);

    TFHEpp::TLWE<tgtP> ct_gap{};
    ClearChapterBitAssign<brP>(ct_gap, ct, weighted_bit);
    return SignMetaPBSWithOffset<brP>(
        ct_gap,
        HalfGapOffsetForChapterBit<typename tgtP::T>(options.p, options.k),
        bkfft, trkeys, cfg);
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

    TFHEpp::TLWE<tgtP> current = ct;
    GapMSBOptions round_options = options;
    for (int r = 0; r < rounds; r++) {
        auto weighted_bit =
            ExtractWeightedChapterBit<brP>(
                current, bkfft, trkeys, cfg, round_options, prune_stats);
        TFHEpp::TLWE<tgtP> cleared{};
        ClearChapterBitAssign<brP>(cleared, current, weighted_bit);
        current = cleared;
        if (round_options.k + 1 < round_options.p)
            round_options.k++;
    }
    return SignMetaPBSWithOffset<brP>(
        current,
        HalfGapOffsetForChapterBit<typename tgtP::T>(options.p, options.k),
        bkfft, trkeys, cfg);
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
// LogicalBitToArithmeticWeight_Lvl02: lightweight LOG_to_ARI.
//
// IKS lvl2→lvl0 + GateBS lvl0→lvl2 (636 CMUXes, not 2048!)
// Output stays at lvl2 for direct subtraction.
// =============================================================
template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
LogicalBitToArithmeticWeight_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& bit_ct,
    typename brP_logari::targetP::T weight_torus,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match LOG_to_ARI BK domain");

    // Center: 0/Q/2 → -Q/4/+Q/4
    TFHEpp::TLWE<inP> centered = bit_ct;
    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    centered[inP::k * inP::n] -=
        static_cast<typename inP::T>(typename inP::T(1) << (digits - 2));

    // IKS: lvl2 → lvl0
    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, centered, iksk);

    // Sign TV: ±half_weight
    auto half_w = static_cast<typename outP::T>(weight_torus >> 1);
    TFHEpp::Polynomial<outP> tv{};
    tv.fill(half_w);

    // GateBS: lvl0 → lvl2 (636 CMUXes)
    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += half_w;
    return out;
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

    auto bit_ct = BitExtract<brP_metapbs>(
        ct, bkfft, trkeys, cfg, bit_options, prune_stats);
    return LogicalBitToArithmeticWeight_Lvl02<iksP, brP_logari>(
        bit_ct, weight, iksk, bkfft_logari);
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
    if (p <= options.kappa)
        return NaiveSignPBS_Lvl01<iksP, brP_base>(ct, p, iksk, bkfft_base);

    HomMSBOptions round_options = options;
    const int chapter_bit = options.kappa;
    auto weighted = ExtractWeightedChapterBitFast<brP_metapbs, iksP, brP_logari>(
        ct, p, chapter_bit, bkfft, trkeys, cfg, iksk, bkfft_logari,
        round_options, prune_stats);

    TFHEpp::TLWE<tgtP> ct_gap{};
    ClearChapterBitAssign<brP_metapbs>(ct_gap, ct, weighted);

    return HomGapMSB<brP_metapbs, brP_logari, brP_base, iksP>(
        ct_gap, p - (options.kappa - 1), bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
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

    if (options.kappa <= 0)
        throw std::invalid_argument("HomMSBOptions.kappa must be positive");
    if (encoding_p <= 0 || remaining_p <= 0)
        throw std::invalid_argument("HomMSB requires positive plaintext precision");
    if (remaining_p <= options.kappa) {
        if (chapter_bit <= 0 || chapter_bit >= encoding_p)
            return NaiveSignPBS_Lvl01<iksP, brP_base>(
                ct, encoding_p, iksk, bkfft_base);
        const auto offset =
            HalfGapOffsetForChapterBit<typename domP::T>(encoding_p, chapter_bit);
        return NaiveSignPBS_Lvl01<iksP, brP_base>(
            ct, offset, iksk, bkfft_base);
    }
    if (chapter_bit <= 0 || chapter_bit >= encoding_p)
        throw std::invalid_argument("HomMSB chapter bit is outside the encoded message");

    auto weighted = ExtractWeightedChapterBitFast<brP_metapbs, iksP, brP_logari>(
        ct, encoding_p, chapter_bit, bkfft, trkeys, cfg, iksk, bkfft_logari,
        options, prune_stats);

    TFHEpp::TLWE<tgtP> ct_gap{};
    ClearChapterBitAssign<brP_metapbs>(ct_gap, ct, weighted);

    return HomGapMSBAtOriginalScale<brP_metapbs, brP_logari, brP_base, iksP>(
        ct_gap, encoding_p, remaining_p - (options.kappa - 1),
        chapter_bit + (options.kappa - 1), bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
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
    return HomGapMSBAtOriginalScale<brP_metapbs, brP_logari, brP_base, iksP>(
        ct, plain_bits, plain_bits, options.kappa, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
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
