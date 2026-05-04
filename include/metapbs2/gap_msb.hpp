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

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "metapbs2/bit_extraction.hpp"
#include "tlwe.hpp"

namespace MetaPBS2 {

struct GapMSBOptions {
    int p;       // Message precision, cfg.t must be 2^p.
    int k;       // Chapter-3/MSB-side bit to clear. Current binary encoding requires 1 <= k < p.
    bool enable_periodic_pruning = true;
    int period = 0;  // Optional override for the first BitExtract blind-rotate period.
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

    const auto bit_ct = BitExtract<brP>(
        ct, bkfft, trkeys, cfg,
        BitExtractOptions{
            .p = options.p,
            .k = options.k,
            .enable_periodic_pruning = options.enable_periodic_pruning,
            .period = options.period,
        },
        prune_stats);

    return LogicalBitToArithmeticWeight<brP>(
        bit_ct, ArithmeticWeightForChapterBit<typename tgtP::T>(options.p, options.k), bkfft);
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

    std::vector<int> redundancies(cfg.K + 1);
    std::vector<int> deltas(cfg.K);
    redundancies[0] = 2 * N / cfg.t;
    for (int rnd = 0; rnd < cfg.K; rnd++) {
        deltas[rnd] = DeltaOffset(redundancies[rnd], cfg.rounds[rnd].beta);
        redundancies[rnd + 1] = redundancies[rnd] * cfg.rounds[rnd].beta;
    }

    TFHEpp::TLWE<domP> cquo0, crem0;
    HomDivRemLWE<domP>(cquo0, crem0, shifted, static_cast<typename domP::T>(2 * N));

    TFHEpp::TRLWE<tgtP> C0;
    BlindRotateGLWEFromQuotient<brP>(
        C0, cquo0, MakeAccumulator<tgtP>(tv), bkfft, N);

    int current_mod = 2 * N;
    TFHEpp::TRLWE<tgtP> prev_GLWE = C0;
    TFHEpp::TLWE<domP> prev_crem = crem0;

    for (int rnd = 0; rnd < cfg.K; rnd++) {
        const auto& round = cfg.rounds[rnd];
        const int delta_k = deltas[rnd];

        TFHEpp::TRLWE<tgtP> Ck_prime;
        HomTruncRepeatShifted<tgtP>(Ck_prime, prev_GLWE,
                                     -round.T, round.T, round.beta, delta_k,
                                     trkeys[rnd]);

        TFHEpp::TLWE<domP> cquo_k, crem_k;
        HomDivRemAtScale<domP>(cquo_k, crem_k, prev_crem, current_mod, round.beta);

        TFHEpp::TRLWE<tgtP> Ck;
        BlindRotateGLWEFromQuotient<brP>(Ck, cquo_k, Ck_prime, bkfft, N);

        current_mod *= round.beta;
        prev_GLWE = Ck;
        prev_crem = crem_k;
    }

    TFHEpp::TLWE<tgtP> out;
    SampleExtractIndex0<tgtP>(out, prev_GLWE);
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

    auto weighted_bit =
        ExtractWeightedChapterBit<brP>(ct, bkfft, trkeys, cfg, options, prune_stats);

    TFHEpp::TLWE<tgtP> ct_gap{};
    ClearChapterBitAssign<brP>(ct_gap, ct, weighted_bit);
    return SignMetaPBSWithOffset<brP>(
        ct_gap,
        HalfGapOffsetForChapterBit<typename tgtP::T>(options.p, options.k),
        bkfft, trkeys, cfg);
}

// HE3DB-style recursive shell: reduce a plaintext precision by repeatedly
// applying GapMSB's clear-bit transform. This version keeps the current
// MetaPBS2 paper-row constraint, so practical encrypted tests use p=log2(2N).
template <class brP>
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
            ExtractWeightedChapterBit<brP>(current, bkfft, trkeys, cfg,
                                           round_options, prune_stats);
        ClearChapterBitAssign<brP>(current, current, weighted_bit);
        if (round_options.k + 1 < round_options.p)
            round_options.k++;
    }
    return SignMetaPBSWithOffset<brP>(
        current,
        HalfGapOffsetForChapterBit<typename tgtP::T>(options.p, options.k),
        bkfft, trkeys, cfg);
}

}  // namespace MetaPBS2
