#pragma once
// =============================================================
// hom_compare.hpp — Homomorphic comparison operators (optimized).
//
// Template parameters:
//   brP_metapbs: lvl2→lvl2 for MetaPBS BitExtract
//   brP_logari:  lvl0→lvl2 for LOG_to_ARI (636 CMUXes)
//   brP_base:    lvl0→lvl1 for base case sign PBS
//   iksP:        lvl2→lvl0 for dimension reduction
// =============================================================

#include <cstdint>
#include <vector>

#include "metapbs2/gap_msb.hpp"

namespace MetaPBS2 {

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomGreaterThan(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
    int plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options = HomMSBOptions{},
    BlindRotatePruneStats* prune_stats = nullptr) {
    using domP = typename brP_metapbs::domainP;
    constexpr int n1 = domP::k * domP::n + 1;
    TFHEpp::TLWE<domP> sub{};
    for (int i = 0; i < n1; i++) sub[i] = c1[i] - c0[i];
    return HomMSB<brP_metapbs, brP_logari, brP_base, iksP>(
        sub, plain_bits + 1, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomGreaterThan(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
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
    return HomGreaterThan<brP_metapbs, brP_logari, brP_base, iksP>(
        c0, c1, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomLessThan(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
    int plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options = HomMSBOptions{},
    BlindRotatePruneStats* prune_stats = nullptr) {
    using domP = typename brP_metapbs::domainP;
    constexpr int n1 = domP::k * domP::n + 1;
    TFHEpp::TLWE<domP> sub{};
    for (int i = 0; i < n1; i++) sub[i] = c0[i] - c1[i];
    return HomMSB<brP_metapbs, brP_logari, brP_base, iksP>(
        sub, plain_bits + 1, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomLessThan(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
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
    return HomLessThan<brP_metapbs, brP_logari, brP_base, iksP>(
        c0, c1, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomGreaterThanEqual(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
    int plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options = HomMSBOptions{},
    BlindRotatePruneStats* prune_stats = nullptr) {
    auto lt = HomLessThan<brP_metapbs, brP_logari, brP_base, iksP>(
        c0, c1, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
    for (size_t i = 0; i < lt.size(); i++) lt[i] = -lt[i];
    return lt;
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomGreaterThanEqual(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
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
    return HomGreaterThanEqual<brP_metapbs, brP_logari, brP_base, iksP>(
        c0, c1, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomLessThanEqual(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
    int plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const HomMSBOptions& options = HomMSBOptions{},
    BlindRotatePruneStats* prune_stats = nullptr) {
    auto gt = HomGreaterThan<brP_metapbs, brP_logari, brP_base, iksP>(
        c0, c1, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
    for (size_t i = 0; i < gt.size(); i++) gt[i] = -gt[i];
    return gt;
}

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP>
HomLessThanEqual(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c0,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& c1,
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
    return HomLessThanEqual<brP_metapbs, brP_logari, brP_base, iksP>(
        c0, c1, plain_bits, bkfft, trkeys, cfg,
        bkfft_logari, bkfft_base, iksk, options, prune_stats);
}

}  // namespace MetaPBS2
