#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "my_ethmsb_pruned_bitextract.hpp"

namespace my_ethmsb_pruned_fused_final {

using namespace my_ethmsb_pruned_gap;

inline constexpr const char* FUSED_FINAL_IMPL_NAME =
    "ethmsb_h3compat_pruned_gap_fused_final";

enum class FusedInputRoundingMode {
    None,
    DeltaQuarter,
    DeltaHalf,
    ThreeDeltaQuarter,
    AutoCalibrated
};

struct FusedFinalStats {
    FusedInputRoundingMode rounding_mode = FusedInputRoundingMode::None;
    Torus rounding_offset = 0;
    int p = 0;
    int bit_index_k = 0;
    bool antisymmetry_ok = false;
};

struct FusedGapStats {
    PruneStats bitextract_stats;
    FusedFinalStats fused_final_stats;
    bool conversion_pbs_removed = true;
    int total_pbs_count = 2;
    bool direct_mask_used = false;
    std::uint64_t W = 0;
    std::uint64_t M = 0;
};

inline Torus mask_for_bits(const int p)
{
    my_ethmsb::require_k(p);
    if (p == 64) return ~Torus{0};
    return (Torus{1} << p) - 1;
}

inline bool upper_half(const Torus x) { return (x & Q_HALF) != 0; }

inline bool fused_final_truth_bit(const int p, const int bit_index_k, Torus c)
{
    require_bit_index(p, bit_index_k);
    c &= mask_for_bits(p);
    const Torus b = (c >> (p - 1 - bit_index_k)) & Torus{1};
    const Torus half_plain = Torus{1} << (p - 1);
    const Torus m = b == 0 ? c : (c ^ half_plain);
    const Torus W = bit_weight_msb_index(p, bit_index_k);
    const Torus m_gap = (m - b * W) & mask_for_bits(p);
    const Torus y =
        static_cast<Torus>(Wide{delta_bits(p)} * Wide{m_gap} +
                           Wide{gap_offset_for_bit_index(p, bit_index_k)});
    return upper_half(y);
}

inline bool direct_guard_mask_periodic_possible(const int p,
                                                const int bit_index_k)
{
    return guard_value_for_bit_index(p, bit_index_k) == Q_HALF;
}

inline bool fused_truth_antisymmetry_ok(const int p, const int bit_index_k)
{
    require_bit_index(p, bit_index_k);
    const Torus half_plain = Torus{1} << (p - 1);
    const Torus limit = std::min<Torus>(half_plain, Torus{1} << 16);
    for (Torus c = 0; c < limit; ++c) {
        const bool a = fused_final_truth_bit(p, bit_index_k, c);
        const bool b = fused_final_truth_bit(p, bit_index_k, c ^ half_plain);
        if (a == b) return false;
    }
    return true;
}

inline Torus rounding_offset_for_mode(const int p,
                                      const FusedInputRoundingMode mode)
{
    const Torus d = delta_bits(p);
    switch (mode) {
        case FusedInputRoundingMode::None:
            return 0;
        case FusedInputRoundingMode::DeltaQuarter:
            return d / 4;
        case FusedInputRoundingMode::DeltaHalf:
            return d / 2;
        case FusedInputRoundingMode::ThreeDeltaQuarter:
            return (d / 4) * 3;
        case FusedInputRoundingMode::AutoCalibrated:
            return d / 2;
    }
    throw std::invalid_argument("unknown fused final rounding mode");
}

inline TFHEpp::Polynomial<P2> make_fused_final_poly(const int p,
                                                    const int bit_index_k,
                                                    const Torus out_value)
{
    if ((out_value & Torus{1}) != 0)
        throw std::invalid_argument("fused final out_value must be even");
    const Torus half = out_value / 2;
    TFHEpp::Polynomial<P2> poly;
    for (std::uint64_t j = 0; j < P2::n; ++j) {
        const Torus c = representative_plain_floor_for_br_index(j, p);
        poly[j] = fused_final_truth_bit(p, bit_index_k, c) ? half
                                                           : Torus{0} - half;
    }
    return poly;
}

inline void pbs_fused_clear_bit_msb_l2_to_l2(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& combined_ct, const int p,
    const int bit_index_k, const Torus out_value,
    const FusedInputRoundingMode rounding_mode,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    FusedFinalStats* stats)
{
    require_bit_index(p, bit_index_k);
    TFHEpp::TLWE<P2> in = combined_ct;
    const Torus rounding = rounding_offset_for_mode(p, rounding_mode);
    my_ethmsb::add_const_inplace<P2>(in, rounding);

    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20);

    const TFHEpp::Polynomial<P2> testvector =
        make_fused_final_poly(p, bit_index_k, out_value);
    alignas(64) TFHEpp::TLWE<P2> tmp;
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(tmp, in0, bkfft02,
                                                testvector);
    out = tmp;
    my_ethmsb::add_const_inplace<P2>(out, out_value / 2);

    if (stats != nullptr) {
        stats->rounding_mode = rounding_mode;
        stats->rounding_offset = rounding;
        stats->p = p;
        stats->bit_index_k = bit_index_k;
        stats->antisymmetry_ok = fused_truth_antisymmetry_ok(p, bit_index_k);
    }
}

inline void gapmsb_pruned_bitk_fused_final_l2_to_l2(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& ct, const int p,
    const int bit_index_k, const Torus out_value,
    const FusedInputRoundingMode rounding_mode,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02, FusedGapStats* stats)
{
    require_bit_index(p, bit_index_k);
    alignas(64) TFHEpp::TLWE<P2> bit_qhalf;
    PruneStats bit_stats;
    bit_extract_pruned_qhalf_l2_to_l2(bit_qhalf, ct, p, bit_index_k, iksk20,
                                      bkfft02, &bit_stats);

    alignas(64) TFHEpp::TLWE<P2> combined;
    my_ethmsb::add<P2>(combined, ct, bit_qhalf);

    FusedFinalStats final_stats;
    pbs_fused_clear_bit_msb_l2_to_l2(out, combined, p, bit_index_k, out_value,
                                     rounding_mode, iksk20, bkfft02,
                                     &final_stats);

    if (stats != nullptr) {
        stats->bitextract_stats = bit_stats;
        stats->fused_final_stats = final_stats;
        stats->conversion_pbs_removed = true;
        stats->total_pbs_count = 2;
        stats->direct_mask_used = false;
        stats->W = bit_weight_msb_index(p, bit_index_k);
        stats->M = bit_period_plain(p, bit_index_k);
    }
}

inline void homcomp_lt_pruned_gap_fused_final(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const FusedInputRoundingMode rounding_mode,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    FusedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, a, b);
    gapmsb_pruned_bitk_fused_final_l2_to_l2(out, diff, t + 1, bit_index_k,
                                            out_value, rounding_mode, iksk20,
                                            bkfft02, stats);
}

inline void homcomp_gt_pruned_gap_fused_final(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const FusedInputRoundingMode rounding_mode,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    FusedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, b, a);
    gapmsb_pruned_bitk_fused_final_l2_to_l2(out, diff, t + 1, bit_index_k,
                                            out_value, rounding_mode, iksk20,
                                            bkfft02, stats);
}

}  // namespace my_ethmsb_pruned_fused_final
