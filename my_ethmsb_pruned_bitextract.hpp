#pragma once

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "detwfa.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace my_ethmsb_pruned_gap {

using P0 = my_ethmsb_params::my_h3_lvl0param;
using P2 = my_ethmsb_params::my_h3_lvl2param;
using KS20 = my_ethmsb_params::my_h3_lvl20param;
using BR02 = my_ethmsb_params::my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

static_assert(std::numeric_limits<Torus>::digits == 64);
static_assert(std::is_same_v<P2::T, Torus>);
static_assert(BR02::Addends == 1,
              "pruned BitExtract currently requires unbundled BR02 keys");

inline constexpr Torus BOOL_ONE = Torus{1} << 61;
inline constexpr Torus Q_HALF = Torus{1} << 63;
inline constexpr Torus Q_QUARTER = Torus{1} << 62;
inline constexpr std::uint64_t BR_CYCLE = 2ULL * P2::n;
inline constexpr int BR_INDEX_BITS = P2::nbit + 1;
inline constexpr Torus DEFAULT_SWITCH_BAND_WIDTH =
    (Torus{1} << (64 - BR_INDEX_BITS)) / 8;
inline constexpr const char* SINGLE_ROUND_IMPL_NAME =
    "ethmsb_h3compat_pruned_gap_single_round";

struct PruneStats {
    std::uint64_t total_terms = 0;
    std::uint64_t executed_terms = 0;
    std::uint64_t zero_skipped = 0;
    std::uint64_t periodic_skipped = 0;
    std::string disabled_reason;
    std::uint64_t M_k = 0;
    std::uint64_t w_k = 0;
    std::uint64_t rotation_period_index = 0;
    std::uint64_t rotation_weight_index = 0;
    int p = 0;
    int bit_index_k = 0;
    bool pruning_enabled = false;
    bool bit_exact_in_br_index = false;
};

struct PrunedGapStats {
    PruneStats bitextract_stats;
    int conversion_pbs_count = 1;
    int final_pbs_count = 1;
    int total_pbs_count = 3;
    std::string gap_offset_hex;
    std::string guard_value_hex;
    std::uint64_t w_k = 0;
    std::uint64_t M_k = 0;
    bool pruning_enabled = false;
    std::uint64_t periodic_skipped = 0;
    std::uint64_t zero_skipped = 0;
    std::uint64_t executed_terms = 0;
};

struct BitIndexInfo {
    std::uint64_t w_k = 0;
    std::uint64_t M_k = 0;
    std::uint64_t rotation_weight_index = 0;
    std::uint64_t rotation_period_index = 0;
    bool exact_in_br_index = false;
    bool supports_period_pruning = false;
    std::string disabled_reason;
};

struct PrunedSweepCandidate {
    int p = 0;
    int bit_index_k = 0;
    double boundary_failure_rate = 1.0;
    double random_failure_rate = 1.0;
    double random_upper95 = 1.0;
    double observed_prune_ratio = 0.0;
    bool selected = false;
};

struct PrunedSweepSummary {
    std::vector<PrunedSweepCandidate> candidates;
};

inline std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

inline Torus centered_distance(const Torus a, const Torus b)
{
    return my_ethmsb::torus_abs_centered(a - b);
}

inline Torus delta_bits(const int p) { return my_ethmsb::delta_bits(p); }

inline Torus q_half() { return Q_HALF; }

inline Torus q_quarter() { return Q_QUARTER; }

inline bool decrypt_arith_bit_by_distance(const Torus phase,
                                           const Torus out_value)
{
    return centered_distance(phase, out_value) < centered_distance(phase, 0);
}

inline bool closer_to_value(const Torus phase, const Torus value)
{
    return decrypt_arith_bit_by_distance(phase, value);
}

inline bool closer_to_qhalf(const Torus phase)
{
    return closer_to_value(phase, Q_HALF);
}

inline Torus phase_l2(const TFHEpp::TLWE<P2>& ct,
                      const TFHEpp::Key<P2>& key)
{
    return TFHEpp::tlweSymPhase<P2>(ct, key);
}

inline TFHEpp::TLWE<P2> trivial_encrypt_phase_l2(const Torus phase)
{
    return my_ethmsb_params::trivial_tlwe<P2>(phase);
}

inline Torus encode_unsigned_phase(const Torus m, const int p)
{
    return static_cast<Torus>(Wide{m} * Wide{my_ethmsb::delta_bits(p)});
}

inline void require_bit_index(const int p, const int bit_index_k)
{
    my_ethmsb::require_k(p);
    if (bit_index_k <= 0 || bit_index_k >= p)
        throw std::invalid_argument(
            "bit_index_k must satisfy 1 <= k <= p-1; k=0 is forbidden");
    if (p > 63)
        throw std::invalid_argument("p > 63 would overflow 64-bit shifts");
}

inline std::uint64_t bit_weight_msb_index(const int p, const int bit_index_k)
{
    require_bit_index(p, bit_index_k);
    const int shift = p - 1 - bit_index_k;
    if (shift < 0 || shift >= 63)
        throw std::invalid_argument("bit weight shift is outside uint64_t");
    return std::uint64_t{1} << shift;
}

inline std::uint64_t bit_period_msb_index(const int p, const int bit_index_k)
{
    require_bit_index(p, bit_index_k);
    const int shift = p - bit_index_k;
    if (shift < 0 || shift >= 63)
        throw std::invalid_argument("bit period shift is outside uint64_t");
    return std::uint64_t{1} << shift;
}

inline std::uint64_t bit_period_plain(const int p, const int bit_index_k)
{
    return bit_period_msb_index(p, bit_index_k);
}

inline Torus guard_value_for_bit_index(const int p, const int bit_index_k)
{
    const std::uint64_t w = bit_weight_msb_index(p, bit_index_k);
    return static_cast<Torus>(Wide{my_ethmsb::delta_bits(p)} * Wide{w});
}

inline Torus gap_offset_for_bit_index(const int p, const int bit_index_k)
{
    const std::uint64_t w = bit_weight_msb_index(p, bit_index_k);
    return static_cast<Torus>((Wide{w + 1} *
                               Wide{my_ethmsb::delta_bits(p)}) /
                              Wide{2});
}

inline BitIndexInfo bit_index_info(const int p, const int bit_index_k)
{
    BitIndexInfo info;
    info.w_k = bit_weight_msb_index(p, bit_index_k);
    info.M_k = bit_period_msb_index(p, bit_index_k);

    if (p <= BR_INDEX_BITS) {
        const int shift = BR_INDEX_BITS - p;
        if ((shift >= 64) || (info.w_k > (std::numeric_limits<std::uint64_t>::max() >> shift)) ||
            (info.M_k > (std::numeric_limits<std::uint64_t>::max() >> shift)))
            throw std::overflow_error("BR index period overflow");
        info.rotation_weight_index = info.w_k << shift;
        info.rotation_period_index = info.M_k << shift;
        info.exact_in_br_index = true;
    }
    else {
        const int shift = p - BR_INDEX_BITS;
        if (shift >= 63) {
            info.disabled_reason = "p_minus_BR_INDEX_BITS_too_large";
            return info;
        }
        const std::uint64_t divisor = std::uint64_t{1} << shift;
        if ((info.w_k % divisor) == 0 && (info.M_k % divisor) == 0) {
            info.rotation_weight_index = info.w_k / divisor;
            info.rotation_period_index = info.M_k / divisor;
            info.exact_in_br_index = true;
        }
        else {
            info.disabled_reason = "bit_not_representable_in_2N_BR_index";
            return info;
        }
    }

    if (info.rotation_period_index == 0) {
        info.disabled_reason = "zero_rotation_period";
    }
    else if ((BR_CYCLE % info.rotation_period_index) != 0) {
        info.disabled_reason = "rotation_period_not_dividing_2N";
    }
    else {
        info.supports_period_pruning = true;
    }
    return info;
}

inline std::uint64_t representative_plain_for_br_index(const std::uint64_t j,
                                                       const int p)
{
    my_ethmsb::require_k(p);
    if (p <= BR_INDEX_BITS) {
        const int shift = BR_INDEX_BITS - p;
        if (shift == 0) return j & ((std::uint64_t{1} << p) - 1);
        const std::uint64_t rounded =
            (j + (std::uint64_t{1} << (shift - 1))) >> shift;
        return rounded & ((std::uint64_t{1} << p) - 1);
    }
    const int shift = p - BR_INDEX_BITS;
    if (shift >= 63)
        throw std::invalid_argument("p is too large for representative index");
    return j << shift;
}

inline std::uint64_t representative_plain_floor_for_br_index(
    const std::uint64_t j, const int p)
{
    my_ethmsb::require_k(p);
    if (p <= BR_INDEX_BITS) {
        const int shift = BR_INDEX_BITS - p;
        return (j >> shift) & ((std::uint64_t{1} << p) - 1);
    }
    const int shift = p - BR_INDEX_BITS;
    if (shift >= 63)
        throw std::invalid_argument("p is too large for representative index");
    return j << shift;
}

inline bool bit_value_for_br_index(const std::uint64_t j, const int p,
                                   const int bit_index_k)
{
    require_bit_index(p, bit_index_k);
    const std::uint64_t m = representative_plain_for_br_index(j, p);
    const int shift = p - 1 - bit_index_k;
    return ((m >> shift) & std::uint64_t{1}) != 0;
}

inline bool bit_value_floor_for_br_index(const std::uint64_t j, const int p,
                                         const int bit_index_k)
{
    require_bit_index(p, bit_index_k);
    const std::uint64_t m = representative_plain_floor_for_br_index(j, p);
    const int shift = p - 1 - bit_index_k;
    return ((m >> shift) & std::uint64_t{1}) != 0;
}

enum class BitExtractGridMode {
    FloorGridRepresentative,
    NearestGridRepresentative
};

inline TFHEpp::Polynomial<P2> make_bitextract_qhalf_poly(
    const int p, const int bit_index_k,
    const BitExtractGridMode mode = BitExtractGridMode::FloorGridRepresentative)
{
    (void)bit_index_info(p, bit_index_k);
    TFHEpp::Polynomial<P2> poly;
    for (std::uint64_t j = 0; j < P2::n; ++j)
        poly[j] = (mode == BitExtractGridMode::NearestGridRepresentative
                       ? bit_value_for_br_index(j, p, bit_index_k)
                       : bit_value_floor_for_br_index(j, p, bit_index_k))
                      ? Q_HALF
                      : 0;
    return poly;
}

inline TFHEpp::Polynomial<P2> make_bitextract_qhalf_floor_poly(
    const int p, const int bit_index_k)
{
    (void)bit_index_info(p, bit_index_k);
    TFHEpp::Polynomial<P2> poly;
    for (std::uint64_t j = 0; j < P2::n; ++j)
        poly[j] =
            bit_value_floor_for_br_index(j, p, bit_index_k) ? Q_HALF : 0;
    return poly;
}

inline TFHEpp::Polynomial<P2> make_direct_guard_mask_poly_experimental(
    const int p, const int bit_index_k, const Torus guard_value)
{
    require_bit_index(p, bit_index_k);
    if (guard_value != Q_HALF)
        throw std::logic_error(
            "direct periodic 0/guard_value LUT is impossible under standard "
            "negacyclic BR unless guard_value == Q/2");
    return make_bitextract_qhalf_poly(p, bit_index_k);
}

inline bool polynomial_rotation_equal(const TFHEpp::Polynomial<P2>& poly,
                                      const std::uint64_t a_mod_2n)
{
    TFHEpp::Polynomial<P2> rotated;
    TFHEpp::PolynomialMulByXai<P2>(rotated, poly,
                                   static_cast<P2::T>(a_mod_2n % BR_CYCLE));
    return rotated == poly;
}

inline void fill_stats_header(PruneStats& stats, const int p,
                              const int bit_index_k, const BitIndexInfo& info)
{
    stats = {};
    stats.p = p;
    stats.bit_index_k = bit_index_k;
    stats.w_k = info.w_k;
    stats.M_k = info.M_k;
    stats.rotation_weight_index = info.rotation_weight_index;
    stats.rotation_period_index = info.rotation_period_index;
    stats.bit_exact_in_br_index = info.exact_in_br_index;
    stats.pruning_enabled = info.supports_period_pruning;
    stats.disabled_reason = info.disabled_reason;
}

inline void blind_rotate_pruned_l0_to_l2(
    TFHEpp::TRLWE<P2>& acc, const TFHEpp::TLWE<P0>& in0,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft,
    const TFHEpp::Polynomial<P2>& testvector,
    const BitIndexInfo& info, const bool enable_periodic_pruning,
    PruneStats* stats)
{
    TFHEpp::ModswitchTLWE<P0> moded;
    TFHEpp::BRModSwitch<BR02, 1>(moded, in0);

    acc = {};
    TFHEpp::PolynomialMulByXai<P2>(acc[P2::k], testvector,
                                   moded[P0::k * P0::n] % BR_CYCLE);

    if (stats != nullptr) stats->total_terms = P0::k * P0::n;
    for (int i = 0; i < P0::k * P0::n; ++i) {
        const std::uint64_t a_mod =
            static_cast<std::uint64_t>(moded[i]) % BR_CYCLE;
        if (a_mod == 0) {
            if (stats != nullptr) ++stats->zero_skipped;
            continue;
        }
        if (enable_periodic_pruning && info.supports_period_pruning &&
            (a_mod % info.rotation_period_index) == 0) {
            if (stats != nullptr) ++stats->periodic_skipped;
            continue;
        }
        if (stats != nullptr) ++stats->executed_terms;
        TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<BR02>(
            acc, bkfft[i], static_cast<int>(a_mod));
    }
}

// CMUX pruning is a semantic-preserving optimization: if the public
// mod-switched rotation keeps the self-negating LUT invariant, skipping that
// CMUX leaves this PBS plaintext function unchanged. The single-round
// BitExtract LUT itself remains boundary-sensitive near switching points.

inline void gate_bootstrapping_pruned_l0_to_l2(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P0>& in0,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft,
    const TFHEpp::Polynomial<P2>& testvector,
    const BitIndexInfo& info, const bool enable_periodic_pruning,
    PruneStats* stats)
{
    alignas(64) TFHEpp::TRLWE<P2> acc;
    blind_rotate_pruned_l0_to_l2(acc, in0, bkfft, testvector, info,
                                 enable_periodic_pruning, stats);
    TFHEpp::SampleExtractIndex<P2>(out, acc, 0);
}

inline void bit_extract_unpruned_reference_qhalf_l2_to_l2(
    TFHEpp::TLWE<P2>& out_qhalf, const TFHEpp::TLWE<P2>& in, const int p,
    const int bit_index_k, const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
{
    (void)bit_index_info(p, bit_index_k);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20);
    const TFHEpp::Polynomial<P2> testvector =
        make_bitextract_qhalf_poly(p, bit_index_k);
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out_qhalf, in0, bkfft02,
                                                testvector);
}

inline void bit_extract_no_periodic_qhalf_l2_to_l2(
    TFHEpp::TLWE<P2>& out_qhalf, const TFHEpp::TLWE<P2>& in, const int p,
    const int bit_index_k, const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02, PruneStats* stats)
{
    const BitIndexInfo info = bit_index_info(p, bit_index_k);
    if (stats != nullptr) fill_stats_header(*stats, p, bit_index_k, info);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20);
    const TFHEpp::Polynomial<P2> testvector =
        make_bitextract_qhalf_poly(p, bit_index_k);
    gate_bootstrapping_pruned_l0_to_l2(out_qhalf, in0, bkfft02, testvector,
                                       info, false, stats);
}

inline void bit_extract_pruned_qhalf_l2_to_l2(
    TFHEpp::TLWE<P2>& out_qhalf, const TFHEpp::TLWE<P2>& in, const int p,
    const int bit_index_k, const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02, PruneStats* stats)
{
    const BitIndexInfo info = bit_index_info(p, bit_index_k);
    if (stats != nullptr) fill_stats_header(*stats, p, bit_index_k, info);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20);
    const TFHEpp::Polynomial<P2> testvector =
        make_bitextract_qhalf_poly(p, bit_index_k);
    gate_bootstrapping_pruned_l0_to_l2(out_qhalf, in0, bkfft02, testvector,
                                       info, true, stats);
}

inline void bit_extract_qhalf_with_input_offset_l2_to_l2(
    TFHEpp::TLWE<P2>& out_qhalf, const TFHEpp::TLWE<P2>& in, const int p,
    const int bit_index_k, const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02, PruneStats* stats)
{
    const BitIndexInfo info = bit_index_info(p, bit_index_k);
    if (stats != nullptr) fill_stats_header(*stats, p, bit_index_k, info);
    alignas(64) TFHEpp::TLWE<P2> shifted = in;
    my_ethmsb::add_const_inplace<P2>(shifted, my_ethmsb::delta_bits(p) / 2);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, shifted, iksk20);
    const TFHEpp::Polynomial<P2> testvector =
        make_bitextract_qhalf_floor_poly(p, bit_index_k);
    gate_bootstrapping_pruned_l0_to_l2(out_qhalf, in0, bkfft02, testvector,
                                       info, true, stats);
}

inline void bit_extract_qhalf_with_input_offset_no_periodic_l2_to_l2(
    TFHEpp::TLWE<P2>& out_qhalf, const TFHEpp::TLWE<P2>& in, const int p,
    const int bit_index_k, const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02, PruneStats* stats)
{
    const BitIndexInfo info = bit_index_info(p, bit_index_k);
    if (stats != nullptr) fill_stats_header(*stats, p, bit_index_k, info);
    alignas(64) TFHEpp::TLWE<P2> shifted = in;
    my_ethmsb::add_const_inplace<P2>(shifted, my_ethmsb::delta_bits(p) / 2);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, shifted, iksk20);
    const TFHEpp::Polynomial<P2> testvector =
        make_bitextract_qhalf_floor_poly(p, bit_index_k);
    gate_bootstrapping_pruned_l0_to_l2(out_qhalf, in0, bkfft02, testvector,
                                       info, false, stats);
}

inline void bit_extract_unpruned_reference_qhalf_with_input_offset_l2_to_l2(
    TFHEpp::TLWE<P2>& out_qhalf, const TFHEpp::TLWE<P2>& in, const int p,
    const int bit_index_k, const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
{
    (void)bit_index_info(p, bit_index_k);
    alignas(64) TFHEpp::TLWE<P2> shifted = in;
    my_ethmsb::add_const_inplace<P2>(shifted, my_ethmsb::delta_bits(p) / 2);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, shifted, iksk20);
    const TFHEpp::Polynomial<P2> testvector =
        make_bitextract_qhalf_floor_poly(p, bit_index_k);
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out_qhalf, in0, bkfft02,
                                                testvector);
}

inline void pbs_msb_value_l2_to_l2(
    TFHEpp::TLWE<P2>& out, TFHEpp::TLWE<P2> in, const int input_bits,
    const Torus offset_l2, const Torus out_value_l2,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
{
    my_ethmsb::require_k(input_bits);
    if ((out_value_l2 & Torus{1}) != 0)
        throw std::invalid_argument(
            "pbs_msb_value_l2_to_l2 requires even out_value");

    my_ethmsb::add_const_inplace<P2>(in, offset_l2);
    alignas(64) TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20);

    const Torus half = out_value_l2 / 2;
    TFHEpp::Polynomial<P2> testvector;
    testvector.fill(Torus{0} - half);

    alignas(64) TFHEpp::TLWE<P2> tmp;
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(tmp, in0, bkfft02,
                                                testvector);
    out = tmp;
    my_ethmsb::add_const_inplace<P2>(out, half);
}

inline void bool_qhalf_to_value_l2_to_l2(
    TFHEpp::TLWE<P2>& out_value_ct, const TFHEpp::TLWE<P2>& bit_qhalf_ct,
    const Torus guard_value,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
{
    if ((guard_value & Torus{1}) != 0)
        throw std::invalid_argument("guard_value must be even");
    pbs_msb_value_l2_to_l2(out_value_ct, bit_qhalf_ct, 1, Q_QUARTER,
                           guard_value, iksk20, bkfft02);
}

inline void gapmsb_pruned_bitk_l2_to_l2(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& ct, const int p,
    const int bit_index_k, const Torus out_value,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02, PrunedGapStats* stats)
{
    require_bit_index(p, bit_index_k);
    const std::uint64_t w = bit_weight_msb_index(p, bit_index_k);
    const std::uint64_t M = bit_period_msb_index(p, bit_index_k);
    const Torus guard_value = guard_value_for_bit_index(p, bit_index_k);
    const Torus gap_offset = gap_offset_for_bit_index(p, bit_index_k);

    alignas(64) TFHEpp::TLWE<P2> bit_qhalf;
    PruneStats bit_stats;
    bit_extract_pruned_qhalf_l2_to_l2(bit_qhalf, ct, p, bit_index_k, iksk20,
                                      bkfft02, &bit_stats);

    alignas(64) TFHEpp::TLWE<P2> guard_mask;
    bool_qhalf_to_value_l2_to_l2(guard_mask, bit_qhalf, guard_value, iksk20,
                                 bkfft02);

    alignas(64) TFHEpp::TLWE<P2> ct_gap;
    my_ethmsb::sub<P2>(ct_gap, ct, guard_mask);

    pbs_msb_value_l2_to_l2(out, ct_gap, p, gap_offset, out_value, iksk20,
                           bkfft02);

    if (stats != nullptr) {
        stats->bitextract_stats = bit_stats;
        stats->conversion_pbs_count = 1;
        stats->final_pbs_count = 1;
        stats->total_pbs_count = 3;
        stats->gap_offset_hex = hex64(gap_offset);
        stats->guard_value_hex = hex64(guard_value);
        stats->w_k = w;
        stats->M_k = M;
        stats->pruning_enabled = bit_stats.pruning_enabled;
        stats->periodic_skipped = bit_stats.periodic_skipped;
        stats->zero_skipped = bit_stats.zero_skipped;
        stats->executed_terms = bit_stats.executed_terms;
    }
}

// Gap offset only increases the final MSB margin after BitExtract and
// BoolToValue conversion are correct. This single-round variant is therefore a
// probabilistic comparator and must be profiled with empirical failure rates.

inline std::optional<int> choose_pruned_gap_bit_index(
    const int p, const double target_failure_rate,
    const PrunedSweepSummary* sweep_summary)
{
    my_ethmsb::require_k(p);
    if (sweep_summary == nullptr) return std::nullopt;
    std::optional<PrunedSweepCandidate> best;
    for (const PrunedSweepCandidate& c : sweep_summary->candidates) {
        if (c.p != p || !c.selected) continue;
        if (c.boundary_failure_rate > target_failure_rate) continue;
        if (c.random_upper95 > target_failure_rate) continue;
        if (!best.has_value() ||
            c.observed_prune_ratio > best->observed_prune_ratio)
            best = c;
    }
    if (!best.has_value()) return std::nullopt;
    return best->bit_index_k;
}

inline int default_zero_bit_index_for_p(const int p)
{
    (void)p;
    throw std::runtime_error(
        "no correctness-supported default k for this p under single-round "
        "pruned variant");
}

inline int parse_zero_bit_choice(const std::string& choice, const int p)
{
    if (choice == "default") return default_zero_bit_index_for_p(p);
    if (choice == "p-5") return std::max(1, p - 5);
    if (choice == "p-4") return std::max(1, p - 4);
    if (choice == "p-3") return std::max(1, p - 3);
    if (choice == "p-2") return std::max(1, p - 2);
    if (choice == "p-1") return p - 1;
    return std::stoi(choice);
}

inline void homcomp_lt_pruned_gap(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    PrunedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, a, b);
    gapmsb_pruned_bitk_l2_to_l2(out, diff, t + 1, bit_index_k, out_value,
                                iksk20, bkfft02, stats);
}

inline void homcomp_gt_pruned_gap(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    PrunedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, b, a);
    gapmsb_pruned_bitk_l2_to_l2(out, diff, t + 1, bit_index_k, out_value,
                                iksk20, bkfft02, stats);
}

inline void homcomp_eq_pruned_gap(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    PrunedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    alignas(64) TFHEpp::TLWE<P2> is_lt;
    alignas(64) TFHEpp::TLWE<P2> is_gt;
    PrunedGapStats lt_stats;
    homcomp_lt_pruned_gap(is_lt, a, b, t, bit_index_k, iksk20, bkfft02,
                          stats != nullptr ? &lt_stats : nullptr, out_value);
    homcomp_gt_pruned_gap(is_gt, a, b, t, bit_index_k, iksk20, bkfft02,
                          nullptr, out_value);
    alignas(64) TFHEpp::TLWE<P2> neq;
    my_ethmsb::add<P2>(neq, is_lt, is_gt);
    my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(out_value), neq);
    if (stats != nullptr) {
        *stats = lt_stats;
        stats->total_pbs_count = 6;
    }
}

inline void homcomp_lt_pruned_gap_single_round(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    PrunedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    homcomp_lt_pruned_gap(out, a, b, t, bit_index_k, iksk20, bkfft02, stats,
                          out_value);
}

inline void homcomp_gt_pruned_gap_single_round(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    PrunedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    homcomp_gt_pruned_gap(out, a, b, t, bit_index_k, iksk20, bkfft02, stats,
                          out_value);
}

inline void homcomp_eq_pruned_gap_single_round(
    TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
    const TFHEpp::TLWE<P2>& b, const int t, const int bit_index_k,
    const TFHEpp::KeySwitchingKey<KS20>& iksk20,
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02,
    PrunedGapStats* stats = nullptr, const Torus out_value = BOOL_ONE)
{
    homcomp_eq_pruned_gap(out, a, b, t, bit_index_k, iksk20, bkfft02, stats,
                          out_value);
}

}  // namespace my_ethmsb_pruned_gap
