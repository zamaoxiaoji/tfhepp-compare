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

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
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
    bool use_lightweight_gap_pbs = true;
};

template <typename T>
struct PrecisionState {
    int p = 0;      // Current effective message precision.
    T delta = 0;    // Current server-side torus scale, Q / 2^p.
    int round = 0;  // Recursive/MetaPBS round that produced this state.
};

struct ResidualInterval {
    double lo = 0.0;
    double hi = 0.0;
};

template <typename T>
struct RecursivePrecisionState {
    int p = 0;
    T delta = 0;
    ResidualInterval rho{};
    int round = 0;
    std::string semantic;
};

template <typename T>
struct GapRoundParams {
    int round = 0;
    int p_in = 0;
    int p_out = 0;
    int s = 0;
    int k = 0;
    std::uint64_t w = 0;
    T delta_in = 0;
    T delta_out = 0;
    T A = 0;
    T offset = 0;
    ResidualInterval rho_in{};
    ResidualInterval rho_out{};
    double residual_margin = 0.0;
};

struct SimRoundState {
    int round = 0;
    int p_in = 0;
    int p_out = 0;
    int s = 0;
    std::uint64_t q = 0;
    std::uint64_t r = 0;
    int bit = 0;
    double rho_in = 0.0;
    double rho_out = 0.0;
};

struct SimResult {
    std::uint64_t m_original = 0;
    int p_original = 0;
    int p_out = 0;
    std::uint64_t m_work = 0;
    double rho_work = 0.0;
    std::vector<SimRoundState> rounds;
    bool msb_ok = false;
    double min_rho = 0.0;
    double max_rho = 0.0;
    double min_final_margin_units = 0.0;
    std::string schedule;
};

template <class P>
struct RecursiveGapRoundOutput {
    int round = 0;
    int p_in = 0;
    int p_out = 0;
    int s = 0;
    ResidualInterval rho_out{};
    TFHEpp::TLWE<P> ct_out{};
};

template <class P>
struct RecursiveGapReduceResult {
    TFHEpp::TLWE<P> ct_out{};
    int p_out = 0;
    typename P::T delta_out = 0;
    ResidualInterval rho_out{};
    std::vector<GapRoundParams<typename P::T>> rounds;
    std::vector<RecursiveGapRoundOutput<P>> round_outputs;
    std::string path_kind;
    std::string schedule;
    std::string semantic;
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

inline std::string RecursiveGapScaleScheduleString(int p_original,
                                                   int p_target = 12) {
    if (p_original <= 0 || p_target <= 0)
        throw std::invalid_argument("recursive schedule requires positive precisions");
    std::ostringstream os;
    os << p_original;
    int p = p_original;
    while (p > p_target) {
        const int s = std::min(5, p - p_target);
        p -= s;
        os << "->" << p;
    }
    if (p_original < p_target) os << "->zeroextend" << p_target;
    return os.str();
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
inline T GapScaleWindowWeight(int p_cur, int s) {
    static_assert(std::is_unsigned_v<T>, "GapScaleWindowWeight expects unsigned torus");
    if (s <= 0 || s > 5 || p_cur <= s)
        throw std::invalid_argument("invalid GapScale window");
    const unsigned __int128 delta = TorusScaleForPrecision<T>(p_cur);
    const unsigned __int128 w = static_cast<unsigned __int128>(1) << (s - 1);
    return static_cast<T>(w * delta);
}

template <typename T>
inline T GapScaleWindowOffset(int p_cur, int s) {
    static_assert(std::is_unsigned_v<T>, "GapScaleWindowOffset expects unsigned torus");
    if (s <= 0 || s > 5 || p_cur <= s)
        throw std::invalid_argument("invalid GapScale window");
    const unsigned __int128 delta = TorusScaleForPrecision<T>(p_cur);
    const unsigned __int128 w = static_cast<unsigned __int128>(1) << (s - 1);
    return static_cast<T>(((w + 1) * delta + 1) >> 1);
}

inline ResidualInterval UpdateGapScaleResidual(ResidualInterval rho_in,
                                               int s) {
    if (s <= 0 || s > 5)
        throw std::invalid_argument("invalid residual update window");
    const double denom = static_cast<double>(std::uint64_t{1} << s);
    const double w = static_cast<double>(std::uint64_t{1} << (s - 1));
    ResidualInterval out{
        .lo = (rho_in.lo + (w + 1.0) / 2.0) / denom,
        .hi = ((w - 1.0) + rho_in.hi + (w + 1.0) / 2.0) / denom,
    };
    if (out.lo < 0.0 || !(out.hi < 1.0))
        throw std::logic_error("GapScale residual interval escaped [0,1)");
    return out;
}

template <typename T>
inline GapRoundParams<T> MakeGapScaleRoundParams(
    int round, int p_in, int s, ResidualInterval rho_in) {
    static_assert(std::is_unsigned_v<T>, "MakeGapScaleRoundParams expects unsigned torus");
    if (p_in <= 1 || s <= 0 || s > 5 || p_in <= s)
        throw std::invalid_argument("invalid GapScale round precision");
    const int p_out = p_in - s;
    const auto rho_out = UpdateGapScaleResidual(rho_in, s);
    const auto w = std::uint64_t{1} << (s - 1);
    return GapRoundParams<T>{
        .round = round,
        .p_in = p_in,
        .p_out = p_out,
        .s = s,
        .k = p_in - s,
        .w = w,
        .delta_in = TorusScaleForPrecision<T>(p_in),
        .delta_out = TorusScaleForPrecision<T>(p_out),
        .A = GapScaleWindowWeight<T>(p_in, s),
        .offset = GapScaleWindowOffset<T>(p_in, s),
        .rho_in = rho_in,
        .rho_out = rho_out,
        .residual_margin = std::min(rho_out.lo, 1.0 - rho_out.hi),
    };
}

inline double FinalGapEffectiveMarginUnits(int p_work, int k,
                                           ResidualInterval rho) {
    const int lsb_k = LSBIndexFromChapterBit(p_work, k);
    const auto w = static_cast<double>(std::uint64_t{1} << lsb_k);
    const double offset_units = (w + 1.0 - rho.lo - rho.hi) / 2.0;
    return std::min(offset_units + rho.lo,
                    (w + 1.0) - offset_units - rho.hi);
}

inline double FinalGapOffsetUnits(int p_work, int k, ResidualInterval rho) {
    const int lsb_k = LSBIndexFromChapterBit(p_work, k);
    const auto w = static_cast<double>(std::uint64_t{1} << lsb_k);
    return (w + 1.0 - rho.lo - rho.hi) / 2.0;
}

template <typename T>
inline T FinalGapOffsetForResidual(int p_work, int k, ResidualInterval rho) {
    static_assert(std::is_unsigned_v<T>,
                  "FinalGapOffsetForResidual expects unsigned torus");
    const long double units =
        static_cast<long double>(FinalGapOffsetUnits(p_work, k, rho));
    if (!(units > 0.0L))
        throw std::logic_error("recursive residual leaves invalid final offset");
    const long double delta =
        static_cast<long double>(TorusScaleForPrecision<T>(p_work));
    return static_cast<T>(std::llround(units * delta));
}

inline SimResult SimulateRecursiveGapReduce(
    std::uint64_t m_original, int p_original, int p_target = 12,
    int k_final = 7) {
    if (p_original <= 0 || p_original > 63)
        throw std::invalid_argument("simulator precision must be in [1,63]");
    if (p_target <= 0 || p_target > 63)
        throw std::invalid_argument("simulator target precision must be in [1,63]");
    if (p_original < 64 && m_original >= (std::uint64_t{1} << p_original))
        throw std::invalid_argument("simulator message outside precision domain");

    SimResult result{};
    result.m_original = m_original;
    result.p_original = p_original;
    result.schedule = RecursiveGapScaleScheduleString(p_original, p_target);
    std::uint64_t m = m_original;
    double rho = 0.0;
    int p_cur = p_original;
    int round = 0;
    while (p_cur > p_target) {
        const int s = std::min(5, p_cur - p_target);
        const std::uint64_t mask = (std::uint64_t{1} << s) - 1;
        const std::uint64_t w = std::uint64_t{1} << (s - 1);
        const std::uint64_t q = m >> s;
        const std::uint64_t r = m & mask;
        const int bit = r >= w ? 1 : 0;
        const std::uint64_t r_clear = r - (bit ? w : 0);
        const double rho_next =
            (static_cast<double>(r_clear) + rho +
             (static_cast<double>(w) + 1.0) / 2.0) /
            static_cast<double>(std::uint64_t{1} << s);
        if (rho_next < 0.0 || !(rho_next < 1.0))
            throw std::logic_error("simulated residual escaped [0,1)");
        result.rounds.push_back(SimRoundState{
            .round = round,
            .p_in = p_cur,
            .p_out = p_cur - s,
            .s = s,
            .q = q,
            .r = r,
            .bit = bit,
            .rho_in = rho,
            .rho_out = rho_next,
        });
        m = q;
        rho = rho_next;
        p_cur -= s;
        round++;
    }
    if (p_cur < p_target) {
        m <<= (p_target - p_cur);
        p_cur = p_target;
    }
    result.p_out = p_cur;
    result.m_work = m;
    result.rho_work = rho;
    result.min_rho = rho;
    result.max_rho = rho;
    result.min_final_margin_units =
        FinalGapEffectiveMarginUnits(p_target, k_final,
                                     ResidualInterval{.lo = rho, .hi = rho});
    const bool original_msb =
        ((m_original >> (p_original - 1)) & std::uint64_t{1}) != 0;
    const bool reduced_msb =
        static_cast<long double>(m) + static_cast<long double>(rho) >=
        std::ldexp(1.0L, p_target - 1);
    result.msb_ok = original_msb == reduced_msb;
    return result;
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
inline void RecordLvl2ToLvl0GatePBSStats(
    BlindRotatePruneStats* stats,
    int target_period = 0) {
    if (!stats) return;
    using midP = typename iksP::targetP;
    using outP = typename brP_base::targetP;
    stats->pbs_calls++;
    stats->pbs_count_gapmsb++;
    stats->key_switch_count++;
    stats->total += midP::k * midP::n;
    stats->cmux_calls += midP::k * midP::n;
    stats->periods.push_back(target_period > 0 ? target_period : 2 * outP::n);
    stats->total_by_pbs.push_back(midP::k * midP::n);
    stats->cmux_by_pbs.push_back(midP::k * midP::n);
    stats->skipped_by_pbs.push_back(0);
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

template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
NaiveSignPBS_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    typename iksP::domainP::T offset,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match lvl0->lvl2 sign BK domain");
    static_assert(std::is_same_v<inP, outP>,
                  "lvl0->lvl2 sign PBS currently returns to the input lvl2 key");

    TFHEpp::TLWE<inP> shifted = ct;
    shifted[inP::k * inP::n] += offset;

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    TFHEpp::Polynomial<outP> tv{};
    tv.fill(static_cast<typename outP::T>(outP::μ));

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    return out;
}

// HE3DB-style lightweight bit extraction for the native GapMSB window.
//
// A ciphertext Enc_{Q/2^p}(m) is multiplied by 2^{p-(ell+1)} so that the
// requested LSB-side bit ell becomes the sign bit of the current window.  The
// following lvl2->lvl0 IKS + lvl0->lvl2 gate bootstrap outputs BoolHalf
// arithmetic encoding 0/Q/2.  This preserves the Chapter-3 clear-one-bit
// semantics while avoiding a full lvl2 MetaPBS bit-extraction round.
template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
BitExtractBoolHalf_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    int encoding_p,
    int chapter_bit,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match lightweight bit-extract BK domain");
    static_assert(std::is_same_v<inP, outP>,
                  "lightweight bit extract currently returns to the input lvl2 key");

    if (encoding_p <= 1)
        throw std::invalid_argument("BitExtractBoolHalf_Lvl02 requires p >= 2");
    if (chapter_bit <= 0 || chapter_bit >= encoding_p)
        throw std::invalid_argument("BitExtractBoolHalf_Lvl02 requires 1 <= k < p");

    const int lsb_k = LSBIndexFromChapterBit(encoding_p, chapter_bit);
    const int window_p = lsb_k + 1;
    const int shift_bits = encoding_p - window_p;
    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    if (shift_bits < 0 || shift_bits >= digits)
        throw std::invalid_argument("BitExtractBoolHalf_Lvl02 shift is out of range");

    TFHEpp::TLWE<inP> shifted{};
    for (std::size_t i = 0; i < ct.size(); i++)
        shifted[i] = static_cast<typename inP::T>(ct[i] << shift_bits);
    shifted[inP::k * inP::n] +=
        StandardMSBOffsetForPrecision<typename inP::T>(window_p);

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    const auto half_bool = static_cast<typename outP::T>(outP::μ << 1);
    TFHEpp::Polynomial<outP> tv{};
    tv.fill(typename outP::T(0) - half_bool);

    if (stats) {
        RecordLvl2ToLvl0GatePBSStats<iksP, brP_logari>(stats, 2 * outP::n);
        stats->pbs_count_bit_extract++;
    }

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += half_bool;
    return out;
}

template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
ExtractWeightedChapterBitCentered_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    int encoding_p,
    int chapter_bit,
    typename brP_logari::targetP::T weight_torus,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match centered weighted-bit BK domain");
    static_assert(std::is_same_v<inP, outP>,
                  "centered weighted-bit extraction currently returns to lvl2");

    if ((weight_torus & typename outP::T(1)) != 0)
        throw std::invalid_argument(
            "centered weighted-bit extraction requires an even torus weight");
    const auto half_w = static_cast<typename outP::T>(weight_torus >> 1);
    if (half_w == 0)
        throw std::invalid_argument("centered weighted-bit extraction requires nonzero weight");

    if (encoding_p <= 1)
        throw std::invalid_argument("ExtractWeightedChapterBitCentered_Lvl02 requires p >= 2");
    if (chapter_bit <= 0 || chapter_bit >= encoding_p)
        throw std::invalid_argument("ExtractWeightedChapterBitCentered_Lvl02 requires 1 <= k < p");

    const int lsb_k = LSBIndexFromChapterBit(encoding_p, chapter_bit);
    const int window_p = lsb_k + 1;
    const int shift_bits = encoding_p - window_p;
    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    if (shift_bits < 0 || shift_bits >= digits)
        throw std::invalid_argument("ExtractWeightedChapterBitCentered_Lvl02 shift is out of range");

    TFHEpp::TLWE<inP> shifted{};
    for (std::size_t i = 0; i < ct.size(); i++)
        shifted[i] = static_cast<typename inP::T>(ct[i] << shift_bits);
    shifted[inP::k * inP::n] +=
        StandardMSBOffsetForPrecision<typename inP::T>(window_p);

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, shifted, iksk);

    TFHEpp::Polynomial<outP> tv{};
    tv.fill(typename outP::T(0) - half_w);

    if (stats) {
        RecordLvl2ToLvl0GatePBSStats<iksP, brP_logari>(stats, 2 * outP::n);
        stats->pbs_count_bit_extract++;
    }

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += half_w;
    return out;
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
    BlindRotatePruneStats* stats = nullptr,
    bool recursive_round = false) {
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
        if (recursive_round) stats->pbs_count_recursive_bool_to_weight++;
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

// Extract b = 1 iff (m mod 2^s) + rho >= 2^(s-1) from a ciphertext whose
// current phase is Enc_{Q/2^p_cur}(m + rho).  This is the recursive GapScale
// low-window predicate; it never builds a p_cur-bit LUT.
template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
BitExtractLowWindowBoolPruned(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    int p_cur,
    int s,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match low-window bit BK domain");
    static_assert(std::is_same_v<inP, outP>,
                  "low-window bit extraction must return to the input level");

    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    if (s <= 0 || s > 5 || p_cur <= s || p_cur >= digits)
        throw std::invalid_argument("BitExtractLowWindowBoolPruned invalid precision/window");
    const int shift_bits = p_cur - s;
    if (shift_bits < 0 || shift_bits >= digits)
        throw std::invalid_argument("BitExtractLowWindowBoolPruned shift out of range");

    TFHEpp::TLWE<inP> window{};
    for (std::size_t i = 0; i < ct.size(); i++)
        window[i] = static_cast<typename inP::T>(ct[i] << shift_bits);
    window[inP::k * inP::n] +=
        StandardMSBOffsetForPrecision<typename inP::T>(s);

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, window, iksk);

    const auto half_bool =
        static_cast<typename outP::T>(static_cast<typename outP::T>(outP::μ) << 1);
    TFHEpp::Polynomial<outP> tv{};
    tv.fill(typename outP::T(0) - half_bool);

    if (stats) {
        RecordLvl2ToLvl0GatePBSStats<iksP, brP_logari>(stats, 2 * outP::n);
        stats->pbs_count_bit_extract++;
        stats->pbs_count_recursive_bit_extract++;
    }

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += half_bool;
    return out;
}

template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
ExtractWeightedLowWindowPBS_Lvl02(
    const TFHEpp::TLWE<typename iksP::domainP>& ct,
    int p_cur,
    int s,
    typename brP_logari::targetP::T weight_torus,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using midP = typename iksP::targetP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<midP, typename brP_logari::domainP>,
                  "IKS target must match weighted low-window BK domain");
    static_assert(std::is_same_v<inP, outP>,
                  "weighted low-window extraction must return to the input level");

    constexpr int digits = std::numeric_limits<typename inP::T>::digits;
    if (s <= 0 || s > 6 || p_cur <= s || p_cur >= digits)
        throw std::invalid_argument("ExtractWeightedLowWindowPBS_Lvl02 invalid precision/window");
    if ((weight_torus & typename outP::T(1)) != 0)
        throw std::invalid_argument("weighted low-window extraction requires even weight");
    const auto half_w = static_cast<typename outP::T>(weight_torus >> 1);
    if (half_w == 0)
        throw std::invalid_argument("weighted low-window extraction requires nonzero weight");

    const int shift_bits = p_cur - s;
    TFHEpp::TLWE<inP> window{};
    for (std::size_t i = 0; i < ct.size(); i++)
        window[i] = static_cast<typename inP::T>(ct[i] << shift_bits);
    window[inP::k * inP::n] +=
        StandardMSBOffsetForPrecision<typename inP::T>(s);

    TFHEpp::TLWE<midP> tlwe_mid{};
    TFHEpp::IdentityKeySwitch<iksP>(tlwe_mid, window, iksk);

    TFHEpp::Polynomial<outP> tv{};
    tv.fill(typename outP::T(0) - half_w);

    if (stats) {
        RecordLvl2ToLvl0GatePBSStats<iksP, brP_logari>(stats, 2 * outP::n);
        stats->pbs_count_bit_extract++;
    }

    TFHEpp::TLWE<outP> out{};
    TFHEpp::GateBootstrappingTLWE2TLWE<brP_logari>(
        out, tlwe_mid, bkfft_logari, tv);
    out[outP::k * outP::n] += half_w;
    return out;
}

template <class iksP, class brP_logari>
TFHEpp::TLWE<typename brP_logari::targetP>
GapScaleRound(
    const TFHEpp::TLWE<typename iksP::domainP>& ct_cur,
    const GapRoundParams<typename brP_logari::targetP::T>& params,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<inP, outP>,
                  "GapScaleRound currently keeps the ciphertext at lvl2");

    const auto bit_ct = BitExtractLowWindowBoolPruned<iksP, brP_logari>(
        ct_cur, params.p_in, params.s, iksk, bkfft_logari, stats);
    const auto weight_ct = BoolToWeightPBS_Lvl02<iksP, brP_logari>(
        bit_ct, params.A, iksk, bkfft_logari, stats, true);

    TFHEpp::TLWE<outP> ct_next{};
    for (std::size_t i = 0; i < ct_next.size(); i++)
        ct_next[i] = ct_cur[i] - weight_ct[i];
    ct_next[outP::k * outP::n] += params.offset;
    return ct_next;
}

template <class iksP, class brP_logari>
RecursiveGapReduceResult<typename brP_logari::targetP>
RecursiveGapReduceToNative(
    const TFHEpp::TLWE<typename iksP::domainP>& ct_in,
    int p_original,
    int p_native,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    BlindRotatePruneStats* stats = nullptr) {
    using inP = typename iksP::domainP;
    using outP = typename brP_logari::targetP;
    static_assert(std::is_same_v<inP, outP>,
                  "RecursiveGapReduceToNative must keep the lvl2 ciphertext key");
    if (p_original <= 0 || p_native <= 0)
        throw std::invalid_argument("RecursiveGapReduceToNative requires positive precisions");
    if (p_original > 33)
        throw std::invalid_argument("RecursiveGapReduceToNative supports p_original <= 33");

    RecursiveGapReduceResult<outP> result{};
    result.ct_out = ct_in;
    result.schedule = RecursiveGapScaleScheduleString(p_original, p_native);
    if (p_original < p_native) {
        result.p_out = p_native;
        result.delta_out = TorusScaleForPrecision<typename outP::T>(p_native);
        result.rho_out = ResidualInterval{.lo = 0.0, .hi = 0.0};
        result.path_kind = "ZeroExtendToNative";
        result.semantic = "exact zero-extension, m_work=m_original<<(p_native-p_original)";
        return result;
    }
    if (p_original == p_native) {
        result.p_out = p_native;
        result.delta_out = TorusScaleForPrecision<typename outP::T>(p_native);
        result.rho_out = ResidualInterval{.lo = 0.0, .hi = 0.0};
        result.path_kind = "NativeExact";
        result.semantic = "exact native precision";
        return result;
    }

    result.path_kind = "RecursiveGapScale";
    ResidualInterval rho{.lo = 0.0, .hi = 0.0};
    int p_cur = p_original;
    int round = 0;
    TFHEpp::TLWE<outP> current = ct_in;
    while (p_cur > p_native) {
        const int s = std::min(5, p_cur - p_native);
        auto params = MakeGapScaleRoundParams<typename outP::T>(
            round, p_cur, s, rho);
        current = GapScaleRound<iksP, brP_logari>(
            current, params, iksk, bkfft_logari, stats);
        result.rounds.push_back(params);
        result.round_outputs.push_back(RecursiveGapRoundOutput<outP>{
            .round = round,
            .p_in = params.p_in,
            .p_out = params.p_out,
            .s = params.s,
            .rho_out = params.rho_out,
            .ct_out = current,
        });
        rho = params.rho_out;
        p_cur = params.p_out;
        round++;
    }
    result.ct_out = current;
    result.p_out = p_cur;
    result.delta_out = TorusScaleForPrecision<typename outP::T>(p_cur);
    result.rho_out = rho;
    result.semantic = "Enc_{Delta_next}(q + rho_next), residual tracked";
    return result;
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
    const int s = LSBIndexFromChapterBit(p_metapbs, chapter_bit) + 1;
    const auto weight =
        ArithmeticWeightForChapterBit<typename tgtP::T>(encoding_p, chapter_bit);
    if (options.weighted_mode == WeightedBitMode::DirectOnly)
        throw std::invalid_argument(
            "DirectOnly weighted-bit extraction is disabled on the lightweight "
            "GapMSB path; use BoolHalf extraction plus BoolToWeightPBS_Lvl02");

    (void)bkfft;
    (void)trkeys;
    (void)cfg;
    return ExtractWeightedLowWindowPBS_Lvl02<iksP, brP_logari>(
        ct, p_metapbs, s, weight, iksk, bkfft_logari, prune_stats);
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

    if (plain_bits <= options.kappa) {
        if (prune_stats) {
            RecordLvl2ToLvl0GatePBSStats<iksP, brP_base>(prune_stats);
            prune_stats->pbs_count_final_msb++;
        }
        const auto offset = StandardMSBOffsetForPrecision<typename inP::T>(plain_bits);
        return NaiveSignPBS_Lvl01<iksP, brP_base>(
            ct, offset, iksk, bkfft_base);
    }

    const int p_native = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    TFHEpp::TLWE<inP> ct_work = ct;
    ResidualInterval rho_work{.lo = 0.0, .hi = 0.0};
    if (plain_bits > p_native) {
        auto reduced = RecursiveGapReduceToNative<iksP, brP_logari>(
            ct_work, plain_bits, p_native, iksk, bkfft_logari, prune_stats);
        ct_work = reduced.ct_out;
        rho_work = reduced.rho_out;
    }

    // For plain_bits < p_native the physical torus phase is unchanged:
    // m * Q/2^p == (m << (p_native - p)) * Q/2^p_native.  For high precision
    // inputs, RecursiveGapReduceToNative has already changed the ciphertext
    // phase by ct_next = ct_cur - bit_weight + offset in every round and left a
    // tracked residual rho_work in [0,1).  The final Chapter-3 GapMSB is always
    // defined at the native/current precision.
    const int p_work = p_native;
    const int chapter_bit = p_work - options.kappa;
    if (chapter_bit <= 0 || chapter_bit >= p_work)
        throw std::invalid_argument("HomMSBOptions.kappa selects an invalid native GapMSB bit");
    if (FinalGapEffectiveMarginUnits(p_work, chapter_bit, rho_work) <= 0.0)
        throw std::logic_error("recursive residual leaves no final GapMSB margin");
    TFHEpp::TLWE<typename brP_base::targetP> out{};
    if (options.use_lightweight_gap_pbs) {
        auto weighted_bit = ExtractWeightedChapterBitFast<
            brP_metapbs, iksP, brP_logari>(
                ct_work, p_work, chapter_bit, bkfft, trkeys, cfg, iksk,
                bkfft_logari, options, prune_stats);

        TFHEpp::TLWE<inP> ct_gap{};
        ClearChapterBitAssign<brP_metapbs>(ct_gap, ct_work, weighted_bit);

        if (prune_stats) {
            RecordLvl2ToLvl0GatePBSStats<iksP, brP_base>(prune_stats);
            prune_stats->pbs_count_final_msb++;
        }
        out = NaiveSignPBS_Lvl01<iksP, brP_base>(
            ct_gap,
            FinalGapOffsetForResidual<typename inP::T>(
                p_work, chapter_bit, rho_work),
            iksk, bkfft_base);
    } else {
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

        // GapMSB already performed the aligned sign decision and returns a lvl2
        // sign ciphertext. This base PBS is only a key/level-compatible refresh
        // to the public HomMSB lvl1 output type.
        out = NaiveSignPBS_Lvl01<iksP, brP_base>(
            lvl2_sign, typename inP::T(0), iksk, bkfft_base);
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
