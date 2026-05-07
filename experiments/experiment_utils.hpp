#pragma once

#include "algorithms.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace OursExperiments {

struct br_lvl22param {
    using domainP = TFHEpp::lvl2param;
    using targetP = TFHEpp::lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

using brP_meta = br_lvl22param;
using brP_logari = TFHEpp::lvl02param;
using brP_base = TFHEpp::lvl01param;
using iksP_t = TFHEpp::lvl20param;
using PIn = TFHEpp::lvl2param;
using POut = TFHEpp::lvl1param;

struct Args {
    int argc;
    char** argv;

    bool Has(const std::string& key) const {
        for (int i = 1; i < argc; i++)
            if (argv[i] == key) return true;
        return false;
    }

    std::string Get(const std::string& key, const std::string& fallback) const {
        for (int i = 1; i + 1 < argc; i++)
            if (argv[i] == key) return argv[i + 1];
        return fallback;
    }

    int GetInt(const std::string& key, int fallback) const {
        return std::stoi(Get(key, std::to_string(fallback)));
    }

    std::uint64_t GetU64(const std::string& key, std::uint64_t fallback) const {
        return static_cast<std::uint64_t>(
            std::stoull(Get(key, std::to_string(fallback))));
    }
};

inline std::vector<std::string> Split(const std::string& text, char delim = ',') {
    std::vector<std::string> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

inline int DefaultK(int p) {
    return p >= 6 ? p - 5 : std::max(1, p - 2);
}

inline int ParseK(const std::string& text, int p) {
    if (text == "default") return DefaultK(p);
    if (text.rfind("p-", 0) == 0)
        return p - std::stoi(text.substr(2));
    if (text.rfind("p_minus_", 0) == 0)
        return p - std::stoi(text.substr(8));
    return std::stoi(text);
}

inline std::vector<std::pair<int, std::string>> ParseKList(
    const std::string& text, int p) {
    std::vector<std::pair<int, std::string>> values;
    for (const auto& item : Split(text)) values.push_back({ParseK(item, p), item});
    return values;
}

inline int NativePrecision(const MetaPBS2::Algorithm1Config& cfg) {
    return MetaPBS2::MessagePrecisionFromPowerOfTwoModulus(cfg.t);
}

inline int WorkPrecisionForComparison(int L,
                                      const MetaPBS2::Algorithm1Config& cfg) {
    const int p_original = L + 1;
    const int p_native = NativePrecision(cfg);
    (void)p_original;
    return p_native;
}

inline std::string PrecisionPathKind(int L,
                                     const MetaPBS2::Algorithm1Config& cfg) {
    const int p_original = L + 1;
    const int p_native = NativePrecision(cfg);
    if (p_original < p_native) return "ZeroExtendToNative";
    if (p_original == p_native) return "NativeExact";
    const auto schedule =
        MetaPBS2::MakeHE3DBStylePrecisionSchedule(p_original, p_native);
    return schedule.p_final < p_native
               ? "HE3DBStylePrecisionReducePlusZeroExtend"
               : "HE3DBStylePrecisionReduce";
}

inline std::string PrecisionScheduleString(
    int L, const MetaPBS2::Algorithm1Config& cfg) {
    const int p_original = L + 1;
    const int p_native = NativePrecision(cfg);
    if (p_original <= p_native)
        return p_original < p_native
                   ? std::to_string(p_original) + "->zeroextend" + std::to_string(p_native)
                   : std::to_string(p_original);
    const auto schedule =
        MetaPBS2::MakeHE3DBStylePrecisionSchedule(p_original, p_native);
    return MetaPBS2::ScheduleString(schedule, p_native);
}

inline std::vector<std::string> ParseOps(const std::string& text) {
    if (text == "all") return {"lt", "le", "gt", "ge", "eq", "ne"};
    return Split(text);
}

inline bool ExpectedCompare(std::uint64_t a, std::uint64_t b,
                            const std::string& op) {
    if (op == "lt") return a < b;
    if (op == "le") return a <= b;
    if (op == "gt") return a > b;
    if (op == "ge") return a >= b;
    if (op == "eq") return a == b;
    if (op == "ne") return a != b;
    throw std::invalid_argument("unsupported op: " + op);
}

template <class P>
typename P::T EncodeMessage(std::uint64_t message, int precision_bits) {
    return static_cast<typename P::T>(message)
           << (std::numeric_limits<typename P::T>::digits - precision_bits);
}

template <class P>
bool DecodeBool(const TFHEpp::TLWE<P>& ct, const TFHEpp::SecretKey& sk) {
    const auto phase = TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
    return static_cast<std::make_signed_t<typename P::T>>(phase) < 0;
}

inline double MsSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

struct LatencyStats {
    double avg = 0;
    double median = 0;
    double min = 0;
    double max = 0;
    double stddev = 0;
};

inline LatencyStats ComputeLatencyStats(std::vector<double> values) {
    LatencyStats s{};
    if (values.empty()) return s;
    std::sort(values.begin(), values.end());
    s.min = values.front();
    s.max = values.back();
    s.median = values.size() % 2 == 0
                   ? (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2.0
                   : values[values.size() / 2];
    s.avg = std::accumulate(values.begin(), values.end(), 0.0) /
            static_cast<double>(values.size());
    double var = 0;
    for (double v : values) var += (v - s.avg) * (v - s.avg);
    s.stddev = std::sqrt(var / static_cast<double>(values.size()));
    return s;
}

struct GapMSBCounters {
    std::uint64_t pbs_count = 0;
    std::uint64_t key_switch_count = 0;
    std::uint64_t pbs_count_reducer = 0;
    std::uint64_t pbs_count_gapmsb = 0;
    std::uint64_t pbs_count_bit_extract = 0;
    std::uint64_t pbs_count_bool_to_weight = 0;
    std::uint64_t pbs_count_final_msb = 0;
    std::uint64_t total_cmux = 0;
    std::uint64_t skipped_cmux = 0;
    std::uint64_t first_round_total_cmux = 0;
    std::uint64_t first_round_skipped_cmux = 0;
    int first_round_slot_period = 0;
};

inline GapMSBCounters ToCounters(const MetaPBS2::BlindRotatePruneStats& stats) {
    GapMSBCounters c{};
    c.pbs_count = stats.pbs_calls;
    c.key_switch_count = stats.key_switch_count == 0 ? stats.pbs_calls : stats.key_switch_count;
    c.pbs_count_reducer = stats.pbs_count_reducer;
    c.pbs_count_gapmsb = stats.pbs_count_gapmsb;
    c.pbs_count_bit_extract = stats.pbs_count_bit_extract;
    c.pbs_count_bool_to_weight = stats.pbs_count_bool_to_weight;
    c.pbs_count_final_msb = stats.pbs_count_final_msb;
    c.total_cmux = stats.total;
    c.skipped_cmux = stats.skipped;
    if (!stats.periods.empty()) c.first_round_slot_period = stats.periods.front();
    if (!stats.total_by_pbs.empty()) c.first_round_total_cmux = stats.total_by_pbs.front();
    if (!stats.skipped_by_pbs.empty()) c.first_round_skipped_cmux = stats.skipped_by_pbs.front();
    return c;
}

inline double Percent(std::uint64_t num, std::uint64_t den) {
    return den == 0 ? 0.0 : 100.0 * static_cast<double>(num) / static_cast<double>(den);
}

inline int ExactBitExtractPeriod(int p, int k,
                                 const MetaPBS2::Algorithm1Config& cfg) {
    const int lsb_k = MetaPBS2::LSBIndexFromChapterBit(p, k);
    auto tv = MetaPBS2::BuildKthBitTestVector<PIn>(lsb_k, cfg.t);
    return MetaPBS2::ExactNegacyclicPeriod<PIn>(tv);
}

struct OursRuntime {
    MetaPBS2::Algorithm1Config cfg;
    TFHEpp::SecretKey sk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_meta>> bk_meta;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_logari>> bk_logari;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_base>> bk_base;
    std::unique_ptr<TFHEpp::KeySwitchingKey<iksP_t>> iksk;
    std::vector<MetaPBS2::TruncRepeatKey<PIn>> trkeys;

    OursRuntime() : cfg(PaperReview::Chapter3MetaPBSConfig()) {
        bk_meta = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_meta>>();
        bk_logari = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_logari>>();
        bk_base = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_base>>();
        iksk = std::make_unique<TFHEpp::KeySwitchingKey<iksP_t>>();
        TFHEpp::bkfftgen<brP_meta>(*bk_meta, sk);
        TFHEpp::bkfftgen<brP_logari>(*bk_logari, sk);
        TFHEpp::bkfftgen<brP_base>(*bk_base, sk);
        TFHEpp::ikskgen<iksP_t>(*iksk, sk);
        for (const auto& round : cfg.rounds)
            trkeys.push_back(MetaPBS2::GenerateTruncRepeatKey<PIn>(
                sk.key.get<PIn>(), round.beta));
    }
};

struct CompareResult {
    bool got = false;
    double ms = 0;
    MetaPBS2::BlindRotatePruneStats stats;
};

inline PaperReview::ComparePredicate RelationPredicate(const std::string& op) {
    if (op == "lt" || op == "ge") return PaperReview::ComparePredicate::LessThan;
    if (op == "gt" || op == "le") return PaperReview::ComparePredicate::GreaterThan;
    throw std::invalid_argument("not a relation op: " + op);
}

inline bool RunOneRelation(
    OursRuntime& rt,
    const TFHEpp::TLWE<PIn>& a_ct,
    const TFHEpp::TLWE<PIn>& b_ct,
    int L,
    int k,
    const std::string& op,
    double& ms,
    MetaPBS2::BlindRotatePruneStats& stats) {
    MetaPBS2::HomMSBOptions options;
    const int p_work = WorkPrecisionForComparison(L, rt.cfg);
    options.kappa = std::max(1, p_work - k);
    options.enable_periodic_pruning = true;
    const auto start = std::chrono::steady_clock::now();
    const auto out = PaperReview::Chapter3HomCompare<
        brP_meta, brP_logari, brP_base, iksP_t>(
            a_ct, b_ct, L, RelationPredicate(op), *rt.bk_meta, rt.trkeys,
            rt.cfg, *rt.bk_logari, *rt.bk_base, *rt.iksk, options, &stats);
    ms += MsSince(start);
    return DecodeBool<POut>(out, rt.sk);
}

inline CompareResult RunOursCompare(
    OursRuntime& rt,
    std::uint64_t a,
    std::uint64_t b,
    int L,
    int k,
    const std::string& op) {
    const int p = L + 1;
    TFHEpp::TLWE<PIn> a_ct{}, b_ct{};
    TFHEpp::tlweSymEncrypt<PIn>(
        a_ct, EncodeMessage<PIn>(a, p), PIn::α, rt.sk.key.get<PIn>());
    TFHEpp::tlweSymEncrypt<PIn>(
        b_ct, EncodeMessage<PIn>(b, p), PIn::α, rt.sk.key.get<PIn>());

    CompareResult result{};
    if (op == "lt" || op == "gt") {
        result.got = RunOneRelation(rt, a_ct, b_ct, L, k, op, result.ms, result.stats);
    } else if (op == "le") {
        result.got = !RunOneRelation(rt, a_ct, b_ct, L, k, "gt", result.ms, result.stats);
    } else if (op == "ge") {
        result.got = !RunOneRelation(rt, a_ct, b_ct, L, k, "lt", result.ms, result.stats);
    } else if (op == "eq" || op == "ne") {
        const bool lt = RunOneRelation(rt, a_ct, b_ct, L, k, "lt", result.ms, result.stats);
        const bool gt = RunOneRelation(rt, a_ct, b_ct, L, k, "gt", result.ms, result.stats);
        const bool neq = lt || gt;
        result.got = (op == "ne") ? neq : !neq;
    } else {
        throw std::invalid_argument("unsupported op: " + op);
    }
    return result;
}

inline void PrintExperimentHeader(const std::string& name, int L, int p, int k,
                                  int trials, std::uint64_t seed) {
    std::cout << "# experiment: " << name << "\n"
              << "# scheme: ours\n"
              << "# N: " << PIn::n << "\n"
              << "# n: " << (PIn::k * PIn::n) << "\n"
              << "# L: " << L << "\n"
              << "# p: " << p << "\n"
              << "# k: " << k << "\n"
              << "# trials: " << trials << "\n"
              << "# seed: " << seed << "\n";
}

inline std::uint64_t MaxValueForL(int L) {
    if (L <= 0 || L >= 63) throw std::invalid_argument("L must be in [1,62]");
    return (std::uint64_t{1} << L) - 1;
}

}  // namespace OursExperiments
