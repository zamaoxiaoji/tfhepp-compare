#include "native_baseline_experiments.hpp"

#include "baselines/arcedb_native/native_compare.hpp"
#include "baselines/he3db_native/native_compare.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

namespace PaperReview {
namespace {

using Clock = std::chrono::steady_clock;

template <typename Fn>
double TimeMs(Fn&& fn) {
    const auto start = Clock::now();
    fn();
    const auto stop = Clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

std::vector<NativeBaselineSystem> SystemsToRun(NativeBaselineSystem system) {
    if (system == NativeBaselineSystem::All)
        return {NativeBaselineSystem::HE3DB, NativeBaselineSystem::ArcEDB};
    if (system == NativeBaselineSystem::Ours)
        throw std::invalid_argument("native baseline runner only handles he3db/arcedb/all");
    return {system};
}

std::uint32_t BitsForDomain(std::size_t domain) {
    if (domain == 0) throw std::invalid_argument("domain must be nonzero");
    std::uint32_t bits = 1;
    std::size_t capacity = 2;
    while (capacity < domain) {
        capacity <<= 1;
        bits++;
    }
    return bits;
}

void AddHE3DBEvalKeys(HE3DBNative::TFHEEvalKey& ek,
                      HE3DBNative::TFHESecretKey& sk,
                      bool need_lvl2) {
    ek.emplacebkfft<HE3DBNative::Lvl01>(sk);
    ek.emplaceiksk<HE3DBNative::Lvl10>(sk);
    if (need_lvl2)
        throw std::invalid_argument("this vendored HE3DB baseline uses the lvl1 native comparison subset");
}

void AddArcEDBEvalKeys(ArcEDBNative::TFHEEvalKey& ek,
                       ArcEDBNative::TFHESecretKey& sk) {
    ek.emplacebkfft<ArcEDBNative::Lvl01>(sk);
    ek.emplaceiksk<ArcEDBNative::Lvl10>(sk);
}

std::vector<int> HE3DBEqualMasks(
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    HE3DBNative::TFHESecretKey& sk,
    HE3DBNative::TFHEEvalKey& ek) {
    std::vector<int> masks(values.size(), 0);
    auto target_ct = HE3DBNative::EncryptInt<HE3DBNative::Lvl1>(target, bits, sk);
    for (std::size_t i = 0; i < values.size(); i++) {
        auto lhs = HE3DBNative::EncryptInt<HE3DBNative::Lvl1>(values[i], bits, sk);
        HE3DBNative::TLWELvl1 res;
        HE3DBNative::equal<HE3DBNative::Lvl1>(
            lhs, target_ct, res, bits, ek, HE3DBNative::LOGIC);
        masks[i] = HE3DBNative::DecryptLogic(res, sk);
    }
    return masks;
}

std::vector<int> ArcEDBEqualMasks(
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    ArcEDBNative::TFHESecretKey& sk,
    ArcEDBNative::TFHEEvalKey& ek) {
    std::vector<int> masks(values.size(), 0);
    std::vector<ArcEDBNative::TRGSWLvl1> target_ct;
    ArcEDBNative::exponent_encrypt_rgsw<ArcEDBNative::Lvl1>(
        target, bits, target_ct, sk, true);
    for (std::size_t i = 0; i < values.size(); i++) {
        std::vector<ArcEDBNative::TRLWELvl1> lhs;
        ArcEDBNative::exponent_encrypt<ArcEDBNative::Lvl1>(
            values[i], bits, lhs, sk);
        ArcEDBNative::TLWELvl1 res;
        ArcEDBNative::equality_tfhepp(lhs, target_ct, target_ct.size(), res, ek, sk);
        masks[i] = ArcEDBNative::DecryptLogic(res, sk);
    }
    return masks;
}

std::vector<int> HE3DBLessMasks(
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    HE3DBNative::TFHESecretKey& sk,
    HE3DBNative::TFHEEvalKey& ek) {
    std::vector<int> masks(values.size(), 0);
    auto target_ct = HE3DBNative::EncryptInt<HE3DBNative::Lvl1>(target, bits, sk);
    for (std::size_t i = 0; i < values.size(); i++) {
        auto lhs = HE3DBNative::EncryptInt<HE3DBNative::Lvl1>(values[i], bits, sk);
        HE3DBNative::TLWELvl1 res;
        HE3DBNative::less_than<HE3DBNative::Lvl1>(
            lhs, target_ct, res, bits, ek, HE3DBNative::LOGIC);
        masks[i] = HE3DBNative::DecryptLogic(res, sk);
    }
    return masks;
}

std::vector<int> HE3DBGreaterMasks(
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    HE3DBNative::TFHESecretKey& sk,
    HE3DBNative::TFHEEvalKey& ek) {
    std::vector<int> masks(values.size(), 0);
    auto target_ct = HE3DBNative::EncryptInt<HE3DBNative::Lvl1>(target, bits, sk);
    for (std::size_t i = 0; i < values.size(); i++) {
        auto lhs = HE3DBNative::EncryptInt<HE3DBNative::Lvl1>(values[i], bits, sk);
        HE3DBNative::TLWELvl1 res;
        HE3DBNative::less_than<HE3DBNative::Lvl1>(
            target_ct, lhs, res, bits, ek, HE3DBNative::LOGIC);
        masks[i] = HE3DBNative::DecryptLogic(res, sk);
    }
    return masks;
}

std::vector<int> ArcEDBLessMasks(
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    ArcEDBNative::TFHESecretKey& sk,
    ArcEDBNative::TFHEEvalKey& ek) {
    std::vector<int> masks(values.size(), 0);
    std::vector<ArcEDBNative::TRGSWLvl1> target_ct;
    ArcEDBNative::exponent_encrypt_rgsw<ArcEDBNative::Lvl1>(
        target, bits, target_ct, sk, true);
    for (std::size_t i = 0; i < values.size(); i++) {
        std::vector<ArcEDBNative::TRLWELvl1> lhs;
        ArcEDBNative::exponent_encrypt<ArcEDBNative::Lvl1>(
            values[i], bits, lhs, sk);
        ArcEDBNative::TLWELvl1 res;
        ArcEDBNative::less_than_tfhepp(lhs, target_ct, target_ct.size(), res, ek, sk);
        masks[i] = ArcEDBNative::DecryptLogic(res, sk);
    }
    return masks;
}

std::vector<int> ArcEDBGreaterMasks(
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    ArcEDBNative::TFHESecretKey& sk,
    ArcEDBNative::TFHEEvalKey& ek) {
    std::vector<int> masks(values.size(), 0);
    std::vector<ArcEDBNative::TRGSWLvl1> target_ct;
    ArcEDBNative::exponent_encrypt_rgsw<ArcEDBNative::Lvl1>(
        target, bits, target_ct, sk, true);
    for (std::size_t i = 0; i < values.size(); i++) {
        std::vector<ArcEDBNative::TRLWELvl1> lhs;
        ArcEDBNative::exponent_encrypt<ArcEDBNative::Lvl1>(
            values[i], bits, lhs, sk);
        ArcEDBNative::TLWELvl1 res;
        ArcEDBNative::greater_than_tfhepp(lhs, target_ct, target_ct.size(), res, ek, sk);
        masks[i] = ArcEDBNative::DecryptLogic(res, sk);
    }
    return masks;
}

std::vector<int> NativeEqualMasks(
    NativeBaselineSystem system,
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    HE3DBNative::TFHESecretKey* he_sk,
    HE3DBNative::TFHEEvalKey* he_ek,
    ArcEDBNative::TFHESecretKey* arc_sk,
    ArcEDBNative::TFHEEvalKey* arc_ek) {
    if (system == NativeBaselineSystem::HE3DB)
        return HE3DBEqualMasks(values, target, bits, *he_sk, *he_ek);
    if (system == NativeBaselineSystem::ArcEDB)
        return ArcEDBEqualMasks(values, target, bits, *arc_sk, *arc_ek);
    throw std::invalid_argument("unsupported native baseline system");
}

std::vector<int> NativeLessMasks(
    NativeBaselineSystem system,
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    HE3DBNative::TFHESecretKey* he_sk,
    HE3DBNative::TFHEEvalKey* he_ek,
    ArcEDBNative::TFHESecretKey* arc_sk,
    ArcEDBNative::TFHEEvalKey* arc_ek) {
    if (system == NativeBaselineSystem::HE3DB)
        return HE3DBLessMasks(values, target, bits, *he_sk, *he_ek);
    if (system == NativeBaselineSystem::ArcEDB)
        return ArcEDBLessMasks(values, target, bits, *arc_sk, *arc_ek);
    throw std::invalid_argument("unsupported native baseline system");
}

std::vector<int> NativeGreaterMasks(
    NativeBaselineSystem system,
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    std::uint32_t bits,
    HE3DBNative::TFHESecretKey* he_sk,
    HE3DBNative::TFHEEvalKey* he_ek,
    ArcEDBNative::TFHESecretKey* arc_sk,
    ArcEDBNative::TFHEEvalKey* arc_ek) {
    if (system == NativeBaselineSystem::HE3DB)
        return HE3DBGreaterMasks(values, target, bits, *he_sk, *he_ek);
    if (system == NativeBaselineSystem::ArcEDB)
        return ArcEDBGreaterMasks(values, target, bits, *arc_sk, *arc_ek);
    throw std::invalid_argument("unsupported native baseline system");
}

NativeBaselineResult MakeBaseResult(
    const std::string& chapter,
    const std::string& op,
    const std::string& query,
    NativeBaselineSystem system,
    std::size_t rows,
    std::size_t groups,
    std::size_t key_domain) {
    NativeBaselineResult r;
    r.chapter = chapter;
    r.op = op;
    r.query = query;
    r.system = NativeBaselineSystemName(system);
    r.rows = rows;
    r.groups = groups;
    r.key_domain = key_domain;
    return r;
}

void Finalize(NativeBaselineResult& r) {
    r.total_query_time_ms = r.filter_time_ms + r.aggregation_time_ms;
    r.plain_result = std::accumulate(r.plain_groups.begin(), r.plain_groups.end(), 0.0);
    r.encrypted_result = std::accumulate(
        r.encrypted_groups.begin(), r.encrypted_groups.end(), 0.0);
    double sum_err = 0.0;
    std::size_t correct = 0;
    for (std::size_t i = 0; i < r.plain_groups.size(); i++) {
        const double err = std::abs(r.plain_groups[i] - r.encrypted_groups[i]);
        sum_err += err;
        r.max_abs_error = std::max(r.max_abs_error, err);
        if (err <= 1e-6) correct++;
    }
    if (!r.plain_groups.empty()) {
        r.avg_abs_error = sum_err / r.plain_groups.size();
        r.accuracy = static_cast<double>(correct) / r.plain_groups.size();
    }
}

struct GroupData {
    std::vector<std::uint32_t> keys;
    std::vector<double> values;
    std::vector<double> plain;
};

GroupData MakeGroupData(std::size_t rows, std::size_t groups) {
    GroupData data;
    data.keys.resize(rows);
    data.values.resize(rows);
    data.plain.assign(groups, 0.0);
    for (std::size_t i = 0; i < rows; i++) {
        data.keys[i] = static_cast<std::uint32_t>(i % groups);
        data.values[i] = static_cast<double>(i + 1);
        data.plain[data.keys[i]] += data.values[i];
    }
    return data;
}

struct JoinData {
    std::vector<std::uint32_t> left_key;
    std::vector<double> right_payload;
    std::vector<double> expected;
};

JoinData MakeJoinData(std::size_t rows, std::size_t key_domain, std::uint64_t seed) {
    JoinData data;
    data.left_key.resize(rows);
    data.right_payload.resize(key_domain);
    data.expected.resize(rows);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> key_dist(0, static_cast<int>(key_domain - 1));
    std::uniform_int_distribution<int> payload_dist(1, 1000);
    for (std::size_t k = 0; k < key_domain; k++)
        data.right_payload[k] = static_cast<double>(payload_dist(rng));
    for (std::size_t i = 0; i < rows; i++) {
        data.left_key[i] = static_cast<std::uint32_t>(key_dist(rng));
        data.expected[i] = data.right_payload[data.left_key[i]];
    }
    return data;
}

NativeBaselineResult RunGroupByOne(
    NativeBaselineSystem system,
    const NativeBaselineOptions& options) {
    const auto data = MakeGroupData(options.rows, options.groups);
    const std::uint32_t bits = BitsForDomain(options.groups);
    auto result = MakeBaseResult(
        "4", "group_by", "chap4_group_by", system,
        options.rows, options.groups, 0);
    result.enumerated_groups = options.groups;
    result.equality_calls = options.rows * options.groups;
    result.plain_groups = data.plain;
    result.encrypted_groups.assign(options.groups, 0.0);

    HE3DBNative::TFHESecretKey he_sk;
    HE3DBNative::TFHEEvalKey he_ek;
    ArcEDBNative::TFHESecretKey arc_sk;
    ArcEDBNative::TFHEEvalKey arc_ek;
    if (system == NativeBaselineSystem::HE3DB) AddHE3DBEvalKeys(he_ek, he_sk, false);
    if (system == NativeBaselineSystem::ArcEDB) AddArcEDBEvalKeys(arc_ek, arc_sk);

    result.aggregation_time_ms = TimeMs([&] {
        for (std::size_t g = 0; g < options.groups; g++) {
            const auto masks = NativeEqualMasks(
                system, data.keys, static_cast<std::uint32_t>(g), bits,
                &he_sk, &he_ek, &arc_sk, &arc_ek);
            for (std::size_t i = 0; i < options.rows; i++)
                result.encrypted_groups[g] += masks[i] * data.values[i];
        }
    });
    Finalize(result);
    return result;
}

NativeBaselineResult RunJoinOne(
    NativeBaselineSystem system,
    const NativeBaselineOptions& options) {
    const auto data = MakeJoinData(options.rows, options.key_domain, options.seed);
    const std::uint32_t bits = BitsForDomain(options.key_domain);
    auto result = MakeBaseResult(
        "4", "join", "chap4_join", system,
        options.rows, 0, options.key_domain);
    result.enumerated_keys = options.key_domain;
    result.equality_calls = options.rows * options.key_domain;
    result.plain_groups = data.expected;
    result.encrypted_groups.assign(options.rows, 0.0);

    HE3DBNative::TFHESecretKey he_sk;
    HE3DBNative::TFHEEvalKey he_ek;
    ArcEDBNative::TFHESecretKey arc_sk;
    ArcEDBNative::TFHEEvalKey arc_ek;
    if (system == NativeBaselineSystem::HE3DB) AddHE3DBEvalKeys(he_ek, he_sk, false);
    if (system == NativeBaselineSystem::ArcEDB) AddArcEDBEvalKeys(arc_ek, arc_sk);

    result.aggregation_time_ms = TimeMs([&] {
        for (std::size_t k = 0; k < options.key_domain; k++) {
            const auto masks = NativeEqualMasks(
                system, data.left_key, static_cast<std::uint32_t>(k), bits,
                &he_sk, &he_ek, &arc_sk, &arc_ek);
            for (std::size_t i = 0; i < options.rows; i++)
                result.encrypted_groups[i] += masks[i] * data.right_payload[k];
        }
    });
    Finalize(result);
    return result;
}

NativeBaselineResult RunQ3One(
    NativeBaselineSystem system,
    const NativeBaselineOptions& options) {
    const std::size_t key_domain = options.key_domain;
    const std::size_t priority_domain = options.priority_domain;
    const auto join_data = MakeJoinData(options.rows, key_domain, options.seed);
    auto result = MakeBaseResult(
        "5", "tpch_q3", "q3", system,
        options.rows, key_domain * priority_domain, key_domain);
    result.join_layers = 2;
    result.enumerated_keys = key_domain * result.join_layers;
    result.enumerated_groups = key_domain * priority_domain;
    result.equality_calls =
        options.rows * result.enumerated_keys +
        options.rows * result.enumerated_groups +
        key_domain * 2 + options.rows;
    result.plain_groups.assign(key_domain * priority_domain, 0.0);
    result.encrypted_groups.assign(result.plain_groups.size(), 0.0);

    std::vector<std::uint32_t> order_priority(key_domain);
    std::vector<std::uint32_t> customer_segment(key_domain);
    std::vector<std::uint32_t> order_date(key_domain);
    std::vector<std::uint32_t> line_shipdate(options.rows);
    std::vector<double> revenue(options.rows);
    for (std::size_t k = 0; k < key_domain; k++) {
        order_priority[k] = static_cast<std::uint32_t>(k % priority_domain);
        customer_segment[k] = static_cast<std::uint32_t>(k % 2);
        order_date[k] = static_cast<std::uint32_t>(10 + k);
    }
    for (std::size_t i = 0; i < options.rows; i++) {
        revenue[i] = 30.0 + static_cast<double>(i % 17);
        line_shipdate[i] = static_cast<std::uint32_t>(20 + i);
        const auto key = join_data.left_key[i];
        const auto pr = order_priority[key];
        const bool order_pass =
            customer_segment[key] == 0 && order_date[key] < 10 + key_domain &&
            line_shipdate[i] > 19;
        if (order_pass)
            result.plain_groups[key + key_domain * pr] += revenue[i];
    }

    HE3DBNative::TFHESecretKey he_sk;
    HE3DBNative::TFHEEvalKey he_ek;
    ArcEDBNative::TFHESecretKey arc_sk;
    ArcEDBNative::TFHEEvalKey arc_ek;
    if (system == NativeBaselineSystem::HE3DB) AddHE3DBEvalKeys(he_ek, he_sk, false);
    if (system == NativeBaselineSystem::ArcEDB) AddArcEDBEvalKeys(arc_ek, arc_sk);
    const std::uint32_t key_bits = BitsForDomain(key_domain);
    const std::uint32_t group_bits = BitsForDomain(key_domain * priority_domain);

    std::vector<int> customer_segment_mask;
    std::vector<int> order_date_mask;
    std::vector<int> line_ship_mask;
    result.filter_time_ms = TimeMs([&] {
        customer_segment_mask = NativeEqualMasks(
            system, customer_segment, 0, 1, &he_sk, &he_ek, &arc_sk, &arc_ek);
        order_date_mask = NativeLessMasks(
            system, order_date, static_cast<std::uint32_t>(10 + key_domain),
            BitsForDomain(10 + key_domain + 1), &he_sk, &he_ek, &arc_sk, &arc_ek);
        line_ship_mask = NativeGreaterMasks(
            system, line_shipdate, 19, BitsForDomain(20 + options.rows + 1),
            &he_sk, &he_ek, &arc_sk, &arc_ek);
    });

    result.aggregation_time_ms = TimeMs([&] {
        std::vector<double> line_pass(options.rows, 0.0);
        std::vector<std::uint32_t> line_group(options.rows, 0);
        for (std::size_t k = 0; k < key_domain; k++) {
            const auto masks = NativeEqualMasks(
                system, join_data.left_key, static_cast<std::uint32_t>(k),
                key_bits, &he_sk, &he_ek, &arc_sk, &arc_ek);
            for (std::size_t i = 0; i < options.rows; i++) {
                if (masks[i]) {
                    line_pass[i] +=
                        customer_segment_mask[k] * order_date_mask[k] *
                        line_ship_mask[i];
                    line_group[i] = static_cast<std::uint32_t>(
                        k + key_domain * order_priority[k]);
                }
            }
        }
        for (std::size_t g = 0; g < result.plain_groups.size(); g++) {
            const auto masks = NativeEqualMasks(
                system, line_group, static_cast<std::uint32_t>(g),
                group_bits, &he_sk, &he_ek, &arc_sk, &arc_ek);
            for (std::size_t i = 0; i < options.rows; i++)
                result.encrypted_groups[g] += masks[i] * line_pass[i] * revenue[i];
        }
    });
    Finalize(result);
    return result;
}

NativeBaselineResult RunQ5One(
    NativeBaselineSystem system,
    const NativeBaselineOptions& options) {
    const std::size_t key_domain = options.key_domain;
    const std::size_t nation_domain = options.nation_domain;
    const auto join_data = MakeJoinData(options.rows, key_domain, options.seed);
    auto result = MakeBaseResult(
        "5", "tpch_q5", "q5", system,
        options.rows, nation_domain, key_domain);
    result.join_layers = 5;
    result.enumerated_keys = key_domain * result.join_layers;
    result.enumerated_groups = nation_domain;
    result.equality_calls =
        options.rows * result.enumerated_keys +
        options.rows * result.enumerated_groups +
        key_domain * 2;
    result.plain_groups.assign(nation_domain, 0.0);
    result.encrypted_groups.assign(nation_domain, 0.0);

    std::vector<std::uint32_t> key_nation(key_domain);
    std::vector<std::uint32_t> nation_region(nation_domain);
    std::vector<std::uint32_t> order_date(key_domain);
    std::vector<double> revenue(options.rows);
    for (std::size_t n = 0; n < nation_domain; n++)
        nation_region[n] = static_cast<std::uint32_t>(n % 2);
    for (std::size_t k = 0; k < key_domain; k++) {
        key_nation[k] = static_cast<std::uint32_t>(k % nation_domain);
        order_date[k] = static_cast<std::uint32_t>(10 + k);
    }
    for (std::size_t i = 0; i < options.rows; i++) {
        revenue[i] = 40.0 + static_cast<double>(i % 19);
        const auto key = join_data.left_key[i];
        const bool pass =
            order_date[key] < 10 + key_domain && nation_region[key_nation[key]] == 0;
        if (pass)
            result.plain_groups[key_nation[key]] += revenue[i];
    }

    HE3DBNative::TFHESecretKey he_sk;
    HE3DBNative::TFHEEvalKey he_ek;
    ArcEDBNative::TFHESecretKey arc_sk;
    ArcEDBNative::TFHEEvalKey arc_ek;
    if (system == NativeBaselineSystem::HE3DB) AddHE3DBEvalKeys(he_ek, he_sk, false);
    if (system == NativeBaselineSystem::ArcEDB) AddArcEDBEvalKeys(arc_ek, arc_sk);
    const std::uint32_t key_bits = BitsForDomain(key_domain);
    const std::uint32_t nation_bits = BitsForDomain(nation_domain);

    std::vector<int> order_date_mask;
    std::vector<int> nation_region_mask;
    result.filter_time_ms = TimeMs([&] {
        order_date_mask = NativeLessMasks(
            system, order_date, static_cast<std::uint32_t>(10 + key_domain),
            BitsForDomain(10 + key_domain + 1), &he_sk, &he_ek, &arc_sk, &arc_ek);
        nation_region_mask = NativeEqualMasks(
            system, nation_region, 0, 1, &he_sk, &he_ek, &arc_sk, &arc_ek);
    });

    result.aggregation_time_ms = TimeMs([&] {
        std::vector<double> line_pass(options.rows, 0.0);
        std::vector<std::uint32_t> line_nation(options.rows, 0);
        for (std::size_t k = 0; k < key_domain; k++) {
            const auto masks = NativeEqualMasks(
                system, join_data.left_key, static_cast<std::uint32_t>(k),
                key_bits, &he_sk, &he_ek, &arc_sk, &arc_ek);
            for (std::size_t i = 0; i < options.rows; i++) {
                if (masks[i]) {
                    const auto nation = key_nation[k];
                    line_pass[i] += order_date_mask[k] * nation_region_mask[nation];
                    line_nation[i] = key_nation[k];
                }
            }
        }
        for (std::size_t n = 0; n < nation_domain; n++) {
            const auto masks = NativeEqualMasks(
                system, line_nation, static_cast<std::uint32_t>(n),
                nation_bits, &he_sk, &he_ek, &arc_sk, &arc_ek);
            for (std::size_t i = 0; i < options.rows; i++)
                result.encrypted_groups[n] += masks[i] * line_pass[i] * revenue[i];
        }
    });
    Finalize(result);
    return result;
}

}  // namespace

NativeBaselineSystem ParseNativeBaselineSystem(const std::string& value) {
    if (value == "ours") return NativeBaselineSystem::Ours;
    if (value == "he3db") return NativeBaselineSystem::HE3DB;
    if (value == "arcedb") return NativeBaselineSystem::ArcEDB;
    if (value == "all") return NativeBaselineSystem::All;
    throw std::invalid_argument("unknown system: " + value);
}

std::string NativeBaselineSystemName(NativeBaselineSystem system) {
    switch (system) {
    case NativeBaselineSystem::Ours:
        return "ours";
    case NativeBaselineSystem::HE3DB:
        return "he3db";
    case NativeBaselineSystem::ArcEDB:
        return "arcedb";
    case NativeBaselineSystem::All:
        return "all";
    }
    throw std::invalid_argument("unknown native baseline system");
}

void PrintNativeBaselineResult(
    const NativeBaselineResult& result,
    std::ostream& os) {
    os << "metric,value\n";
    os << "chapter," << result.chapter << "\n";
    os << "operator," << result.op << "\n";
    os << "query," << result.query << "\n";
    os << "system," << result.system << "\n";
    os << "rows," << result.rows << "\n";
    os << "groups," << result.groups << "\n";
    os << "key_domain," << result.key_domain << "\n";
    os << "enumerated_groups," << result.enumerated_groups << "\n";
    os << "enumerated_keys," << result.enumerated_keys << "\n";
    os << "equality_calls," << result.equality_calls << "\n";
    os << "join_layers," << result.join_layers << "\n";
    os << "filter_time_ms," << result.filter_time_ms << "\n";
    os << "aggregation_time_ms," << result.aggregation_time_ms << "\n";
    os << "total_query_time_ms," << result.total_query_time_ms << "\n";
    os << "latency_seconds," << result.total_query_time_ms / 1000.0 << "\n";
    os << std::fixed << std::setprecision(6);
    os << "plain_result," << result.plain_result << "\n";
    os << "encrypted_result," << result.encrypted_result << "\n";
    os << "accuracy," << result.accuracy << "\n";
    os << "avg_abs_error," << result.avg_abs_error << "\n";
    os << "max_abs_error," << result.max_abs_error << "\n";
    os << "time_source,measured_native_operator\n";
    if (!result.plain_groups.empty()) {
        os << "group,plain,encrypted,abs_error\n";
        for (std::size_t i = 0; i < result.plain_groups.size(); i++)
            os << i << "," << result.plain_groups[i] << ","
               << result.encrypted_groups[i] << ","
               << std::abs(result.plain_groups[i] - result.encrypted_groups[i])
               << "\n";
    }
}

void WriteNativeBaselineResultIfRequested(
    const NativeBaselineResult& result,
    const std::string& output_path) {
    if (output_path.empty()) return;
    std::ofstream out(output_path, std::ios::app);
    if (!out) throw std::runtime_error("failed to open output file: " + output_path);
    PrintNativeBaselineResult(result, out);
}

std::vector<NativeBaselineResult> RunChap4GroupByBaselines(
    const NativeBaselineOptions& options) {
    if (options.rows == 0 || options.groups == 0)
        throw std::invalid_argument("rows and groups must be nonzero");
    std::vector<NativeBaselineResult> out;
    for (auto system : SystemsToRun(options.system))
        out.push_back(RunGroupByOne(system, options));
    return out;
}

std::vector<NativeBaselineResult> RunChap4JoinBaselines(
    const NativeBaselineOptions& options) {
    if (options.rows == 0 || options.key_domain == 0)
        throw std::invalid_argument("rows and key domain must be nonzero");
    std::vector<NativeBaselineResult> out;
    for (auto system : SystemsToRun(options.system))
        out.push_back(RunJoinOne(system, options));
    return out;
}

std::vector<NativeBaselineResult> RunChap5Q3Baselines(
    const NativeBaselineOptions& options) {
    if (options.rows == 0 || options.key_domain == 0 || options.priority_domain == 0)
        throw std::invalid_argument("q3 rows/key-domain/priority-domain must be nonzero");
    std::vector<NativeBaselineResult> out;
    for (auto system : SystemsToRun(options.system))
        out.push_back(RunQ3One(system, options));
    return out;
}

std::vector<NativeBaselineResult> RunChap5Q5Baselines(
    const NativeBaselineOptions& options) {
    if (options.rows == 0 || options.key_domain == 0 || options.nation_domain == 0)
        throw std::invalid_argument("q5 rows/key-domain/nation-domain must be nonzero");
    std::vector<NativeBaselineResult> out;
    for (auto system : SystemsToRun(options.system))
        out.push_back(RunQ5One(system, options));
    return out;
}

}  // namespace PaperReview
