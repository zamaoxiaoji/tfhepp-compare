#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace PaperReview {

enum class NativeBaselineSystem {
    Ours,
    HE3DB,
    ArcEDB,
    All,
};

struct NativeBaselineOptions {
    NativeBaselineSystem system = NativeBaselineSystem::Ours;
    std::size_t rows = 0;
    std::size_t groups = 0;
    std::size_t key_domain = 0;
    std::size_t priority_domain = 2;
    std::size_t nation_domain = 4;
    std::uint64_t seed = 42;
    std::string output_path;
};

struct NativeBaselineResult {
    std::string chapter;
    std::string op;
    std::string query;
    std::string system;
    std::size_t rows = 0;
    std::size_t groups = 0;
    std::size_t key_domain = 0;
    std::size_t enumerated_groups = 0;
    std::size_t enumerated_keys = 0;
    std::size_t equality_calls = 0;
    std::size_t join_layers = 0;
    double filter_time_ms = 0.0;
    double aggregation_time_ms = 0.0;
    double total_query_time_ms = 0.0;
    double plain_result = 0.0;
    double encrypted_result = 0.0;
    double accuracy = 1.0;
    double avg_abs_error = 0.0;
    double max_abs_error = 0.0;
    std::vector<double> plain_groups;
    std::vector<double> encrypted_groups;
};

NativeBaselineSystem ParseNativeBaselineSystem(const std::string& value);
std::string NativeBaselineSystemName(NativeBaselineSystem system);

void PrintNativeBaselineResult(
    const NativeBaselineResult& result,
    std::ostream& os);
void WriteNativeBaselineResultIfRequested(
    const NativeBaselineResult& result,
    const std::string& output_path);

std::vector<NativeBaselineResult> RunChap4GroupByBaselines(
    const NativeBaselineOptions& options);
std::vector<NativeBaselineResult> RunChap4JoinBaselines(
    const NativeBaselineOptions& options);
std::vector<NativeBaselineResult> RunChap5Q3Baselines(
    const NativeBaselineOptions& options);
std::vector<NativeBaselineResult> RunChap5Q5Baselines(
    const NativeBaselineOptions& options);

}  // namespace PaperReview
