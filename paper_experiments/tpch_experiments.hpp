#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace PaperReview {

struct ExperimentOptions {
    std::size_t rows = 16;
    std::uint64_t seed = 42;
    std::uint32_t ckks_depth_override = 0;
    std::string output_path;
    std::string system = "ours";

    int compare_kappa = 9;
    int date_bits = 8;
    int quantity_bits = 7;
    int discount_bits = 5;
    int key_bits = 4;

    std::uint32_t q6_shipdate_data_min = 100;
    std::uint32_t q6_shipdate_data_max = 200;
    std::uint32_t q6_shipdate_lower = 130;
    std::uint32_t q6_shipdate_upper = 170;
    std::uint32_t q6_discount_lower = 8;
    std::uint32_t q6_discount_upper = 10;
    std::uint32_t q6_quantity_upper = 32;

    std::uint32_t q14_shipdate_lower = 130;
    std::uint32_t q14_shipdate_upper = 170;
    std::uint32_t q14_shipdate_data_min = 100;
    std::uint32_t q14_shipdate_data_max = 200;
    std::size_t q14_type_domain = 4;
    std::uint32_t q14_promo_type = 0;

    std::size_t q3_key_domain = 4;
    std::size_t q3_priority_domain = 2;
    std::size_t q3_segment_domain = 2;
    std::uint32_t q3_customer_segment = 0;
    std::uint32_t q3_orderdate_base = 120;
    std::uint32_t q3_orderdate_step = 10;
    std::uint32_t q3_orderdate_upper = 150;
    std::uint32_t q3_shipdate_data_min = 100;
    std::uint32_t q3_shipdate_data_max = 200;
    std::uint32_t q3_shipdate_lower = 150;

    std::size_t q5_key_domain = 4;
    std::size_t q5_nation_domain = 4;
    std::size_t q5_region_domain = 2;
    std::uint32_t q5_region = 0;
    std::uint32_t q5_orderdate_base = 120;
    std::uint32_t q5_orderdate_step = 10;
    std::uint32_t q5_orderdate_lower = 130;
    std::uint32_t q5_orderdate_upper = 170;
};

struct StageTiming {
    std::string category;
    std::string name;
    double milliseconds = 0.0;
};

struct ExperimentResult {
    std::string query_name;
    std::size_t rows = 0;
    std::size_t slots = 0;
    std::string mask_mode;
    double plain_scalar = 0.0;
    double encrypted_scalar = 0.0;
    double abs_error = 0.0;
    double filter_time_ms = 0.0;
    double aggregation_time_ms = 0.0;
    double total_query_time_ms = 0.0;
    std::size_t predicate_errors = 0;
    std::string system = "ours";
    std::vector<double> plain_groups;
    std::vector<double> encrypted_groups;
    std::vector<StageTiming> timings;
    std::vector<std::pair<std::string, std::string>> parameters;
};

ExperimentOptions ParseExperimentOptions(int argc, char** argv);
void PrintExperimentResult(const ExperimentResult& result, std::ostream& os);
void WriteExperimentResultIfRequested(
    const ExperimentResult& result,
    const ExperimentOptions& options);

ExperimentResult RunTpchQ6Experiment(const ExperimentOptions& options);
ExperimentResult RunTpchQ14Experiment(const ExperimentOptions& options);
ExperimentResult RunTpchQ3Experiment(const ExperimentOptions& options);
ExperimentResult RunTpchQ5Experiment(const ExperimentOptions& options);

}  // namespace PaperReview
