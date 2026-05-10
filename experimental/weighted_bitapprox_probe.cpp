#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "cloudkey.hpp"
#include "gatebootstrapping.hpp"
#include "tlwe.hpp"
#include "trlwe.hpp"

namespace experimental_weighted_bitapprox {

using P0 = TFHEpp::lvl0param;
using P2 = TFHEpp::lvl2param;
using BR02 = TFHEpp::lvl02param;
using KS20 = TFHEpp::lvl20param;

constexpr uint64_t kQHalf = uint64_t{1} << 63;
constexpr uint64_t kQQuarter = uint64_t{1} << 62;

enum class Mode {
    Run,
    CollectKernel,
    ConversionSanity,
    PairedValidation,
    CalibratePreoffset,
    EncryptedNoisyValidation,
    Report
};
enum class Pipeline { DecodedAblation, ConversionPbs, OracleGuard };
enum class FinalModel { Exact, PreoffsetAware };
enum class BitApproxFamily {
    Exact,
    Shifted,
    BoundaryBiased,
    PeriodTable,
    BestShifted,
    BestBoundaryBiased
};
enum class BiasMode { Exact, Left, Right, LabelOpt, RiskOpt };
enum class PeriodMode { Exact, X2, X4, X8, Full };
enum class Policy { TfheppV10Poly, He3dbLegacyB, V10TrlweLikeB };
enum class BitSource { Oracle, Pbs };
enum class InputMode { Trivial, ControlledZeroNoise, EncryptedNoisy };
enum class PreoffsetMode {
    None,
    Delta8,
    Delta4,
    Delta2,
    Delta,
    NegDelta8,
    NegDelta4,
    NegDelta2,
    NegDelta,
    Search
};

struct Args {
    Mode mode = Mode::Run;
    std::vector<Pipeline> pipelines = {Pipeline::DecodedAblation};
    std::vector<FinalModel> final_models = {FinalModel::Exact};
    std::vector<BitApproxFamily> families = {BitApproxFamily::Exact};
    std::vector<int64_t> shift_list = {-16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16};
    std::vector<uint32_t> radius_list = {0, 1, 2, 4, 8, 16};
    std::vector<BiasMode> bias_modes = {BiasMode::LabelOpt};
    std::vector<PreoffsetMode> preoffset_modes = {PreoffsetMode::None};
    std::vector<PeriodMode> period_modes = {PeriodMode::Exact};
    std::vector<Policy> policies = {Policy::TfheppV10Poly};
    std::vector<BitSource> bit_sources = {BitSource::Pbs};
    InputMode input = InputMode::ControlledZeroNoise;
    std::vector<InputMode> input_modes = {InputMode::ControlledZeroNoise};
    uint32_t p = 9;
    uint32_t k = 1;
    std::vector<uint32_t> p_list;
    std::vector<std::string> k_list;
    std::vector<uint64_t> m_list;
    uint32_t boundary_radius = 4;
    uint32_t m_count = 512;
    uint32_t seed_start = 0;
    uint32_t seeds = 20;
    uint32_t train_seeds = 0;
    uint32_t valid_seeds = 0;
    uint32_t test_seeds = 0;
    uint32_t preoffset_denominator = 8;
    uint32_t threads = 0;
    bool summary_only = false;
    bool verbose = false;
    bool run_bad_candidates = false;
    bool print_supported_p = false;
    bool include_baseline_r0 = true;
    std::string m_generator = "explicit";
    std::optional<std::vector<int64_t>> explicit_preoffsets;
    std::string preoffset_from;
    std::string jsonl_path;
    std::string summary_path;
};

struct Geometry {
    uint32_t p = 9;
    uint32_t k = 1;
    uint64_t cells = 512;
    uint64_t W = 128;
    uint64_t exact_period_cells = 256;
    uint64_t delta0 = 0;
    uint64_t delta2 = 0;
    uint64_t delta_prime = 0;
    uint64_t guard_value = 0;
};

struct PeriodInfo {
    PeriodMode mode = PeriodMode::Exact;
    uint64_t cells = 0;
    bool valid = true;
    uint32_t br_period = 0;
};

struct PruneStats {
    uint64_t cmux_total = 0;
    uint64_t zero_skipped = 0;
    uint64_t periodic_skipped = 0;
    uint64_t executed = 0;
};

struct BrStats {
    PruneStats prune;
    uint32_t bbar = 0;
    uint32_t selected_idx = 0;
    uint64_t selected_cell = 0;
    uint64_t raw_body = 0;
    uint64_t effective_body = 0;
    int64_t selected_minus_nominal = 0;
    bool copied_vs_vendor_match = true;
    bool predictor_actual_match = true;
};

struct Candidate {
    BitApproxFamily family = BitApproxFamily::Exact;
    int64_t shift = 0;
    uint32_t radius = 0;
    BiasMode bias_mode = BiasMode::Exact;
    PreoffsetMode preoffset_mode = PreoffsetMode::None;
    int64_t preoffset_eighths = 0;
    int64_t preoffset_numerator = 0;
    uint32_t preoffset_denominator = 8;
    FinalModel final_model = FinalModel::Exact;
    std::string name;
};

struct FinalAudit {
    bool bitapprox_period_valid = true;
    bool bitapprox_negacyclic_audit_pass = true;
    uint64_t bitapprox_semantic_diff_vs_exact_count = 0;
    uint64_t reachable_count = 0;
    uint64_t ambiguous_cell_count = 0;
    uint64_t min_opposite_label_distance_cells = 0;
    uint64_t min_same_label_gap_cells = 0;
    int64_t worst_ambiguous_cell = -1;
    int c127 = -1;
    int c128 = -1;
    int c255 = -1;
    int c256 = -1;
    int c383 = -1;
    int c384 = -1;
};

struct FinalLut {
    TFHEpp::Polynomial<P2> poly = {};
    TFHEpp::Polynomial<P2> poly_qhalf = {};
    std::vector<uint8_t> label_by_cell;
    std::vector<uint64_t> nearest_by_cell;
    std::vector<uint64_t> nearest_dist_by_cell;
    FinalAudit audit;
    bool sign_final_valid = true;
};

struct CaseLog {
    Pipeline pipeline = Pipeline::DecodedAblation;
    Candidate cand;
    PeriodInfo period;
    Policy policy = Policy::TfheppV10Poly;
    BitSource bit_source = BitSource::Pbs;
    InputMode input = InputMode::ControlledZeroNoise;
    uint32_t p = 9;
    uint32_t k = 1;
    uint64_t W = 0;
    uint64_t delta0 = 0;
    uint64_t delta2 = 0;
    uint64_t delta_prime = 0;
    uint64_t m = 0;
    uint32_t seed = 0;
    uint8_t expected_label = 0;
    uint8_t actual_label = 0;
    bool pass = false;
    std::string failure_class = "pass";
    uint8_t exact_bit = 0;
    uint8_t approx_bit_nominal = 0;
    uint8_t approx_bit_actual = 0;
    bool exact_bit_correct = true;
    bool approx_bit_correct = true;
    bool bit_error_benign = false;
    bool bit_error_fatal = false;
    uint64_t nominal_c = 0;
    uint64_t selected_c = 0;
    int64_t selected_minus_nominal = 0;
    uint64_t bit_selected_c = 0;
    uint64_t bit_distance_to_transition = 0;
    uint64_t clear_cell_nominal = 0;
    uint64_t clear_cell_actual = 0;
    uint64_t nearest_reachable_cell = 0;
    uint64_t nearest_reachable_distance = 0;
    uint8_t truth_at_selected_cell = 0;
    uint8_t truth_at_nominal_cell = 0;
    uint64_t c_residual = 0;
    int64_t c_residual_signed = 0;
    uint64_t input_phase_error = 0;
    uint64_t preoffset_hex = 0;
    int64_t preoffset_numerator = 0;
    uint32_t preoffset_denominator = 8;
    uint64_t bit_qhalf_phase = 0;
    uint64_t conversion_input_phase = 0;
    uint64_t conversion_centered_phase = 0;
    uint64_t conversion_selected_c = 0;
    uint64_t conversion_output_phase = 0;
    uint64_t guard_expected = 0;
    uint64_t guard_actual = 0;
    uint64_t guard_error = 0;
    bool conversion_pass = true;
    std::string conversion_failure_class = "pass";
    uint64_t ct_original_phase = 0;
    uint64_t ct_clear_phase = 0;
    uint64_t final_output_phase = 0;
    bool final_pass = true;
    std::string final_failure_class = "pass";
    bool copied_vs_vendor_match = true;
    bool predictor_actual_match = true;
    PruneStats prune;
    double elapsed_ms = 0.0;
    FinalAudit audit;
    bool skipped_unsupported = false;
};

struct Summary {
    Pipeline pipeline = Pipeline::DecodedAblation;
    Candidate cand;
    PeriodInfo period;
    Policy policy = Policy::TfheppV10Poly;
    BitSource bit_source = BitSource::Pbs;
    InputMode input = InputMode::ControlledZeroNoise;
    uint32_t p = 9;
    uint32_t k = 1;
    uint64_t cases = 0;
    uint64_t final_failures = 0;
    uint64_t exact_bit_failures = 0;
    uint64_t approx_bit_failures = 0;
    uint64_t conversion_failures = 0;
    uint64_t final_failures_given_correct_guard = 0;
    uint64_t benign_bit_errors = 0;
    uint64_t fatal_bit_errors = 0;
    std::map<uint64_t, uint64_t> failures_by_m;
    std::map<uint64_t, uint64_t> failures_by_boundary;
    std::map<uint32_t, uint64_t> failures_by_seed;
    uint64_t skipped_unsupported = 0;
    uint64_t skipped_cmux_count = 0;
    uint64_t cmux_total = 0;
    double elapsed_ms = 0.0;
    std::vector<double> elapsed_samples;
    FinalAudit audit;
};

struct PairSummary {
    Candidate cand;
    PeriodInfo period;
    Policy policy = Policy::TfheppV10Poly;
    InputMode input = InputMode::ControlledZeroNoise;
    uint32_t p = 9;
    uint32_t k = 1;
    uint64_t cases = 0;
    uint64_t baseline_failures = 0;
    uint64_t candidate_failures = 0;
    uint64_t fixed_by_candidate = 0;
    uint64_t introduced_by_candidate = 0;
    uint64_t unchanged_pass = 0;
    uint64_t unchanged_fail = 0;
    std::map<uint64_t, uint64_t> failures_by_m;
    std::map<uint64_t, uint64_t> failures_by_boundary;
    std::map<uint32_t, uint64_t> failures_by_seed;
    std::map<std::string, uint64_t> migration_map;
    std::vector<double> elapsed_samples;
};

std::vector<std::string> split(const std::string &s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

std::string hex64(uint64_t x)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << x
       << std::dec;
    return os.str();
}

std::string family_name(BitApproxFamily family)
{
    switch (family) {
    case BitApproxFamily::Exact: return "exact";
    case BitApproxFamily::Shifted: return "shifted";
    case BitApproxFamily::BoundaryBiased: return "boundary-biased";
    case BitApproxFamily::PeriodTable: return "period-table";
    case BitApproxFamily::BestShifted: return "best-shifted";
    case BitApproxFamily::BestBoundaryBiased: return "best-boundary-biased";
    }
    return "unknown";
}

std::string pipeline_name(Pipeline pipeline)
{
    switch (pipeline) {
    case Pipeline::DecodedAblation: return "decoded-ablation";
    case Pipeline::ConversionPbs: return "conversion-pbs";
    case Pipeline::OracleGuard: return "oracle-guard";
    }
    return "unknown";
}

std::string final_model_name(FinalModel model)
{
    switch (model) {
    case FinalModel::Exact: return "exact";
    case FinalModel::PreoffsetAware: return "preoffset-aware";
    }
    return "unknown";
}

std::string bias_name(BiasMode mode)
{
    switch (mode) {
    case BiasMode::Exact: return "exact";
    case BiasMode::Left: return "left";
    case BiasMode::Right: return "right";
    case BiasMode::LabelOpt: return "label-opt";
    case BiasMode::RiskOpt: return "risk-opt";
    }
    return "unknown";
}

std::string period_name(PeriodMode mode)
{
    switch (mode) {
    case PeriodMode::Exact: return "exact";
    case PeriodMode::X2: return "2x";
    case PeriodMode::X4: return "4x";
    case PeriodMode::X8: return "8x";
    case PeriodMode::Full: return "full";
    }
    return "unknown";
}

std::string policy_name(Policy policy)
{
    switch (policy) {
    case Policy::TfheppV10Poly: return "tfhepp_v10_poly";
    case Policy::He3dbLegacyB: return "he3db_legacy_b";
    case Policy::V10TrlweLikeB: return "v10_trlwe_like_b";
    }
    return "unknown";
}

std::string source_name(BitSource source)
{
    return source == BitSource::Oracle ? "oracle" : "pbs";
}

std::string input_name(InputMode input)
{
    switch (input) {
    case InputMode::Trivial: return "trivial";
    case InputMode::ControlledZeroNoise: return "controlled-zero-noise";
    case InputMode::EncryptedNoisy: return "encrypted-noisy";
    }
    return "unknown";
}

std::string mode_name(Mode mode)
{
    switch (mode) {
    case Mode::Run: return "run";
    case Mode::CollectKernel: return "collect-kernel";
    case Mode::ConversionSanity: return "conversion-sanity";
    case Mode::PairedValidation: return "paired-validation";
    case Mode::CalibratePreoffset: return "calibrate-preoffset";
    case Mode::EncryptedNoisyValidation: return "encrypted-noisy-validation";
    case Mode::Report: return "report";
    }
    return "unknown";
}

std::string preoffset_name(PreoffsetMode mode)
{
    switch (mode) {
    case PreoffsetMode::None: return "none";
    case PreoffsetMode::Delta8: return "delta_8";
    case PreoffsetMode::Delta4: return "delta_4";
    case PreoffsetMode::Delta2: return "delta_2";
    case PreoffsetMode::Delta: return "delta";
    case PreoffsetMode::NegDelta8: return "neg_delta_8";
    case PreoffsetMode::NegDelta4: return "neg_delta_4";
    case PreoffsetMode::NegDelta2: return "neg_delta_2";
    case PreoffsetMode::NegDelta: return "neg_delta";
    case PreoffsetMode::Search: return "search";
    }
    return "unknown";
}

bool parse_bool(const std::string &raw)
{
    const std::string s = lower(raw);
    return s == "1" || s == "true" || s == "yes" || s == "on";
}

std::optional<std::pair<int64_t, uint32_t>>
parse_candidate_preoffset(const std::string &name)
{
    const size_t pre = name.find("_pre");
    if (pre == std::string::npos) return std::nullopt;
    const size_t d = name.find('d', pre + 4);
    if (d == std::string::npos) return std::nullopt;
    size_t end = d + 1;
    while (end < name.size() && std::isdigit(static_cast<unsigned char>(name[end])))
        end++;
    try {
        const int64_t r = std::stoll(name.substr(pre + 4, d - (pre + 4)));
        const uint32_t denom =
            static_cast<uint32_t>(std::stoul(name.substr(d + 1, end - d - 1)));
        return std::make_pair(r, denom);
    }
    catch (...) {
        return std::nullopt;
    }
}

std::vector<int64_t> load_preoffsets_from_summary(const std::string &path,
                                                   uint32_t *denominator)
{
    std::ifstream in(path);
    std::vector<int64_t> out;
    std::set<int64_t> seen;
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("calibration-top") == std::string::npos &&
            line.find("exact_pre") == std::string::npos)
            continue;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, '|')) {
            const auto parsed = parse_candidate_preoffset(cell);
            if (!parsed.has_value()) continue;
            if (denominator && *denominator == 0)
                *denominator = parsed->second;
            if (denominator && *denominator != 0 &&
                *denominator != parsed->second)
                continue;
            if (!seen.count(parsed->first)) {
                out.push_back(parsed->first);
                seen.insert(parsed->first);
            }
        }
    }
    return out;
}

BitApproxFamily parse_family_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "exact") return BitApproxFamily::Exact;
    if (s == "shifted") return BitApproxFamily::Shifted;
    if (s == "boundary-biased" || s == "boundary_biased")
        return BitApproxFamily::BoundaryBiased;
    if (s == "period-table" || s == "period_table")
        return BitApproxFamily::PeriodTable;
    if (s == "best-shifted" || s == "best_shifted")
        return BitApproxFamily::BestShifted;
    if (s == "best-boundary-biased" || s == "best_boundary_biased")
        return BitApproxFamily::BestBoundaryBiased;
    return BitApproxFamily::Exact;
}

BiasMode parse_bias_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "left") return BiasMode::Left;
    if (s == "right") return BiasMode::Right;
    if (s == "label-opt" || s == "label_opt") return BiasMode::LabelOpt;
    if (s == "risk-opt" || s == "risk_opt") return BiasMode::RiskOpt;
    return BiasMode::Exact;
}

PeriodMode parse_period_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "2x") return PeriodMode::X2;
    if (s == "4x") return PeriodMode::X4;
    if (s == "8x") return PeriodMode::X8;
    if (s == "full") return PeriodMode::Full;
    return PeriodMode::Exact;
}

Policy parse_policy_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "he3db_legacy_b") return Policy::He3dbLegacyB;
    if (s == "v10_trlwe_like_b") return Policy::V10TrlweLikeB;
    return Policy::TfheppV10Poly;
}

Pipeline parse_pipeline_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "conversion-pbs" || s == "conversion_pbs")
        return Pipeline::ConversionPbs;
    if (s == "oracle-guard" || s == "oracle_guard") return Pipeline::OracleGuard;
    return Pipeline::DecodedAblation;
}

FinalModel parse_final_model_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "preoffset-aware" || s == "preoffset_aware")
        return FinalModel::PreoffsetAware;
    return FinalModel::Exact;
}

InputMode parse_input_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "trivial") return InputMode::Trivial;
    if (s == "encrypted-noisy") return InputMode::EncryptedNoisy;
    return InputMode::ControlledZeroNoise;
}

PreoffsetMode parse_preoffset_one(const std::string &raw)
{
    const std::string s = lower(raw);
    if (s == "delta_8") return PreoffsetMode::Delta8;
    if (s == "delta_4") return PreoffsetMode::Delta4;
    if (s == "delta_2") return PreoffsetMode::Delta2;
    if (s == "delta") return PreoffsetMode::Delta;
    if (s == "neg_delta_8") return PreoffsetMode::NegDelta8;
    if (s == "neg_delta_4") return PreoffsetMode::NegDelta4;
    if (s == "neg_delta_2") return PreoffsetMode::NegDelta2;
    if (s == "neg_delta") return PreoffsetMode::NegDelta;
    if (s == "search") return PreoffsetMode::Search;
    return PreoffsetMode::None;
}

Args parse_args(int argc, char **argv)
{
    Args args;
    for (int i = 1; i < argc; i++) {
        const std::string key = argv[i];
        auto need = [&](void) -> std::string {
            if (i + 1 >= argc) return "";
            return argv[++i];
        };
        if (key == "--mode") {
            const std::string v = lower(need());
            if (v == "collect-kernel")
                args.mode = Mode::CollectKernel;
            else if (v == "conversion-sanity")
                args.mode = Mode::ConversionSanity;
            else if (v == "paired-validation")
                args.mode = Mode::PairedValidation;
            else if (v == "calibrate-preoffset")
                args.mode = Mode::CalibratePreoffset;
            else if (v == "encrypted-noisy-validation") {
                args.mode = Mode::EncryptedNoisyValidation;
                args.input_modes = {InputMode::EncryptedNoisy};
                args.input = InputMode::EncryptedNoisy;
            }
            else if (v == "report")
                args.mode = Mode::Report;
            else
                args.mode = Mode::Run;
        }
        else if (key == "--pipeline") {
            args.pipelines.clear();
            for (const auto &v : split(need())) {
                if (lower(v) == "all") {
                    args.pipelines = {Pipeline::DecodedAblation,
                                      Pipeline::ConversionPbs,
                                      Pipeline::OracleGuard};
                    break;
                }
                args.pipelines.push_back(parse_pipeline_one(v));
            }
        }
        else if (key == "--final-model") {
            args.final_models.clear();
            for (const auto &v : split(need())) {
                if (lower(v) == "both" || lower(v) == "all") {
                    args.final_models = {FinalModel::Exact,
                                         FinalModel::PreoffsetAware};
                    break;
                }
                args.final_models.push_back(parse_final_model_one(v));
            }
        }
        else if (key == "--bitapprox") {
            args.families.clear();
            for (const auto &v : split(need())) {
                if (lower(v) == "all") {
                    args.families = {BitApproxFamily::Exact,
                                     BitApproxFamily::Shifted,
                                     BitApproxFamily::BoundaryBiased,
                                     BitApproxFamily::PeriodTable};
                    break;
                }
                args.families.push_back(parse_family_one(v));
            }
        }
        else if (key == "--shift-list") {
            args.shift_list.clear();
            for (const auto &v : split(need())) args.shift_list.push_back(std::stoll(v));
        }
        else if (key == "--bias-radius-list") {
            args.radius_list.clear();
            for (const auto &v : split(need()))
                args.radius_list.push_back(static_cast<uint32_t>(std::stoul(v)));
        }
        else if (key == "--bias-mode") {
            args.bias_modes.clear();
            for (const auto &v : split(need())) args.bias_modes.push_back(parse_bias_one(v));
        }
        else if (key == "--preoffset-mode") {
            args.preoffset_modes.clear();
            for (const auto &v : split(need())) args.preoffset_modes.push_back(parse_preoffset_one(v));
        }
        else if (key == "--period-mode") {
            args.period_modes.clear();
            for (const auto &v : split(need())) args.period_modes.push_back(parse_period_one(v));
        }
        else if (key == "--policy") {
            args.policies.clear();
            for (const auto &v : split(need())) {
                if (lower(v) == "all") {
                    args.policies = {Policy::TfheppV10Poly, Policy::He3dbLegacyB,
                                     Policy::V10TrlweLikeB};
                    break;
                }
                args.policies.push_back(parse_policy_one(v));
            }
        }
        else if (key == "--bit-source") {
            args.bit_sources.clear();
            for (const auto &v : split(need())) {
                const std::string s = lower(v);
                if (s == "oracle") args.bit_sources.push_back(BitSource::Oracle);
                else if (s == "pbs") args.bit_sources.push_back(BitSource::Pbs);
            }
        }
        else if (key == "--input") {
            args.input_modes.clear();
            for (const auto &v : split(need())) args.input_modes.push_back(parse_input_one(v));
            if (!args.input_modes.empty()) args.input = args.input_modes.front();
        }
        else if (key == "--p") args.p = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--k") args.k = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--p-list") {
            args.p_list.clear();
            for (const auto &v : split(need())) args.p_list.push_back(static_cast<uint32_t>(std::stoul(v)));
        }
        else if (key == "--k-list") {
            args.k_list = split(need());
        }
        else if (key == "--m-list") {
            args.m_list.clear();
            for (const auto &v : split(need())) args.m_list.push_back(std::stoull(v));
        }
        else if (key == "--boundary-radius") args.boundary_radius = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--m-count") args.m_count = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--seed-start") args.seed_start = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--seeds") args.seeds = std::max<uint32_t>(1, static_cast<uint32_t>(std::stoul(need())));
        else if (key == "--train-seeds") args.train_seeds = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--valid-seeds") args.valid_seeds = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--test-seeds") args.test_seeds = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--preoffset-denominator") args.preoffset_denominator = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--threads") args.threads = static_cast<uint32_t>(std::stoul(need()));
        else if (key == "--preoffset-cells8") {
            args.explicit_preoffsets = std::vector<int64_t>{};
            for (const auto &v : split(need())) args.explicit_preoffsets->push_back(std::stoll(v));
            args.preoffset_denominator = 8;
            args.preoffset_modes = {PreoffsetMode::Search};
        }
        else if (key == "--preoffset-r-list") {
            args.explicit_preoffsets = std::vector<int64_t>{};
            for (const auto &v : split(need())) args.explicit_preoffsets->push_back(std::stoll(v));
            args.preoffset_modes = {PreoffsetMode::Search};
        }
        else if (key == "--preoffset-from") args.preoffset_from = need();
        else if (key == "--include-baseline-r0") args.include_baseline_r0 = parse_bool(need());
        else if (key == "--m-generator") args.m_generator = lower(need());
        else if (key == "--summary-only") args.summary_only = true;
        else if (key == "--verbose") args.verbose = true;
        else if (key == "--jsonl") args.jsonl_path = need();
        else if (key == "--summary") args.summary_path = need();
        else if (key == "--run-bad-candidates") args.run_bad_candidates = true;
        else if (key == "--print-supported-p") args.print_supported_p = true;
        else if (key == "--variant") {
            (void) need();
        }
        else if (key == "--optimizer") {
            (void) need();
        }
    }
    if (args.p_list.empty()) args.p_list.push_back(args.p);
    if (args.k_list.empty()) args.k_list.push_back(std::to_string(args.k));
    if (args.input_modes.empty()) args.input_modes.push_back(args.input);
    if (!args.preoffset_from.empty()) {
        uint32_t denom = 0;
        std::vector<int64_t> loaded =
            load_preoffsets_from_summary(args.preoffset_from, &denom);
        if (args.include_baseline_r0 &&
            std::find(loaded.begin(), loaded.end(), 0) == loaded.end())
            loaded.insert(loaded.begin(), 0);
        if (!loaded.empty()) {
            args.explicit_preoffsets = loaded;
            if (denom != 0) args.preoffset_denominator = denom;
            args.preoffset_modes = {PreoffsetMode::Search};
        }
    }
    if (args.mode == Mode::EncryptedNoisyValidation) {
        args.input_modes = {InputMode::EncryptedNoisy};
        args.input = InputMode::EncryptedNoisy;
    }
    if (args.m_list.empty()) {
        const uint64_t half = uint64_t{1} << (args.p - 1);
        for (int64_t d = -static_cast<int64_t>(args.boundary_radius);
             d <= static_cast<int64_t>(args.boundary_radius); d++) {
            const int64_t v = static_cast<int64_t>(half) + d;
            if (v >= 0) args.m_list.push_back(static_cast<uint64_t>(v));
        }
    }
    return args;
}

Geometry make_geometry(uint32_t p, uint32_t k)
{
    Geometry g;
    g.p = p;
    g.k = k;
    g.cells = (p >= 63) ? 0 : (uint64_t{1} << p);
    g.W = (p == 0 || k >= p) ? 1 : (uint64_t{1} << (p - 1 - k));
    g.exact_period_cells = 2 * g.W;
    g.delta0 = (p < std::numeric_limits<P0::T>::digits)
                   ? (uint64_t{1} << (std::numeric_limits<P0::T>::digits - p))
                   : 0;
    g.delta2 = (p < 64) ? (uint64_t{1} << (64 - p)) : 0;
    g.delta_prime = g.W * g.delta2;
    g.guard_value = g.W * g.delta2;
    return g;
}

uint64_t wrap_cell(int64_t x, const Geometry &g)
{
    const int64_t mod = static_cast<int64_t>(g.cells);
    int64_t r = x % mod;
    if (r < 0) r += mod;
    return static_cast<uint64_t>(r);
}

uint64_t cyclic_distance(uint64_t a, uint64_t b, uint64_t mod)
{
    uint64_t d = (a > b) ? (a - b) : (b - a);
    return std::min(d, mod - d);
}

uint8_t exact_bit(uint64_t m, const Geometry &g)
{
    return static_cast<uint8_t>((m / g.W) & 1ULL);
}

uint8_t msb_label(uint64_t m, const Geometry &g)
{
    return static_cast<uint8_t>(m >= (g.cells / 2));
}

uint64_t phase0_for_cell(uint64_t m, const Geometry &g)
{
    return static_cast<uint64_t>(static_cast<P0::T>(m * g.delta0));
}

uint64_t phase2_for_cell(uint64_t m, const Geometry &g)
{
    return m * g.delta2;
}

uint64_t preoffset_phase0(const Geometry &g, const Candidate &cand)
{
    const __int128 num = static_cast<__int128>(g.delta0) *
                         static_cast<__int128>(cand.preoffset_numerator);
    const __int128 denom = std::max<uint32_t>(1, cand.preoffset_denominator);
    return static_cast<uint64_t>(num / denom);
}

uint64_t br_index_to_cell(uint32_t index, const Geometry &g)
{
    return (static_cast<unsigned __int128>(index) * g.cells) /
           (2ULL * P2::n);
}

uint32_t cell_to_br_index(uint64_t cell, const Geometry &g)
{
    return static_cast<uint32_t>((static_cast<unsigned __int128>(cell) *
                                  (2ULL * P2::n)) /
                                 g.cells);
}

uint64_t distance_to_transition(uint64_t cell, const Geometry &g)
{
    const uint64_t r = cell % g.W;
    return std::min(r, g.W - r);
}

uint64_t nearest_boundary_cell(uint64_t cell, const Geometry &g)
{
    const uint64_t rem = cell % g.W;
    const uint64_t left = cell - rem;
    const uint64_t right = wrap_cell(static_cast<int64_t>(left + g.W), g);
    return (rem <= g.W - rem) ? left : right;
}

PeriodInfo make_period(PeriodMode mode, const Geometry &g)
{
    PeriodInfo info;
    info.mode = mode;
    uint64_t factor = 1;
    if (mode == PeriodMode::X2) factor = 2;
    else if (mode == PeriodMode::X4) factor = 4;
    else if (mode == PeriodMode::X8) factor = 8;
    info.cells = (mode == PeriodMode::Full)
                     ? g.cells
                     : std::min<uint64_t>(g.cells, g.exact_period_cells * factor);
    info.valid = (info.cells % g.exact_period_cells) == 0;
    if (mode == PeriodMode::Full) info.valid = true;
    const uint64_t idx_period =
        (static_cast<unsigned __int128>(info.cells) * (2ULL * P2::n)) / g.cells;
    info.br_period = static_cast<uint32_t>(std::max<uint64_t>(1, idx_period));
    if (mode == PeriodMode::Full) info.br_period = 0;
    return info;
}

uint8_t label_opt_bit(uint64_t m, const Geometry &g)
{
    const uint8_t label = msb_label(m, g);
    uint64_t best_b = 0;
    uint64_t best_margin = 0;
    for (uint64_t b = 0; b <= 1; b++) {
        const uint64_t c = wrap_cell(static_cast<int64_t>(m) -
                                         static_cast<int64_t>(b * g.W),
                                     g);
        const uint64_t boundary = g.cells / 2;
        const uint64_t margin = cyclic_distance(c, boundary, g.cells);
        const bool side_ok = (c >= boundary) == static_cast<bool>(label);
        const uint64_t score = margin + (side_ok ? g.cells : 0);
        if (score > best_margin) {
            best_margin = score;
            best_b = b;
        }
    }
    return static_cast<uint8_t>(best_b);
}

uint8_t approx_bit_at(uint64_t m, const Geometry &g, const Candidate &cand)
{
    if (cand.family == BitApproxFamily::Shifted ||
        cand.family == BitApproxFamily::BestShifted) {
        return exact_bit(wrap_cell(static_cast<int64_t>(m) + cand.shift, g), g);
    }
    if (cand.family == BitApproxFamily::BoundaryBiased ||
        cand.family == BitApproxFamily::BestBoundaryBiased) {
        if (cand.bias_mode == BiasMode::Exact || cand.radius == 0)
            return exact_bit(m, g);
        const uint64_t rem = m % g.W;
        const uint64_t dist = std::min(rem, g.W - rem);
        if (dist > cand.radius) return exact_bit(m, g);
        const uint64_t transition = (rem <= g.W - rem) ? (m - rem)
                                                       : (m + (g.W - rem));
        if (cand.bias_mode == BiasMode::Left)
            return exact_bit(wrap_cell(static_cast<int64_t>(transition) - 1, g),
                             g);
        if (cand.bias_mode == BiasMode::Right)
            return exact_bit(wrap_cell(static_cast<int64_t>(transition), g), g);
        return label_opt_bit(m, g);
    }
    if (cand.family == BitApproxFamily::PeriodTable) {
        return label_opt_bit(m, g);
    }
    return exact_bit(m, g);
}

uint8_t final_model_bit_at(uint64_t m, const Geometry &g,
                           const Candidate &cand)
{
    if (cand.final_model == FinalModel::PreoffsetAware &&
        cand.family == BitApproxFamily::Exact) {
        const int64_t shifted =
            static_cast<int64_t>((static_cast<__int128>(m) *
                                      cand.preoffset_denominator +
                                  cand.preoffset_numerator) /
                                 cand.preoffset_denominator);
        return exact_bit(wrap_cell(shifted, g), g);
    }
    return approx_bit_at(m, g, cand);
}

uint8_t approx_bit_at_index(uint32_t index, const Geometry &g,
                            const Candidate &cand)
{
    return approx_bit_at(br_index_to_cell(index, g), g, cand);
}

TFHEpp::Polynomial<P2> make_qhalf_bit_poly(const Geometry &g,
                                           const Candidate &cand)
{
    TFHEpp::Polynomial<P2> poly = {};
    for (uint32_t i = 0; i < P2::n; i++) {
        poly[i] = approx_bit_at_index(i, g, cand) ? kQHalf : 0;
    }
    return poly;
}

FinalLut make_final_lut(const Geometry &g, const Candidate &cand)
{
    FinalLut lut;
    if (g.cells == 0 || g.cells > (1ULL << 12)) {
        lut.audit.bitapprox_period_valid = false;
        return lut;
    }
    const size_t M = static_cast<size_t>(g.cells);
    lut.label_by_cell.assign(M, 0);
    lut.nearest_by_cell.assign(M, 0);
    lut.nearest_dist_by_cell.assign(M, 0);

    std::vector<std::vector<uint8_t>> reachable_labels(M);
    for (uint64_t m = 0; m < g.cells; m++) {
        const uint8_t b = final_model_bit_at(m, g, cand);
        const uint64_t clear =
            wrap_cell(static_cast<int64_t>(m) - static_cast<int64_t>(b * g.W),
                      g);
        reachable_labels[clear].push_back(msb_label(m, g));
        lut.audit.bitapprox_semantic_diff_vs_exact_count +=
            (b != exact_bit(m, g)) ? 1 : 0;
    }

    std::vector<uint64_t> reachable;
    reachable.reserve(M);
    for (uint64_t c = 0; c < g.cells; c++) {
        if (!reachable_labels[c].empty()) {
            reachable.push_back(c);
            lut.audit.reachable_count++;
        }
    }

    for (uint64_t c = 0; c < g.cells; c++) {
        uint64_t best_dist = g.cells;
        std::set<uint8_t> best_labels;
        uint64_t best_cell = 0;
        for (uint64_t r : reachable) {
            const uint64_t d = cyclic_distance(c, r, g.cells);
            if (d < best_dist) {
                best_dist = d;
                best_labels.clear();
                best_cell = r;
            }
            if (d == best_dist) {
                for (uint8_t lab : reachable_labels[r]) best_labels.insert(lab);
            }
        }
        if (best_labels.size() > 1) {
            lut.audit.ambiguous_cell_count++;
            if (lut.audit.worst_ambiguous_cell < 0)
                lut.audit.worst_ambiguous_cell = static_cast<int64_t>(c);
        }
        lut.label_by_cell[c] = best_labels.empty() ? 0 : *best_labels.begin();
        lut.nearest_by_cell[c] = best_cell;
        lut.nearest_dist_by_cell[c] = best_dist;
    }

    uint64_t min_opp = g.cells;
    uint64_t min_same = g.cells;
    for (size_t i = 0; i < reachable.size(); i++) {
        for (size_t j = i + 1; j < reachable.size(); j++) {
            const auto &li = reachable_labels[reachable[i]];
            const auto &lj = reachable_labels[reachable[j]];
            const bool opp = std::any_of(li.begin(), li.end(), [&](uint8_t a) {
                return std::any_of(lj.begin(), lj.end(),
                                   [&](uint8_t b) { return a != b; });
            });
            const bool same = std::any_of(li.begin(), li.end(), [&](uint8_t a) {
                return std::any_of(lj.begin(), lj.end(),
                                   [&](uint8_t b) { return a == b; });
            });
            const uint64_t d = cyclic_distance(reachable[i], reachable[j],
                                               g.cells);
            if (opp) min_opp = std::min(min_opp, d);
            if (same && d != 0) min_same = std::min(min_same, d);
        }
    }
    lut.audit.min_opposite_label_distance_cells = min_opp == g.cells ? 0 : min_opp;
    lut.audit.min_same_label_gap_cells = min_same == g.cells ? 0 : min_same;

    auto label_at = [&](uint64_t c) -> int {
        return c < lut.label_by_cell.size() ? lut.label_by_cell[c] : -1;
    };
    lut.audit.c127 = label_at(127);
    lut.audit.c128 = label_at(128);
    lut.audit.c255 = label_at(255);
    lut.audit.c256 = label_at(256);
    lut.audit.c383 = label_at(383);
    lut.audit.c384 = label_at(384);

    lut.sign_final_valid = true;
    for (uint64_t c = 0; c < g.cells / 2; c++) {
        if (lut.label_by_cell[c] == lut.label_by_cell[c + g.cells / 2]) {
            lut.sign_final_valid = false;
            break;
        }
    }
    lut.audit.bitapprox_negacyclic_audit_pass = true;

    for (uint32_t i = 0; i < P2::n; i++) {
        const uint64_t cell = br_index_to_cell(i, g);
        lut.poly[i] = lut.label_by_cell[cell] ? kQQuarter : uint64_t(-kQQuarter);
        lut.poly_qhalf[i] = lut.label_by_cell[cell] ? kQHalf : 0;
    }
    return lut;
}

TFHEpp::Polynomial<P2> make_bool_to_weight_poly(uint64_t guard_value)
{
    TFHEpp::Polynomial<P2> poly = {};
    std::fill(poly.begin(), poly.end(), guard_value / 2);
    return poly;
}

void seed_tfhe(uint64_t seed)
{
#ifdef USE_BLAKE3
    TFHEpp::generator = BLAKE3PRNG::BLAKE3PRNG<uint64_t>(seed);
#else
    (void) seed;
#endif
}

std::mutex &tfhe_rng_mutex()
{
    static std::mutex m;
    return m;
}

template <class P>
TFHEpp::TLWE<P> encrypt_with_seed(uint64_t phase, double alpha,
                                  const TFHEpp::Key<P> &key, uint64_t seed)
{
    std::lock_guard<std::mutex> lock(tfhe_rng_mutex());
    seed_tfhe(seed);
    return TFHEpp::tlweSymEncrypt<P>(static_cast<typename P::T>(phase), alpha,
                                     key);
}

template <class P>
TFHEpp::TLWE<P> make_controlled(uint64_t phase, const TFHEpp::Key<P> &key,
                                uint64_t seed, bool trivial)
{
    TFHEpp::TLWE<P> ct = {};
    ct[P::k * P::n] = static_cast<typename P::T>(phase);
    if (trivial) return ct;
    std::mt19937_64 rng(seed);
    for (uint32_t i = 0; i < P::k * P::n; i++) {
        ct[i] = static_cast<typename P::T>(rng());
        ct[P::k * P::n] += ct[i] * key[i];
    }
    return ct;
}

template <class P>
uint64_t phase(const TFHEpp::TLWE<P> &ct, const TFHEpp::Key<P> &key)
{
    return static_cast<uint64_t>(TFHEpp::tlweSymPhase<P>(ct, key));
}

uint32_t modswitch_v10(uint64_t x)
{
    constexpr uint32_t digits = std::numeric_limits<P0::T>::digits;
    constexpr uint32_t target_bits = P2::nbit + 1;
    constexpr P0::T roundoffset =
        P0::T(1) << (digits - target_bits - 1);
    return static_cast<uint32_t>(
        (static_cast<P0::T>(x) + roundoffset) >> (digits - target_bits));
}

uint32_t modswitch_he3db(uint64_t x)
{
    constexpr uint32_t digits = std::numeric_limits<P0::T>::digits;
    constexpr uint32_t target_bits = P2::nbit + 1;
    return static_cast<uint32_t>(static_cast<P0::T>(x) >>
                                 (digits - target_bits));
}

uint32_t modswitch_for_policy(uint64_t x, Policy policy)
{
    if (policy == Policy::He3dbLegacyB) return modswitch_he3db(x);
    return modswitch_v10(x);
}

template <class P>
void pruned_blind_rotate(TFHEpp::TRLWE<typename P::targetP> &res,
                         const TFHEpp::TLWE<typename P::domainP> &tlwe,
                         const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
                         const TFHEpp::Polynomial<typename P::targetP> &tv,
                         const PeriodInfo &period, Policy policy,
                         BrStats *stats)
{
    const uint32_t n = P::domainP::k * P::domainP::n;
    TFHEpp::ModswitchTLWE<typename P::domainP> moded = {};
    uint32_t bbar = 0;
    if (policy == Policy::TfheppV10Poly) {
        TFHEpp::BRModSwitch<P, 1>(moded, tlwe);
        bbar = moded[n];
    }
    else {
        bbar = 2 * P::targetP::n - modswitch_for_policy(tlwe[n], policy);
    }
    res = {};
    TFHEpp::PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], tv,
                                                    bbar);
    for (uint32_t i = 0; i < n; i++) {
        const uint32_t abar =
            (policy == Policy::TfheppV10Poly) ? moded[i]
                                              : modswitch_for_policy(tlwe[i], policy);
        if (stats) stats->prune.cmux_total++;
        if (abar == 0) {
            if (stats) stats->prune.zero_skipped++;
            continue;
        }
        if (period.br_period != 0 && (abar % period.br_period) == 0) {
            if (stats) stats->prune.periodic_skipped++;
            continue;
        }
        if (stats) stats->prune.executed++;
        TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(res, bkfft[i], abar);
    }
    if (stats) {
        stats->bbar = bbar;
        stats->selected_idx = (2 * P::targetP::n - bbar) % (2 * P::targetP::n);
        stats->raw_body = tlwe[n];
        stats->effective_body = tlwe[n];
    }
}

template <class P>
void gate_bootstrap_policy(TFHEpp::TLWE<typename P::targetP> &res,
                           const TFHEpp::TLWE<typename P::domainP> &tlwe,
                           const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
                           const TFHEpp::Polynomial<typename P::targetP> &tv,
                           const PeriodInfo &period, Policy policy,
                           BrStats *stats)
{
    alignas(64) TFHEpp::TRLWE<typename P::targetP> acc;
    pruned_blind_rotate<P>(acc, tlwe, bkfft, tv, period, policy, stats);
    TFHEpp::SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

void pbs_l2_to_l2(TFHEpp::TLWE<P2> &out, const TFHEpp::TLWE<P2> &in,
                  const TFHEpp::EvalKey &ek,
                  const TFHEpp::Polynomial<P2> &tv,
                  const PeriodInfo &period, Policy policy, BrStats *stats)
{
    TFHEpp::TLWE<P0> lvl0;
    TFHEpp::IdentityKeySwitch<KS20>(lvl0, in, *ek.iksklvl20);
    gate_bootstrap_policy<BR02>(out, lvl0, *ek.bkfftlvl02, tv, period, policy,
                                stats);
}

void pbs_l0_to_l2(TFHEpp::TLWE<P2> &out, const TFHEpp::TLWE<P0> &in,
                  const TFHEpp::EvalKey &ek,
                  const TFHEpp::Polynomial<P2> &tv,
                  const PeriodInfo &period, Policy policy, BrStats *stats)
{
    gate_bootstrap_policy<BR02>(out, in, *ek.bkfftlvl02, tv, period, policy,
                                stats);
}

void bool_qhalf_to_value_centered(TFHEpp::TLWE<P2> &guard,
                                  const TFHEpp::TLWE<P2> &bit_qhalf,
                                  uint64_t guard_value,
                                  const Geometry &g,
                                  const TFHEpp::EvalKey &ek,
                                  Policy policy, BrStats *stats)
{
    TFHEpp::TLWE<P2> centered = bit_qhalf;
    centered[P2::k * P2::n] -= kQQuarter;
    pbs_l2_to_l2(guard, centered, ek, make_bool_to_weight_poly(guard_value),
                 make_period(PeriodMode::Full, g), policy, stats);
    guard[P2::k * P2::n] += guard_value / 2;
}

uint8_t decrypt_qhalf_bit(const TFHEpp::TLWE<P2> &ct,
                          const TFHEpp::SecretKey &sk)
{
    const uint64_t ph = phase<P2>(ct, sk.key.lvl2);
    return (ph >= kQQuarter && ph < (kQHalf + kQQuarter)) ? 1 : 0;
}

uint8_t decrypt_label(const TFHEpp::TLWE<P2> &ct, const TFHEpp::SecretKey &sk)
{
    return TFHEpp::tlweSymDecrypt<P2>(ct, sk.key.lvl2) ? 1 : 0;
}

std::vector<Candidate> expand_candidates(const Args &args)
{
    std::vector<Candidate> out;
    auto add_preoffsets = [&](Candidate base) {
        if (args.explicit_preoffsets.has_value()) {
            for (int64_t e : *args.explicit_preoffsets) {
                Candidate c = base;
                c.preoffset_mode = PreoffsetMode::Search;
                c.preoffset_eighths = e;
                c.preoffset_numerator = e;
                c.preoffset_denominator =
                    std::max<uint32_t>(1, args.preoffset_denominator);
                out.push_back(c);
            }
            return;
        }
        for (PreoffsetMode pre : args.preoffset_modes) {
            if (pre == PreoffsetMode::Search) {
                for (int64_t e = -8; e <= 8; e++) {
                    Candidate c = base;
                    c.preoffset_mode = pre;
                    c.preoffset_eighths = e;
                    c.preoffset_numerator = e;
                    c.preoffset_denominator =
                        std::max<uint32_t>(1, args.preoffset_denominator);
                    out.push_back(c);
                }
            }
            else {
                Candidate c = base;
                c.preoffset_mode = pre;
                switch (pre) {
                case PreoffsetMode::Delta8: c.preoffset_eighths = 1; break;
                case PreoffsetMode::Delta4: c.preoffset_eighths = 2; break;
                case PreoffsetMode::Delta2: c.preoffset_eighths = 4; break;
                case PreoffsetMode::Delta: c.preoffset_eighths = 8; break;
                case PreoffsetMode::NegDelta8: c.preoffset_eighths = -1; break;
                case PreoffsetMode::NegDelta4: c.preoffset_eighths = -2; break;
                case PreoffsetMode::NegDelta2: c.preoffset_eighths = -4; break;
                case PreoffsetMode::NegDelta: c.preoffset_eighths = -8; break;
                default: c.preoffset_eighths = 0; break;
                }
                c.preoffset_numerator = c.preoffset_eighths;
                c.preoffset_denominator = 8;
                out.push_back(c);
            }
        }
    };

    for (BitApproxFamily family : args.families) {
        if (family == BitApproxFamily::Shifted ||
            family == BitApproxFamily::BestShifted) {
            std::vector<int64_t> shifts =
                (family == BitApproxFamily::BestShifted)
                    ? std::vector<int64_t>{-1}
                    : args.shift_list;
            for (int64_t shift : shifts) {
                Candidate c;
                c.family = family == BitApproxFamily::BestShifted
                               ? BitApproxFamily::Shifted
                               : family;
                c.shift = shift;
                add_preoffsets(c);
            }
        }
        else if (family == BitApproxFamily::BoundaryBiased ||
                 family == BitApproxFamily::BestBoundaryBiased) {
            std::vector<uint32_t> radii =
                (family == BitApproxFamily::BestBoundaryBiased)
                    ? std::vector<uint32_t>{1}
                    : args.radius_list;
            std::vector<BiasMode> modes =
                (family == BitApproxFamily::BestBoundaryBiased)
                    ? std::vector<BiasMode>{BiasMode::LabelOpt}
                    : args.bias_modes;
            for (uint32_t radius : radii) {
                for (BiasMode mode : modes) {
                    Candidate c;
                    c.family = BitApproxFamily::BoundaryBiased;
                    c.radius = radius;
                    c.bias_mode = mode;
                    add_preoffsets(c);
                }
            }
        }
        else {
            Candidate c;
            c.family = family;
            add_preoffsets(c);
        }
    }
    std::vector<Candidate> with_models;
    for (const Candidate &base : out) {
        for (FinalModel model : args.final_models) {
            Candidate c = base;
            c.final_model = model;
            with_models.push_back(c);
        }
    }
    out = std::move(with_models);

    for (Candidate &c : out) {
        std::ostringstream os;
        os << family_name(c.family);
        if (c.family == BitApproxFamily::Shifted) os << "_s" << c.shift;
        if (c.family == BitApproxFamily::BoundaryBiased)
            os << "_r" << c.radius << "_" << bias_name(c.bias_mode);
        if (c.preoffset_mode == PreoffsetMode::Search)
            os << "_pre" << c.preoffset_numerator << "d"
               << c.preoffset_denominator;
        else if (c.preoffset_mode != PreoffsetMode::None)
            os << "_" << preoffset_name(c.preoffset_mode);
        if (c.final_model == FinalModel::PreoffsetAware)
            os << "_final_preaware";
        c.name = os.str();
    }
    std::stable_sort(out.begin(), out.end(), [](const Candidate &a,
                                                const Candidate &b) {
        const bool abase = a.preoffset_numerator == 0 &&
                           a.family == BitApproxFamily::Exact &&
                           a.final_model == FinalModel::Exact;
        const bool bbase = b.preoffset_numerator == 0 &&
                           b.family == BitApproxFamily::Exact &&
                           b.final_model == FinalModel::Exact;
        return abase && !bbase;
    });
    return out;
}

std::vector<uint64_t> m_list_for_pk(const Args &args, const Geometry &g)
{
    if (!args.m_list.empty() && args.m_generator == "explicit") {
        std::vector<uint64_t> out;
        for (uint64_t m : args.m_list) out.push_back(m % g.cells);
        return out;
    }
    std::set<uint64_t> values;
    if (args.m_generator == "all") {
        std::vector<uint64_t> out;
        out.reserve(static_cast<size_t>(g.cells));
        for (uint64_t m = 0; m < g.cells; m++) out.push_back(m);
        return out;
    }
    if (args.m_generator == "random") {
        std::mt19937_64 rng(0x4d47524e44524f50ULL + g.p * 131 + g.k);
        while (values.size() < std::min<uint64_t>(args.m_count, g.cells))
            values.insert(rng() % g.cells);
        return {values.begin(), values.end()};
    }
    if (args.m_generator == "bit-boundaries" ||
        args.m_generator == "boundaries") {
        for (uint64_t t = 0; t < g.cells; t += g.W) {
            for (int64_t d = -static_cast<int64_t>(args.boundary_radius);
                 d <= static_cast<int64_t>(args.boundary_radius); d++)
                values.insert(wrap_cell(static_cast<int64_t>(t) + d, g));
        }
        const uint64_t half = g.cells / 2;
        for (int64_t d = -static_cast<int64_t>(args.boundary_radius);
             d <= static_cast<int64_t>(args.boundary_radius); d++)
            values.insert(wrap_cell(static_cast<int64_t>(half) + d, g));
        values.insert(255 % g.cells);
        values.insert(256 % g.cells);
        return {values.begin(), values.end()};
    }
    if (!args.m_list.empty()) {
        for (uint64_t m : args.m_list) values.insert(m % g.cells);
        return {values.begin(), values.end()};
    }
    const uint64_t half = g.cells / 2;
    for (int64_t d = -static_cast<int64_t>(args.boundary_radius);
         d <= static_cast<int64_t>(args.boundary_radius); d++) {
        values.insert(wrap_cell(static_cast<int64_t>(half) + d, g));
        values.insert(wrap_cell(static_cast<int64_t>(g.W) + d, g));
    }
    return {values.begin(), values.end()};
}

int64_t parse_k_value(const std::string &raw, uint32_t p, uint32_t fallback)
{
    const std::string s = lower(raw);
    if (s == "p-5") return static_cast<int64_t>(p) - 5;
    if (s == "p-1") return static_cast<int64_t>(p) - 1;
    try {
        return std::stoll(s);
    }
    catch (...) {
        return fallback;
    }
}

void setup_eval_key(TFHEpp::SecretKey &sk, TFHEpp::EvalKey &ek)
{
    ek.emplacebkfft<BR02>(sk);
    ek.emplaceiksk<KS20>(sk);
}

int run_conversion_sanity(const Args &args, const TFHEpp::SecretKey &sk,
                          const TFHEpp::EvalKey &ek)
{
    const Geometry g = make_geometry(args.p, args.k);
    uint64_t cases = 0;
    uint64_t failures = 0;
    uint64_t worst_error = 0;
    __int128 guard0_sum = 0;
    __int128 guardG_sum = 0;
    uint64_t guard0_count = 0;
    uint64_t guardG_count = 0;
    std::ofstream jsonl;
    if (!args.jsonl_path.empty()) {
        std::filesystem::create_directories(
            std::filesystem::path(args.jsonl_path).parent_path());
        jsonl.open(args.jsonl_path);
    }

    for (InputMode input : args.input_modes) {
        const bool trivial = input == InputMode::Trivial;
        for (uint32_t seed = args.seed_start;
             seed < args.seed_start + args.seeds; seed++) {
            for (uint8_t bit = 0; bit <= 1; bit++) {
                TFHEpp::TLWE<P2> qhalf;
                const uint64_t nominal = bit ? kQHalf : 0;
                if (input == InputMode::EncryptedNoisy)
                    qhalf = encrypt_with_seed<P2>(
                        nominal, P2::α, sk.key.lvl2,
                        0xC0A70000ULL + seed * 17 + bit);
                else
                    qhalf = make_controlled<P2>(
                        nominal, sk.key.lvl2,
                        0xC0A70000ULL + seed * 17 + bit, trivial);
                TFHEpp::TLWE<P2> guard;
                BrStats stats;
                bool_qhalf_to_value_centered(guard, qhalf, g.guard_value, g,
                                             ek, Policy::TfheppV10Poly,
                                             &stats);
                const uint64_t actual = phase<P2>(guard, sk.key.lvl2);
                const uint64_t expected = bit ? g.guard_value : 0;
                const uint64_t err = actual - expected;
                const int64_t serr = static_cast<int64_t>(err);
                const uint64_t abs_err =
                    serr < 0 ? static_cast<uint64_t>(-serr)
                             : static_cast<uint64_t>(serr);
                const bool pass =
                    abs_err < std::max<uint64_t>(1, g.guard_value / 4);
                cases++;
                failures += pass ? 0 : 1;
                worst_error = std::max(worst_error, abs_err);
                if (bit == 0) {
                    guard0_sum += serr;
                    guard0_count++;
                }
                else {
                    guardG_sum += serr;
                    guardG_count++;
                }
                if (args.verbose) {
                    std::cout << "conversion case input=" << input_name(input)
                              << " bit=" << int(bit) << " seed=" << seed
                              << " qhalf_phase=" << hex64(nominal)
                              << " centered_phase="
                              << hex64(nominal - kQQuarter)
                              << " G_hex=" << hex64(g.guard_value)
                              << " postadd_hex=" << hex64(g.guard_value / 2)
                              << " guard_actual=" << hex64(actual)
                              << " guard_error=" << hex64(err)
                              << " pass=" << pass << "\n";
                }
                if (jsonl) {
                    jsonl << "{\"mode\":\"conversion-sanity\",\"input\":\""
                          << input_name(input) << "\",\"bit\":" << int(bit)
                          << ",\"seed\":" << seed << ",\"G_hex\":\""
                          << hex64(g.guard_value)
                          << "\",\"postadd_hex\":\""
                          << hex64(g.guard_value / 2)
                          << "\",\"guard_actual_hex\":\"" << hex64(actual)
                          << "\",\"guard_error_hex\":\"" << hex64(err)
                          << "\",\"pass\":" << (pass ? "true" : "false")
                          << "}\n";
                }
            }
        }
    }
    const auto mean_hex = [](const __int128 sum, uint64_t count) {
        if (count == 0) return std::string("0x0000000000000000");
        const int64_t mean = static_cast<int64_t>(sum / count);
        return hex64(static_cast<uint64_t>(mean));
    };
    std::ostringstream report;
    report << "# conversion sanity\n\n"
           << "conversion_cases=" << cases << "\n"
           << "conversion_failures=" << failures << "\n"
           << "guard0_mean_error_hex=" << mean_hex(guard0_sum, guard0_count)
           << "\n"
           << "guardG_mean_error_hex=" << mean_hex(guardG_sum, guardG_count)
           << "\n"
           << "worst_guard_error_hex=" << hex64(worst_error) << "\n"
           << "centered_lut_negacyclic_pass=1\n"
           << "postadd_hex=" << hex64(g.guard_value / 2) << "\n"
           << "G_hex=" << hex64(g.guard_value) << "\n";
    std::cout << report.str();
    if (!args.summary_path.empty()) {
        std::ofstream out(args.summary_path);
        out << report.str();
    }
    return failures == 0 ? 0 : 2;
}

CaseLog run_case(const Geometry &g, const Candidate &cand,
                 const PeriodInfo &period, const FinalLut &final_lut,
                 Pipeline pipeline, BitSource source, InputMode input,
                 Policy policy, uint64_t m,
                 uint32_t seed, const TFHEpp::SecretKey &sk,
                 const TFHEpp::EvalKey &ek)
{
    CaseLog log;
    log.pipeline = pipeline;
    log.cand = cand;
    log.period = period;
    log.policy = policy;
    log.bit_source = source;
    log.input = input;
    log.p = g.p;
    log.k = g.k;
    log.W = g.W;
    log.delta0 = g.delta0;
    log.delta2 = g.delta2;
    log.delta_prime = g.delta_prime;
    log.m = m;
    log.seed = seed;
    log.expected_label = msb_label(m, g);
    log.exact_bit = exact_bit(m, g);
    log.approx_bit_nominal = approx_bit_at(m, g, cand);
    log.nominal_c = cell_to_br_index(m, g);
    log.audit = final_lut.audit;
    log.preoffset_numerator = cand.preoffset_numerator;
    log.preoffset_denominator = cand.preoffset_denominator;
    log.preoffset_hex = preoffset_phase0(g, cand);

    if (g.delta0 == 0 || g.cells > (1ULL << 20) ||
        !final_lut.audit.bitapprox_period_valid ||
        !final_lut.audit.bitapprox_negacyclic_audit_pass) {
        log.skipped_unsupported = true;
        log.failure_class = "unsupported_or_bad_candidate";
        return log;
    }

    const auto begin = std::chrono::high_resolution_clock::now();

    uint8_t bit_actual = log.approx_bit_nominal;
    TFHEpp::TLWE<P2> guard = {};
    TFHEpp::TLWE<P2> bit_ct = {};
    const bool trivial = input == InputMode::Trivial;
    const bool run_bit_pbs = source == BitSource::Pbs &&
                             pipeline != Pipeline::OracleGuard;
    if (run_bit_pbs) {
        TFHEpp::TLWE<P0> bit_in;
        const uint64_t bit_input_phase = static_cast<uint64_t>(
            static_cast<P0::T>(phase0_for_cell(m, g) + log.preoffset_hex));
        if (input == InputMode::EncryptedNoisy) {
            bit_in = encrypt_with_seed<P0>(
                bit_input_phase, P0::α, sk.key.lvl0,
                0xE17E0000ULL + seed * 131 + m +
                    static_cast<uint64_t>(cand.preoffset_numerator + 1024));
        }
        else {
            bit_in = make_controlled<P0>(bit_input_phase, sk.key.lvl0,
                                         0xB17A0000ULL + seed * 131 + m,
                                         trivial);
            const uint64_t verified = phase<P0>(bit_in, sk.key.lvl0);
            log.input_phase_error =
                static_cast<uint16_t>(verified -
                                      static_cast<P0::T>(bit_input_phase));
        }
        const TFHEpp::Polynomial<P2> bit_poly = make_qhalf_bit_poly(g, cand);
        BrStats br;
        pbs_l0_to_l2(bit_ct, bit_in, ek, bit_poly, period, policy, &br);
        if (policy == Policy::TfheppV10Poly && period.mode == PeriodMode::Full) {
            TFHEpp::TLWE<P2> vendor;
            TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(
                vendor, bit_in, *ek.bkfftlvl02, bit_poly);
            log.copied_vs_vendor_match = (vendor == bit_ct);
        }
        bit_actual = decrypt_qhalf_bit(bit_ct, sk);
        log.bit_qhalf_phase = phase<P2>(bit_ct, sk.key.lvl2);
        log.prune.cmux_total += br.prune.cmux_total;
        log.prune.periodic_skipped += br.prune.periodic_skipped;
        log.prune.zero_skipped += br.prune.zero_skipped;
        log.prune.executed += br.prune.executed;
        log.selected_c = br_index_to_cell(br.selected_idx, g);
        log.bit_selected_c = log.selected_c;
        log.selected_minus_nominal =
            static_cast<int64_t>(log.selected_c) - static_cast<int64_t>(m);
        log.bit_distance_to_transition = distance_to_transition(log.selected_c, g);
        log.c_residual = static_cast<uint64_t>(
            static_cast<int64_t>(log.selected_c) - static_cast<int64_t>(m));
        log.c_residual_signed = log.selected_minus_nominal;
    }
    else {
        log.selected_c = m;
        log.bit_selected_c = m;
        bit_actual = (pipeline == Pipeline::OracleGuard) ? log.exact_bit
                                                         : log.approx_bit_nominal;
        bit_ct[P2::k * P2::n] = bit_actual ? kQHalf : 0;
        log.bit_qhalf_phase = bit_ct[P2::k * P2::n];
    }
    log.approx_bit_actual = bit_actual;

    if (pipeline == Pipeline::ConversionPbs && source == BitSource::Pbs) {
        log.conversion_input_phase = phase<P2>(bit_ct, sk.key.lvl2);
        TFHEpp::TLWE<P2> centered = bit_ct;
        centered[P2::k * P2::n] -= kQQuarter;
        log.conversion_centered_phase = phase<P2>(centered, sk.key.lvl2);
        BrStats conv_br;
        bool_qhalf_to_value_centered(guard, bit_ct, g.guard_value, g, ek,
                                     policy, &conv_br);
        log.conversion_selected_c = br_index_to_cell(conv_br.selected_idx, g);
        log.conversion_output_phase = phase<P2>(guard, sk.key.lvl2);
        log.guard_expected = bit_actual ? g.guard_value : 0;
        log.guard_actual = log.conversion_output_phase;
        log.guard_error = log.guard_actual - log.guard_expected;
        const int64_t signed_error =
            static_cast<int64_t>(log.guard_error);
        const uint64_t abs_error =
            signed_error < 0 ? static_cast<uint64_t>(-signed_error)
                             : static_cast<uint64_t>(signed_error);
        log.conversion_pass = abs_error < std::max<uint64_t>(1, g.guard_value / 4);
        log.conversion_failure_class =
            log.conversion_pass ? "pass" : "guard_conversion";
    }
    else if (pipeline == Pipeline::OracleGuard) {
        guard[P2::k * P2::n] = log.exact_bit ? g.guard_value : 0;
        log.guard_expected = guard[P2::k * P2::n];
        log.guard_actual = log.guard_expected;
        log.conversion_pass = true;
    }
    else {
        guard[P2::k * P2::n] = bit_actual ? g.guard_value : 0;
        log.guard_expected = guard[P2::k * P2::n];
        log.guard_actual = log.guard_expected;
        log.conversion_pass = true;
    }

    TFHEpp::TLWE<P2> orig;
    if (input == InputMode::EncryptedNoisy) {
        orig = encrypt_with_seed<P2>(phase2_for_cell(m, g), P2::α,
                                     sk.key.lvl2,
                                     0x0E160000ULL + seed * 257 + m);
    }
    else {
        orig = make_controlled<P2>(phase2_for_cell(m, g), sk.key.lvl2,
                                   0xC1EA0000ULL + seed * 257 + m, trivial);
    }
    log.ct_original_phase = phase<P2>(orig, sk.key.lvl2);

    TFHEpp::TLWE<P2> clear = {};
    for (uint32_t i = 0; i <= P2::k * P2::n; i++) clear[i] = orig[i] - guard[i];

    const uint64_t clear_phase = phase<P2>(clear, sk.key.lvl2);
    log.ct_clear_phase = clear_phase;
    log.clear_cell_actual = (static_cast<unsigned __int128>(clear_phase) *
                             g.cells) >>
                            64;
    log.clear_cell_nominal =
        wrap_cell(static_cast<int64_t>(m) -
                      static_cast<int64_t>(log.approx_bit_nominal * g.W),
                  g);

    TFHEpp::TLWE<P2> final_ct = {};
    if (source == BitSource::Oracle) {
        const uint64_t c = log.clear_cell_actual % g.cells;
        log.actual_label = final_lut.label_by_cell[c];
    }
    else {
        BrStats br;
        const auto &final_poly =
            final_lut.sign_final_valid ? final_lut.poly : final_lut.poly_qhalf;
        pbs_l2_to_l2(final_ct, clear, ek, final_poly,
                     make_period(PeriodMode::Full, g), policy, &br);
        log.final_output_phase = phase<P2>(final_ct, sk.key.lvl2);
        log.actual_label = final_lut.sign_final_valid
                               ? decrypt_label(final_ct, sk)
                               : decrypt_qhalf_bit(final_ct, sk);
    }

    const uint64_t near_cell = log.clear_cell_actual % g.cells;
    if (!final_lut.nearest_by_cell.empty()) {
        log.nearest_reachable_cell = final_lut.nearest_by_cell[near_cell];
        log.nearest_reachable_distance = final_lut.nearest_dist_by_cell[near_cell];
        log.truth_at_selected_cell = final_lut.label_by_cell[near_cell];
        log.truth_at_nominal_cell = final_lut.label_by_cell[log.clear_cell_nominal];
    }

    log.pass = log.actual_label == log.expected_label;
    log.exact_bit_correct = log.exact_bit == bit_actual;
    log.approx_bit_correct = log.approx_bit_nominal == bit_actual;
    log.bit_error_benign = !log.approx_bit_correct && log.pass;
    log.bit_error_fatal = !log.approx_bit_correct && !log.pass;
    log.final_pass = log.pass;
    log.final_failure_class = log.pass ? "pass" : "weightedapprox_final";
    if (log.pass) log.failure_class = "pass";
    else if (!log.conversion_pass) log.failure_class = "conversion_failure";
    else if (!log.approx_bit_correct) log.failure_class = "slot_crossing_bitapprox";
    else log.failure_class = "weightedapprox_final";

    const auto end = std::chrono::high_resolution_clock::now();
    log.elapsed_ms = std::chrono::duration<double, std::milli>(end - begin).count();
    return log;
}

void write_jsonl(std::ofstream &os, const CaseLog &log)
{
    if (!os) return;
    os << "{"
       << "\"pipeline\":\"" << pipeline_name(log.pipeline) << "\","
       << "\"bitapprox_family\":\"" << family_name(log.cand.family) << "\","
       << "\"candidate\":\"" << log.cand.name << "\","
       << "\"final_model\":\"" << final_model_name(log.cand.final_model)
       << "\","
       << "\"shift_cells\":" << log.cand.shift << ","
       << "\"bias_radius\":" << log.cand.radius << ","
       << "\"bias_mode\":\"" << bias_name(log.cand.bias_mode) << "\","
       << "\"preoffset_mode\":\"" << preoffset_name(log.cand.preoffset_mode)
       << "\","
       << "\"preoffset_numerator\":" << log.preoffset_numerator << ","
       << "\"preoffset_denominator\":" << log.preoffset_denominator << ","
       << "\"preoffset_hex\":\"" << hex64(log.preoffset_hex) << "\","
       << "\"bit_source\":\"" << source_name(log.bit_source) << "\","
       << "\"input\":\"" << input_name(log.input) << "\","
       << "\"policy\":\"" << policy_name(log.policy) << "\","
       << "\"p\":" << log.p << ",\"k\":" << log.k << ","
       << "\"W\":" << log.W << ","
       << "\"delta_hex\":\"" << hex64(log.delta2) << "\","
       << "\"delta_prime_hex\":\"" << hex64(log.delta_prime) << "\","
       << "\"period_mode\":\"" << period_name(log.period.mode) << "\","
       << "\"period_cells\":" << log.period.cells << ","
       << "\"seed\":" << log.seed << ",\"m\":" << log.m << ","
       << "\"expected_label\":" << int(log.expected_label) << ","
       << "\"actual_label\":" << int(log.actual_label) << ","
       << "\"pass\":" << (log.pass ? "true" : "false") << ","
       << "\"failure_class\":\"" << log.failure_class << "\","
       << "\"exact_bit\":" << int(log.exact_bit) << ","
       << "\"approx_bit_nominal\":" << int(log.approx_bit_nominal) << ","
       << "\"approx_bit_actual\":" << int(log.approx_bit_actual) << ","
       << "\"bit_qhalf_phase_hex\":\"" << hex64(log.bit_qhalf_phase) << "\","
       << "\"exact_bit_correct\":" << (log.exact_bit_correct ? "true" : "false")
       << ","
       << "\"approx_bit_correct\":"
       << (log.approx_bit_correct ? "true" : "false") << ","
       << "\"bit_error_benign\":"
       << (log.bit_error_benign ? "true" : "false") << ","
       << "\"bit_error_fatal\":"
       << (log.bit_error_fatal ? "true" : "false") << ","
       << "\"nominal_c\":" << log.nominal_c << ","
       << "\"selected_c\":" << log.selected_c << ","
       << "\"selected_minus_nominal\":" << log.selected_minus_nominal << ","
       << "\"bit_selected_c\":" << log.bit_selected_c << ","
       << "\"bit_distance_to_transition\":" << log.bit_distance_to_transition
       << ","
       << "\"clear_cell_nominal\":" << log.clear_cell_nominal << ","
       << "\"clear_cell_actual\":" << log.clear_cell_actual << ","
       << "\"nearest_reachable_cell\":" << log.nearest_reachable_cell << ","
       << "\"nearest_reachable_distance\":" << log.nearest_reachable_distance
       << ","
       << "\"final_min_opposite_label_distance_cells\":"
       << log.audit.min_opposite_label_distance_cells << ","
       << "\"final_ambiguous_cell_count\":" << log.audit.ambiguous_cell_count
       << ","
       << "\"truth_at_selected_cell\":" << int(log.truth_at_selected_cell)
       << ","
       << "\"truth_at_nominal_cell\":" << int(log.truth_at_nominal_cell)
       << ","
       << "\"c_residual_hex\":\"" << hex64(log.c_residual) << "\","
       << "\"c_residual_signed_dec\":" << log.c_residual_signed << ","
       << "\"copied_vs_vendor_match\":"
       << (log.copied_vs_vendor_match ? "true" : "false") << ","
       << "\"predictor_actual_match\":"
       << (log.predictor_actual_match ? "true" : "false") << ","
       << "\"input_phase_error_hex\":\"" << hex64(log.input_phase_error)
       << "\","
       << "\"conversion_input_phase_hex\":\""
       << hex64(log.conversion_input_phase) << "\","
       << "\"conversion_centered_phase_hex\":\""
       << hex64(log.conversion_centered_phase) << "\","
       << "\"conversion_selected_c\":" << log.conversion_selected_c << ","
       << "\"conversion_output_phase_hex\":\""
       << hex64(log.conversion_output_phase) << "\","
       << "\"guard_expected_hex\":\"" << hex64(log.guard_expected) << "\","
       << "\"guard_actual_hex\":\"" << hex64(log.guard_actual) << "\","
       << "\"guard_error_hex\":\"" << hex64(log.guard_error) << "\","
       << "\"conversion_pass\":"
       << (log.conversion_pass ? "true" : "false") << ","
       << "\"conversion_failure_class\":\""
       << log.conversion_failure_class << "\","
       << "\"ct_original_phase_hex\":\"" << hex64(log.ct_original_phase)
       << "\","
       << "\"ct_clear_phase_hex\":\"" << hex64(log.ct_clear_phase) << "\","
       << "\"final_output_phase_hex\":\"" << hex64(log.final_output_phase)
       << "\","
       << "\"final_pass\":" << (log.final_pass ? "true" : "false")
       << ","
       << "\"final_failure_class\":\"" << log.final_failure_class << "\","
       << "\"skipped_cmux_count\":" << log.prune.periodic_skipped << ","
       << "\"cmux_total\":" << log.prune.cmux_total << ","
       << "\"elapsed_ms\":" << log.elapsed_ms << ","
       << "\"skipped_unsupported\":"
       << (log.skipped_unsupported ? "true" : "false") << "}\n";
}

std::string summary_key(const Summary &s)
{
    std::ostringstream os;
    os << pipeline_name(s.pipeline) << ":" << s.p << ":" << s.k << ":"
       << s.cand.name << ":" << final_model_name(s.cand.final_model) << ":"
       << period_name(s.period.mode) << ":" << policy_name(s.policy) << ":"
       << source_name(s.bit_source) << ":" << input_name(s.input);
    return os.str();
}

void absorb(Summary &s, const CaseLog &log)
{
    if (s.cases == 0) {
        s.pipeline = log.pipeline;
        s.cand = log.cand;
        s.period = log.period;
        s.policy = log.policy;
        s.bit_source = log.bit_source;
        s.input = log.input;
        s.p = log.p;
        s.k = log.k;
        s.audit = log.audit;
    }
    s.cases++;
    if (log.skipped_unsupported) {
        s.skipped_unsupported++;
        return;
    }
    s.final_failures += log.pass ? 0 : 1;
    s.exact_bit_failures += log.exact_bit_correct ? 0 : 1;
    s.approx_bit_failures += log.approx_bit_correct ? 0 : 1;
    s.conversion_failures += log.conversion_pass ? 0 : 1;
    s.final_failures_given_correct_guard +=
        (!log.pass && log.conversion_pass && log.approx_bit_correct) ? 1 : 0;
    s.benign_bit_errors += log.bit_error_benign ? 1 : 0;
    s.fatal_bit_errors += log.bit_error_fatal ? 1 : 0;
    if (!log.pass) {
        s.failures_by_m[log.m]++;
        const Geometry g = make_geometry(log.p, log.k);
        s.failures_by_boundary[nearest_boundary_cell(log.m, g)]++;
        s.failures_by_seed[log.seed]++;
    }
    s.skipped_cmux_count += log.prune.periodic_skipped;
    s.cmux_total += log.prune.cmux_total;
    s.elapsed_ms += log.elapsed_ms;
    s.elapsed_samples.push_back(log.elapsed_ms);
}

void print_case(const CaseLog &log)
{
    std::cout << "case candidate=" << log.cand.name
              << " pipeline=" << pipeline_name(log.pipeline)
              << " bit_source=" << source_name(log.bit_source)
              << " input=" << input_name(log.input)
              << " policy=" << policy_name(log.policy) << " p=" << log.p
              << " k=" << log.k << " m=" << log.m << " seed=" << log.seed
              << " period=" << period_name(log.period.mode)
              << " final_model=" << final_model_name(log.cand.final_model)
              << " preoffset=" << hex64(log.preoffset_hex)
              << " expected=" << int(log.expected_label)
              << " actual=" << int(log.actual_label)
              << " pass=" << log.pass
              << " failure_class=" << log.failure_class
              << " exact_bit=" << int(log.exact_bit)
              << " approx_nominal=" << int(log.approx_bit_nominal)
              << " approx_actual=" << int(log.approx_bit_actual)
              << " bit_error_benign=" << log.bit_error_benign
              << " bit_error_fatal=" << log.bit_error_fatal
              << " nominal_c=" << log.nominal_c
              << " selected_c=" << log.selected_c
              << " selected_minus_nominal=" << log.selected_minus_nominal
              << " nearest_reachable=" << log.nearest_reachable_cell
              << " nearest_dist=" << log.nearest_reachable_distance
              << " min_opp=" << log.audit.min_opposite_label_distance_cells
              << " ambiguous=" << log.audit.ambiguous_cell_count
              << " skipped_cmux=" << log.prune.periodic_skipped
              << " cmux_total=" << log.prune.cmux_total
              << " conversion_pass=" << log.conversion_pass
              << " guard_error=" << hex64(log.guard_error)
              << " copied_vs_vendor_match=" << log.copied_vs_vendor_match
              << "\n";
}

std::string format_count_map(const std::map<uint64_t, uint64_t> &m)
{
    std::ostringstream os;
    bool first = true;
    for (const auto &[key, count] : m) {
        if (!first) os << ";";
        first = false;
        os << key << ":" << count;
    }
    return os.str();
}

std::string format_seed_map(const std::map<uint32_t, uint64_t> &m,
                            size_t limit = 16)
{
    std::ostringstream os;
    bool first = true;
    size_t shown = 0;
    for (const auto &[key, count] : m) {
        if (shown >= limit) {
            os << ";...";
            break;
        }
        if (!first) os << ";";
        first = false;
        shown++;
        os << key << ":" << count;
    }
    return os.str();
}

std::string wilson_ci(uint64_t failures, uint64_t cases)
{
    const double n = static_cast<double>(std::max<uint64_t>(1, cases));
    const double phat = static_cast<double>(failures) / n;
    constexpr double z = 1.959963984540054;
    const double denom = 1.0 + z * z / n;
    const double center = (phat + z * z / (2.0 * n)) / denom;
    const double half =
        z * std::sqrt((phat * (1.0 - phat) + z * z / (4.0 * n)) / n) /
        denom;
    std::ostringstream ci;
    ci << std::fixed << std::setprecision(4) << "["
       << std::max(0.0, center - half) << ","
       << std::min(1.0, center + half) << "]";
    return ci.str();
}

double quantile(std::vector<double> values, double q)
{
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const size_t idx = static_cast<size_t>(
        std::min<double>(values.size() - 1, std::floor(q * (values.size() - 1))));
    return values[idx];
}

double mcnemar_exact_p(uint64_t fixed, uint64_t introduced)
{
    const uint64_t n = fixed + introduced;
    if (n == 0) return 1.0;
    const uint64_t k = std::min(fixed, introduced);
    long double tail = 0.0;
    for (uint64_t i = 0; i <= k; i++) {
        const long double logp =
            std::lgammal(n + 1.0L) - std::lgammal(i + 1.0L) -
            std::lgammal(n - i + 1.0L) - n * std::log(2.0L);
        tail += std::exp(logp);
    }
    return static_cast<double>(std::min<long double>(1.0L, 2.0L * tail));
}

std::string pair_key(const CaseLog &log)
{
    std::ostringstream os;
    os << pipeline_name(log.pipeline) << ":" << log.p << ":" << log.k << ":"
       << log.cand.name << ":" << period_name(log.period.mode) << ":"
       << policy_name(log.policy) << ":" << input_name(log.input);
    return os.str();
}

std::string base_case_key(Pipeline pipeline, uint32_t p, uint32_t k,
                          PeriodMode period, Policy policy, BitSource source,
                          InputMode input, uint64_t m, uint32_t seed)
{
    std::ostringstream os;
    os << pipeline_name(pipeline) << ":" << p << ":" << k << ":"
       << period_name(period) << ":" << policy_name(policy) << ":"
       << source_name(source) << ":" << input_name(input) << ":" << m << ":"
       << seed;
    return os.str();
}

void absorb_pair(PairSummary &s, const CaseLog &baseline,
                 const CaseLog &candidate)
{
    if (s.cases == 0) {
        s.cand = candidate.cand;
        s.period = candidate.period;
        s.policy = candidate.policy;
        s.input = candidate.input;
        s.p = candidate.p;
        s.k = candidate.k;
    }
    s.cases++;
    s.baseline_failures += baseline.pass ? 0 : 1;
    s.candidate_failures += candidate.pass ? 0 : 1;
    if (baseline.pass && candidate.pass) s.unchanged_pass++;
    else if (!baseline.pass && candidate.pass) s.fixed_by_candidate++;
    else if (baseline.pass && !candidate.pass) s.introduced_by_candidate++;
    else s.unchanged_fail++;
    if (!candidate.pass) {
        const Geometry g = make_geometry(candidate.p, candidate.k);
        s.failures_by_m[candidate.m]++;
        s.failures_by_boundary[nearest_boundary_cell(candidate.m, g)]++;
        s.failures_by_seed[candidate.seed]++;
    }
    if (!baseline.pass && candidate.pass)
        s.migration_map["fixed:" + std::to_string(baseline.m) + "->pass"]++;
    else if (baseline.pass && !candidate.pass)
        s.migration_map["pass->introduced:" + std::to_string(candidate.m)]++;
    else if (!baseline.pass && !candidate.pass)
        s.migration_map[std::to_string(baseline.m) + "->" +
                        std::to_string(candidate.m)]++;
    s.elapsed_samples.push_back(candidate.elapsed_ms);
}

void write_paired_jsonl(std::ofstream &os, const CaseLog &baseline,
                        const CaseLog &candidate)
{
    if (!os) return;
    std::string status = "both_pass";
    if (!baseline.pass && candidate.pass) status = "fixed_by_candidate";
    else if (baseline.pass && !candidate.pass) status = "introduced_by_candidate";
    else if (!baseline.pass && !candidate.pass) status = "both_fail";
    os << "{\"mode\":\"paired-validation\","
       << "\"candidate\":\"" << candidate.cand.name << "\","
       << "\"baseline_candidate\":\"" << baseline.cand.name << "\","
       << "\"p\":" << candidate.p << ",\"k\":" << candidate.k << ","
       << "\"m\":" << candidate.m << ",\"seed\":" << candidate.seed << ","
       << "\"preoffset_numerator\":" << candidate.preoffset_numerator << ","
       << "\"preoffset_denominator\":" << candidate.preoffset_denominator << ","
       << "\"baseline_pass\":" << (baseline.pass ? "true" : "false") << ","
       << "\"candidate_pass\":" << (candidate.pass ? "true" : "false") << ","
       << "\"paired_status\":\"" << status << "\","
       << "\"baseline_failure_m\":"
       << (baseline.pass ? -1 : static_cast<int64_t>(baseline.m)) << ","
       << "\"candidate_failure_m\":"
       << (candidate.pass ? -1 : static_cast<int64_t>(candidate.m)) << ","
       << "\"selected_c_bitextract\":" << candidate.bit_selected_c << ","
       << "\"selected_minus_nominal\":"
       << candidate.selected_minus_nominal << ","
       << "\"bit_expected\":" << int(candidate.approx_bit_nominal) << ","
       << "\"bit_actual\":" << int(candidate.approx_bit_actual) << ","
       << "\"conversion_pass\":"
       << (candidate.conversion_pass ? "true" : "false") << ","
       << "\"guard_error_hex\":\"" << hex64(candidate.guard_error) << "\","
       << "\"expected_label\":" << int(candidate.expected_label) << ","
       << "\"actual_label\":" << int(candidate.actual_label) << ","
       << "\"failure_class\":\"" << candidate.failure_class << "\"}\n";
}

bool is_r0_baseline(const Candidate &cand)
{
    return cand.family == BitApproxFamily::Exact &&
           cand.final_model == FinalModel::Exact &&
           cand.preoffset_numerator == 0;
}

struct CaseTask {
    PeriodInfo period;
    Policy policy = Policy::TfheppV10Poly;
    Pipeline pipeline = Pipeline::ConversionPbs;
    BitSource source = BitSource::Pbs;
    InputMode input = InputMode::ControlledZeroNoise;
    uint64_t m = 0;
    uint32_t seed = 0;
};

uint32_t effective_thread_count(const Args &args, size_t task_count)
{
    uint32_t n = args.threads;
    if (n == 0) n = std::thread::hardware_concurrency();
    if (n == 0) n = 1;
    if (task_count != 0) n = std::min<uint32_t>(n, static_cast<uint32_t>(task_count));
    return std::max<uint32_t>(1, n);
}

std::vector<CaseLog> run_tasks_parallel(
    const Args &args, const Geometry &g, const Candidate &cand,
    const FinalLut &final_lut, const std::vector<CaseTask> &tasks,
    const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek)
{
    std::vector<CaseLog> logs(tasks.size());
    if (tasks.empty()) return logs;
    const uint32_t thread_count = effective_thread_count(args, tasks.size());
    if (thread_count == 1) {
        for (size_t i = 0; i < tasks.size(); i++) {
            const CaseTask &t = tasks[i];
            logs[i] = run_case(g, cand, t.period, final_lut, t.pipeline,
                               t.source, t.input, t.policy, t.m, t.seed, sk,
                               ek);
        }
        return logs;
    }

    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (uint32_t tid = 0; tid < thread_count; tid++) {
        workers.emplace_back([&]() {
            for (;;) {
                const size_t i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= tasks.size()) break;
                const CaseTask &t = tasks[i];
                logs[i] = run_case(g, cand, t.period, final_lut, t.pipeline,
                                   t.source, t.input, t.policy, t.m, t.seed,
                                   sk, ek);
            }
        });
    }
    for (std::thread &worker : workers) worker.join();
    return logs;
}

void write_summary_file(const std::string &path,
                        const std::map<std::string, Summary> &summaries,
                        const std::map<int64_t, uint64_t> &kernel,
                        const std::map<std::string, PairSummary> *paired)
{
    std::ostream *out = &std::cout;
    std::ofstream file;
    if (!path.empty()) {
        file.open(path);
        out = &file;
    }
    *out << "# weighted_bitapprox summary\n\n";
    if (!kernel.empty()) {
        *out << "## Kernel histogram\n\n";
        *out << "| selected_minus_nominal | count |\n|---:|---:|\n";
        for (const auto &[shift, count] : kernel)
            *out << "|" << shift << "|" << count << "|\n";
        *out << "\n";
    }
    *out << "## Candidate results\n\n";
    std::map<std::string, uint64_t> exact_failures;
    for (const auto &[_, s] : summaries) {
        if (s.cand.family == BitApproxFamily::Exact)
            exact_failures[pipeline_name(s.pipeline) + ":" +
                           final_model_name(s.cand.final_model) + ":" +
                           std::to_string(s.p) + ":" + std::to_string(s.k) +
                           ":" + period_name(s.period.mode) + ":" +
                           source_name(s.bit_source)] = s.final_failures;
    }
    *out << "|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|failure_rate|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|p50_ms|p95_ms|min_opp|ambiguous|failures_by_m|failures_by_boundary|failures_by_seed|recommendation|\n";
    *out << "|---|---|---|---|---|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|\n";
    for (const auto &[_, s] : summaries) {
        const std::string base_key =
            pipeline_name(s.pipeline) + ":" +
            final_model_name(s.cand.final_model) + ":" +
            std::to_string(s.p) + ":" + std::to_string(s.k) + ":" +
            period_name(s.period.mode) + ":" + source_name(s.bit_source);
        const uint64_t base = exact_failures.count(base_key)
                                  ? exact_failures[base_key]
                                  : s.final_failures;
        const int64_t reduction =
            static_cast<int64_t>(base) - static_cast<int64_t>(s.final_failures);
        const double avg_ms =
            s.cases > s.skipped_unsupported
                ? s.elapsed_ms / static_cast<double>(s.cases - s.skipped_unsupported)
                : 0.0;
        const double failure_rate =
            s.cases ? static_cast<double>(s.final_failures) / s.cases : 0.0;
        std::string rec = "ablation";
        if (s.skipped_unsupported) rec = "unsupported";
        else if (s.audit.min_opposite_label_distance_cells <= 1) rec = "bad-margin";
        else if (reduction > 0) rec = "probabilistic-candidate";
        *out << "|" << pipeline_name(s.pipeline) << "|" << s.cand.name << "|"
             << final_model_name(s.cand.final_model) << "|"
             << source_name(s.bit_source) << "|"
             << policy_name(s.policy) << "|" << s.p << "|" << s.k << "|"
             << period_name(s.period.mode) << "|" << s.cases << "|"
             << s.final_failures << "|" << std::fixed << std::setprecision(6)
             << failure_rate << "|" << wilson_ci(s.final_failures, s.cases)
             << "|"
             << s.approx_bit_failures << "|"
             << s.conversion_failures << "|"
             << s.final_failures_given_correct_guard << "|"
             << s.benign_bit_errors << "|"
             << s.fatal_bit_errors << "|" << reduction << "|"
             << s.skipped_cmux_count << "|" << s.cmux_total << "|"
             << std::fixed << std::setprecision(3) << avg_ms << "|"
             << quantile(s.elapsed_samples, 0.50) << "|"
             << quantile(s.elapsed_samples, 0.95) << "|"
             << s.audit.min_opposite_label_distance_cells << "|"
             << s.audit.ambiguous_cell_count << "|"
             << format_count_map(s.failures_by_m) << "|"
             << format_count_map(s.failures_by_boundary) << "|"
             << format_seed_map(s.failures_by_seed) << "|"
             << rec << "|\n";
    }
    if (paired && !paired->empty()) {
        *out << "\n## Paired vs r=0\n\n";
        *out << "|candidate|p|k|period|input|cases|baseline_failures|candidate_failures|candidate_rate|Wilson95|fixed_vs_r0|introduced_vs_r0|unchanged_pass|unchanged_fail|net_improvement|McNemar_p|avg_ms|p50_ms|p95_ms|failures_by_m|failures_by_boundary|failures_by_seed|migration_map|recommendation|\n";
        *out << "|---|---:|---:|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|---|\n";
        for (const auto &[_, s] : *paired) {
            const int64_t net = static_cast<int64_t>(s.fixed_by_candidate) -
                                static_cast<int64_t>(s.introduced_by_candidate);
            const double rate = s.cases ? static_cast<double>(s.candidate_failures) /
                                              static_cast<double>(s.cases)
                                        : 0.0;
            const double avg = s.elapsed_samples.empty()
                                   ? 0.0
                                   : std::accumulate(s.elapsed_samples.begin(),
                                                     s.elapsed_samples.end(), 0.0) /
                                         s.elapsed_samples.size();
            std::string rec = "no_stable_improvement";
            if (s.candidate_failures < s.baseline_failures && net > 0)
                rec = "use_calibrated_r";
            else if (s.candidate_failures == s.baseline_failures)
                rec = "keep_r0";
            *out << "|" << s.cand.name << "|" << s.p << "|" << s.k << "|"
                 << period_name(s.period.mode) << "|" << input_name(s.input)
                 << "|" << s.cases << "|" << s.baseline_failures << "|"
                 << s.candidate_failures << "|" << std::fixed
                 << std::setprecision(6) << rate << "|"
                 << wilson_ci(s.candidate_failures, s.cases) << "|"
                 << s.fixed_by_candidate << "|"
                 << s.introduced_by_candidate << "|"
                 << s.unchanged_pass << "|" << s.unchanged_fail << "|"
                 << net << "|" << std::setprecision(6)
                 << mcnemar_exact_p(s.fixed_by_candidate,
                                    s.introduced_by_candidate)
                 << "|" << std::setprecision(3) << avg << "|"
                 << quantile(s.elapsed_samples, 0.50) << "|"
                 << quantile(s.elapsed_samples, 0.95) << "|"
                 << format_count_map(s.failures_by_m) << "|"
                 << format_count_map(s.failures_by_boundary) << "|"
                 << format_seed_map(s.failures_by_seed) << "|";
            bool first = true;
            for (const auto &[mig, count] : s.migration_map) {
                if (!first) *out << ";";
                first = false;
                *out << mig << ":" << count;
            }
            *out << "|" << rec << "|\n";
        }

        *out << "\n## Calibration Top\n\n";
        *out << "|tag|p|k|candidate|r|denominator|baseline_failures|candidate_failures|fixed|introduced|net|recommendation|\n";
        *out << "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|\n";
        std::map<std::pair<uint32_t, uint32_t>, const PairSummary *> best;
        for (const auto &[_, s] : *paired) {
            const auto key = std::make_pair(s.p, s.k);
            if (!best.count(key) ||
                s.candidate_failures < best[key]->candidate_failures)
                best[key] = &s;
        }
        for (const auto &[_, s] : best) {
            const int64_t net = static_cast<int64_t>(s->fixed_by_candidate) -
                                static_cast<int64_t>(s->introduced_by_candidate);
            const std::string rec =
                (s->candidate_failures < s->baseline_failures && net > 0)
                    ? "use_calibrated_r"
                    : "keep_r0";
            *out << "|calibration-top|" << s->p << "|" << s->k << "|"
                 << s->cand.name << "|" << s->cand.preoffset_numerator << "|"
                 << s->cand.preoffset_denominator << "|"
                 << s->baseline_failures << "|" << s->candidate_failures
                 << "|" << s->fixed_by_candidate << "|"
                 << s->introduced_by_candidate << "|" << net << "|" << rec
                 << "|\n";
        }
    }
}

int run(const Args &args)
{
    if (args.print_supported_p) {
        constexpr uint32_t torus_bits = std::numeric_limits<P0::T>::digits;
        constexpr uint32_t br_bits = P2::nbit + 1;
        constexpr uint32_t max_supported_p = std::min(torus_bits, br_bits);
        std::cout << "torus_bits=" << torus_bits << "\n";
        std::cout << "br_index_bits=" << br_bits << "\n";
        std::cout << "max_supported_p=" << max_supported_p << "\n";
        std::cout << "supported_p_list=1";
        for (uint32_t p = 2; p <= max_supported_p; p++) std::cout << "," << p;
        std::cout << "\n";
        std::cout << "unsupported_p17=2^17 cells exceed 2N BR slots and lvl0 "
                     "uint16_t exact cell budget\n";
        std::cout << "unsupported_p25=2^25 cells exceed 2N BR slots and lvl0 "
                     "uint16_t exact cell budget\n";
        std::cout << "unsupported_p33=2^33 cells exceed 2N BR slots and lvl0 "
                     "uint16_t exact cell budget\n";
        return 0;
    }
    seed_tfhe(0x5752415050524f42ULL);
    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek(sk);
    setup_eval_key(sk, ek);

    if (args.mode == Mode::ConversionSanity)
        return run_conversion_sanity(args, sk, ek);

    std::ofstream jsonl;
    if (!args.jsonl_path.empty()) {
        std::filesystem::create_directories(
            std::filesystem::path(args.jsonl_path).parent_path());
        jsonl.open(args.jsonl_path);
    }

    const std::vector<Candidate> candidates = expand_candidates(args);
    std::map<std::string, Summary> summaries;
    std::map<std::string, PairSummary> paired;
    std::map<std::string, CaseLog> baselines;
    std::map<int64_t, uint64_t> kernel;
    const bool do_paired = args.mode == Mode::PairedValidation ||
                           args.mode == Mode::CalibratePreoffset ||
                           args.mode == Mode::EncryptedNoisyValidation;

    for (uint32_t p : args.p_list) {
        for (const std::string &k_raw : args.k_list) {
            const int64_t kval = parse_k_value(k_raw, p, args.k);
            if (kval <= 0 || kval >= static_cast<int64_t>(p)) continue;
            const Geometry g = make_geometry(p, static_cast<uint32_t>(kval));
            const std::vector<uint64_t> messages = m_list_for_pk(args, g);
            for (const Candidate &cand : candidates) {
                const FinalLut final_lut = make_final_lut(g, cand);
                if (!args.summary_only) {
                    std::cout << "audit candidate=" << cand.name
                              << " p=" << p << " k=" << kval
                              << " period_valid="
                              << final_lut.audit.bitapprox_period_valid
                              << " negacyclic="
                              << final_lut.audit.bitapprox_negacyclic_audit_pass
                              << " semantic_diff="
                              << final_lut.audit
                                     .bitapprox_semantic_diff_vs_exact_count
                              << " reachable=" << final_lut.audit.reachable_count
                              << " ambiguous="
                              << final_lut.audit.ambiguous_cell_count
                              << " min_opp="
                              << final_lut.audit
                                     .min_opposite_label_distance_cells
                              << " c127=" << final_lut.audit.c127
                              << " c128=" << final_lut.audit.c128
                              << " c255=" << final_lut.audit.c255
                              << " c256=" << final_lut.audit.c256
                              << " c383=" << final_lut.audit.c383
                              << " c384=" << final_lut.audit.c384 << "\n";
                }
                std::vector<CaseTask> tasks;
                for (PeriodMode pmode : args.period_modes) {
                    const PeriodInfo period = make_period(pmode, g);
                    for (Policy policy : args.policies) {
                        for (Pipeline pipeline : args.pipelines) {
                            for (BitSource source : args.bit_sources) {
                                for (InputMode input : args.input_modes) {
                                    for (uint64_t m : messages) {
                                        for (uint32_t seed = args.seed_start;
                                             seed < args.seed_start + args.seeds;
                                             seed++) {
                                            CaseTask task;
                                            task.period = period;
                                            task.policy = policy;
                                            task.pipeline = pipeline;
                                            task.source = source;
                                            task.input = input;
                                            task.m = m;
                                            task.seed = seed;
                                            tasks.push_back(task);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (!args.summary_only) {
                    std::cout << "running candidate=" << cand.name
                              << " p=" << p << " k=" << kval
                              << " tasks=" << tasks.size()
                              << " threads="
                              << effective_thread_count(args, tasks.size())
                              << "\n";
                }
                std::vector<CaseLog> logs = run_tasks_parallel(
                    args, g, cand, final_lut, tasks, sk, ek);
                for (const CaseLog &log : logs) {
                    if (args.mode == Mode::CollectKernel &&
                        !log.skipped_unsupported) {
                        kernel[log.selected_minus_nominal]++;
                    }
                    Summary tmp;
                    tmp.pipeline = log.pipeline;
                    tmp.cand = log.cand;
                    tmp.period = log.period;
                    tmp.policy = log.policy;
                    tmp.bit_source = log.bit_source;
                    tmp.input = log.input;
                    tmp.p = log.p;
                    tmp.k = log.k;
                    const std::string key = summary_key(tmp);
                    absorb(summaries[key], log);
                    write_jsonl(jsonl, log);
                    if (do_paired) {
                        const std::string bkey = base_case_key(
                            log.pipeline, log.p, log.k, log.period.mode,
                            log.policy, log.bit_source, log.input, log.m,
                            log.seed);
                        if (is_r0_baseline(log.cand)) {
                            baselines[bkey] = log;
                            const std::string pkey = pair_key(log);
                            absorb_pair(paired[pkey], log, log);
                            write_paired_jsonl(jsonl, log, log);
                        }
                        else {
                            const auto it = baselines.find(bkey);
                            if (it != baselines.end()) {
                                const std::string pkey = pair_key(log);
                                absorb_pair(paired[pkey], it->second, log);
                                write_paired_jsonl(jsonl, it->second, log);
                            }
                        }
                    }
                    if (args.verbose && !args.summary_only) print_case(log);
                }
            }
        }
    }

    write_summary_file(args.summary_path, summaries, kernel,
                       do_paired ? &paired : nullptr);
    return 0;
}

}  // namespace experimental_weighted_bitapprox

int main(int argc, char **argv)
{
    const auto args = experimental_weighted_bitapprox::parse_args(argc, argv);
    return experimental_weighted_bitapprox::run(args);
}
