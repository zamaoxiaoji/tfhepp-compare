/**
 * @file bitextract_qhalf_boundary_probe.cpp
 * @brief Boundary diagnostics for window-local BitExtractPrunedQHalf.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../comparison/bitextract_qhalf.h"

using namespace tfhepp_compare;
namespace be = tfhepp_compare::bitextract_qhalf;

namespace
{
    struct AlphaSpec {
        std::string name = "default";
        double alpha = Lvl1::α;
        bool noiseless = false;
    };

    struct Options {
        std::vector<uint32_t> p_list = {9, 12, 16, 24, 32};
        std::vector<uint32_t> k_list = {1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<std::string> modes = {"margin_theory",
                                          "uncentered_no_prune",
                                          "uncentered_pruned",
                                          "centered_no_prune",
                                          "centered_pruned"};
        std::vector<AlphaSpec> alphas = {AlphaSpec{}};
        std::vector<std::string> iks_sweep = {"default"};
        std::vector<int> boundary_bands = {-1};
        int boundary_trials = 256;
        int random_trials = 1000;
        int threads = 1;
        bool append = false;
        std::string jsonl_path =
            "experimental/results/bitextract_qhalf_boundary_summary.jsonl";
    };

    struct EvalResult {
        std::string mode;
        std::string alpha_name = "default";
        std::string iks_name = "default";
        uint32_t p = 0;
        uint32_t k = 0;
        uint64_t message = 0;
        uint64_t W = 0;
        Lvl1::T delta_p = 0;
        int boundary_band = -1;
        double coverage_percent = 100.0;
        uint32_t expected_bit = 0;
        Lvl1::T expected_phase = 0;
        Lvl1::T actual_phase = 0;
        uint32_t decoded_bit = 0;
        int64_t centered_error = 0;
        uint64_t abs_error = 0;
        be::Trace trace;
        bool pass = false;
        be::FailureClass failure_class = be::FailureClass::None;
        double elapsed_ms = 0.0;
    };

    struct Agg {
        std::string mode;
        std::string alpha_name = "default";
        std::string iks_name = "default";
        uint32_t p = 0;
        uint32_t k = 0;
        int boundary_band = -1;
        double coverage_percent = 100.0;
        uint64_t pass_count = 0;
        uint64_t fail_count = 0;
        std::vector<double> abs_pre_pbs_error_delta_units;
        std::string first_failure;
    };

    std::string json_escape(const std::string &s)
    {
        std::ostringstream out;
        for (char c : s) {
            if (c == '"' || c == '\\') out << '\\';
            else if (c == '\n')
                out << "\\n";
            else
                out << c;
        }
        return out.str();
    }

    std::vector<std::string> split_csv_string(const std::string &s)
    {
        std::vector<std::string> out;
        size_t start = 0;
        while (start <= s.size()) {
            const size_t comma = s.find(',', start);
            std::string item =
                s.substr(start, comma == std::string::npos
                                    ? std::string::npos
                                    : comma - start);
            if (item == "plaintext") item = "plaintext_lut_oracle";
            if (item == "qexact") item = "qexact_br_oracle";
            if (item == "qexact_centered") item = "qexact_centered_oracle";
            if (item == "centered_arbitrary_no_noise")
                item = "arbitrary_no_noise_centered";
            if (item == "encrypted_no_prune") item = "uncentered_no_prune";
            if (item == "encrypted_pruned") item = "uncentered_pruned";
            if (!item.empty()) out.push_back(item);
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        return out;
    }

    std::vector<uint32_t> split_csv_u32(const std::string &s)
    {
        std::vector<uint32_t> out;
        for (const std::string &item : split_csv_string(s))
            out.push_back(static_cast<uint32_t>(std::stoul(item)));
        return out;
    }

    std::vector<int> split_csv_i32(const std::string &s)
    {
        std::vector<int> out;
        for (const std::string &item : split_csv_string(s))
            out.push_back(std::stoi(item));
        return out;
    }

    AlphaSpec parse_alpha_one(const std::string &s)
    {
        if (s == "default") return AlphaSpec{};
        if (s == "zero" || s == "noiseless") return {s, 0.0, true};
        const size_t caret = s.find("^-");
        if (s.rfind("2^-", 0) == 0 && caret != std::string::npos) {
            const int exponent = std::stoi(s.substr(caret + 2));
            return {s, std::pow(2.0, -exponent), false};
        }
        return {s, std::stod(s), false};
    }

    std::vector<AlphaSpec> parse_alpha_sweep(const std::string &s)
    {
        std::vector<AlphaSpec> out;
        for (const std::string &item : split_csv_string(s))
            out.push_back(parse_alpha_one(item));
        return out.empty() ? std::vector<AlphaSpec>{AlphaSpec{}} : out;
    }

    Options parse_options(int argc, char **argv)
    {
        Options opts;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto need_value = [&](const char *name) -> const char * {
                if (i + 1 >= argc) {
                    std::cerr << "missing value for " << name << "\n";
                    std::exit(2);
                }
                return argv[++i];
            };
            if (arg == "--p-list")
                opts.p_list = split_csv_u32(need_value("--p-list"));
            else if (arg == "--k-list")
                opts.k_list = split_csv_u32(need_value("--k-list"));
            else if (arg == "--modes")
                opts.modes = split_csv_string(need_value("--modes"));
            else if (arg == "--alpha-sweep")
                opts.alphas = parse_alpha_sweep(need_value("--alpha-sweep"));
            else if (arg == "--iks-sweep")
                opts.iks_sweep = split_csv_string(need_value("--iks-sweep"));
            else if (arg == "--exclude-bit-boundary-band")
                opts.boundary_bands =
                    split_csv_i32(need_value("--exclude-bit-boundary-band"));
            else if (arg == "--boundary-trials")
                opts.boundary_trials =
                    std::stoi(need_value("--boundary-trials"));
            else if (arg == "--random-trials")
                opts.random_trials = std::stoi(need_value("--random-trials"));
            else if (arg == "--threads")
                opts.threads = std::stoi(need_value("--threads"));
            else if (arg == "--append")
                opts.append = true;
            else if (arg == "--jsonl")
                opts.jsonl_path = need_value("--jsonl");
            else
                std::cerr << "unknown option: " << arg << "\n";
        }
        if (opts.boundary_trials < 0) opts.boundary_trials = 0;
        if (opts.random_trials < 0) opts.random_trials = 0;
        if (opts.threads < 1) opts.threads = 1;
        if (opts.boundary_bands.empty()) opts.boundary_bands = {-1};
        return opts;
    }

    bool wants_mode(const Options &opts, const std::string &mode)
    {
        return std::find(opts.modes.begin(), opts.modes.end(), mode) !=
               opts.modes.end();
    }

    bool is_encrypted_mode(const std::string &mode)
    {
        return mode == "uncentered_no_prune" ||
               mode == "uncentered_pruned" ||
               mode == "centered_no_prune" || mode == "centered_pruned";
    }

    uint64_t message_modulus(uint32_t p)
    {
        return p >= 63 ? std::numeric_limits<uint64_t>::max() : (1ULL << p);
    }

    void add_clamped(std::set<uint64_t> &messages, int64_t v, uint64_t modulus)
    {
        if (v < 0) return;
        const uint64_t u = static_cast<uint64_t>(v);
        if (u < modulus) messages.insert(u);
    }

    void add_window(std::set<uint64_t> &messages, int64_t center,
                    int radius, uint64_t modulus)
    {
        for (int d = -radius; d <= radius; ++d)
            add_clamped(messages, center + d, modulus);
    }

    std::vector<uint64_t> build_messages(uint32_t p, uint32_t k,
                                         int boundary_trials,
                                         int random_trials)
    {
        std::set<uint64_t> messages;
        const uint64_t modulus = message_modulus(p);
        if (p <= 12) {
            for (uint64_t m = 0; m < modulus; ++m) messages.insert(m);
        }

        add_window(messages, static_cast<int64_t>(modulus / 2), 128, modulus);

        const uint64_t W = be::WindowWeight(p, k);
        if (W > 0) {
            const uint64_t transitions = (modulus + W - 1) / W;
            const int radius = std::max(4, std::min(boundary_trials, 16));
            for (uint64_t r = 0; r <= transitions; ++r) {
                if (r > std::numeric_limits<int64_t>::max() / W) break;
                add_window(messages, static_cast<int64_t>(r * W), radius,
                           modulus);
                if (p > 12 && r > 32 && r + 32 < transitions) {
                    r = transitions > 32 ? transitions - 32 : r;
                }
            }
        }

        if (p == 9 && k == 1) messages.insert(255);

        std::mt19937_64 rng((uint64_t(p) << 32) ^ (uint64_t(k) << 16) ^
                            0x9e3779b97f4a7c15ULL);
        std::uniform_int_distribution<uint64_t> dist(0, modulus - 1);
        for (int i = 0; i < random_trials; ++i) messages.insert(dist(rng));
        return {messages.begin(), messages.end()};
    }

    bool include_for_boundary_band(uint64_t message, uint64_t W, int band)
    {
        if (band < 0 || W == 0) return true;
        const uint64_t rem = message % W;
        const uint64_t dist = std::min(rem, W - rem);
        return dist > static_cast<uint64_t>(band);
    }

    std::vector<uint64_t> filter_messages(const std::vector<uint64_t> &messages,
                                          uint64_t W, int band)
    {
        std::vector<uint64_t> out;
        for (uint64_t m : messages)
            if (include_for_boundary_band(m, W, band)) out.push_back(m);
        return out;
    }

    Lvl1::T qhalf()
    {
        return Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 1);
    }

    uint64_t torus_abs_centered(Lvl1::T phase, Lvl1::T expected)
    {
        const Lvl1::T diff = phase - expected;
        if ((diff & qhalf()) != 0) return uint64_t(Lvl1::T(0) - diff);
        return uint64_t(diff);
    }

    int64_t torus_centered_error(Lvl1::T phase, Lvl1::T expected)
    {
        const Lvl1::T diff = phase - expected;
        if ((diff & qhalf()) == 0) return static_cast<int64_t>(diff);
        return -static_cast<int64_t>(Lvl1::T(0) - diff);
    }

    uint32_t decode_qhalf_bit(Lvl1::T phase)
    {
        return torus_abs_centered(phase, qhalf()) <
                       torus_abs_centered(phase, Lvl1::T(0))
                   ? 1U
                   : 0U;
    }

    TLWELvl1 trivial_lvl1(Lvl1::T phase)
    {
        TLWELvl1 ct = {};
        ct[Lvl1::k * Lvl1::n] = phase;
        return ct;
    }

    TLWELvl1 noiseless_encrypt_lvl1(Lvl1::T phase,
                                    const TFHEpp::Key<Lvl1> &key)
    {
        std::uniform_int_distribution<Lvl1::T> dist(
            0, std::numeric_limits<Lvl1::T>::max());
        TLWELvl1 res = {};
        res[Lvl1::k * Lvl1::n] = phase;
        for (int i = 0; i < Lvl1::k * Lvl1::n; ++i) {
            res[i] = dist(TFHEpp::generator);
            res[Lvl1::k * Lvl1::n] += res[i] * key[i];
        }
        return res;
    }

    TLWELvl1 encrypt_lvl1(Lvl1::T phase, const AlphaSpec &alpha,
                          const TFHEpp::Key<Lvl1> &key)
    {
        if (alpha.noiseless) return noiseless_encrypt_lvl1(phase, key);
        return TFHEpp::tlweSymEncrypt<Lvl1>(phase, alpha.alpha, key);
    }

    Lvl1::T decrypt_phase_lvl1(const TLWELvl1 &ct,
                               const TFHEpp::Key<Lvl1> &key)
    {
        return TFHEpp::tlweSymPhase<Lvl1>(ct, key);
    }

    Lvl0::T decrypt_phase_lvl0(const TLWELvl0 &ct,
                               const TFHEpp::Key<Lvl0> &key)
    {
        return TFHEpp::tlweSymPhase<Lvl0>(ct, key);
    }

    Lvl1::T plaintext_lut_oracle(uint64_t message, uint32_t p, uint32_t k,
                                 bool centered, be::Trace *trace)
    {
        return centered
                   ? be::QExactOraclePhase_Lvl1_CenteredCell(message, p, k,
                                                             trace)
                   : be::QExactOraclePhase_Lvl1(message, p, k, trace);
    }

    void fill_pre_pbs_trace(be::Trace &trace, const TLWELvl1 &input,
                            uint64_t message, uint32_t p, uint32_t k,
                            bool centered, const TFHESecretKey &sk,
                            const TFHEEvalKey &ek)
    {
        const Lvl1::T center = centered ? be::CenterOffset(p) : Lvl1::T(0);
        const Lvl1::T expected_center =
            static_cast<Lvl1::T>(message) * be::PhaseDelta(p) + center;
        TLWELvl1 centered_input = input;
        centered_input[Lvl1::k * Lvl1::n] += center;

        TLWELvl0 after_iks;
        TFHEpp::IdentityKeySwitch<Lvl10>(after_iks, centered_input,
                                         *ek.iksklvl10);

        const Lvl1::T phase_before =
            decrypt_phase_lvl1(input, sk.key.lvl1);
        const Lvl1::T phase_centered =
            decrypt_phase_lvl1(centered_input, sk.key.lvl1);
        const Lvl0::T phase_after_iks =
            decrypt_phase_lvl0(after_iks, sk.key.lvl0);

        trace.bitextract_centering_enabled = centered;
        trace.center_offset = center;
        trace.original_phase = phase_before;
        trace.centered_phase_for_bitextract = phase_centered;
        trace.expected_nearest_message = message;
        trace.distance_to_nearest_bit_boundary_in_delta_units =
            be::DistanceToNearestBitBoundaryDeltaUnits(message, p, k, centered);
        trace.distance_to_nearest_bit_boundary_in_torus =
            be::DistanceToNearestBitBoundaryTorus(message, p, k, centered);
        trace.phase_before_iks = phase_before;
        trace.phase_after_centering = phase_centered;
        trace.phase_after_iks_before_offset = phase_after_iks;
        trace.phase_after_centering_and_iks = phase_after_iks;
        trace.signed_error_to_expected_center =
            torus_centered_error(static_cast<Lvl1::T>(phase_after_iks),
                                 expected_center);
        trace.signed_error_in_delta_units =
            static_cast<double>(trace.signed_error_to_expected_center) /
            static_cast<double>(be::PhaseDelta(p));
    }

    EvalResult make_result(const std::string &mode, const AlphaSpec &alpha,
                           uint32_t p, uint32_t k, uint64_t message,
                           int boundary_band, double coverage_percent,
                           Lvl1::T actual, const be::Trace &trace,
                           double elapsed_ms,
                           be::FailureClass failure_class_hint)
    {
        EvalResult r;
        r.mode = mode;
        r.alpha_name = alpha.name;
        r.p = p;
        r.k = k;
        r.message = message;
        r.W = be::WindowWeight(p, k);
        r.delta_p = be::PhaseDelta(p);
        r.boundary_band = boundary_band;
        r.coverage_percent = coverage_percent;
        r.expected_phase = be::ExpectedQHalfBit(message, p, k);
        r.expected_bit = r.expected_phase == Lvl1::T(0) ? 0U : 1U;
        r.actual_phase = actual;
        r.decoded_bit = decode_qhalf_bit(actual);
        r.centered_error = torus_centered_error(actual, r.expected_phase);
        r.abs_error = torus_abs_centered(actual, r.expected_phase);
        r.trace = trace;
        r.pass = r.decoded_bit == r.expected_bit && r.abs_error < (qhalf() >> 1);
        r.failure_class =
            r.pass ? be::FailureClass::None : failure_class_hint;
        r.elapsed_ms = elapsed_ms;
        return r;
    }

    double percentile(std::vector<double> v, double q)
    {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        const size_t idx =
            std::min(v.size() - 1,
                     static_cast<size_t>(std::ceil(q * v.size())) - 1);
        return v[idx];
    }

    std::string first_failure_string(const EvalResult &r)
    {
        std::ostringstream out;
        out << "p=" << r.p << ",k=" << r.k << ",m=" << r.message
            << ",expected_bit=" << r.expected_bit
            << ",decoded_bit=" << r.decoded_bit
            << ",pre_pbs_error_delta="
            << r.trace.signed_error_in_delta_units
            << ",failure_class=" << be::FailureClassName(r.failure_class);
        return out.str();
    }

    void write_case_jsonl(std::ofstream &os, const EvalResult &r)
    {
        os << "{"
           << "\"kind\":\"bitextract_qhalf_boundary\","
           << "\"mode\":\"" << json_escape(r.mode) << "\","
           << "\"conversion_mode\":\"bitextract_qhalf\","
           << "\"alpha_name\":\"" << json_escape(r.alpha_name) << "\","
           << "\"iks_name\":\"" << json_escape(r.iks_name) << "\","
           << "\"p\":" << r.p << ","
           << "\"k\":" << r.k << ","
           << "\"message\":" << r.message << ","
           << "\"W\":" << r.W << ","
           << "\"delta_p\":" << r.delta_p << ","
           << "\"boundary_band\":" << r.boundary_band << ","
           << "\"coverage_percent\":" << std::setprecision(9)
           << r.coverage_percent << ","
           << "\"expected_bit\":" << r.expected_bit << ","
           << "\"expected_phase\":" << r.expected_phase << ","
           << "\"actual_phase\":" << r.actual_phase << ","
           << "\"decoded_bit\":" << r.decoded_bit << ","
           << "\"centered_error\":" << r.centered_error << ","
           << "\"abs_error\":" << r.abs_error << ","
           << "\"bitextract_centering_enabled\":"
           << (r.trace.bitextract_centering_enabled ? "true" : "false") << ","
           << "\"center_offset\":" << r.trace.center_offset << ","
           << "\"original_phase\":" << r.trace.original_phase << ","
           << "\"centered_phase_for_bitextract\":"
           << r.trace.centered_phase_for_bitextract << ","
           << "\"expected_nearest_message\":"
           << r.trace.expected_nearest_message << ","
           << "\"distance_to_nearest_bit_boundary_in_torus\":"
           << r.trace.distance_to_nearest_bit_boundary_in_torus << ","
           << "\"distance_to_nearest_bit_boundary_in_delta_units\":"
           << r.trace.distance_to_nearest_bit_boundary_in_delta_units << ","
           << "\"phase_before_iks\":" << r.trace.phase_before_iks << ","
           << "\"phase_after_iks_before_offset\":"
           << r.trace.phase_after_iks_before_offset << ","
           << "\"phase_after_centering\":"
           << r.trace.phase_after_centering << ","
           << "\"phase_after_centering_and_iks\":"
           << r.trace.phase_after_centering_and_iks << ","
           << "\"signed_error_to_expected_center\":"
           << r.trace.signed_error_to_expected_center << ","
           << "\"signed_error_in_delta_units\":"
           << r.trace.signed_error_in_delta_units << ","
           << "\"br_index\":" << r.trace.br_index << ","
           << "\"modswitch_index\":" << r.trace.modswitch_index << ","
           << "\"lut_index\":" << r.trace.lut_index << ","
           << "\"skip_count\":" << r.trace.skip_count << ","
           << "\"cmux_count\":" << r.trace.cmux_count << ","
           << "\"period\":" << r.trace.period << ","
           << "\"offset_used\":" << r.trace.offset_used << ","
           << "\"representative_mode\":\""
           << json_escape(r.trace.representative_mode) << "\","
           << "\"elapsed_ms\":" << r.elapsed_ms << ","
           << "\"status\":\"" << (r.pass ? "pass" : "fail") << "\","
           << "\"failure_class\":\""
           << be::FailureClassName(r.failure_class) << "\""
           << "}\n";
    }

    void write_summary_jsonl(std::ofstream &os, const Agg &a)
    {
        const double p50 = percentile(a.abs_pre_pbs_error_delta_units, 0.50);
        const double p90 = percentile(a.abs_pre_pbs_error_delta_units, 0.90);
        const double p99 = percentile(a.abs_pre_pbs_error_delta_units, 0.99);
        const double mx =
            a.abs_pre_pbs_error_delta_units.empty()
                ? 0.0
                : *std::max_element(a.abs_pre_pbs_error_delta_units.begin(),
                                    a.abs_pre_pbs_error_delta_units.end());
        const char *classification =
            p99 > 0.5 ? "BITEXTRACT_NOISE_EXCEEDS_CELL_MARGIN"
                      : (a.fail_count == 0 ? "PASS" : "ENCRYPTED_NOISE_BUG");
        os << "{"
           << "\"kind\":\"bitextract_qhalf_boundary_summary\","
           << "\"mode\":\"" << json_escape(a.mode) << "\","
           << "\"alpha_name\":\"" << json_escape(a.alpha_name) << "\","
           << "\"iks_name\":\"" << json_escape(a.iks_name) << "\","
           << "\"p\":" << a.p << ","
           << "\"k\":" << a.k << ","
           << "\"boundary_band\":" << a.boundary_band << ","
           << "\"coverage_percent\":" << std::setprecision(9)
           << a.coverage_percent << ","
           << "\"pass_count\":" << a.pass_count << ","
           << "\"fail_count\":" << a.fail_count << ","
           << "\"p50_abs_pre_pbs_error_delta_units\":" << p50 << ","
           << "\"p90_abs_pre_pbs_error_delta_units\":" << p90 << ","
           << "\"p99_abs_pre_pbs_error_delta_units\":" << p99 << ","
           << "\"max_abs_pre_pbs_error_delta_units\":" << mx << ","
           << "\"margin_delta_units\":0.5,"
           << "\"first_failure\":\"" << json_escape(a.first_failure) << "\","
           << "\"classification\":\"" << classification << "\""
           << "}\n";
    }

    void write_margin_theory(std::ofstream &os, uint32_t p, uint32_t k,
                             const std::vector<uint64_t> &messages,
                             int boundary_band, double coverage_percent)
    {
        double min_uncentered = std::numeric_limits<double>::infinity();
        double min_centered = std::numeric_limits<double>::infinity();
        for (uint64_t m : messages) {
            min_uncentered = std::min(
                min_uncentered,
                be::DistanceToNearestBitBoundaryDeltaUnits(m, p, k, false));
            min_centered = std::min(
                min_centered,
                be::DistanceToNearestBitBoundaryDeltaUnits(m, p, k, true));
        }
        if (messages.empty()) {
            min_uncentered = 0.0;
            min_centered = 0.0;
        }
        os << "{"
           << "\"kind\":\"margin_theory\","
           << "\"mode\":\"margin_theory\","
           << "\"p\":" << p << ","
           << "\"k\":" << k << ","
           << "\"W\":" << be::WindowWeight(p, k) << ","
           << "\"delta_p\":" << be::PhaseDelta(p) << ","
           << "\"boundary_band\":" << boundary_band << ","
           << "\"coverage_percent\":" << std::setprecision(9)
           << coverage_percent << ","
           << "\"min_uncentered_margin\":" << min_uncentered << ","
           << "\"min_centered_margin\":" << min_centered << ","
           << "\"min_uncentered_margin_torus\":"
           << min_uncentered * static_cast<double>(be::PhaseDelta(p)) << ","
           << "\"min_centered_margin_torus\":"
           << min_centered * static_cast<double>(be::PhaseDelta(p)) << ","
           << "\"classification\":\"MARGIN_THEORY\""
           << "}\n";
    }

    void print_failure_trace(const EvalResult &r)
    {
        std::cout << "FIRST_FAILURE\n"
                  << "  mode=" << r.mode << "\n"
                  << "  p=" << r.p << " k=" << r.k
                  << " m=" << r.message << "\n"
                  << "  centered="
                  << (r.trace.bitextract_centering_enabled ? "true" : "false")
                  << " center_offset=" << r.trace.center_offset << "\n"
                  << "  Delta_p=" << r.delta_p
                  << " distance_delta="
                  << r.trace.distance_to_nearest_bit_boundary_in_delta_units
                  << " distance_torus="
                  << r.trace.distance_to_nearest_bit_boundary_in_torus << "\n"
                  << "  phase_before_iks=" << r.trace.phase_before_iks
                  << " phase_after_iks_before_offset="
                  << r.trace.phase_after_iks_before_offset << "\n"
                  << "  pre_pbs_error_delta="
                  << r.trace.signed_error_in_delta_units << "\n"
                  << "  expected_bit=" << r.expected_bit
                  << " decoded_bit_qhalf=" << r.decoded_bit
                  << " actual_phase=" << r.actual_phase << "\n"
                  << "  no-prune/pruned failure_class="
                  << be::FailureClassName(r.failure_class) << "\n";
    }
} // namespace

int main(int argc, char **argv)
{
    const Options opts = parse_options(argc, argv);
#ifdef _OPENMP
    omp_set_num_threads(opts.threads);
#endif

    bool need_eval_key = wants_mode(opts, "arbitrary_no_noise") ||
                         wants_mode(opts, "arbitrary_no_noise_centered");
    for (const std::string &mode : opts.modes)
        need_eval_key = need_eval_key || is_encrypted_mode(mode);

    std::filesystem::create_directories(
        std::filesystem::path(opts.jsonl_path).parent_path());
    std::ofstream jsonl(opts.jsonl_path,
                        opts.append ? std::ios::app : std::ios::trunc);
    if (!jsonl) {
        std::cerr << "failed to open jsonl path: " << opts.jsonl_path << "\n";
        return 2;
    }

    std::cout << "BitExtractPrunedQHalf boundary probe\n";
    std::cout << "  UNSAFE_PERFORMANCE_ONLY=true\n";
    std::cout << "  bit_index_convention=MSB_FIRST_WINDOW_LOCAL\n";
    std::cout << "  bit_k(m)=floor(m / W) mod 2\n";
    std::cout << "  W=2^(p-1-k)\n";
    std::cout << "  centered_cell_formula=ct_for_bitextract + Delta_p/2\n";
    std::cout << "  modes=";
    for (size_t i = 0; i < opts.modes.size(); ++i)
        std::cout << (i ? "," : "") << opts.modes[i];
    std::cout << "\n";
    for (const std::string &iks : opts.iks_sweep) {
        if (iks != "default")
            std::cout << "  iks_variant=" << iks
                      << " status=unsupported_in_this_local_probe\n";
    }

    std::optional<TFHESecretKey> sk;
    std::optional<TFHEEvalKey> ek;
    if (need_eval_key) {
        std::cout << "  keygen=Lvl01 bkfft + default Lvl10 iksk\n";
        sk.emplace();
        ek.emplace();
        ek->emplacebkfft<Lvl01>(*sk);
        ek->emplaceiksk<Lvl10>(*sk);
    }

    std::map<std::tuple<std::string, std::string, uint32_t, uint32_t, int>, Agg>
        aggs;
    std::optional<EvalResult> first_failure;

    for (const uint32_t p : opts.p_list) {
        for (const uint32_t k : opts.k_list) {
            if (!be::IsValidWindowLocalK(k) || k >= p) {
                std::cout << "skip p=" << p << " k=" << k
                          << " reason=SCALE_OR_INDEXING_BUG\n";
                continue;
            }

            const uint64_t W = be::WindowWeight(p, k);
            const auto base_messages =
                build_messages(p, k, opts.boundary_trials, opts.random_trials);

            for (const int band : opts.boundary_bands) {
                const auto messages = filter_messages(base_messages, W, band);
                const double coverage =
                    base_messages.empty()
                        ? 0.0
                        : 100.0 * static_cast<double>(messages.size()) /
                              static_cast<double>(base_messages.size());
                std::cout << "  p=" << p << " k=" << k << " W=" << W
                          << " period=" << be::WindowLocalPeriodLvl1(k)
                          << " band=" << band
                          << " messages=" << messages.size()
                          << " coverage=" << std::fixed << std::setprecision(2)
                          << coverage << "%\n";

                if (wants_mode(opts, "margin_theory"))
                    write_margin_theory(jsonl, p, k, messages, band, coverage);

                for (const AlphaSpec &alpha : opts.alphas) {
                    const bool need_encrypted =
                        wants_mode(opts, "uncentered_no_prune") ||
                        wants_mode(opts, "uncentered_pruned") ||
                        wants_mode(opts, "centered_no_prune") ||
                        wants_mode(opts, "centered_pruned");

                    for (const uint64_t message : messages) {
                        auto handle = [&](const EvalResult &r) {
                            write_case_jsonl(jsonl, r);
                            auto key = std::make_tuple(r.mode, r.alpha_name,
                                                       r.p, r.k,
                                                       r.boundary_band);
                            Agg &a = aggs[key];
                            a.mode = r.mode;
                            a.alpha_name = r.alpha_name;
                            a.p = r.p;
                            a.k = r.k;
                            a.boundary_band = r.boundary_band;
                            a.coverage_percent = r.coverage_percent;
                            if (r.pass)
                                a.pass_count++;
                            else {
                                a.fail_count++;
                                if (a.first_failure.empty())
                                    a.first_failure = first_failure_string(r);
                                if (!first_failure) first_failure = r;
                            }
                            a.abs_pre_pbs_error_delta_units.push_back(
                                std::abs(r.trace.signed_error_in_delta_units));
                        };

                        if (wants_mode(opts, "plaintext_lut_oracle")) {
                            be::Trace trace;
                            const Lvl1::T phase = plaintext_lut_oracle(
                                message, p, k, false, &trace);
                            EvalResult r = make_result(
                                "plaintext_lut_oracle", alpha, p, k, message,
                                band, coverage, phase, trace, 0.0,
                                be::FailureClass::PlaintextLutBug);
                            handle(r);
                        }

                        if (wants_mode(opts, "qexact_br_oracle")) {
                            be::Trace trace;
                            const Lvl1::T phase =
                                be::QExactOraclePhase_Lvl1(message, p, k,
                                                           &trace);
                            EvalResult r = make_result(
                                "qexact_br_oracle", alpha, p, k, message, band,
                                coverage, phase, trace, 0.0,
                                be::FailureClass::ModswitchRepresentativeBug);
                            handle(r);
                        }

                        if (wants_mode(opts, "qexact_centered_oracle")) {
                            be::Trace trace;
                            const Lvl1::T phase =
                                be::QExactOraclePhase_Lvl1_CenteredCell(
                                    message, p, k, &trace);
                            EvalResult r = make_result(
                                "qexact_centered_oracle", alpha, p, k,
                                message, band, coverage, phase, trace, 0.0,
                                be::FailureClass::ModswitchRepresentativeBug);
                            handle(r);
                        }

                        if (wants_mode(opts, "arbitrary_no_noise")) {
                            be::Trace trace;
                            TLWELvl1 out;
                            TLWELvl1 in = trivial_lvl1(
                                static_cast<Lvl1::T>(message) *
                                be::PhaseDelta(p));
                            fill_pre_pbs_trace(trace, in, message, p, k, false,
                                               *sk, *ek);
                            const auto t0 = std::chrono::steady_clock::now();
                            be::BitExtractQHalf_NoPrune_Lvl1(out, in, k, *ek,
                                                             &trace);
                            const auto t1 = std::chrono::steady_clock::now();
                            EvalResult r = make_result(
                                "arbitrary_no_noise", alpha, p, k, message,
                                band, coverage,
                                decrypt_phase_lvl1(out, sk->key.lvl1), trace,
                                std::chrono::duration<double, std::milli>(t1 -
                                                                          t0)
                                    .count(),
                                be::FailureClass::CenterOffsetBug);
                            handle(r);
                        }

                        if (wants_mode(opts, "arbitrary_no_noise_centered")) {
                            be::Trace trace;
                            TLWELvl1 out;
                            TLWELvl1 in = trivial_lvl1(
                                static_cast<Lvl1::T>(message) *
                                be::PhaseDelta(p));
                            fill_pre_pbs_trace(trace, in, message, p, k, true,
                                               *sk, *ek);
                            const auto t0 = std::chrono::steady_clock::now();
                            be::BitExtractQHalf_NoPrune_Lvl1_CenteredCell(
                                out, in, p, k, *ek, &trace);
                            const auto t1 = std::chrono::steady_clock::now();
                            EvalResult r = make_result(
                                "arbitrary_no_noise_centered", alpha, p, k,
                                message, band, coverage,
                                decrypt_phase_lvl1(out, sk->key.lvl1), trace,
                                std::chrono::duration<double, std::milli>(t1 -
                                                                          t0)
                                    .count(),
                                be::FailureClass::CenterOffsetBug);
                            handle(r);
                        }

                        if (!need_encrypted) continue;

                        TLWELvl1 encrypted_input = encrypt_lvl1(
                            static_cast<Lvl1::T>(message) * be::PhaseDelta(p),
                            alpha, sk->key.lvl1);

                        std::optional<EvalResult> uncentered_ref;
                        std::optional<EvalResult> centered_ref;

                        auto run_encrypted =
                            [&](const std::string &mode, bool centered,
                                bool pruned,
                                be::FailureClass failure_hint) -> EvalResult {
                            be::Trace trace;
                            fill_pre_pbs_trace(trace, encrypted_input, message,
                                               p, k, centered, *sk, *ek);
                            TLWELvl1 out;
                            const auto t0 = std::chrono::steady_clock::now();
                            if (centered && pruned)
                                be::BitExtractQHalf_Pruned_Lvl1_CenteredCell(
                                    out, encrypted_input, p, k, *ek, &trace);
                            else if (centered)
                                be::BitExtractQHalf_NoPrune_Lvl1_CenteredCell(
                                    out, encrypted_input, p, k, *ek, &trace);
                            else if (pruned)
                                be::BitExtractQHalf_Pruned_Lvl1(
                                    out, encrypted_input, k, *ek, &trace);
                            else
                                be::BitExtractQHalf_NoPrune_Lvl1(
                                    out, encrypted_input, k, *ek, &trace);
                            const auto t1 = std::chrono::steady_clock::now();
                            EvalResult r = make_result(
                                mode, alpha, p, k, message, band, coverage,
                                decrypt_phase_lvl1(out, sk->key.lvl1), trace,
                                std::chrono::duration<double, std::milli>(t1 -
                                                                          t0)
                                    .count(),
                                failure_hint);
                            if (!r.pass &&
                                std::abs(trace.signed_error_in_delta_units) >
                                    0.5)
                                r.failure_class =
                                    be::FailureClass::
                                        BitextractNoiseExceedsCellMargin;
                            return r;
                        };

                        if (wants_mode(opts, "uncentered_no_prune")) {
                            uncentered_ref = run_encrypted(
                                "uncentered_no_prune", false, false,
                                be::FailureClass::EncryptedNoiseBug);
                            handle(*uncentered_ref);
                        }
                        if (wants_mode(opts, "uncentered_pruned")) {
                            EvalResult r = run_encrypted(
                                "uncentered_pruned", false, true,
                                be::FailureClass::EncryptedNoiseBug);
                            if (uncentered_ref && uncentered_ref->pass &&
                                !r.pass)
                                r.failure_class =
                                    be::FailureClass::PeriodicSkipBug;
                            handle(r);
                        }
                        if (wants_mode(opts, "centered_no_prune")) {
                            centered_ref = run_encrypted(
                                "centered_no_prune", true, false,
                                be::FailureClass::EncryptedNoiseBug);
                            handle(*centered_ref);
                        }
                        if (wants_mode(opts, "centered_pruned")) {
                            EvalResult r = run_encrypted(
                                "centered_pruned", true, true,
                                be::FailureClass::EncryptedNoiseBug);
                            if (centered_ref && centered_ref->pass && !r.pass)
                                r.failure_class =
                                    be::FailureClass::PeriodicSkipBug;
                            handle(r);
                        }
                    }
                }
            }
        }
    }

    for (const auto &[_, agg] : aggs) write_summary_jsonl(jsonl, agg);

    std::cout << "\nSummary:\n";
    for (const auto &[_, agg] : aggs) {
        const double p99 = percentile(agg.abs_pre_pbs_error_delta_units, 0.99);
        std::cout << "  " << agg.mode << " alpha=" << agg.alpha_name
                  << " p=" << agg.p << " k=" << agg.k
                  << " band=" << agg.boundary_band
                  << " pass=" << agg.pass_count
                  << " fail=" << agg.fail_count
                  << " p99_pre_pbs_delta=" << std::setprecision(3) << p99
                  << "\n";
    }
    if (first_failure)
        print_failure_trace(*first_failure);
    else
        std::cout << "FIRST_FAILURE none\n";
    std::cout << "jsonl=" << opts.jsonl_path << "\n";
    return first_failure ? 1 : 0;
}
