/**
 * @file micro_pbs_full_pipeline_test.cpp
 * @brief Full 3-PBS pipeline probe with locked unsafe Micro-PBS conversion.
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../comparison/bitextract_qhalf.h"
#include "../comparison/micro_evalkey.h"
#include "../comparison/micro_pbs.h"

using namespace tfhepp_compare;
namespace be = tfhepp_compare::bitextract_qhalf;
namespace mp = tfhepp_compare::micro_pbs;

namespace
{
    using BestMicroCandidate = mp::micro_n64_N1024_l2_b8;

    struct Options {
        std::vector<std::string> conversions = {"boolarith_existing",
                                                "micro_lut_lvl01_compat",
                                                "micro_unsafe_best"};
        std::vector<uint32_t> p_list = {9, 12, 16, 24, 32};
        std::vector<uint32_t> k_list = {5, 6, 7, 8};
        int boundary_trials = 256;
        int random_trials = 5000;
        int threads = 1;
        bool append = false;
        bool centered_bitextract = false;
        std::string jsonl_path =
            "experimental/results/micro_pbs_full_pipeline_summary.jsonl";
    };

    struct Row {
        std::string conversion;
        bool centered_bitextract = false;
        uint32_t p = 0;
        uint32_t k = 0;
        uint64_t message = 0;
        uint64_t W = 0;
        Lvl1::T delta_p = 0;
        Lvl1::T guard_value = 0;
        Lvl1::T A = 0;
        Lvl1::T gap_offset = 0;
        uint32_t expected_bit = 0;
        uint32_t decoded_bit_qhalf = 0;
        Lvl1::T input_phase = 0;
        Lvl1::T bit_qhalf_phase_raw = 0;
        Lvl1::T expected_guard = 0;
        Lvl1::T decoded_guard = 0;
        int64_t guard_error = 0;
        Lvl1::T gap_raw_phase = 0;
        Lvl1::T expected_gap = 0;
        int64_t gap_error = 0;
        Lvl1::T final_output_raw_phase = 0;
        uint32_t final_output = 0;
        uint32_t expected_output = 0;
        double t_encrypt_ms = 0.0;
        double t_bitextract_ms = 0.0;
        double t_micro_pre_ks_ms = 0.0;
        double t_micro_pbs_ms = 0.0;
        double t_micro_post_ks_ms = 0.0;
        double t_conversion_total_ms = 0.0;
        double t_clear_sub_ms = 0.0;
        double t_final_gap_ms = 0.0;
        double t_total_ms = 0.0;
        std::string status = "pass";
        std::string failure_class = "NONE";
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
            const std::string item =
                s.substr(start, comma == std::string::npos
                                    ? std::string::npos
                                    : comma - start);
            if (!item.empty()) out.push_back(item);
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        return out;
    }

    std::vector<uint32_t> split_csv_u32(const std::string &s)
    {
        std::vector<uint32_t> out;
        size_t start = 0;
        while (start <= s.size()) {
            const size_t comma = s.find(',', start);
            const std::string item =
                s.substr(start, comma == std::string::npos
                                    ? std::string::npos
                                    : comma - start);
            if (!item.empty()) out.push_back(static_cast<uint32_t>(std::stoul(item)));
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        return out;
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
            if (arg == "--conversion")
                opts.conversions = split_csv_string(need_value("--conversion"));
            else if (arg == "--p-list")
                opts.p_list = split_csv_u32(need_value("--p-list"));
            else if (arg == "--k-list")
                opts.k_list = split_csv_u32(need_value("--k-list"));
            else if (arg == "--boundary-trials")
                opts.boundary_trials =
                    std::stoi(need_value("--boundary-trials"));
            else if (arg == "--random-trials")
                opts.random_trials = std::stoi(need_value("--random-trials"));
            else if (arg == "--threads")
                opts.threads = std::stoi(need_value("--threads"));
            else if (arg == "--append")
                opts.append = true;
            else if (arg == "--centered-bitextract") {
                const std::string value = need_value("--centered-bitextract");
                opts.centered_bitextract =
                    value == "true" || value == "1" || value == "yes";
            }
            else if (arg == "--jsonl")
                opts.jsonl_path = need_value("--jsonl");
            else
                std::cerr << "unknown option: " << arg << "\n";
        }
        if (opts.boundary_trials < 0) opts.boundary_trials = 0;
        if (opts.random_trials < 0) opts.random_trials = 0;
        if (opts.threads < 1) opts.threads = 1;
        return opts;
    }

    bool wants_conversion(const Options &opts, const std::string &name)
    {
        return std::find(opts.conversions.begin(), opts.conversions.end(),
                         name) != opts.conversions.end();
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
        add_window(messages, static_cast<int64_t>(modulus / 2),
                   std::min(boundary_trials, 128), modulus);

        const uint64_t W = be::WindowWeight(p, k);
        if (W > 0) {
            const uint64_t transitions = (modulus + W - 1) / W;
            for (uint64_t r = 0; r <= transitions; ++r) {
                if (r > std::numeric_limits<int64_t>::max() / W) break;
                add_window(messages, static_cast<int64_t>(r * W),
                           std::min(boundary_trials, 16), modulus);
            }
        }

        std::mt19937_64 rng((uint64_t(p) << 32) ^ (uint64_t(k) << 16) ^
                            0x243f6a8885a308d3ULL);
        std::uniform_int_distribution<uint64_t> dist(0, modulus - 1);
        for (int i = 0; i < random_trials; ++i) messages.insert(dist(rng));
        return {messages.begin(), messages.end()};
    }

    Lvl1::T phase_delta(uint32_t p)
    {
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        return p >= digits ? Lvl1::T(1) : (Lvl1::T(1) << (digits - p));
    }

    Lvl1::T encode_message(uint64_t message, uint32_t p)
    {
        return static_cast<Lvl1::T>(message) * phase_delta(p);
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

    Lvl1::T decrypt_phase(const TLWELvl1 &ct, const TFHEpp::Key<Lvl1> &key)
    {
        return TFHEpp::tlweSymPhase<Lvl1>(ct, key);
    }

    uint32_t expected_final(uint64_t message, uint32_t p, uint32_t k)
    {
        const uint64_t modulus = message_modulus(p);
        const uint64_t W = be::WindowWeight(p, k);
        const uint64_t bit = W == 0 ? 0 : ((message / W) & 1ULL);
        const uint64_t gap = message - bit * W;
        const uint64_t shifted = (gap + ((W + 1) / 2)) % modulus;
        return shifted >= (modulus / 2) ? 1U : 0U;
    }

    void write_jsonl(std::ofstream &os, const Row &r)
    {
        os << "{"
           << "\"kind\":\"micro_pbs_full_pipeline\","
           << "\"conversion_mode\":\"" << json_escape(r.conversion) << "\","
           << "\"centered_bitextract\":"
           << (r.centered_bitextract ? "true" : "false") << ","
           << "\"param_name\":\""
           << (r.conversion == "micro_unsafe_best" ? BestMicroCandidate::name
                                                    : r.conversion)
           << "\","
           << "\"p\":" << r.p << ","
           << "\"k\":" << r.k << ","
           << "\"message\":" << r.message << ","
           << "\"W\":" << r.W << ","
           << "\"delta_p\":" << r.delta_p << ","
           << "\"guard_value\":" << r.guard_value << ","
           << "\"A\":" << r.A << ","
           << "\"gap_offset\":" << r.gap_offset << ","
           << "\"expected_bit_k\":" << r.expected_bit << ","
           << "\"decoded_bit_qhalf\":" << r.decoded_bit_qhalf << ","
           << "\"input_phase\":" << r.input_phase << ","
           << "\"bit_qhalf_phase_raw\":" << r.bit_qhalf_phase_raw << ","
           << "\"expected_guard\":" << r.expected_guard << ","
           << "\"decoded_guard\":" << r.decoded_guard << ","
           << "\"guard_error\":" << r.guard_error << ","
           << "\"gap_raw_phase\":" << r.gap_raw_phase << ","
           << "\"expected_gap\":" << r.expected_gap << ","
           << "\"gap_error\":" << r.gap_error << ","
           << "\"final_output_raw_phase\":" << r.final_output_raw_phase << ","
           << "\"final_output\":" << r.final_output << ","
           << "\"expected_output\":" << r.expected_output << ","
           << "\"t_encrypt_ms\":" << std::setprecision(9) << r.t_encrypt_ms
           << ","
           << "\"t_bitextract_ms\":" << r.t_bitextract_ms << ","
           << "\"t_micro_pre_ks_ms\":" << r.t_micro_pre_ks_ms << ","
           << "\"t_micro_pbs_ms\":" << r.t_micro_pbs_ms << ","
           << "\"t_micro_post_ks_ms\":" << r.t_micro_post_ks_ms << ","
           << "\"t_conversion_total_ms\":" << r.t_conversion_total_ms << ","
           << "\"t_clear_sub_ms\":" << r.t_clear_sub_ms << ","
           << "\"t_final_gap_ms\":" << r.t_final_gap_ms << ","
           << "\"t_total_ms\":" << r.t_total_ms << ","
           << "\"status\":\"" << r.status << "\","
           << "\"failure_class\":\"" << r.failure_class << "\""
           << "}\n";
    }

    void print_failure_trace(const Row &r)
    {
        std::cout << "FIRST_FAILURE\n"
                  << "  conversion=" << r.conversion << "\n"
                  << "  p=" << r.p << " k=" << r.k
                  << " m=" << r.message << "\n"
                  << "  W=" << r.W << " delta_p=" << r.delta_p
                  << " guard_value=" << r.guard_value << " A=" << r.A
                  << " gap_offset=" << r.gap_offset << "\n"
                  << "  expected_bit_k=" << r.expected_bit
                  << " decoded_bit_qhalf=" << r.decoded_bit_qhalf
                  << " bit_qhalf_phase_raw=" << r.bit_qhalf_phase_raw
                  << "\n"
                  << "  expected_guard=" << r.expected_guard
                  << " decoded_guard=" << r.decoded_guard
                  << " guard_error=" << r.guard_error << "\n"
                  << "  decoded_gap=" << r.gap_raw_phase
                  << " expected_gap=" << r.expected_gap
                  << " gap_error=" << r.gap_error << "\n"
                  << "  final_output=" << r.final_output
                  << " expected_output=" << r.expected_output
                  << " final_output_raw_phase="
                  << r.final_output_raw_phase << "\n"
                  << "  input_phase=" << r.input_phase
                  << " guard_raw_phase=" << r.decoded_guard
                  << " gap_raw_phase=" << r.gap_raw_phase << "\n"
                  << "  failure_class=" << r.failure_class << "\n";
    }

    template <class Pack>
    Row run_case(const std::string &conversion, uint32_t p, uint32_t k,
                 uint64_t message, const TFHESecretKey &sk,
                 const TFHEEvalKey &ek, const Pack *micro_pack,
                 bool centered_bitextract)
    {
        Row r;
        r.conversion = conversion;
        r.centered_bitextract = centered_bitextract;
        r.p = p;
        r.k = k;
        r.message = message;
        r.W = be::WindowWeight(p, k);
        r.delta_p = phase_delta(p);
        r.guard_value = Lvl1::T(1)
                        << (std::numeric_limits<Lvl1::T>::digits - 1 - k);
        r.A = r.guard_value >> 1;
        r.gap_offset = static_cast<Lvl1::T>(((r.W + 1) * uint64_t(r.delta_p)) /
                                            2);
        r.expected_bit = r.W == 0 ? 0 : ((message / r.W) & 1ULL);
        r.expected_guard = r.expected_bit ? r.guard_value : Lvl1::T(0);
        const uint64_t gap_message = message - uint64_t(r.expected_bit) * r.W;
        r.expected_gap = encode_message(gap_message, p);
        r.expected_output = expected_final(message, p, k);

        const auto total0 = std::chrono::steady_clock::now();
        const auto e0 = std::chrono::steady_clock::now();
        TLWELvl1 ct =
            TFHEpp::tlweSymEncrypt<Lvl1>(encode_message(message, p), Lvl1::α,
                                         sk.key.lvl1);
        const auto e1 = std::chrono::steady_clock::now();
        r.t_encrypt_ms = std::chrono::duration<double, std::milli>(e1 - e0).count();
        r.input_phase = decrypt_phase(ct, sk.key.lvl1);

        TLWELvl1 bit_qhalf;
        const auto b0 = std::chrono::steady_clock::now();
        if (centered_bitextract)
            be::BitExtractQHalf_Pruned_Lvl1_CenteredCell(bit_qhalf, ct, p, k,
                                                         ek);
        else
            be::BitExtractQHalf_Pruned_Lvl1(bit_qhalf, ct, k, ek);
        const auto b1 = std::chrono::steady_clock::now();
        r.t_bitextract_ms =
            std::chrono::duration<double, std::milli>(b1 - b0).count();
        r.bit_qhalf_phase_raw = decrypt_phase(bit_qhalf, sk.key.lvl1);
        r.decoded_bit_qhalf = decode_qhalf_bit(r.bit_qhalf_phase_raw);

        TLWELvl1 guard;
        const auto c0 = std::chrono::steady_clock::now();
        if (conversion == "boolarith_existing") {
            mp::LegacyBool2WeightBaseline_Lvl1(guard, bit_qhalf, r.guard_value,
                                               ek);
        }
        else if (conversion == "micro_lut_lvl01_compat") {
            mp::MicroPBSQHalfToGuardValue_Lvl1Compat(guard, bit_qhalf,
                                                     r.guard_value, ek);
        }
        else if (conversion == "micro_unsafe_best") {
            const mp::MicroUnsafeBestPlan plan =
                mp::SelectMicroUnsafeBestPlan(k);
            if (!plan.supported) {
                r.status = "skipped";
                r.failure_class = plan.classification;
                return r;
            }
            mp::MicroTiming timing;
            mp::MicroUnsafeKsPreQHalfToGuardValue_Lvl1<BestMicroCandidate>(
                guard, bit_qhalf, r.guard_value, plan.offset_value,
                *micro_pack, &timing);
            r.t_micro_pre_ks_ms = timing.t_pre_ks_ms;
            r.t_micro_pbs_ms = timing.t_micro_pbs_ms;
            r.t_micro_post_ks_ms = timing.t_post_ks_ms;
        }
        else {
            r.status = "skipped";
            r.failure_class = "UNKNOWN_CONVERSION_MODE";
            return r;
        }
        const auto c1 = std::chrono::steady_clock::now();
        r.t_conversion_total_ms =
            std::chrono::duration<double, std::milli>(c1 - c0).count();
        if (conversion != "micro_unsafe_best")
            r.t_micro_pbs_ms = r.t_conversion_total_ms;

        r.decoded_guard = decrypt_phase(guard, sk.key.lvl1);
        r.guard_error =
            torus_centered_error(r.decoded_guard, r.expected_guard);

        TLWELvl1 gap;
        const auto s0 = std::chrono::steady_clock::now();
        for (int i = 0; i <= Lvl1::k * Lvl1::n; ++i) gap[i] = ct[i] - guard[i];
        const auto s1 = std::chrono::steady_clock::now();
        r.t_clear_sub_ms =
            std::chrono::duration<double, std::milli>(s1 - s0).count();
        r.gap_raw_phase = decrypt_phase(gap, sk.key.lvl1);
        r.gap_error = torus_centered_error(r.gap_raw_phase, r.expected_gap);

        TLWELvl1 shifted = gap;
        shifted[Lvl1::k * Lvl1::n] += r.gap_offset;
        TLWELvl1 out;
        const auto f0 = std::chrono::steady_clock::now();
        MSBGateBootstrapping(out, shifted, ek, LOGIC);
        const auto f1 = std::chrono::steady_clock::now();
        r.t_final_gap_ms =
            std::chrono::duration<double, std::milli>(f1 - f0).count();
        r.final_output_raw_phase = decrypt_phase(out, sk.key.lvl1);
        r.final_output = TFHEpp::tlweSymDecrypt<Lvl1>(out, sk.key.lvl1) ? 1U : 0U;

        if (r.decoded_bit_qhalf != r.expected_bit) {
            r.status = "fail";
            r.failure_class = "BITEXTRACT_BOUNDARY_FAILURE";
        }
        else if (torus_abs_centered(r.decoded_guard, r.expected_guard) >= r.A) {
            r.status = "fail";
            r.failure_class = "MICRO_CONVERSION_FAILURE";
        }
        else if (r.final_output != r.expected_output) {
            r.status = "fail";
            r.failure_class = "FINAL_GAP_FAILURE";
        }

        const auto total1 = std::chrono::steady_clock::now();
        r.t_total_ms =
            std::chrono::duration<double, std::milli>(total1 - total0).count();
        return r;
    }
} // namespace

int main(int argc, char **argv)
{
    const Options opts = parse_options(argc, argv);
#ifdef _OPENMP
    omp_set_num_threads(opts.threads);
#endif

    std::filesystem::create_directories(
        std::filesystem::path(opts.jsonl_path).parent_path());
    std::ofstream jsonl(opts.jsonl_path,
                        opts.append ? std::ios::app : std::ios::trunc);
    if (!jsonl) {
        std::cerr << "failed to open jsonl path: " << opts.jsonl_path << "\n";
        return 2;
    }

    std::cout << "Micro-PBS full 3-PBS pipeline probe\n";
    std::cout << "  UNSAFE_PERFORMANCE_ONLY=true\n";
    std::cout << "  bit_index_convention=MSB_FIRST_WINDOW_LOCAL\n";
    std::cout << "  GLOBAL_K_USED_IN_MICRO_CONVERSION=false\n";
    std::cout << "  centered_bitextract="
              << (opts.centered_bitextract ? "true" : "false") << "\n";
    std::cout << "  micro_unsafe_best="
              << BestMicroCandidate::name
              << " k<=6:Q/32 k<=8:Q/16 k>8:unsupported\n";

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);

    std::optional<mp::MicroEvalKeyPack<BestMicroCandidate>> micro_pack;
    if (wants_conversion(opts, "micro_unsafe_best")) {
        std::cout << "  keygen=micro unsafe best pack\n";
        micro_pack = mp::GenerateMicroEvalKeyPack<BestMicroCandidate>(
            sk, /*with_direct=*/false);
    }

    std::map<std::string, uint64_t> pass_count;
    std::map<std::string, uint64_t> fail_count;
    std::map<std::string, uint64_t> skipped_count;
    std::optional<Row> first_failure;

    for (const uint32_t p : opts.p_list) {
        for (const uint32_t k : opts.k_list) {
            if (!be::IsValidWindowLocalK(k) || k >= p) continue;
            const auto messages =
                build_messages(p, k, opts.boundary_trials, opts.random_trials);
            std::cout << "  p=" << p << " k=" << k
                      << " W=" << be::WindowWeight(p, k)
                      << " messages=" << messages.size() << "\n";

            for (const std::string &conversion : opts.conversions) {
                for (const uint64_t message : messages) {
                    const Row row =
                        run_case(conversion, p, k, message, sk, ek,
                                 micro_pack ? &*micro_pack : nullptr,
                                 opts.centered_bitextract);
                    write_jsonl(jsonl, row);
                    if (row.status == "pass")
                        pass_count[conversion]++;
                    else if (row.status == "skipped")
                        skipped_count[conversion]++;
                    else {
                        fail_count[conversion]++;
                        if (!first_failure) first_failure = row;
                    }
                }
            }
        }
    }

    std::cout << "\nSummary:\n";
    for (const std::string &conversion : opts.conversions) {
        std::cout << "  " << conversion << " pass="
                  << pass_count[conversion] << " fail="
                  << fail_count[conversion] << " skipped="
                  << skipped_count[conversion] << "\n";
    }
    if (first_failure)
        print_failure_trace(*first_failure);
    else
        std::cout << "FIRST_FAILURE none\n";
    std::cout << "jsonl=" << opts.jsonl_path << "\n";
    return first_failure ? 1 : 0;
}
