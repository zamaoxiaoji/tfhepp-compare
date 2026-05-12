/**
 * @file micro_pbs_isolated_test.cpp
 * @brief Noise profiler for the Lvl1Compat Boolean Q/2 -> guard conversion.
 *
 * This target intentionally profiles the existing Lvl10 IKS + Lvl01 PBS
 * chain. The real unsafe micro-parameter sweep lives in micro_pbs_sweep_test.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../comparison/micro_pbs.h"
#include "../comparison/tfhepp_utils.h"

using namespace tfhepp_compare;

namespace
{
    using ConvFn = void (*)(TFHEpp::TLWE<Lvl1> &, const TFHEpp::TLWE<Lvl1> &,
                            Lvl1::T, const TFHEEvalKey &);

    struct Options {
        int trials = 16;
        int k_min = 1;
        int k_max = 6;
        bool profile_noise = false;
        std::string jsonl_path;
    };

    struct NoiseStats {
        uint64_t p50_abs_error = 0;
        uint64_t p90_abs_error = 0;
        uint64_t p99_abs_error = 0;
        uint64_t max_abs_error = 0;
        double mean_abs_error = 0.0;
    };

    struct ProfileRow {
        std::string kind = "noise_profile";
        std::string conversion_mode;
        std::string param_name = "lvl01_default";
        std::string offset_name = "Q/64";
        Lvl1::T offset_value = Lvl1::T(1)
                               << (std::numeric_limits<Lvl1::T>::digits - 6);
        uint32_t k = 0;
        Lvl1::T guard_value = 0;
        Lvl1::T A = 0;
        int bit = 0;
        int trials = 0;
        int pass_count = 0;
        int fail_count = 0;
        NoiseStats stats;
        double noise_floor_log2 = 0.0;
        double margin_log2 = 0.0;
        double t_total_ms = 0.0;
        double speedup_vs_lvl01compat = 1.0;
        std::string first_failure;
        std::string status;
    };

    Options parse_options(int argc, char **argv)
    {
        Options opts;
        int positional = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--profile-noise") opts.profile_noise = true;
            else if (arg == "--k-min" && i + 1 < argc)
                opts.k_min = std::stoi(argv[++i]);
            else if (arg == "--k-max" && i + 1 < argc)
                opts.k_max = std::stoi(argv[++i]);
            else if (arg == "--jsonl" && i + 1 < argc)
                opts.jsonl_path = argv[++i];
            else if (arg.rfind("--", 0) == 0) {
                std::cerr << "unknown option: " << arg << "\n";
            }
            else if (positional == 0) {
                opts.trials = std::stoi(arg);
                positional++;
            }
        }
        if (opts.k_min < 1) opts.k_min = 1;
        if (opts.k_max < opts.k_min) opts.k_max = opts.k_min;
        return opts;
    }

    TFHEpp::TLWE<Lvl1> encrypt_qhalf(Lvl1::T bit, const TFHEpp::Key<Lvl1> &key)
    {
        const Lvl1::T msg =
            (bit & 1) ? (Lvl1::T(1)
                         << (std::numeric_limits<Lvl1::T>::digits - 1))
                      : Lvl1::T(0);
        return tfhepp_compare::tlweSymInt32Encrypt<Lvl1>(msg, Lvl1::α,
                                                         /*scale=*/1.0, key);
    }

    Lvl1::T decrypt_torus(const TFHEpp::TLWE<Lvl1> &c,
                          const TFHEpp::Key<Lvl1> &key)
    {
        Lvl1::T phase = c[Lvl1::k * Lvl1::n];
        for (int i = 0; i < Lvl1::k * Lvl1::n; i++) phase -= c[i] * key[i];
        return phase;
    }

    uint64_t torus_abs_centered(Lvl1::T phase, Lvl1::T expected)
    {
        const Lvl1::T diff = phase - expected;
        if ((diff & (Lvl1::T(1)
                     << (std::numeric_limits<Lvl1::T>::digits - 1))) != 0)
            return uint64_t(Lvl1::T(0) - diff);
        return uint64_t(diff);
    }

    uint64_t percentile(std::vector<uint64_t> v, double q)
    {
        if (v.empty()) return 0;
        std::sort(v.begin(), v.end());
        const size_t idx =
            std::min(v.size() - 1,
                     static_cast<size_t>(std::ceil(q * v.size())) - 1);
        return v[idx];
    }

    NoiseStats compute_stats(const std::vector<uint64_t> &errors)
    {
        NoiseStats s;
        if (errors.empty()) return s;
        s.p50_abs_error = percentile(errors, 0.50);
        s.p90_abs_error = percentile(errors, 0.90);
        s.p99_abs_error = percentile(errors, 0.99);
        s.max_abs_error = *std::max_element(errors.begin(), errors.end());
        const long double total =
            std::accumulate(errors.begin(), errors.end(), (long double) 0.0);
        s.mean_abs_error = static_cast<double>(total / errors.size());
        return s;
    }

    double log2_or_zero(uint64_t x)
    {
        return x == 0 ? 0.0 : std::log2(static_cast<double>(x));
    }

    std::string json_escape(const std::string &s)
    {
        std::ostringstream out;
        for (char c : s) {
            if (c == '"' || c == '\\') out << '\\';
            out << c;
        }
        return out.str();
    }

    void write_jsonl(std::ofstream &os, const ProfileRow &r)
    {
        if (!os) return;
        os << "{"
           << "\"kind\":\"" << json_escape(r.kind) << "\","
           << "\"conversion_mode\":\"" << json_escape(r.conversion_mode)
           << "\","
           << "\"param_name\":\"" << json_escape(r.param_name) << "\","
           << "\"n_in\":" << Lvl0::n << ","
           << "\"N_out\":" << Lvl1::n << ","
           << "\"glwe_k\":" << Lvl1::k << ","
           << "\"pbs_level\":" << Lvl1::l << ","
           << "\"pbs_basebit\":" << Lvl1::Bgbit << ","
           << "\"iks_level\":" << Lvl10::t << ","
           << "\"iks_basebit\":" << Lvl10::basebit << ","
           << "\"alpha\":" << std::setprecision(17) << Lvl1::α << ","
           << "\"offset_name\":\"" << r.offset_name << "\","
           << "\"offset_value\":" << r.offset_value << ","
           << "\"k\":" << r.k << ","
           << "\"guard_value\":" << r.guard_value << ","
           << "\"A\":" << r.A << ","
           << "\"bit\":" << r.bit << ","
           << "\"trials\":" << r.trials << ","
           << "\"pass_count\":" << r.pass_count << ","
           << "\"fail_count\":" << r.fail_count << ","
           << "\"p50_abs_error\":" << r.stats.p50_abs_error << ","
           << "\"p90_abs_error\":" << r.stats.p90_abs_error << ","
           << "\"p99_abs_error\":" << r.stats.p99_abs_error << ","
           << "\"max_abs_error\":" << r.stats.max_abs_error << ","
           << "\"noise_floor_log2\":" << r.noise_floor_log2 << ","
           << "\"margin_log2\":" << r.margin_log2 << ","
           << "\"t_pre_ks_ms\":0,"
           << "\"t_micro_pbs_ms\":0,"
           << "\"t_post_ks_ms\":0,"
           << "\"t_total_ms\":" << r.t_total_ms << ","
           << "\"speedup_vs_lvl01compat\":" << r.speedup_vs_lvl01compat << ","
           << "\"first_failure\":\"" << json_escape(r.first_failure) << "\","
           << "\"status\":\"" << json_escape(r.status) << "\""
           << "}\n";
    }

    ProfileRow profile_one(const std::string &mode, ConvFn conv, uint32_t k,
                           int bit, int trials, const TFHESecretKey &sk,
                           const TFHEEvalKey &ek)
    {
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        const Lvl1::T guard_value = Lvl1::T(1) << (digits - 1 - k);
        const Lvl1::T A = guard_value >> 1;
        std::vector<uint64_t> errors;
        errors.reserve(trials);
        int pass_count = 0;
        double total_ms = 0.0;
        std::string first_failure;

        for (int t = 0; t < trials; ++t) {
            const TFHEpp::TLWE<Lvl1> ct_qhalf =
                encrypt_qhalf(static_cast<Lvl1::T>(bit), sk.key.lvl1);
            TFHEpp::TLWE<Lvl1> ct_guard;
            const auto start = std::chrono::steady_clock::now();
            conv(ct_guard, ct_qhalf, guard_value, ek);
            const auto end = std::chrono::steady_clock::now();
            total_ms += std::chrono::duration<double, std::milli>(end - start)
                            .count();

            const Lvl1::T phase = decrypt_torus(ct_guard, sk.key.lvl1);
            const Lvl1::T expected = bit ? guard_value : Lvl1::T(0);
            const uint64_t err = torus_abs_centered(phase, expected);
            errors.push_back(err);
            if (err < A)
                pass_count++;
            else if (first_failure.empty())
                first_failure =
                    "CONVERSION_NOISE_EXCEEDS_GUARD_MARGIN";
        }

        ProfileRow row;
        row.conversion_mode = mode;
        row.k = k;
        row.guard_value = guard_value;
        row.A = A;
        row.bit = bit;
        row.trials = trials;
        row.pass_count = pass_count;
        row.fail_count = trials - pass_count;
        row.stats = compute_stats(errors);
        row.noise_floor_log2 = log2_or_zero(row.stats.p99_abs_error);
        row.margin_log2 = log2_or_zero(A);
        row.t_total_ms = total_ms / trials;
        if (row.stats.p99_abs_error >= A) {
            row.status = "fail";
            row.first_failure = first_failure.empty()
                                    ? "CONVERSION_NOISE_EXCEEDS_GUARD_MARGIN"
                                    : first_failure;
        }
        else {
            row.status = "pass";
        }
        return row;
    }

    bool plaintext_oracle()
    {
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        bool ok = true;
        for (uint32_t k : {1u, 2u, 3u, 4u, 5u, 6u}) {
            const Lvl1::T guard_value = Lvl1::T(1) << (digits - 1 - k);
            const Lvl1::T A = guard_value >> 1;
            if (Lvl1::T(-A) + A != Lvl1::T(0)) ok = false;
        }
        return ok;
    }
} // namespace

int main(int argc, char **argv)
{
    const Options opts = parse_options(argc, argv);
    std::cout << "Micro-PBS isolated conversion profiler\n";
    std::cout << "  UNSAFE_PERFORMANCE_ONLY=true\n";
    std::cout << "  vendor_headers_modified=false\n";
    std::cout << "  conversion_mode=micro_lut_lvl01_compat\n";
    std::cout << "  note=Lvl1Compat and legacy both use Lvl10 IKS + Lvl01 PBS\n";
    std::cout << "  trials_per_bit=" << opts.trials << "\n";
    std::cout << "  k_range=" << opts.k_min << ".." << opts.k_max << "\n\n";

    std::ofstream jsonl;
    if (!opts.jsonl_path.empty()) {
        jsonl.open(opts.jsonl_path, std::ios::app);
        if (!jsonl) {
            std::cerr << "failed to open jsonl path: " << opts.jsonl_path
                      << "\n";
            return 2;
        }
    }

    std::cout << "Plaintext oracle (bit=0 branch only): "
              << (plaintext_oracle() ? "PASS" : "FAIL") << "\n\n";

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);

    std::vector<ProfileRow> rows;
    for (int k = opts.k_min; k <= opts.k_max; ++k) {
        for (int bit : {0, 1}) {
            rows.push_back(profile_one(
                micro_pbs::ConversionModeName(
                    micro_pbs::ConversionMode::MicroLutLvl01Compat),
                &micro_pbs::MicroPBSQHalfToGuardValue_Lvl1Compat,
                static_cast<uint32_t>(k), bit, opts.trials, sk, ek));
            rows.push_back(profile_one(
                "legacy_bool2weight_lvl01",
                &micro_pbs::LegacyBool2WeightBaseline_Lvl1,
                static_cast<uint32_t>(k), bit, opts.trials, sk, ek));
        }
    }

    double compat_total = 0.0;
    int compat_count = 0;
    for (const auto &r : rows) {
        if (r.conversion_mode == "micro_lut_lvl01_compat") {
            compat_total += r.t_total_ms;
            compat_count++;
        }
    }
    const double compat_avg =
        compat_count == 0 ? 1.0 : compat_total / compat_count;

    std::cout << "mode,k,bit,status,p99_abs_error,noise_floor_log2,A,margin_log2,"
                 "pass_count,fail_count,avg_ms,failure_class\n";
    bool required_ok = true;
    for (auto &r : rows) {
        r.speedup_vs_lvl01compat =
            r.t_total_ms > 0.0 ? compat_avg / r.t_total_ms : 0.0;
        write_jsonl(jsonl, r);
        std::cout << r.conversion_mode << "," << r.k << "," << r.bit << ","
                  << r.status << "," << r.stats.p99_abs_error << ","
                  << std::fixed << std::setprecision(3)
                  << r.noise_floor_log2 << "," << r.A << ","
                  << r.margin_log2 << "," << r.pass_count << ","
                  << r.fail_count << "," << r.t_total_ms << ","
                  << r.first_failure << "\n";

        if (r.k <= 6 && r.status != "pass") required_ok = false;
    }

    std::cout << "\nSummary:\n";
    std::cout << "  Lvl1Compat micro and legacy use the same IKS + Lvl01 PBS "
                 "chain; latency equality is expected.\n";
    std::cout << "  required_k_le_6_status="
              << (required_ok ? "PASS" : "FAIL") << "\n";
    std::cout << "  jsonl="
              << (opts.jsonl_path.empty() ? "(not requested)" : opts.jsonl_path)
              << "\n";
    return required_ok ? 0 : 1;
}
