/**
 * @file micro_pbs_sweep_test.cpp
 * @brief Unsafe micro-parameter sweep for Boolean Q/2 -> guard conversion.
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
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../comparison/micro_evalkey.h"
#include "../comparison/micro_pbs.h"
#include "../comparison/tfhepp_utils.h"

using namespace tfhepp_compare;

namespace
{
    struct Options {
        std::string mode = "isolated";
        std::string params = "all";
        std::string offsets = "all";
        std::string conversion = "micro_unsafe_best";
        int trials = 64;
        int runs = 20;
        int k_min = 1;
        int k_max = 12;
        int threads = 1;
        int boundary_trials = 64;
        int random_trials = 256;
        std::string jsonl_path = "experimental/results/micro_pbs_summary.jsonl";
    };

    struct OffsetSpec {
        std::string name;
        Lvl1::T value = 0;
    };

    struct NoiseStats {
        uint64_t p50_abs_error = 0;
        uint64_t p90_abs_error = 0;
        uint64_t p99_abs_error = 0;
        uint64_t max_abs_error = 0;
        double mean_abs_error = 0.0;
    };

    struct SweepRow {
        std::string kind = "isolated_sweep";
        std::string conversion_mode;
        std::string param_name;
        uint32_t n_in = 0;
        uint32_t N_out = 0;
        uint32_t glwe_k = 0;
        uint32_t pbs_level = 0;
        uint32_t pbs_basebit = 0;
        uint32_t iks_level = 0;
        uint32_t iks_basebit = 0;
        double alpha = 0.0;
        std::string offset_name;
        Lvl1::T offset_value = 0;
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
        double t_pre_ks_ms = 0.0;
        double t_micro_pbs_ms = 0.0;
        double t_post_ks_ms = 0.0;
        double t_total_ms = 0.0;
        double speedup_vs_lvl01compat = 0.0;
        std::string first_failure;
        std::string status;
    };

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
            if (arg == "--mode")
                opts.mode = need_value("--mode");
            else if (arg == "--params")
                opts.params = need_value("--params");
            else if (arg == "--offsets")
                opts.offsets = need_value("--offsets");
            else if (arg == "--conversion")
                opts.conversion = need_value("--conversion");
            else if (arg == "--trials")
                opts.trials = std::stoi(need_value("--trials"));
            else if (arg == "--runs")
                opts.runs = std::stoi(need_value("--runs"));
            else if (arg == "--k-min")
                opts.k_min = std::stoi(need_value("--k-min"));
            else if (arg == "--k-max")
                opts.k_max = std::stoi(need_value("--k-max"));
            else if (arg == "--threads")
                opts.threads = std::stoi(need_value("--threads"));
            else if (arg == "--jsonl")
                opts.jsonl_path = need_value("--jsonl");
            else if (arg == "--boundary-trials")
                opts.boundary_trials = std::stoi(need_value("--boundary-trials"));
            else if (arg == "--random-trials")
                opts.random_trials = std::stoi(need_value("--random-trials"));
            else
                std::cerr << "unknown option: " << arg << "\n";
        }
        if (opts.trials < 1) opts.trials = 1;
        if (opts.runs < 1) opts.runs = 1;
        if (opts.k_min < 1) opts.k_min = 1;
        if (opts.k_max < opts.k_min) opts.k_max = opts.k_min;
        if (opts.threads < 1) opts.threads = 1;
        return opts;
    }

    bool wants_token(const std::string &list, const std::string &token)
    {
        if (list == "all") return true;
        size_t start = 0;
        while (start <= list.size()) {
            const size_t comma = list.find(',', start);
            const std::string item =
                list.substr(start, comma == std::string::npos
                                       ? std::string::npos
                                       : comma - start);
            if (item == token) return true;
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        return false;
    }

    std::vector<OffsetSpec> all_offsets()
    {
        constexpr uint32_t d = std::numeric_limits<Lvl1::T>::digits;
        return {{"Q/4", Lvl1::T(1) << (d - 2)},
                {"Q/8", Lvl1::T(1) << (d - 3)},
                {"Q/16", Lvl1::T(1) << (d - 4)},
                {"Q/32", Lvl1::T(1) << (d - 5)},
                {"Q/64", Lvl1::T(1) << (d - 6)},
                {"Q/128", Lvl1::T(1) << (d - 7)},
                {"0", Lvl1::T(0)}};
    }

    std::vector<OffsetSpec> selected_offsets(const std::string &spec)
    {
        std::vector<OffsetSpec> out;
        for (const auto &o : all_offsets())
            if (wants_token(spec, o.name)) out.push_back(o);
        return out.empty() ? all_offsets() : out;
    }

    TFHEpp::TLWE<Lvl1> encrypt_qhalf(int bit, const TFHEpp::Key<Lvl1> &key)
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

    void write_jsonl(std::ofstream &os, const SweepRow &r)
    {
        if (!os) return;
        os << "{"
           << "\"kind\":\"" << json_escape(r.kind) << "\","
           << "\"conversion_mode\":\"" << json_escape(r.conversion_mode)
           << "\","
           << "\"param_name\":\"" << json_escape(r.param_name) << "\","
           << "\"n_in\":" << r.n_in << ","
           << "\"N_out\":" << r.N_out << ","
           << "\"glwe_k\":" << r.glwe_k << ","
           << "\"pbs_level\":" << r.pbs_level << ","
           << "\"pbs_basebit\":" << r.pbs_basebit << ","
           << "\"iks_level\":" << r.iks_level << ","
           << "\"iks_basebit\":" << r.iks_basebit << ","
           << "\"alpha\":" << std::setprecision(17) << r.alpha << ","
           << "\"offset_name\":\"" << json_escape(r.offset_name) << "\","
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
           << "\"t_pre_ks_ms\":" << r.t_pre_ks_ms << ","
           << "\"t_micro_pbs_ms\":" << r.t_micro_pbs_ms << ","
           << "\"t_post_ks_ms\":" << r.t_post_ks_ms << ","
           << "\"t_total_ms\":" << r.t_total_ms << ","
           << "\"speedup_vs_lvl01compat\":" << r.speedup_vs_lvl01compat << ","
           << "\"first_failure\":\"" << json_escape(r.first_failure) << "\","
           << "\"status\":\"" << json_escape(r.status) << "\""
           << "}\n";
    }

    double measure_lvl01compat_ms(const TFHESecretKey &sk, const TFHEEvalKey &ek,
                                  int runs)
    {
        constexpr uint32_t k = 1;
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        const Lvl1::T guard_value = Lvl1::T(1) << (digits - 1 - k);
        double total_ms = 0.0;
        for (int i = 0; i < runs; ++i) {
            const auto ct = encrypt_qhalf(i & 1, sk.key.lvl1);
            TLWELvl1 out;
            const auto a = std::chrono::steady_clock::now();
            micro_pbs::MicroPBSQHalfToGuardValue_Lvl1Compat(out, ct,
                                                            guard_value, ek);
            const auto b = std::chrono::steady_clock::now();
            total_ms += std::chrono::duration<double, std::milli>(b - a).count();
        }
        return total_ms / runs;
    }

    template <class Candidate>
    SweepRow profile_candidate_bit(
        const std::string &mode, const OffsetSpec &offset, uint32_t k, int bit,
        int trials, int runs, double baseline_ms, const TFHESecretKey &sk,
        const micro_pbs::MicroEvalKeyPack<Candidate> &pack)
    {
        using MicroOut = typename Candidate::MicroOut;
        constexpr uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        const Lvl1::T guard_value = Lvl1::T(1) << (digits - 1 - k);
        const Lvl1::T A = guard_value >> 1;
        const int timed_limit = std::min(trials, runs);

        std::vector<uint64_t> errors;
        errors.reserve(trials);
        micro_pbs::MicroTiming timing;
        int pass_count = 0;
        std::string first_failure;

        for (int t = 0; t < trials; ++t) {
            const auto ct = encrypt_qhalf(bit, sk.key.lvl1);
            TLWELvl1 out;
            micro_pbs::MicroTiming *maybe_timing =
                (t < timed_limit) ? &timing : nullptr;
            if (mode == "micro_unsafe_direct") {
                micro_pbs::MicroUnsafeDirectQHalfToGuardValue_Lvl1<Candidate>(
                    out, ct, guard_value, offset.value, pack, maybe_timing);
            }
            else {
                micro_pbs::MicroUnsafeKsPreQHalfToGuardValue_Lvl1<Candidate>(
                    out, ct, guard_value, offset.value, pack, maybe_timing);
            }

            const Lvl1::T phase = decrypt_torus(out, sk.key.lvl1);
            const Lvl1::T expected = bit ? guard_value : Lvl1::T(0);
            const uint64_t err = torus_abs_centered(phase, expected);
            errors.push_back(err);
            if (err < A)
                pass_count++;
            else if (first_failure.empty())
                first_failure =
                    "CONVERSION_NOISE_EXCEEDS_GUARD_MARGIN";
        }

        SweepRow row;
        row.conversion_mode = mode;
        row.param_name = Candidate::name;
        row.n_in = mode == "micro_unsafe_direct" ? Lvl1::n
                                                  : Candidate::MicroIn::n;
        row.N_out = MicroOut::n;
        row.glwe_k = MicroOut::k;
        row.pbs_level = MicroOut::l;
        row.pbs_basebit = MicroOut::Bgbit;
        row.iks_level = Candidate::iks_level;
        row.iks_basebit = Candidate::iks_basebit;
        row.alpha = MicroOut::α;
        row.offset_name = offset.name;
        row.offset_value = offset.value;
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
        row.t_pre_ks_ms = timing.t_pre_ks_ms / timed_limit;
        row.t_micro_pbs_ms = timing.t_micro_pbs_ms / timed_limit;
        row.t_post_ks_ms = timing.t_post_ks_ms / timed_limit;
        row.t_total_ms =
            row.t_pre_ks_ms + row.t_micro_pbs_ms + row.t_post_ks_ms;
        row.speedup_vs_lvl01compat =
            row.t_total_ms > 0.0 ? baseline_ms / row.t_total_ms : 0.0;
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

    template <class Candidate>
    void run_candidate(const Options &opts, const std::vector<OffsetSpec> &offsets,
                       double baseline_ms, const TFHESecretKey &sk,
                       std::ofstream &jsonl, std::vector<SweepRow> &rows)
    {
        if (!wants_token(opts.params, Candidate::name)) return;
        std::cout << "\n[keygen] " << Candidate::name << " n_in="
                  << Candidate::MicroIn::n << " N_out="
                  << Candidate::MicroOut::n << " l="
                  << Candidate::MicroOut::l << " Bgbit="
                  << Candidate::MicroOut::Bgbit << " iks=("
                  << Candidate::iks_level << "," << Candidate::iks_basebit
                  << ")\n";
        auto pack = micro_pbs::GenerateMicroEvalKeyPack<Candidate>(
            sk, /*with_direct=*/true);

        std::vector<std::string> modes;
        if (opts.conversion == "micro_unsafe_direct")
            modes.push_back("micro_unsafe_direct");
        else if (opts.conversion == "all")
            modes = {"micro_unsafe_direct", "micro_unsafe_ks_pre"};
        else
            modes.push_back("micro_unsafe_ks_pre");

        for (const auto &mode : modes) {
            for (const auto &offset : offsets) {
                int local_fail = 0;
                double local_ms = 0.0;
                int local_rows = 0;
                for (int k = opts.k_min; k <= opts.k_max; ++k) {
                    for (int bit : {0, 1}) {
                        SweepRow row = profile_candidate_bit<Candidate>(
                            mode, offset, static_cast<uint32_t>(k), bit,
                            opts.trials, opts.runs, baseline_ms, sk, pack);
                        if (row.status != "pass") local_fail++;
                        local_ms += row.t_total_ms;
                        local_rows++;
                        write_jsonl(jsonl, row);
                        rows.push_back(row);
                    }
                }
                std::cout << "  " << mode << " offset=" << std::setw(5)
                          << offset.name << " failures=" << local_fail
                          << " avg_ms="
                          << std::fixed << std::setprecision(3)
                          << (local_rows ? local_ms / local_rows : 0.0)
                          << "\n";
            }
        }
    }

    bool row_in_group(const SweepRow &r, const SweepRow &g)
    {
        return r.conversion_mode == g.conversion_mode &&
               r.param_name == g.param_name && r.offset_name == g.offset_name;
    }

    void print_best_for_sets(const std::vector<SweepRow> &rows)
    {
        struct SetSpec {
            const char *name;
            uint32_t max_k;
        };
        const std::vector<SetSpec> sets = {{"set_A_fast", 6},
                                           {"set_B_deeper", 8},
                                           {"set_C_aggressive", 10},
                                           {"set_D_stress", 12}};
        std::cout << "\nSelection summary:\n";
        for (const auto &set : sets) {
            std::optional<SweepRow> best;
            double best_ms = 0.0;
            for (const auto &candidate : rows) {
                if (candidate.k > set.max_k) continue;
                bool all_pass = true;
                double total_ms = 0.0;
                int count = 0;
                for (const auto &r : rows) {
                    if (!row_in_group(r, candidate)) continue;
                    if (r.k > set.max_k) continue;
                    if (r.status != "pass") all_pass = false;
                    total_ms += r.t_total_ms;
                    count++;
                }
                if (!all_pass || count == 0) continue;
                const double avg_ms = total_ms / count;
                if (!best || avg_ms < best_ms) {
                    best = candidate;
                    best_ms = avg_ms;
                }
            }
            if (best) {
                std::cout << "  " << set.name << " k<= " << set.max_k
                          << ": " << best->conversion_mode << " "
                          << best->param_name << " offset="
                          << best->offset_name << " avg_ms=" << std::fixed
                          << std::setprecision(3) << best_ms
                          << " speedup_vs_lvl01compat="
                          << (best_ms > 0.0
                                  ? (best->speedup_vs_lvl01compat *
                                     best->t_total_ms / best_ms)
                                  : 0.0)
                          << "\n";
            }
            else {
                std::cout << "  " << set.name << " k<= " << set.max_k
                          << ": NO_PASSING_CANDIDATE\n";
            }
        }
    }

    int run_full_mode(const Options &opts)
    {
        std::ofstream jsonl(opts.jsonl_path, std::ios::app);
        SweepRow row;
        row.kind = "full_pipeline_status";
        row.conversion_mode = opts.conversion;
        row.param_name = "micro_unsafe_best";
        row.status = "skipped";
        row.first_failure =
            "BITEXTRACT_BOUNDARY_FAILURE: full wrapper debug hook is not "
            "connected in this sweep target; isolated conversion evidence is "
            "reported separately";
        write_jsonl(jsonl, row);
        std::cout << "Full pipeline status: SKIPPED\n";
        std::cout << "  classification=BITEXTRACT_BOUNDARY_FAILURE\n";
        std::cout << "  reason=current work focuses on isolated conversion; "
                     "the previous BitExtract boundary bug must be isolated "
                     "before full-pipeline claims are meaningful\n";
        std::cout << "  boundary_trials=" << opts.boundary_trials
                  << " random_trials=" << opts.random_trials << "\n";
        return 0;
    }
} // namespace

int main(int argc, char **argv)
{
    const Options opts = parse_options(argc, argv);
#ifdef _OPENMP
    omp_set_num_threads(opts.threads);
#endif

    if (opts.mode == "full") return run_full_mode(opts);

    std::ofstream jsonl(opts.jsonl_path, std::ios::app);
    if (!jsonl) {
        std::cerr << "failed to open jsonl path: " << opts.jsonl_path << "\n";
        return 2;
    }

    std::cout << "Micro-PBS unsafe parameter sweep\n";
    std::cout << "  UNSAFE_PERFORMANCE_ONLY=true\n";
    std::cout << "  security_not_evaluated=true\n";
    std::cout << "  vendor_headers_modified=false\n";
    std::cout << "  params=" << opts.params << " offsets=" << opts.offsets
              << " trials=" << opts.trials << " runs=" << opts.runs
              << " k_range=" << opts.k_min << ".." << opts.k_max
              << " threads=" << opts.threads << "\n";

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    const double baseline_ms = measure_lvl01compat_ms(sk, ek, opts.runs);
    std::cout << "  measured_lvl01compat_ms=" << std::fixed
              << std::setprecision(3) << baseline_ms << "\n";

    const auto offsets = selected_offsets(opts.offsets);
    std::vector<SweepRow> rows;

    run_candidate<micro_pbs::micro_n64_N1024_l2_b8>(
        opts, offsets, baseline_ms, sk, jsonl, rows);
    run_candidate<micro_pbs::micro_n96_N1024_l2_b8>(
        opts, offsets, baseline_ms, sk, jsonl, rows);
    run_candidate<micro_pbs::micro_n128_N1024_l2_b8>(
        opts, offsets, baseline_ms, sk, jsonl, rows);
    run_candidate<micro_pbs::micro_n192_N1024_l2_b8>(
        opts, offsets, baseline_ms, sk, jsonl, rows);
    run_candidate<micro_pbs::micro_n256_N1024_l2_b8>(
        opts, offsets, baseline_ms, sk, jsonl, rows);
    run_candidate<micro_pbs::micro_n64_N1024_ahlvl1_l4_b5>(
        opts, offsets, baseline_ms, sk, jsonl, rows);

    if (rows.empty()) {
        std::cerr << "no candidate rows produced; check --params/--offsets\n";
        return 2;
    }
    print_best_for_sets(rows);
    std::cout << "\njsonl=" << opts.jsonl_path << "\n";
    return 0;
}
