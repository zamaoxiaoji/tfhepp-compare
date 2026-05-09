#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "my_ethmsb_pruned_bitextract.hpp"

namespace {

using namespace my_ethmsb_pruned_gap;
using Clock = std::chrono::steady_clock;

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct Options {
    std::vector<int> bits = {8, 16, 24, 32};
    std::string zero_bit_choice;
    std::string selected_k_json;
    std::vector<std::string> ops = {"lt", "gt", "eq"};
    int trials = 3;
    int warmup = 1;
    std::uint64_t seed = 0;
    std::string csv_path;
};

struct TimingStats {
    double mean = 0;
    double median = 0;
    double p95 = 0;
    double stddev = 0;
    double min = 0;
    double max = 0;
};

Context make_context(const std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    const auto keys = my_ethmsb_params::make_h3compat_keys(rng);
    Context ctx;
    ctx.key2 = keys.key2;
    ctx.key0 = keys.key0;
    ctx.iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    return ctx;
}

std::vector<int> parse_bits(const std::string& s)
{
    if (s == "all") return {8, 16, 24, 32};
    return {std::stoi(s)};
}

std::vector<std::string> parse_ops(const std::string& s)
{
    std::vector<std::string> out;
    std::size_t pos = 0;
    while (pos < s.size()) {
        const std::size_t comma = s.find(',', pos);
        out.push_back(s.substr(pos, comma - pos));
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
}

std::optional<int> selected_k_from_summary(const std::string& path,
                                           const int bits)
{
    if (path.empty()) return std::nullopt;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open selected-k summary: " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    const std::string bits_pat = "\"bits\":" + std::to_string(bits);
    const std::string bits_pat_spaced =
        "\"bits\": " + std::to_string(bits);
    std::size_t pos = text.find(bits_pat);
    if (pos == std::string::npos) pos = text.find(bits_pat_spaced);
    if (pos == std::string::npos) return std::nullopt;
    const std::size_t obj_end = text.find('}', pos);
    const std::string obj =
        text.substr(pos, obj_end == std::string::npos ? std::string::npos
                                                      : obj_end - pos);
    const std::string key = "\"selected_k\"";
    const std::size_t key_pos = obj.find(key);
    if (key_pos == std::string::npos) return std::nullopt;
    const std::size_t colon = obj.find(':', key_pos);
    if (colon == std::string::npos) return std::nullopt;
    std::size_t val = colon + 1;
    while (val < obj.size() && std::isspace(static_cast<unsigned char>(obj[val])))
        ++val;
    if (obj.compare(val, 4, "null") == 0) return std::nullopt;
    return std::stoi(obj.substr(val));
}

Options parse(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--bits")
            opt.bits = parse_bits(need("--bits"));
        else if (a == "--impl") {
            const std::string impl = need("--impl");
            if (impl != SINGLE_ROUND_IMPL_NAME)
                throw std::runtime_error(
                    "unsupported --impl for this benchmark target");
        }
        else if (a == "--zero-bit")
            opt.zero_bit_choice = need("--zero-bit");
        else if (a == "--use-selected-k")
            opt.selected_k_json = need("--use-selected-k");
        else if (a == "--ops")
            opt.ops = parse_ops(need("--ops"));
        else if (a == "--trials")
            opt.trials = std::stoi(need("--trials"));
        else if (a == "--warmup")
            opt.warmup = std::stoi(need("--warmup"));
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    return opt;
}

int choose_k_or_throw(const Options& opt, const int bits)
{
    const int p = bits + 1;
    if (!opt.zero_bit_choice.empty()) {
        const int k = parse_zero_bit_choice(opt.zero_bit_choice, p);
        if (p == 33 && k == p - 5)
            std::cerr << "WARNING: p-5 is margin-unsafe/not BR-index "
                         "supported under current params.\n";
        return k;
    }
    const std::optional<int> selected =
        selected_k_from_summary(opt.selected_k_json, bits);
    if (selected.has_value()) return *selected;
    throw std::runtime_error(
        "no correctness-supported default k for this p under single-round "
        "pruned variant");
}

TFHEpp::TLWE<P2> encrypt_compare_value(const Torus x, const int t,
                                       const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned_phase(x, t + 1),
                                      ctx.key2);
}

bool arith_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_value(phase_l2(ct, ctx.key2), BOOL_ONE);
}

bool eval_op(const Context& ctx, const int bits, const int k,
             const std::string& op, const Torus a0, const Torus b0,
             PrunedGapStats& stats)
{
    const auto a = encrypt_compare_value(a0, bits, ctx);
    const auto b = encrypt_compare_value(b0, bits, ctx);
    TFHEpp::TLWE<P2> out;
    if (op == "lt") {
        homcomp_lt_pruned_gap(out, a, b, bits, k, *ctx.iksk, *ctx.bkfft,
                              &stats);
    }
    else if (op == "gt") {
        homcomp_gt_pruned_gap(out, a, b, bits, k, *ctx.iksk, *ctx.bkfft,
                              &stats);
    }
    else if (op == "eq") {
        homcomp_eq_pruned_gap(out, a, b, bits, k, *ctx.iksk, *ctx.bkfft,
                              &stats);
    }
    else {
        throw std::runtime_error("unknown op: " + op);
    }
    return arith_bit(out, ctx);
}

std::set<std::pair<Torus, Torus>> correctness_pairs(const int bits)
{
    const Torus maxv = (Torus{1} << bits) - 1;
    const Torus mid = Torus{1} << (bits - 1);
    return {{0, 0},
            {0, 1},
            {1, 0},
            {0, maxv},
            {maxv, 0},
            {mid - 1, mid},
            {mid, mid - 1},
            {maxv, maxv}};
}

std::string correctness_gate(const Context& ctx, const int bits, const int k,
                             const std::string& op)
{
    const BitIndexInfo info = bit_index_info(bits + 1, k);
    for (const auto& [a, b] : correctness_pairs(bits)) {
        PrunedGapStats stats;
        const bool actual = eval_op(ctx, bits, k, op, a, b, stats);
        bool expected = false;
        if (op == "lt")
            expected = a < b;
        else if (op == "gt")
            expected = a > b;
        else if (op == "eq")
            expected = a == b;
        if (actual != expected)
            return info.exact_in_br_index ? "gap margin"
                                          : "BitExtract precision";
    }
    return "none";
}

TimingStats summarize(std::vector<double> xs)
{
    TimingStats s;
    if (xs.empty()) return s;
    std::sort(xs.begin(), xs.end());
    s.min = xs.front();
    s.max = xs.back();
    s.median = xs[xs.size() / 2];
    s.p95 = xs[std::min(xs.size() - 1,
                        static_cast<std::size_t>(
                            std::ceil(xs.size() * 0.95) - 1))];
    s.mean = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    double var = 0;
    for (const double x : xs) var += (x - s.mean) * (x - s.mean);
    s.stddev = std::sqrt(var / xs.size());
    return s;
}

void emit_header(std::ostream& out)
{
    out << "impl,p,k,w_k,M_k,rotation_period_index,gap_margin_hex,"
           "switch_band_hex,safe_margin_hex,prune_ratio_observed,op,trials,"
           "mean_ms,median_ms,p95_ms,correctness_failures\n";
}

void emit_row(std::ostream& out, const int bits, const int k,
              const std::string& op, const int trials, const int warmup,
              const TimingStats& timing, const PrunedGapStats& avg_stats,
              const int correctness_failures,
              const std::string& failure_class)
{
    const int p = bits + 1;
    const Torus w = bit_weight_msb_index(p, k);
    const Torus M = bit_period_msb_index(p, k);
    const BitIndexInfo info = bit_index_info(p, k);
    const Torus gap_margin = gap_offset_for_bit_index(p, k);
    const Torus safe_margin =
        gap_margin > DEFAULT_SWITCH_BAND_WIDTH
            ? gap_margin - DEFAULT_SWITCH_BAND_WIDTH
            : 0;
    const double total = static_cast<double>(
        avg_stats.periodic_skipped + avg_stats.zero_skipped +
        avg_stats.executed_terms);
    const double prune_ratio =
        total == 0 ? 0.0
                   : static_cast<double>(avg_stats.periodic_skipped +
                                         avg_stats.zero_skipped) /
                         total;
    (void)bits;
    (void)warmup;
    (void)failure_class;
    out << SINGLE_ROUND_IMPL_NAME << ',' << p << ',' << k
        << ',' << w << ',' << M << ',' << info.rotation_period_index << ','
        << hex64(gap_margin) << ',' << hex64(DEFAULT_SWITCH_BAND_WIDTH)
        << ',' << hex64(safe_margin) << ',' << prune_ratio << ',' << op
        << ',' << trials << ',' << timing.mean << ',' << timing.median
        << ',' << timing.p95 << ',' << correctness_failures << "\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const Context ctx = make_context(opt.seed);

        std::ofstream file;
        std::ostream* out = &std::cout;
        if (!opt.csv_path.empty()) {
            file.open(opt.csv_path);
            out = &file;
        }
        emit_header(*out);

        std::mt19937_64 rng(opt.seed);
        for (const int bits : opt.bits) {
            const int p = bits + 1;
            const int k = choose_k_or_throw(opt, bits);
            require_bit_index(p, k);
            const Torus maxv = (Torus{1} << bits) - 1;
            std::uniform_int_distribution<Torus> dist(0, maxv);
            for (const std::string& op : opt.ops) {
                const std::string failure_class =
                    correctness_gate(ctx, bits, k, op);
                if (failure_class != "none") {
                    emit_row(*out, bits, k, op, opt.trials, opt.warmup, {},
                             {}, 1, failure_class);
                    continue;
                }

                for (int i = 0; i < opt.warmup; ++i) {
                    PrunedGapStats stats;
                    (void)eval_op(ctx, bits, k, op, dist(rng), dist(rng),
                                  stats);
                }

                std::vector<double> samples;
                PrunedGapStats sums;
                for (int i = 0; i < opt.trials; ++i) {
                    const Torus a = dist(rng);
                    const Torus b = dist(rng);
                    PrunedGapStats stats;
                    const auto t0 = Clock::now();
                    (void)eval_op(ctx, bits, k, op, a, b, stats);
                    const auto t1 = Clock::now();
                    samples.push_back(
                        std::chrono::duration<double, std::milli>(t1 - t0)
                            .count());
                    sums.periodic_skipped += stats.periodic_skipped;
                    sums.zero_skipped += stats.zero_skipped;
                    sums.executed_terms += stats.executed_terms;
                    sums.pruning_enabled =
                        sums.pruning_enabled || stats.pruning_enabled;
                }
                if (opt.trials > 0) {
                    sums.periodic_skipped /= opt.trials;
                    sums.zero_skipped /= opt.trials;
                    sums.executed_terms /= opt.trials;
                }
                emit_row(*out, bits, k, op, opt.trials, opt.warmup,
                         summarize(samples), sums, 0, "none");
            }
        }
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: benchmark\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr
            << "  ./build-ethmsb-release/my_ethmsb_pruned_gapmsb_bench_smoke\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
