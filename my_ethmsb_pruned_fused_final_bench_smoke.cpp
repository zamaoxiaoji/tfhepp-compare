#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "my_ethmsb_pruned_fused_final.hpp"

namespace {

using namespace my_ethmsb_pruned_fused_final;
using Clock = std::chrono::steady_clock;

struct Options {
    std::vector<int> bits = {8, 16, 24, 32};
    std::string zero_bit_choice;
    std::string selected_k_json;
    int trials = 3;
    int warmup = 1;
    std::uint64_t seed = 0;
    std::string csv_path;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct Timing {
    double mean = 0;
    double median = 0;
    double p95 = 0;
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
    std::vector<int> out;
    std::size_t pos = 0;
    while (pos < s.size()) {
        const std::size_t comma = s.find(',', pos);
        out.push_back(std::stoi(s.substr(pos, comma - pos)));
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
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
        if (a == "--impl") {
            const std::string impl = need("--impl");
            if (impl != FUSED_FINAL_IMPL_NAME)
                throw std::runtime_error("unsupported --impl: " + impl);
        }
        else if (a == "--bits")
            opt.bits = parse_bits(need("--bits"));
        else if (a == "--zero-bit")
            opt.zero_bit_choice = need("--zero-bit");
        else if (a == "--use-selected-k")
            opt.selected_k_json = need("--use-selected-k");
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
    std::size_t pos = text.find(bits_pat);
    if (pos == std::string::npos)
        pos = text.find("\"bits\": " + std::to_string(bits));
    if (pos == std::string::npos) return std::nullopt;
    const std::size_t key = text.find("\"selected_k\"", pos);
    if (key == std::string::npos) return std::nullopt;
    const std::size_t colon = text.find(':', key);
    std::size_t val = colon + 1;
    while (val < text.size() && std::isspace(static_cast<unsigned char>(text[val])))
        ++val;
    if (text.compare(val, 4, "null") == 0) return std::nullopt;
    return std::stoi(text.substr(val));
}

int choose_k(const Options& opt, const int bits)
{
    const int p = bits + 1;
    if (!opt.zero_bit_choice.empty()) return parse_zero_bit_choice(opt.zero_bit_choice, p);
    const auto selected = selected_k_from_summary(opt.selected_k_json, bits);
    if (selected.has_value()) return *selected;
    throw std::runtime_error(
        "no correctness-supported default k for fused final experimental "
        "branch");
}

TFHEpp::TLWE<P2> encrypt_m(const Torus m, const int p, const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned_phase(m, p), ctx.key2);
}

bool out_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_value(phase_l2(ct, ctx.key2), BOOL_ONE);
}

bool correctness_ok(const Context& ctx, const int bits, const int k)
{
    const int p = bits + 1;
    for (const Torus m : {Torus{0}, Torus{1}, (Torus{1} << (p - 1)) - 1,
                          Torus{1} << (p - 1),
                          (Torus{1} << p) - 1}) {
        TFHEpp::TLWE<P2> out;
        FusedGapStats stats;
        gapmsb_pruned_bitk_fused_final_l2_to_l2(
            out, encrypt_m(m, p, ctx), p, k, BOOL_ONE,
            FusedInputRoundingMode::DeltaHalf, *ctx.iksk, *ctx.bkfft, &stats);
        const bool expected = ((m >> (p - 1)) & 1) != 0;
        if (out_bit(out, ctx) != expected) return false;
    }
    return true;
}

Timing summarize(std::vector<double> xs)
{
    Timing t;
    if (xs.empty()) return t;
    std::sort(xs.begin(), xs.end());
    t.mean = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    t.median = xs[xs.size() / 2];
    t.p95 = xs[std::min(xs.size() - 1,
                        static_cast<std::size_t>(
                            std::ceil(xs.size() * 0.95) - 1))];
    return t;
}

void header(std::ostream& out)
{
    out << "impl,bits,p,k,W,M,pbs_count,op,trials,warmup,mean_ms,median_ms,"
           "p95_ms,correctness_failures,failure_class,observed_prune_ratio,"
           "periodic_skipped_mean,zero_skipped_mean,executed_mean\n";
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
            if (!file) throw std::runtime_error("cannot open csv output");
            out = &file;
        }
        header(*out);

        std::mt19937_64 rng(opt.seed);
        for (const int bits : opt.bits) {
            const int p = bits + 1;
            const int k = choose_k(opt, bits);
            const Torus maxv = (Torus{1} << p) - 1;
            std::uniform_int_distribution<Torus> dist(0, maxv);
            const Torus W = bit_weight_msb_index(p, k);
            const Torus M = bit_period_plain(p, k);
            if (!correctness_ok(ctx, bits, k)) {
                *out << FUSED_FINAL_IMPL_NAME << ',' << bits << ',' << p
                     << ',' << k << ',' << W << ',' << M
                     << ",2,msb,0,0,0,0,0,1,fused_final_oracle_boundary,0,0,0,0\n";
                continue;
            }

            for (int i = 0; i < opt.warmup; ++i) {
                TFHEpp::TLWE<P2> out_ct;
                FusedGapStats stats;
                gapmsb_pruned_bitk_fused_final_l2_to_l2(
                    out_ct, encrypt_m(dist(rng), p, ctx), p, k, BOOL_ONE,
                    FusedInputRoundingMode::DeltaHalf, *ctx.iksk, *ctx.bkfft,
                    &stats);
            }
            std::vector<double> xs;
            FusedGapStats last;
            for (int i = 0; i < opt.trials; ++i) {
                TFHEpp::TLWE<P2> out_ct;
                const auto t0 = Clock::now();
                gapmsb_pruned_bitk_fused_final_l2_to_l2(
                    out_ct, encrypt_m(dist(rng), p, ctx), p, k, BOOL_ONE,
                    FusedInputRoundingMode::DeltaHalf, *ctx.iksk, *ctx.bkfft,
                    &last);
                const auto t1 = Clock::now();
                xs.push_back(
                    std::chrono::duration<double, std::milli>(t1 - t0)
                        .count());
            }
            const Timing t = summarize(xs);
            const double total = static_cast<double>(
                last.bitextract_stats.periodic_skipped +
                last.bitextract_stats.zero_skipped +
                last.bitextract_stats.executed_terms);
            const double ratio =
                total == 0.0
                    ? 0.0
                    : static_cast<double>(
                          last.bitextract_stats.periodic_skipped +
                          last.bitextract_stats.zero_skipped) /
                          total;
            *out << FUSED_FINAL_IMPL_NAME << ',' << bits << ',' << p << ','
                 << k << ',' << W << ',' << M << ",2,msb," << opt.trials
                 << ',' << opt.warmup << ',' << t.mean << ',' << t.median
                 << ',' << t.p95
                 << ",0,none," << ratio << ','
                 << last.bitextract_stats.periodic_skipped << ','
                 << last.bitextract_stats.zero_skipped << ','
                 << last.bitextract_stats.executed_terms << "\n";
        }
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\nstage: fused_final\nwhat_failed: "
                  << ex.what() << "\nEND_NEED_INFO\n";
        return 1;
    }
}
