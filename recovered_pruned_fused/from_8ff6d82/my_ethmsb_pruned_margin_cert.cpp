#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "my_ethmsb_pruned_bitextract.hpp"

namespace {

using namespace my_ethmsb_pruned_gap;

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct Options {
    std::vector<int> ps = {9, 17, 25, 33};
    std::string zero_bit_choice = "p-5";
    std::string mode = "trivial";
    std::uint64_t seed = 0;
    Torus switch_band_width = (Torus{1} << (64 - (P2::nbit + 1))) / 8;
};

struct Summary {
    Torus min_distance = std::numeric_limits<Torus>::max();
    Torus min_safe_margin = std::numeric_limits<Torus>::max();
    int failures = 0;
    std::uint64_t periodic_sum = 0;
    std::uint64_t zero_sum = 0;
    std::uint64_t executed_sum = 0;
    std::uint64_t count = 0;
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
        if (a == "--p") {
            opt.ps = {std::stoi(need("--p"))};
        }
        else if (a == "--zero-bit") {
            opt.zero_bit_choice = need("--zero-bit");
        }
        else if (a == "--mode") {
            opt.mode = need("--mode");
        }
        else if (a == "--seed") {
            opt.seed = std::stoull(need("--seed"));
        }
        else if (a == "--switch-band-width-hex") {
            opt.switch_band_width =
                static_cast<Torus>(std::strtoull(need("--switch-band-width-hex"),
                                                 nullptr, 0));
        }
        else {
            throw std::runtime_error("unknown argument: " + a);
        }
    }
    return opt;
}

std::vector<Torus> candidates(const int p, const int k)
{
    std::set<Torus> vals = {0,
                            1,
                            (Torus{1} << (p - 1)) - 1,
                            Torus{1} << (p - 1),
                            (Torus{1} << (p - 1)) + 1,
                            (Torus{1} << p) - 1};
    const Torus half = Torus{1} << (p - 1);
    const Torus w = bit_weight_msb_index(p, k);
    for (const Torus v :
         {half - 1 - w, half - w, half - 1, half, w - 1, w, w + 1,
          2 * w - 1, 2 * w, 2 * w + 1}) {
        if (v < (Torus{1} << p)) vals.insert(v);
    }
    return {vals.begin(), vals.end()};
}

TFHEpp::TLWE<P2> make_input(const Torus m, const int p, const Context& ctx,
                            const std::string& mode)
{
    const Torus ph = encode_unsigned_phase(m, p);
    if (mode == "trivial") return trivial_encrypt_phase_l2(ph);
    if (mode == "encrypted") return TFHEpp::tlweSymEncrypt<P2>(ph, ctx.key2);
    throw std::runtime_error("unknown --mode");
}

void emit_header()
{
    std::cout
        << "p,k,w_k,M_k,guard_value,gap_offset,phase_before_offset,"
           "phase_after_offset,distance_to_half,measured_switch_band,"
           "safe_margin,pruning_enabled,periodic_skipped,zero_skipped,"
           "executed,expected,actual,pass\n";
}

int run_case(const Context& ctx, const Options& opt, const int p, const int k,
             const Torus m, Summary& summary)
{
    const Torus w = bit_weight_msb_index(p, k);
    const Torus guard = guard_value_for_bit_index(p, k);
    const Torus offset = gap_offset_for_bit_index(p, k);
    const auto ct = make_input(m, p, ctx, opt.mode);

    TFHEpp::TLWE<P2> bit;
    PruneStats prune;
    bit_extract_pruned_qhalf_l2_to_l2(bit, ct, p, k, *ctx.iksk, *ctx.bkfft,
                                      &prune);
    TFHEpp::TLWE<P2> mask;
    bool_qhalf_to_value_l2_to_l2(mask, bit, guard, *ctx.iksk, *ctx.bkfft);
    TFHEpp::TLWE<P2> gap;
    my_ethmsb::sub<P2>(gap, ct, mask);
    const Torus before = phase_l2(gap, ctx.key2);
    const Torus after = before + offset;
    const Torus dist = centered_distance(after, Q_HALF);
    const Torus safe =
        dist > opt.switch_band_width ? dist - opt.switch_band_width : 0;

    TFHEpp::TLWE<P2> out;
    pbs_msb_value_l2_to_l2(out, gap, p, offset, BOOL_ONE, *ctx.iksk,
                           *ctx.bkfft);
    const bool expected = ((m >> (p - 1)) & 1) != 0;
    const bool actual = closer_to_value(phase_l2(out, ctx.key2), BOOL_ONE);
    const bool pass = expected == actual;

    std::cout << p << ',' << k << ',' << w << ',' << bit_period_msb_index(p, k)
              << ',' << hex64(guard) << ',' << hex64(offset) << ','
              << hex64(before) << ',' << hex64(after) << ',' << hex64(dist)
              << ',' << hex64(opt.switch_band_width) << ',' << hex64(safe)
              << ',' << prune.pruning_enabled << ',' << prune.periodic_skipped
              << ',' << prune.zero_skipped << ',' << prune.executed_terms
              << ',' << expected << ',' << actual << ',' << pass << "\n";

    summary.min_distance = std::min(summary.min_distance, dist);
    summary.min_safe_margin = std::min(summary.min_safe_margin, safe);
    summary.periodic_sum += prune.periodic_skipped;
    summary.zero_sum += prune.zero_skipped;
    summary.executed_sum += prune.executed_terms;
    ++summary.count;
    if (!pass) ++summary.failures;
    return pass ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const Context ctx = make_context(opt.seed);
        emit_header();
        for (const int p : opt.ps) {
            const int k = parse_zero_bit_choice(opt.zero_bit_choice, p);
            Summary summary;
            for (const Torus m : candidates(p, k))
                (void)run_case(ctx, opt, p, k, m, summary);
            const Torus w = bit_weight_msb_index(p, k);
            const Torus M = bit_period_msb_index(p, k);
            const Torus guard = guard_value_for_bit_index(p, k);
            const Torus offset = gap_offset_for_bit_index(p, k);
            const BitIndexInfo info = bit_index_info(p, k);
            const std::uint64_t denom = std::max<std::uint64_t>(1, summary.count);
            std::cerr << "PRUNED_MARGIN_SUMMARY p=" << p << " k=" << k
                      << " w_k=" << w << " M_k=" << M
                      << " guard_value=" << hex64(guard)
                      << " gap_offset=" << hex64(offset)
                      << " min_distance_to_half="
                      << hex64(summary.min_distance)
                      << " min_safe_margin="
                      << hex64(summary.min_safe_margin)
                      << " pruning_enabled=" << info.supports_period_pruning
                      << " periodic_skipped_mean="
                      << (summary.periodic_sum / denom)
                      << " periodic_skipped_p95=not_computed"
                      << " failures=" << summary.failures << "\n";
        }
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: margin\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_margin_cert\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
