#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
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

struct Options {
    std::vector<int> bits = {8, 16, 24, 32};
    int p_override = 0;
    std::string bit_index = "all";
    int trials_per_case = 10;
    std::uint64_t seed = 0;
    std::string mode = "encrypted";
    std::string csv_path;
    int max_cases = 0;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct Summary {
    std::uint64_t total = 0;
    std::uint64_t failures = 0;
    std::uint64_t bitextract_failures = 0;
    std::uint64_t conversion_failures = 0;
    std::uint64_t final_gap_failures = 0;
    std::uint64_t wrapper_failures = 0;
};

std::vector<int> parse_bits_arg(const std::string& s)
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
        if (a == "--bits")
            opt.bits = parse_bits_arg(need("--bits"));
        else if (a == "--p")
            opt.p_override = std::stoi(need("--p"));
        else if (a == "--bit-index")
            opt.bit_index = need("--bit-index");
        else if (a == "--trials-per-case")
            opt.trials_per_case = std::stoi(need("--trials-per-case"));
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
        else if (a == "--mode")
            opt.mode = need("--mode");
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else if (a == "--max-cases")
            opt.max_cases = std::stoi(need("--max-cases"));
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (opt.p_override != 0) opt.bits = {opt.p_override - 1};
    return opt;
}

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

TFHEpp::TLWE<P2> make_ct(const Torus m, const int p, const Context& ctx,
                         const std::string& mode)
{
    const Torus phase = encode_unsigned_phase(m, p);
    if (mode == "trivial") return trivial_encrypt_phase_l2(phase);
    if (mode == "encrypted") return TFHEpp::tlweSymEncrypt<P2>(phase, ctx.key2);
    throw std::runtime_error("unknown --mode");
}

bool arith_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx,
               const Torus out_value = BOOL_ONE)
{
    return closer_to_value(phase_l2(ct, ctx.key2), out_value);
}

std::vector<int> k_values(const int p, const std::string& arg)
{
    if (arg != "all") return {std::stoi(arg)};
    std::vector<int> out;
    for (int k = 1; k <= p - 1; ++k) out.push_back(k);
    return out;
}

std::vector<std::pair<Torus, Torus>> base_pairs(const int bits)
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

std::string classify_failure(const bool bitextract_pass,
                             const bool conversion_pass,
                             const bool final_pass_given_bit_correct,
                             const bool overall_pass)
{
    if (!bitextract_pass) return "bitextract_boundary";
    if (!conversion_pass) return "conversion";
    if (!final_pass_given_bit_correct) return "final_gap_margin";
    if (!overall_pass) return "wrapper";
    return "none";
}

void emit_header(std::ostream& out)
{
    out << "bits,p,k,a,b,m,expected_msb,expected_bit_k,bit_actual,"
           "bitextract_pass,conversion_pass,final_pass_given_bit_correct,"
           "overall_pass,failure_class,bit_qhalf_phase_hex,guard_value_hex,"
           "guard_mask_phase_hex,gap_offset_hex,final_phase_hex,"
           "periodic_skipped,zero_skipped,executed,observed_prune_ratio\n";
}

void eval_pair(std::ostream& out, Summary& summary, const Context& ctx,
               const Options& opt, const int bits, const int k, const Torus a0,
               const Torus b0)
{
    const int p = bits + 1;
    const Torus p_mask = (Torus{1} << p) - 1;
    const Torus m = (a0 - b0) & p_mask;
    const bool expected_msb = ((m >> (p - 1)) & 1) != 0;
    const bool expected_bit_k = ((m >> (p - 1 - k)) & 1) != 0;
    const Torus guard_value = guard_value_for_bit_index(p, k);
    const Torus gap_offset = gap_offset_for_bit_index(p, k);
    const Torus w = bit_weight_msb_index(p, k);
    const Torus expected_m_gap = m - (expected_bit_k ? w : 0);
    (void)expected_m_gap;

    const auto a = make_ct(a0, p, ctx, opt.mode);
    const auto b = make_ct(b0, p, ctx, opt.mode);
    TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, a, b);

    TFHEpp::TLWE<P2> bit;
    PruneStats prune;
    bit_extract_pruned_qhalf_l2_to_l2(bit, diff, p, k, *ctx.iksk,
                                      *ctx.bkfft, &prune);
    const Torus bit_phase = phase_l2(bit, ctx.key2);
    const bool bit_actual = closer_to_qhalf(bit_phase);
    const bool bitextract_pass = bit_actual == expected_bit_k;

    TFHEpp::TLWE<P2> guard_mask;
    bool_qhalf_to_value_l2_to_l2(guard_mask, bit, guard_value, *ctx.iksk,
                                 *ctx.bkfft);
    const Torus guard_phase = phase_l2(guard_mask, ctx.key2);
    const bool guard_actual = closer_to_value(guard_phase, guard_value);
    const bool conversion_pass = guard_actual == bit_actual;

    TFHEpp::TLWE<P2> gap;
    my_ethmsb::sub<P2>(gap, diff, guard_mask);
    const Torus guarded_phase = phase_l2(gap, ctx.key2);
    (void)guarded_phase;
    const bool clear_pass_if_bit_correct =
        !bitextract_pass || guard_actual == expected_bit_k;
    (void)clear_pass_if_bit_correct;

    TFHEpp::TLWE<P2> final_out;
    pbs_msb_value_l2_to_l2(final_out, gap, p, gap_offset, BOOL_ONE,
                           *ctx.iksk, *ctx.bkfft);
    const Torus final_phase = phase_l2(final_out, ctx.key2);
    const bool final_actual_msb = closer_to_value(final_phase, BOOL_ONE);
    const bool final_pass_given_bit_correct =
        bitextract_pass && conversion_pass &&
        (final_actual_msb == expected_msb);

    TFHEpp::TLWE<P2> lt;
    TFHEpp::TLWE<P2> gt;
    TFHEpp::TLWE<P2> eq;
    homcomp_lt_pruned_gap(lt, a, b, bits, k, *ctx.iksk, *ctx.bkfft);
    homcomp_gt_pruned_gap(gt, a, b, bits, k, *ctx.iksk, *ctx.bkfft);
    homcomp_eq_pruned_gap(eq, a, b, bits, k, *ctx.iksk, *ctx.bkfft);
    const bool overall_pass = arith_bit(lt, ctx) == (a0 < b0) &&
                              arith_bit(gt, ctx) == (a0 > b0) &&
                              arith_bit(eq, ctx) == (a0 == b0);

    const std::string failure_class = classify_failure(
        bitextract_pass, conversion_pass, final_pass_given_bit_correct,
        overall_pass);
    const double total_terms = static_cast<double>(
        prune.periodic_skipped + prune.zero_skipped + prune.executed_terms);
    const double observed_prune_ratio =
        total_terms == 0.0
            ? 0.0
            : static_cast<double>(prune.periodic_skipped +
                                  prune.zero_skipped) /
                  total_terms;

    out << bits << ',' << p << ',' << k << ',' << a0 << ',' << b0 << ','
        << m << ',' << expected_msb << ',' << expected_bit_k << ','
        << bit_actual << ',' << bitextract_pass << ',' << conversion_pass
        << ',' << final_pass_given_bit_correct << ',' << overall_pass << ','
        << failure_class << ',' << hex64(bit_phase) << ','
        << hex64(guard_value) << ',' << hex64(guard_phase) << ','
        << hex64(gap_offset) << ',' << hex64(final_phase) << ','
        << prune.periodic_skipped << ',' << prune.zero_skipped << ','
        << prune.executed_terms << ',' << observed_prune_ratio << "\n";

    ++summary.total;
    if (failure_class != "none") ++summary.failures;
    if (failure_class == "bitextract_boundary") ++summary.bitextract_failures;
    if (failure_class == "conversion") ++summary.conversion_failures;
    if (failure_class == "final_gap_margin") ++summary.final_gap_failures;
    if (failure_class == "wrapper") ++summary.wrapper_failures;
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

        Summary summary;
        std::mt19937_64 rng(opt.seed);
        for (const int bits : opt.bits) {
            const int p = bits + 1;
            const Torus maxv = (Torus{1} << bits) - 1;
            std::uniform_int_distribution<Torus> dist(0, maxv);
            for (const int k : k_values(p, opt.bit_index)) {
                require_bit_index(p, k);
                for (const auto& [a, b] : base_pairs(bits)) {
                    if (opt.max_cases > 0 &&
                        summary.total >=
                            static_cast<std::uint64_t>(opt.max_cases))
                        goto done;
                    eval_pair(*out, summary, ctx, opt, bits, k, a, b);
                }
                for (int i = 0; i < opt.trials_per_case; ++i) {
                    if (opt.max_cases > 0 &&
                        summary.total >=
                            static_cast<std::uint64_t>(opt.max_cases))
                        goto done;
                    eval_pair(*out, summary, ctx, opt, bits, k, dist(rng),
                              dist(rng));
                }
            }
        }

    done:
        std::cerr << "PRUNED_FAILURE_DECOMP_SUMMARY total=" << summary.total
                  << " failures=" << summary.failures
                  << " bitextract_failures=" << summary.bitextract_failures
                  << " conversion_failures=" << summary.conversion_failures
                  << " final_gap_failures=" << summary.final_gap_failures
                  << " wrapper_failures=" << summary.wrapper_failures
                  << "\n";
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: failure_decomp\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr
            << "  ./build-ethmsb-release/my_ethmsb_pruned_failure_decomp\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
