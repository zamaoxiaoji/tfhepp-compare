#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
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
    std::string k_range = "safe";
    int trials_random = 100;
    int trials_boundary = 5;
    std::uint64_t seed = 0;
    double target_failure_rate = 1e-4;
    std::string csv_path;
    std::string json_path;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct Counts {
    std::uint64_t total = 0;
    std::uint64_t failures = 0;
    std::uint64_t bitextract_failures = 0;
    std::uint64_t conversion_failures = 0;
    std::uint64_t final_gap_failures = 0;
    std::uint64_t periodic_sum = 0;
    std::uint64_t zero_sum = 0;
    std::uint64_t executed_sum = 0;
};

struct Row {
    int bits = 0;
    int p = 0;
    int k = 0;
    std::uint64_t w = 0;
    std::uint64_t M = 0;
    Torus gap_margin = 0;
    Torus switch_band = DEFAULT_SWITCH_BAND_WIDTH;
    Torus safe_margin = 0;
    double expected_prune_ratio = 0.0;
    double observed_prune_ratio_mean = 0.0;
    Counts boundary;
    Counts random;
    double boundary_failure_rate = 0.0;
    double random_failure_rate = 0.0;
    double bitextract_failure_rate = 0.0;
    double conversion_failure_rate = 0.0;
    double final_gap_failure_rate = 0.0;
    double random_ci_low = 0.0;
    double random_ci_high = 1.0;
    bool selected = false;
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
        else if (a == "--all")
            opt.bits = {8, 16, 24, 32};
        else if (a == "--k-range")
            opt.k_range = need("--k-range");
        else if (a == "--trials-random")
            opt.trials_random = std::stoi(need("--trials-random"));
        else if (a == "--trials-boundary")
            opt.trials_boundary = std::stoi(need("--trials-boundary"));
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
        else if (a == "--target-failure-rate")
            opt.target_failure_rate = std::stod(need("--target-failure-rate"));
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else if (a == "--json")
            opt.json_path = need("--json");
        else
            throw std::runtime_error("unknown argument: " + a);
    }
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

TFHEpp::TLWE<P2> encrypt_phase_m(const Torus m, const int p,
                                 const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned_phase(m, p), ctx.key2);
}

std::vector<int> k_values_for(const int p, const std::string& range)
{
    std::set<int> ks;
    if (range == "all" || range == "safe") {
        for (int k = 1; k <= std::min(p - 1, 11); ++k) ks.insert(k);
    }
    else {
        const std::size_t colon = range.find(':');
        if (colon == std::string::npos)
            throw std::runtime_error("unsupported --k-range");
        const int lo = std::stoi(range.substr(0, colon));
        const int hi = std::stoi(range.substr(colon + 1));
        for (int k = lo; k <= hi && k <= p - 1; ++k)
            if (k >= 1) ks.insert(k);
    }
    const int pminus5 = std::max(1, p - 5);
    if (pminus5 >= 1 && pminus5 <= p - 1) ks.insert(pminus5);
    if (range == "safe") {
        for (auto it = ks.begin(); it != ks.end();) {
            const BitIndexInfo info = bit_index_info(p, *it);
            const Torus gap = gap_offset_for_bit_index(p, *it);
            if (!info.exact_in_br_index || !info.supports_period_pruning ||
                gap <= Torus{4} * DEFAULT_SWITCH_BAND_WIDTH)
                it = ks.erase(it);
            else
                ++it;
        }
    }
    return {ks.begin(), ks.end()};
}

std::pair<double, double> wilson95(const std::uint64_t failures,
                                   const std::uint64_t total)
{
    if (total == 0) return {0.0, 1.0};
    constexpr double z = 1.959963984540054;
    const double n = static_cast<double>(total);
    const double phat = static_cast<double>(failures) / n;
    const double denom = 1.0 + z * z / n;
    const double center = (phat + z * z / (2.0 * n)) / denom;
    const double half =
        z * std::sqrt((phat * (1.0 - phat) + z * z / (4.0 * n)) / n) /
        denom;
    return {std::max(0.0, center - half), std::min(1.0, center + half)};
}

std::vector<Torus> boundary_ms(const int p, const int k)
{
    std::set<Torus> vals = {0,
                            1,
                            (Torus{1} << (p - 1)) - 1,
                            Torus{1} << (p - 1),
                            (Torus{1} << (p - 1)) + 1,
                            (Torus{1} << p) - 1};
    const Torus w = bit_weight_msb_index(p, k);
    const Torus half = Torus{1} << (p - 1);
    const Torus top = Torus{1} << p;
    const Torus mid_r = half / w;
    const Torus top_r = top / w;
    std::set<Torus> rs = {0, 1, 2, 3};
    for (int off = -2; off <= 2; ++off) {
        const std::int64_t r = static_cast<std::int64_t>(mid_r) + off;
        if (r >= 0) rs.insert(static_cast<Torus>(r));
    }
    for (int off = -4; off <= -1; ++off) {
        const std::int64_t r = static_cast<std::int64_t>(top_r) + off;
        if (r >= 0) rs.insert(static_cast<Torus>(r));
    }
    for (const Torus r : rs) {
        for (int d = -2; d <= 2; ++d) {
            const std::int64_t v =
                static_cast<std::int64_t>(r * w) + static_cast<std::int64_t>(d);
            if (v >= 0 && static_cast<Torus>(v) < top)
                vals.insert(static_cast<Torus>(v));
        }
    }
    for (const Torus v :
         {half - 1 - w, half - w, half - 1, half, half + 1,
          top - 1 - w, top - 1}) {
        if (v < top) vals.insert(v);
    }
    return {vals.begin(), vals.end()};
}

std::vector<std::pair<Torus, Torus>> adversarial_pairs(const int bits)
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

std::string eval_diff(const Context& ctx, const TFHEpp::TLWE<P2>& diff,
                      const Torus m, const int p, const int k,
                      Counts& counts)
{
    const bool expected_msb = ((m >> (p - 1)) & 1) != 0;
    const bool expected_bit = ((m >> (p - 1 - k)) & 1) != 0;
    const Torus guard_value = guard_value_for_bit_index(p, k);
    const Torus gap_offset = gap_offset_for_bit_index(p, k);

    TFHEpp::TLWE<P2> bit;
    PruneStats prune;
    bit_extract_pruned_qhalf_l2_to_l2(bit, diff, p, k, *ctx.iksk,
                                      *ctx.bkfft, &prune);
    const bool bit_actual = closer_to_qhalf(phase_l2(bit, ctx.key2));
    counts.periodic_sum += prune.periodic_skipped;
    counts.zero_sum += prune.zero_skipped;
    counts.executed_sum += prune.executed_terms;

    ++counts.total;
    if (bit_actual != expected_bit) {
        ++counts.failures;
        ++counts.bitextract_failures;
        return "bitextract_boundary";
    }

    TFHEpp::TLWE<P2> mask;
    bool_qhalf_to_value_l2_to_l2(mask, bit, guard_value, *ctx.iksk,
                                 *ctx.bkfft);
    const bool guard_actual =
        closer_to_value(phase_l2(mask, ctx.key2), guard_value);
    if (guard_actual != expected_bit) {
        ++counts.failures;
        ++counts.conversion_failures;
        return "conversion";
    }

    TFHEpp::TLWE<P2> gap;
    my_ethmsb::sub<P2>(gap, diff, mask);
    TFHEpp::TLWE<P2> out;
    pbs_msb_value_l2_to_l2(out, gap, p, gap_offset, BOOL_ONE, *ctx.iksk,
                           *ctx.bkfft);
    const bool final_actual = closer_to_value(phase_l2(out, ctx.key2), BOOL_ONE);
    if (final_actual != expected_msb) {
        ++counts.failures;
        ++counts.final_gap_failures;
        return "final_gap_margin";
    }
    return "none";
}

void eval_m_case(const Context& ctx, const Torus m, const int p, const int k,
                 Counts& counts)
{
    (void)eval_diff(ctx, encrypt_phase_m(m, p, ctx), m, p, k, counts);
}

void eval_pair_case(const Context& ctx, const int bits, const int k,
                    const Torus a0, const Torus b0, Counts& counts)
{
    const int p = bits + 1;
    const auto a = encrypt_phase_m(a0, p, ctx);
    const auto b = encrypt_phase_m(b0, p, ctx);
    TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, a, b);
    const Torus m = (a0 - b0) & ((Torus{1} << p) - 1);
    (void)eval_diff(ctx, diff, m, p, k, counts);
}

double ratio(const std::uint64_t n, const std::uint64_t d)
{
    return d == 0 ? 0.0 : static_cast<double>(n) / static_cast<double>(d);
}

void finalize_row(Row& row)
{
    row.boundary_failure_rate =
        ratio(row.boundary.failures, row.boundary.total);
    row.random_failure_rate = ratio(row.random.failures, row.random.total);
    const std::uint64_t all_total = row.boundary.total + row.random.total;
    const std::uint64_t bit_fail =
        row.boundary.bitextract_failures + row.random.bitextract_failures;
    const std::uint64_t conv_fail =
        row.boundary.conversion_failures + row.random.conversion_failures;
    const std::uint64_t gap_fail =
        row.boundary.final_gap_failures + row.random.final_gap_failures;
    row.bitextract_failure_rate = ratio(bit_fail, all_total);
    row.conversion_failure_rate = ratio(conv_fail, all_total);
    row.final_gap_failure_rate = ratio(gap_fail, all_total);
    const auto [lo, hi] = wilson95(row.random.failures, row.random.total);
    row.random_ci_low = lo;
    row.random_ci_high = hi;
    const double skipped = static_cast<double>(
        row.boundary.periodic_sum + row.boundary.zero_sum +
        row.random.periodic_sum + row.random.zero_sum);
    const double terms = skipped +
                         static_cast<double>(row.boundary.executed_sum +
                                             row.random.executed_sum);
    row.observed_prune_ratio_mean = terms == 0.0 ? 0.0 : skipped / terms;
}

void choose_selected(std::vector<Row>& rows, const double target)
{
    for (const int bits : {8, 16, 24, 32}) {
        Row* best = nullptr;
        for (Row& r : rows) {
            if (r.bits != bits) continue;
            if (r.boundary_failure_rate > target) continue;
            if (r.random_ci_high > target) continue;
            if (best == nullptr ||
                r.observed_prune_ratio_mean > best->observed_prune_ratio_mean)
                best = &r;
        }
        if (best != nullptr) best->selected = true;
    }
}

void emit_csv(std::ostream& out, const std::vector<Row>& rows)
{
    out << "bits,p,k,w,M,gap_margin_hex,measured_switch_band_hex,"
           "safe_margin_hex,expected_prune_ratio,observed_prune_ratio_mean,"
           "boundary_cases,boundary_failures,boundary_failure_rate,"
           "random_trials,random_failures,random_failure_rate,"
           "bitextract_failure_rate,conversion_failure_rate,"
           "final_gap_failure_rate,selected_default_candidate,"
           "random_wilson_low,random_wilson_high\n";
    for (const Row& r : rows) {
        out << r.bits << ',' << r.p << ',' << r.k << ',' << r.w << ','
            << r.M << ',' << hex64(r.gap_margin) << ','
            << hex64(r.switch_band) << ',' << hex64(r.safe_margin) << ','
            << r.expected_prune_ratio << ',' << r.observed_prune_ratio_mean
            << ',' << r.boundary.total << ',' << r.boundary.failures << ','
            << r.boundary_failure_rate << ',' << r.random.total << ','
            << r.random.failures << ',' << r.random_failure_rate << ','
            << r.bitextract_failure_rate << ',' << r.conversion_failure_rate
            << ',' << r.final_gap_failure_rate << ',' << r.selected << ','
            << r.random_ci_low << ',' << r.random_ci_high << "\n";
    }
}

void emit_json(const std::string& path, const std::vector<Row>& rows,
               const Options& opt)
{
    if (path.empty()) return;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open json output: " + path);
    out << "{\n  \"impl\":\"" << SINGLE_ROUND_IMPL_NAME << "\",\n";
    out << "  \"target_failure_rate\":" << opt.target_failure_rate << ",\n";
    out << "  \"selected_candidates\":[\n";
    bool first = true;
    for (const int bits : opt.bits) {
        const Row* selected = nullptr;
        for (const Row& r : rows)
            if (r.bits == bits && r.selected) selected = &r;
        if (!first) out << ",\n";
        first = false;
        out << "    {\"bits\":" << bits << ",\"p\":" << (bits + 1)
            << ",\"selected_k\":";
        if (selected == nullptr)
            out << "null";
        else
            out << selected->k;
        out << ",\"random_rate\":"
            << (selected == nullptr ? 0.0 : selected->random_failure_rate)
            << ",\"random_upper95\":"
            << (selected == nullptr ? 1.0 : selected->random_ci_high)
            << ",\"boundary_rate\":"
            << (selected == nullptr ? 1.0
                                    : selected->boundary_failure_rate)
            << ",\"prune_ratio\":"
            << (selected == nullptr ? 0.0
                                    : selected->observed_prune_ratio_mean)
            << "}";
    }
    out << "\n  ]\n}\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const Context ctx = make_context(opt.seed);
        std::vector<Row> rows;
        std::mt19937_64 rng(opt.seed);
        for (const int bits : opt.bits) {
            const int p = bits + 1;
            const Torus maxv = (Torus{1} << bits) - 1;
            std::uniform_int_distribution<Torus> dist(0, maxv);
            for (const int k : k_values_for(p, opt.k_range)) {
                Row row;
                row.bits = bits;
                row.p = p;
                row.k = k;
                row.w = bit_weight_msb_index(p, k);
                row.M = bit_period_msb_index(p, k);
                row.gap_margin = gap_offset_for_bit_index(p, k);
                row.safe_margin =
                    row.gap_margin > row.switch_band
                        ? row.gap_margin - row.switch_band
                        : 0;
                const BitIndexInfo info = bit_index_info(p, k);
                row.expected_prune_ratio =
                    info.supports_period_pruning &&
                            info.rotation_period_index != 0
                        ? 1.0 /
                              static_cast<double>(info.rotation_period_index)
                        : 0.0;

                const int boundary_repeats = std::max(1, opt.trials_boundary);
                for (const Torus m : boundary_ms(p, k))
                    for (int r = 0; r < boundary_repeats; ++r)
                        eval_m_case(ctx, m, p, k, row.boundary);
                for (const auto& [a, b] : adversarial_pairs(bits))
                    for (int r = 0; r < boundary_repeats; ++r)
                        eval_pair_case(ctx, bits, k, a, b, row.boundary);
                for (int i = 0; i < opt.trials_random; ++i)
                    eval_pair_case(ctx, bits, k, dist(rng), dist(rng),
                                   row.random);
                finalize_row(row);
                rows.push_back(row);
            }
        }
        choose_selected(rows, opt.target_failure_rate);

        std::ofstream csv_file;
        std::ostream* csv = &std::cout;
        if (!opt.csv_path.empty()) {
            csv_file.open(opt.csv_path);
            if (!csv_file)
                throw std::runtime_error("cannot open csv output: " +
                                         opt.csv_path);
            csv = &csv_file;
        }
        emit_csv(*csv, rows);
        emit_json(opt.json_path, rows, opt);

        for (const int bits : opt.bits) {
            const Row* selected = nullptr;
            for (const Row& r : rows)
                if (r.bits == bits && r.selected) selected = &r;
            std::cerr << "PRUNED_FAILURE_SWEEP_SELECTED bits=" << bits
                      << " k="
                      << (selected == nullptr ? std::string("none")
                                              : std::to_string(selected->k))
                      << "\n";
        }
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: failure_sweep\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr
            << "  ./build-ethmsb-release/my_ethmsb_pruned_failure_sweep\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
