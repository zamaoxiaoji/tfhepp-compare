#include <algorithm>
#include <cstdint>
#include <iostream>
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
    int max_trivial_p = 10;
    bool encrypted_targeted = true;
    bool comparison = true;
    std::uint64_t seed = 0;
};

struct Failure {
    std::string stage = "none";
    std::string failure_class = "none";
    int p = 0;
    int k = 0;
    Torus w = 0;
    Torus M = 0;
    Torus m = 0;
    Torus a = 0;
    Torus b = 0;
    bool expected = false;
    bool actual = false;
    Torus phase_in = 0;
    Torus bit_qhalf_phase = 0;
    Torus guard_value = 0;
    Torus guard_mask_phase = 0;
    Torus gap_offset = 0;
    Torus final_output_phase = 0;
    PruneStats prune;
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
        if (a == "--max-trivial-p")
            opt.max_trivial_p = std::stoi(need("--max-trivial-p"));
        else if (a == "--no-encrypted-targeted")
            opt.encrypted_targeted = false;
        else if (a == "--no-comparison")
            opt.comparison = false;
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    return opt;
}

std::vector<int> k_choices_for_p(const int p)
{
    std::set<int> ks;
    if (p > 1) ks.insert(1);
    for (const int k : {std::max(1, p - 5), p - 4, p - 3, p - 2, p - 1})
        if (k >= 1 && k <= p - 1) ks.insert(k);
    return {ks.begin(), ks.end()};
}

std::vector<Torus> gap_target_cases(const int p, const int k)
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

std::set<std::pair<Torus, Torus>> comparison_pairs(const int t)
{
    const Torus maxv = (Torus{1} << t) - 1;
    const Torus mid = Torus{1} << (t - 1);
    std::set<std::pair<Torus, Torus>> pairs = {
        {0, 0}, {0, 1}, {1, 0}, {0, maxv}, {maxv, 0},
        {mid - 1, mid}, {mid, mid - 1}, {maxv, maxv}};
    const std::vector<Torus> xs = {Torus{0}, Torus{1}, Torus{2},
                                   mid - 2, mid - 1, mid, mid + 1,
                                   maxv - 2, maxv - 1, maxv};
    for (const Torus x : xs) {
        if (x > maxv) continue;
        pairs.insert({x, x});
        if (x + 1 <= maxv) {
            pairs.insert({x, x + 1});
            pairs.insert({x + 1, x});
        }
        pairs.insert({0, x});
        pairs.insert({x, 0});
        pairs.insert({maxv, x});
        pairs.insert({x, maxv});
    }
    return pairs;
}

TFHEpp::TLWE<P2> make_input(const Torus m, const int p, const Context& ctx,
                            const bool encrypted)
{
    const Torus ph = encode_unsigned_phase(m, p);
    if (!encrypted) return trivial_encrypt_phase_l2(ph);
    return TFHEpp::tlweSymEncrypt<P2>(ph, ctx.key2);
}

bool arith_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx,
               const Torus out_value = BOOL_ONE)
{
    return closer_to_value(phase_l2(ct, ctx.key2), out_value);
}

bool qhalf_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_qhalf(phase_l2(ct, ctx.key2));
}

int helper_formula_tests()
{
    int failures = 0;
    for (const int p : {9, 17, 25, 33}) {
        for (const int k : k_choices_for_p(p)) {
            const Torus w = bit_weight_msb_index(p, k);
            const Torus M = bit_period_msb_index(p, k);
            const Torus guard = guard_value_for_bit_index(p, k);
            const Torus offset = gap_offset_for_bit_index(p, k);
            const Torus gap_width = static_cast<Torus>(
                Wide{w + 1} * Wide{my_ethmsb::delta_bits(p)});
            std::cout << "HELPER_FORMULA p=" << p << " k=" << k
                      << " w_k=" << w << " M_k=" << M
                      << " guard_value=" << hex64(guard)
                      << " gap_offset=" << hex64(offset)
                      << " gap_width=" << hex64(gap_width)
                      << " expected_prune_ratio=1/" << M << "\n";
            if (p == 33 && k == p - 5) {
                failures += w != 16;
                failures += M != 32;
                failures += guard != (Torus{1} << 35);
                failures += offset != (Torus{17} << 30);
            }
        }
    }
    return failures;
}

int conversion_tests(const Context& ctx, Failure& first)
{
    int failures = 0;
    for (const int p : {9, 17, 25, 33}) {
        for (const int k : k_choices_for_p(p)) {
            const Torus guard = guard_value_for_bit_index(p, k);
            for (const Torus in_phase : {Torus{0}, Q_HALF}) {
                TFHEpp::TLWE<P2> out;
                bool_qhalf_to_value_l2_to_l2(
                    out, trivial_encrypt_phase_l2(in_phase), guard,
                    *ctx.iksk, *ctx.bkfft);
                const bool expected = in_phase == Q_HALF;
                const bool actual = arith_bit(out, ctx, guard);
                if (expected != actual) {
                    ++failures;
                    if (first.stage == "none") {
                        first.stage = "conversion";
                        first.failure_class = "conversion";
                        first.p = p;
                        first.k = k;
                        first.w = bit_weight_msb_index(p, k);
                        first.M = bit_period_msb_index(p, k);
                        first.expected = expected;
                        first.actual = actual;
                        first.bit_qhalf_phase = in_phase;
                        first.guard_value = guard;
                        first.guard_mask_phase = phase_l2(out, ctx.key2);
                    }
                    return failures;
                }
            }
        }
    }
    return failures;
}

int bitextract_small_tests(const Context& ctx, Failure& first)
{
    int failures = 0;
    for (int p = 2; p <= std::min(10, BR_INDEX_BITS); ++p) {
        for (int k = 1; k <= p - 1; ++k) {
            for (Torus m = 0; m < (Torus{1} << p); ++m) {
                TFHEpp::TLWE<P2> out;
                PruneStats stats;
                bit_extract_pruned_qhalf_l2_to_l2(
                    out, make_input(m, p, ctx, false), p, k, *ctx.iksk,
                    *ctx.bkfft, &stats);
                const bool expected = ((m >> (p - 1 - k)) & 1) != 0;
                const bool actual = qhalf_bit(out, ctx);
                if (expected != actual) {
                    ++failures;
                    if (first.stage == "none") {
                        first.stage = "bitextract_functional";
                        first.failure_class =
                            stats.bit_exact_in_br_index
                                ? "decode"
                                : "BitExtract precision";
                        first.p = p;
                        first.k = k;
                        first.w = stats.w_k;
                        first.M = stats.M_k;
                        first.m = m;
                        first.expected = expected;
                        first.actual = actual;
                        first.phase_in =
                            phase_l2(make_input(m, p, ctx, false), ctx.key2);
                        first.bit_qhalf_phase = phase_l2(out, ctx.key2);
                        first.prune = stats;
                    }
                    return failures;
                }
            }
        }
    }
    return failures;
}

int gapmsb_case(const Context& ctx, const int p, const int k, const Torus m,
                const bool encrypted, Failure& first)
{
    const auto ct = make_input(m, p, ctx, encrypted);
    const bool expected = ((m >> (p - 1)) & 1) != 0;

    TFHEpp::TLWE<P2> bit;
    PruneStats bit_stats;
    bit_extract_pruned_qhalf_l2_to_l2(bit, ct, p, k, *ctx.iksk, *ctx.bkfft,
                                      &bit_stats);
    const bool expected_bit = ((m >> (p - 1 - k)) & 1) != 0;
    const bool actual_bit = qhalf_bit(bit, ctx);

    const Torus guard = guard_value_for_bit_index(p, k);
    TFHEpp::TLWE<P2> mask;
    bool_qhalf_to_value_l2_to_l2(mask, bit, guard, *ctx.iksk, *ctx.bkfft);
    const bool actual_mask = arith_bit(mask, ctx, guard);

    TFHEpp::TLWE<P2> out;
    PrunedGapStats stats;
    gapmsb_pruned_bitk_l2_to_l2(out, ct, p, k, BOOL_ONE, *ctx.iksk,
                                *ctx.bkfft, &stats);
    const bool actual = arith_bit(out, ctx, BOOL_ONE);
    if (expected == actual) return 0;

    if (first.stage == "none") {
        first.stage = encrypted ? "gapmsb_encrypted_targeted"
                                : "gapmsb_trivial_exhaustive";
        if (expected_bit != actual_bit)
            first.failure_class = bit_stats.bit_exact_in_br_index
                                      ? "decode"
                                      : "BitExtract precision";
        else if (expected_bit != actual_mask)
            first.failure_class = "conversion";
        else
            first.failure_class = "gap margin";
        first.p = p;
        first.k = k;
        first.w = bit_weight_msb_index(p, k);
        first.M = bit_period_msb_index(p, k);
        first.m = m;
        first.expected = expected;
        first.actual = actual;
        first.phase_in = phase_l2(ct, ctx.key2);
        first.bit_qhalf_phase = phase_l2(bit, ctx.key2);
        first.guard_value = guard;
        first.guard_mask_phase = phase_l2(mask, ctx.key2);
        first.gap_offset = gap_offset_for_bit_index(p, k);
        first.final_output_phase = phase_l2(out, ctx.key2);
        first.prune = stats.bitextract_stats;
        if (expected_bit != actual_bit && encrypted &&
            bit_stats.bit_exact_in_br_index)
            first.failure_class = "BitExtract switch-band";
    }
    return 1;
}

int gapmsb_trivial_tests(const Context& ctx, const int max_p, Failure& first)
{
    int failures = 0;
    for (int p = 2; p <= std::min(max_p, 10); ++p) {
        for (int k = 1; k <= p - 1; ++k) {
            for (Torus m = 0; m < (Torus{1} << p); ++m) {
                failures += gapmsb_case(ctx, p, k, m, false, first);
                if (failures) return failures;
            }
        }
    }
    return failures;
}

int gapmsb_encrypted_tests(const Context& ctx, Failure& first)
{
    int failures = 0;
    for (const int p : {9, 17, 25, 33}) {
        for (const int k : k_choices_for_p(p)) {
            for (const Torus m : gap_target_cases(p, k)) {
                failures += gapmsb_case(ctx, p, k, m, true, first);
                if (failures) {
                    const BitIndexInfo info = bit_index_info(p, k);
                    if (p > BR_INDEX_BITS || !info.exact_in_br_index) {
                        std::cout << "GAPMSB_EXPECTED_LIMITATION p=" << p
                                  << " k=" << k << " m=" << m
                                  << " failure_class="
                                  << first.failure_class << "\n";
                        return failures;
                    }
                    return failures;
                }
            }
        }
    }
    return failures;
}

TFHEpp::TLWE<P2> encrypt_compare_value(const Torus x, const int t,
                                       const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned_phase(x, t + 1),
                                      ctx.key2);
}

int comparison_tests(const Context& ctx, Failure& first)
{
    int failures = 0;
    for (const int t : {8, 16, 24, 32}) {
        const int p = t + 1;
        for (const int k : {std::max(1, p - 5),
                            std::max(1, p - 4), std::max(1, p - 3)}) {
            for (const auto& [a0, b0] : comparison_pairs(t)) {
                const auto a = encrypt_compare_value(a0, t, ctx);
                const auto b = encrypt_compare_value(b0, t, ctx);
                TFHEpp::TLWE<P2> lt;
                TFHEpp::TLWE<P2> gt;
                TFHEpp::TLWE<P2> eq;
                homcomp_lt_pruned_gap(lt, a, b, t, k, *ctx.iksk,
                                      *ctx.bkfft);
                homcomp_gt_pruned_gap(gt, a, b, t, k, *ctx.iksk,
                                      *ctx.bkfft);
                homcomp_eq_pruned_gap(eq, a, b, t, k, *ctx.iksk,
                                      *ctx.bkfft);
                const bool got_lt = arith_bit(lt, ctx);
                const bool got_gt = arith_bit(gt, ctx);
                const bool got_eq = arith_bit(eq, ctx);
                const bool pass = got_lt == (a0 < b0) &&
                                  got_gt == (a0 > b0) &&
                                  got_eq == (a0 == b0);
                if (!pass) {
                    ++failures;
                    if (first.stage == "none") {
                        std::string failed_op = "lt";
                        bool expected = a0 < b0;
                        bool actual = got_lt;
                        Torus out_phase = phase_l2(lt, ctx.key2);
                        if (got_lt == (a0 < b0) &&
                            got_gt != (a0 > b0)) {
                            failed_op = "gt";
                            expected = a0 > b0;
                            actual = got_gt;
                            out_phase = phase_l2(gt, ctx.key2);
                        }
                        else if (got_lt == (a0 < b0) &&
                                 got_gt == (a0 > b0) &&
                                 got_eq != (a0 == b0)) {
                            failed_op = "eq";
                            expected = a0 == b0;
                            actual = got_eq;
                            out_phase = phase_l2(eq, ctx.key2);
                        }
                        first.stage = "comparison";
                        first.failure_class =
                            bit_index_info(p, k).exact_in_br_index
                                ? "gap margin"
                                : "BitExtract precision";
                        first.p = p;
                        first.k = k;
                        first.w = bit_weight_msb_index(p, k);
                        first.M = bit_period_msb_index(p, k);
                        first.a = a0;
                        first.b = b0;
                        first.expected = expected;
                        first.actual = actual;
                        first.final_output_phase = out_phase;
                        first.failure_class += ":" + failed_op;
                    }
                    return failures;
                }
            }
        }
    }
    return failures;
}

void print_need_info(const Failure& f)
{
    std::cerr << "BEGIN_NEED_INFO\n";
    std::cerr << "stage: " << f.stage << "\n";
    std::cerr << "what_failed: pruned GapMSB correctness failure\n";
    std::cerr << "commands_run:\n";
    std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_gapmsb_tests\n";
    std::cerr << "first_failure:\n";
    std::cerr << "  p=" << f.p << "\n";
    std::cerr << "  bit_index_k=" << f.k << "\n";
    std::cerr << "  w_k=" << f.w << "\n";
    std::cerr << "  M_k=" << f.M << "\n";
    if (f.stage == "comparison")
        std::cerr << "  m_or_a_b=(" << f.a << "," << f.b << ")\n";
    else
        std::cerr << "  m_or_a_b=" << f.m << "\n";
    std::cerr << "  expected=" << f.expected << "\n";
    std::cerr << "  actual=" << f.actual << "\n";
    std::cerr << "  phase_in=" << hex64(f.phase_in) << "\n";
    std::cerr << "  bit_qhalf_phase=" << hex64(f.bit_qhalf_phase) << "\n";
    std::cerr << "  guard_value=" << hex64(f.guard_value) << "\n";
    std::cerr << "  guard_mask_phase=" << hex64(f.guard_mask_phase) << "\n";
    std::cerr << "  gap_offset=" << hex64(f.gap_offset) << "\n";
    std::cerr << "  final_output_phase=" << hex64(f.final_output_phase)
              << "\n";
    std::cerr << "  failure_class=" << f.failure_class << "\n";
    std::cerr << "pruning_stats:\n";
    std::cerr << "  total_terms=" << f.prune.total_terms << "\n";
    std::cerr << "  periodic_skipped=" << f.prune.periodic_skipped << "\n";
    std::cerr << "  zero_skipped=" << f.prune.zero_skipped << "\n";
    std::cerr << "  executed=" << f.prune.executed_terms << "\n";
    std::cerr << "  pruning_enabled=" << f.prune.pruning_enabled << "\n";
    std::cerr << "  disabled_reason=" << f.prune.disabled_reason << "\n";
    std::cerr << "END_NEED_INFO\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const Context ctx = make_context(opt.seed);
        Failure first;
        const int helper_failures = helper_formula_tests();
        const int conversion_failures = conversion_tests(ctx, first);
        const int bitextract_failures = first.stage == "none"
                                            ? bitextract_small_tests(ctx, first)
                                            : 0;
        const int trivial_failures =
            first.stage == "none"
                ? gapmsb_trivial_tests(ctx, opt.max_trivial_p, first)
                : 0;
        const int encrypted_failures =
            (first.stage == "none" && opt.encrypted_targeted)
                ? gapmsb_encrypted_tests(ctx, first)
                : 0;
        const int comparison_failures =
            (first.stage == "none" && opt.comparison)
                ? comparison_tests(ctx, first)
                : 0;

        std::cout << "PRUNED_GAPMSB_TEST_RESULT helper_formula_failures="
                  << helper_failures
                  << " conversion_failures=" << conversion_failures
                  << " bitextract_small_p_failures=" << bitextract_failures
                  << " gapmsb_trivial_exhaustive_failures="
                  << trivial_failures
                  << " gapmsb_encrypted_targeted_failures="
                  << encrypted_failures
                  << " comparison_adversarial_failures="
                  << comparison_failures << "\n";

        const int failures = helper_failures + conversion_failures +
                             bitextract_failures + trivial_failures +
                             encrypted_failures + comparison_failures;
        if (failures != 0) {
            print_need_info(first);
            return 1;
        }
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: gapmsb_correctness\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_gapmsb_tests\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
