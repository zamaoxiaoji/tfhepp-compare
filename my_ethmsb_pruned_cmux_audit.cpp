#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <set>
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

Context make_context(const std::uint64_t seed = 0)
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

std::vector<int> k_choices_for_p(const int p)
{
    std::set<int> ks;
    if (p > 1) ks.insert(1);
    for (const int k : {p - 5, p - 4, p - 3, p - 2, p - 1})
        if (k >= 1 && k <= p - 1) ks.insert(k);
    return {ks.begin(), ks.end()};
}

std::vector<Torus> bitextract_target_cases(const int p, const int k)
{
    std::set<Torus> vals = {0,
                            1,
                            (Torus{1} << (p - 1)) - 1,
                            Torus{1} << (p - 1),
                            (Torus{1} << (p - 1)) + 1,
                            (Torus{1} << p) - 1};
    const Torus w = bit_weight_msb_index(p, k);
    for (const Torus v :
         {w - 1, w, w + 1, 2 * w - 1, 2 * w, 2 * w + 1,
          (Torus{1} << (p - 1)) - w - 1,
          (Torus{1} << (p - 1)) - w, (Torus{1} << (p - 1)) - 1,
          Torus{1} << (p - 1)}) {
        if (v < (Torus{1} << p)) vals.insert(v);
    }
    return {vals.begin(), vals.end()};
}

TFHEpp::TLWE<P2> encrypt_m(const Torus m, const int p, const Context& ctx,
                           const bool encrypted)
{
    const Torus phase = encode_unsigned_phase(m, p);
    if (!encrypted) return trivial_encrypt_phase_l2(phase);
    return TFHEpp::tlweSymEncrypt<P2>(phase, ctx.key2);
}

bool output_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_qhalf(phase_l2(ct, ctx.key2));
}

void print_need_info(const std::string& stage, const std::string& what,
                     const int p, const int k, const Torus m,
                     const bool expected, const bool actual,
                     const Torus phase)
{
    const BitIndexInfo info = bit_index_info(p, k);
    std::cerr << "BEGIN_NEED_INFO\n";
    std::cerr << "stage: " << stage << "\n";
    std::cerr << "what_failed: " << what << "\n";
    std::cerr << "commands_run:\n";
    std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_cmux_audit\n";
    std::cerr << "first_failure:\n";
    std::cerr << "  p=" << p << "\n";
    std::cerr << "  bit_index_k=" << k << "\n";
    std::cerr << "  w_k=" << info.w_k << "\n";
    std::cerr << "  M_k=" << info.M_k << "\n";
    std::cerr << "  rotation_period_index=" << info.rotation_period_index
              << "\n";
    std::cerr << "  m_or_a_b=" << m << "\n";
    std::cerr << "  expected=" << expected << "\n";
    std::cerr << "  actual=" << actual << "\n";
    std::cerr << "  final_output_phase=" << hex64(phase) << "\n";
    std::cerr << "  failure_class="
              << (info.exact_in_br_index ? "decode_or_lut_orientation"
                                          : "BitExtract precision")
              << "\n";
    std::cerr << "END_NEED_INFO\n";
}

int run_periodicity_audit()
{
    int failures = 0;
    for (const int p : {4, 5, 6, 8, 9, 16, 17, 25, 33}) {
        for (const int k : k_choices_for_p(p)) {
            const BitIndexInfo info = bit_index_info(p, k);
            std::cout << "PERIODICITY_AUDIT p=" << p << " k=" << k
                      << " w_k=" << info.w_k << " M_k=" << info.M_k
                      << " twoN_mod_M=" << (BR_CYCLE % info.M_k)
                      << " rotation_period_index="
                      << info.rotation_period_index
                      << " supports_period_pruning="
                      << info.supports_period_pruning
                      << " expected_prune_ratio=1/" << info.M_k
                      << " disabled_reason=" << info.disabled_reason << "\n";

            const auto poly = make_bitextract_qhalf_poly(p, k);
            if (!info.supports_period_pruning) continue;
            const std::uint64_t limit =
                std::min<std::uint64_t>(BR_CYCLE - info.rotation_period_index,
                                        16 * info.rotation_period_index);
            for (std::uint64_t a = 0; a <= limit;
                 a += info.rotation_period_index) {
                if (!polynomial_rotation_equal(poly, a)) {
                    std::cerr << "BEGIN_NEED_INFO\n";
                    std::cerr << "stage: periodicity_audit\n";
                    std::cerr << "what_failed: safe rotation period was not invariant\n";
                    std::cerr << "commands_run:\n";
                    std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_cmux_audit\n";
                    std::cerr << "first_failure:\n";
                    std::cerr << "  p=" << p << "\n";
                    std::cerr << "  bit_index_k=" << k << "\n";
                    std::cerr << "  w_k=" << info.w_k << "\n";
                    std::cerr << "  M_k=" << info.M_k << "\n";
                    std::cerr << "  m_or_a_b=" << a << "\n";
                    std::cerr << "  failure_class=periodicity_audit\n";
                    std::cerr << "END_NEED_INFO\n";
                    return 1;
                }
            }
        }
    }
    return failures;
}

int run_ab_and_functional(const Context& ctx)
{
    int failures = 0;
    int targeted_failures = 0;
    for (int p = 2; p <= std::min(10, BR_INDEX_BITS); ++p) {
        for (int k = 1; k <= p - 1; ++k) {
            for (Torus m = 0; m < (Torus{1} << p); ++m) {
                const auto ct = encrypt_m(m, p, ctx, false);
                TFHEpp::TLWE<P2> official;
                TFHEpp::TLWE<P2> no_periodic;
                TFHEpp::TLWE<P2> pruned;
                PruneStats stats;
                bit_extract_unpruned_reference_qhalf_l2_to_l2(
                    official, ct, p, k, *ctx.iksk, *ctx.bkfft);
                bit_extract_no_periodic_qhalf_l2_to_l2(
                    no_periodic, ct, p, k, *ctx.iksk, *ctx.bkfft, nullptr);
                bit_extract_pruned_qhalf_l2_to_l2(
                    pruned, ct, p, k, *ctx.iksk, *ctx.bkfft, &stats);
                if (official != no_periodic) {
                    std::cerr << "BEGIN_NEED_INFO\n";
                    std::cerr << "stage: pruned_blindrotate\n";
                    std::cerr << "what_failed: no-periodic pruned path differs from official GateBootstrapping path\n";
                    std::cerr << "commands_run:\n";
                    std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_cmux_audit\n";
                    std::cerr << "first_failure:\n";
                    std::cerr << "  p=" << p << "\n";
                    std::cerr << "  bit_index_k=" << k << "\n";
                    std::cerr << "  m_or_a_b=" << m << "\n";
                    std::cerr << "  failure_class=api_or_blindrotate_mismatch\n";
                    std::cerr << "END_NEED_INFO\n";
                    return 1;
                }
                const bool expected = ((m >> (p - 1 - k)) & 1) != 0;
                const bool actual = output_bit(pruned, ctx);
                if (expected != actual) {
                    print_need_info("bitextract_functional",
                                    "small-p BitExtract mismatch", p, k, m,
                                    expected, actual,
                                    phase_l2(pruned, ctx.key2));
                    return 1;
                }
            }
        }
    }

    for (const int p : {9, 17, 25, 33}) {
        for (const int k : k_choices_for_p(p)) {
            const BitIndexInfo info = bit_index_info(p, k);
            for (const Torus m : bitextract_target_cases(p, k)) {
                const auto ct = encrypt_m(m, p, ctx, true);
                TFHEpp::TLWE<P2> out;
                PruneStats stats;
                bit_extract_pruned_qhalf_l2_to_l2(
                    out, ct, p, k, *ctx.iksk, *ctx.bkfft, &stats);
                const bool expected = ((m >> (p - 1 - k)) & 1) != 0;
                const bool actual = output_bit(out, ctx);
                if (expected != actual) {
                    ++targeted_failures;
                    std::cout
                        << "BITEXTRACT_TARGETED_FAILURE p=" << p
                        << " k=" << k << " m=" << m
                        << " expected=" << expected << " actual=" << actual
                        << " phase=" << hex64(phase_l2(out, ctx.key2))
                        << " failure_class="
                        << (info.exact_in_br_index
                                ? "BitExtract switch-band limitation"
                                : "BitExtract precision limitation")
                        << " disabled_reason=" << info.disabled_reason
                        << "\n";
                    break;
                }
            }
        }
    }
    std::cout << "PRUNED_CMUX_AUDIT_RESULT small_p_failures=" << failures
              << " encrypted_targeted_failures=" << targeted_failures
              << "\n";
    return failures;
}

}  // namespace

int main()
{
    try {
        const Context ctx = make_context();
        const int periodic_failures = run_periodicity_audit();
        if (periodic_failures != 0) return 1;
        const int functional_failures = run_ab_and_functional(ctx);
        return functional_failures == 0 ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: cmux_audit\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n";
        std::cerr << "  ./build-ethmsb-release/my_ethmsb_pruned_cmux_audit\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
