#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "my_ethmsb_pruned_fused_final.hpp"

namespace {

using namespace my_ethmsb_pruned_fused_final;

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
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

TFHEpp::TLWE<P2> encrypt_m(const Torus m, const int p, const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned_phase(m, p), ctx.key2);
}

bool arith_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_value(phase_l2(ct, ctx.key2), BOOL_ONE);
}

bool expected_msb(const Torus m, const int p)
{
    return ((m >> (p - 1)) & 1) != 0;
}

bool pure_oracle_ok()
{
    for (int p = 2; p <= 10; ++p) {
        const Torus top = Torus{1} << p;
        for (int k = 1; k <= p - 1; ++k) {
            if (!fused_truth_antisymmetry_ok(p, k)) return false;
            for (Torus m = 0; m < top; ++m) {
                const Torus b = (m >> (p - 1 - k)) & 1;
                const Torus combined = (m + b * (Torus{1} << (p - 1))) &
                                       mask_for_bits(p);
                if (fused_final_truth_bit(p, k, combined) !=
                    expected_msb(m, p)) {
                    std::cerr << "pure_oracle_failure p=" << p << " k=" << k
                              << " m=" << m << " combined=" << combined
                              << "\n";
                    return false;
                }
            }
        }
    }
    return true;
}

std::vector<Torus> targeted_ms(const int p, const int k)
{
    std::set<Torus> vals = {0,
                            1,
                            (Torus{1} << (p - 1)) - 1,
                            Torus{1} << (p - 1),
                            (Torus{1} << (p - 1)) + 1,
                            (Torus{1} << p) - 1};
    const Torus W = bit_weight_msb_index(p, k);
    for (const Torus v : {W - 1, W, W + 1, 2 * W - 1, 2 * W, 2 * W + 1})
        if (v < (Torus{1} << p)) vals.insert(v);
    return {vals.begin(), vals.end()};
}

int encrypted_targeted(const Context& ctx)
{
    int failures = 0;
    for (const int p : {9, 17, 25, 33}) {
        for (int k = 1; k <= std::min(p - 1, 11); ++k) {
            for (const auto mode :
                 {FusedInputRoundingMode::None,
                  FusedInputRoundingMode::DeltaQuarter,
                  FusedInputRoundingMode::DeltaHalf,
                  FusedInputRoundingMode::ThreeDeltaQuarter,
                  FusedInputRoundingMode::AutoCalibrated}) {
                for (const Torus m : targeted_ms(p, k)) {
                    TFHEpp::TLWE<P2> out;
                    FusedGapStats stats;
                    gapmsb_pruned_bitk_fused_final_l2_to_l2(
                        out, encrypt_m(m, p, ctx), p, k, BOOL_ONE, mode,
                        *ctx.iksk, *ctx.bkfft, &stats);
                    const bool actual = arith_bit(out, ctx);
                    const bool expected = expected_msb(m, p);
                    if (actual != expected) {
                        ++failures;
                        std::cerr << "fused_encrypted_failure p=" << p
                                  << " k=" << k << " m=" << m
                                  << " expected=" << expected
                                  << " actual=" << actual
                                  << " phase="
                                  << hex64(phase_l2(out, ctx.key2)) << "\n";
                        return failures;
                    }
                }
            }
        }
    }
    return failures;
}

}  // namespace

int main(int argc, char** argv)
{
    bool run_encrypted = false;
    std::uint64_t seed = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--encrypted-targeted")
            run_encrypted = true;
        else if (a == "--seed" && i + 1 < argc)
            seed = std::stoull(argv[++i]);
        else
            throw std::runtime_error("unknown argument: " + a);
    }

    if (!pure_oracle_ok()) return 1;
    int failures = 0;
    if (run_encrypted) {
        const Context ctx = make_context(seed);
        failures += encrypted_targeted(ctx);
    }
    std::cout << "FUSED_FINAL_TEST_RESULT failures=" << failures
              << " encrypted_targeted=" << run_encrypted << "\n";
    return failures == 0 ? 0 : 1;
}
