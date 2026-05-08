#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"

namespace {

using namespace my_ethmsb_params;
using P0 = my_h3_lvl0param;
using P2 = my_h3_lvl2param;
using KS20 = my_h3_lvl20param;
using BR02 = my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Context {
    TFHEpp::Key<P0> key0;
    TFHEpp::Key<P2> key2;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

Torus encode(const Torus m, const int k)
{
    return static_cast<Torus>(Wide{m} * Wide{my_ethmsb::delta_bits(k)});
}

Torus mask_for(const int k) { return (Torus{1} << k) - 1; }

Torus phase(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return TFHEpp::tlweSymPhase<P2>(ct, ctx.key2);
}

bool decode_arith_phase(const Torus ph, const Torus out_value)
{
    return my_ethmsb::torus_abs_centered(ph - out_value) <
           my_ethmsb::torus_abs_centered(ph);
}

Context make_context()
{
    std::mt19937_64 rng(0);
    const auto keys = make_h3compat_keys(rng);
    Context ctx;
    ctx.key0 = keys.key0;
    ctx.key2 = keys.key2;
    ctx.iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    return ctx;
}

TFHEpp::TLWE<P2> input_ct(const std::string& mode, const Torus m, const int k,
                          const Context& ctx)
{
    if (mode == "trivial") {
        TFHEpp::TLWE<P2> ct = {};
        ct[P2::k * P2::n] = encode(m, k);
        return ct;
    }
    return TFHEpp::tlweSymEncrypt<P2>(encode(m, k), ctx.key2);
}

TFHEpp::TLWE<P2> pbs_trace(const Context& ctx, TFHEpp::TLWE<P2> in,
                           const int k, const int depth,
                           const std::string& stage, const Torus offset,
                           const Torus out_value, int& failures)
{
    const Torus before = phase(in, ctx);
    my_ethmsb::add_const_inplace<P2>(in, offset);
    const Torus after = phase(in, ctx);
    const Torus expected_after = before + offset;
    if (after != expected_after) {
        ++failures;
        std::cout << "OFFSET_FAILURE kind=after_offset_mismatch depth="
                  << depth << " k=" << k << " stage=" << stage
                  << " before=" << hex64(before)
                  << " offset=" << hex64(offset)
                  << " expected_after=" << hex64(expected_after)
                  << " actual_after=" << hex64(after) << "\n";
    }
    TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, *ctx.iksk);
    const Torus half = out_value / 2;
    TFHEpp::Polynomial<P2> tv;
    tv.fill(Torus{0} - half);
    TFHEpp::TLWE<P2> out;
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out, in0, *ctx.bkfft, tv);
    my_ethmsb::add_const_inplace<P2>(out, half);
    const Torus out_phase = phase(out, ctx);
    std::cout << "OFFSET_PBS depth=" << depth << " k=" << k
              << " stage=" << stage << " offset_hex=" << hex64(offset)
              << " phase_before_offset=" << hex64(before)
              << " phase_after_offset=" << hex64(after)
              << " out_value_hex=" << hex64(out_value)
              << " output_phase=" << hex64(out_phase)
              << " decoded=" << decode_arith_phase(out_phase, out_value)
              << "\n";
    return out;
}

TFHEpp::TLWE<P2> trace_ethmsb(const Context& ctx, const TFHEpp::TLWE<P2>& ct,
                              const int k, const int kappa, const Torus m,
                              const Torus out_value, const int depth,
                              const std::string& mode, int& failures)
{
    const Torus plain = m & mask_for(k);
    const Torus input_phase = phase(ct, ctx);
    std::cout << "BEGIN_OFFSET_NODE mode=" << mode << " depth=" << depth
              << " k=" << k << " kappa=" << kappa
              << " input_phase=" << hex64(input_phase)
              << " out_value=" << hex64(out_value) << "\n";

    if (k <= kappa) {
        const Torus base_offset = my_ethmsb::base_offset_for_current_layer(k);
        std::cout << "OFFSET_LAYER mode=" << mode << " depth=" << depth
                  << " k=" << k << " stage=base base_offset_hex="
                  << hex64(base_offset) << "\n";
        auto out = pbs_trace(ctx, ct, k, depth, "base", base_offset, out_value,
                             failures);
        std::cout << "END_OFFSET_NODE mode=" << mode << " depth=" << depth
                  << " k=" << k << "\n";
        return out;
    }

    const int suffix_bits = k - kappa;
    const Torus suffix_plain = plain & mask_for(suffix_bits);
    const Torus local_gap =
        my_ethmsb::gap_offset_for_current_layer(k, kappa);
    const Torus delta_prime =
        my_ethmsb::guard_value_for_parent_scale(k, kappa);
    const Torus guard_weight = my_ethmsb::guard_weight(k, kappa);
    std::cout << "OFFSET_LAYER mode=" << mode << " depth=" << depth
              << " k=" << k
              << " stage=recursive local_gap_offset_hex="
              << hex64(local_gap) << " delta_prime_hex="
              << hex64(delta_prime)
              << " analysis_only=true guard_weight=" << guard_weight << "\n";

    TFHEpp::TLWE<P2> shifted;
    my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);
    const Torus shifted_phase = phase(shifted, ctx);
    const Torus expected_shifted_phase = input_phase << kappa;
    if (shifted_phase != expected_shifted_phase) {
        ++failures;
        std::cout << "OFFSET_FAILURE kind=shifted_contains_unexpected_value"
                  << " mode=" << mode << " depth=" << depth << " k=" << k
                  << " expected_shifted_phase="
                  << hex64(expected_shifted_phase)
                  << " actual_shifted_phase=" << hex64(shifted_phase)
                  << "\n";
    }

    TFHEpp::TLWE<P2> guard =
        trace_ethmsb(ctx, shifted, suffix_bits, kappa, suffix_plain,
                     delta_prime, depth + 1, mode, failures);
    const Torus guard_phase = phase(guard, ctx);
    const bool expected_guard =
        ((suffix_plain >> (suffix_bits - 1)) & Torus{1}) != 0;
    const bool decoded_guard = decode_arith_phase(guard_phase, delta_prime);
    if (decoded_guard != expected_guard) {
        ++failures;
        std::cout << "OFFSET_FAILURE kind=guard_decode_mismatch mode=" << mode
                  << " depth=" << depth << " k=" << k
                  << " expected=" << expected_guard
                  << " actual=" << decoded_guard
                  << " guard_phase=" << hex64(guard_phase)
                  << " guard_out_value=" << hex64(delta_prime) << "\n";
    }

    TFHEpp::TLWE<P2> guarded;
    my_ethmsb::sub<P2>(guarded, ct, guard);
    const Torus guarded_phase = phase(guarded, ctx);
    const Torus expected_guarded_phase = input_phase - guard_phase;
    if (guarded_phase != expected_guarded_phase) {
        ++failures;
        std::cout << "OFFSET_FAILURE kind=guarded_not_input_minus_guard"
                  << " mode=" << mode << " depth=" << depth << " k=" << k
                  << " input_phase=" << hex64(input_phase)
                  << " guard_phase=" << hex64(guard_phase)
                  << " expected_guarded=" << hex64(expected_guarded_phase)
                  << " actual_guarded=" << hex64(guarded_phase) << "\n";
    }

    auto out =
        pbs_trace(ctx, guarded, k, depth, "final", local_gap, out_value,
                  failures);
    std::cout << "END_OFFSET_NODE mode=" << mode << " depth=" << depth
              << " k=" << k << "\n";
    return out;
}

int formula_tests()
{
    int failures = 0;
    const std::vector<std::pair<int, Torus>> gap_expected = {
        {33, 0x0200000040000000ULL}, {28, 0x0200000800000000ULL},
        {23, 0x0200010000000000ULL}, {18, 0x0200200000000000ULL},
        {13, 0x0204000000000000ULL}, {8, 0x0280000000000000ULL},
    };
    for (const auto& [k, expected] : gap_expected) {
        const Torus actual = my_ethmsb::gap_offset_for_current_layer(k, 5);
        std::cout << "OFFSET_FORMULA k=" << k << " kappa=5 actual="
                  << hex64(actual) << " expected=" << hex64(expected)
                  << " pass=" << (actual == expected) << "\n";
        if (actual != expected) ++failures;
    }
    const Torus base = my_ethmsb::base_offset_for_current_layer(3);
    std::cout << "OFFSET_FORMULA k=3 stage=base actual=" << hex64(base)
              << " expected=0x1000000000000000 pass="
              << (base == 0x1000000000000000ULL) << "\n";
    if (base != 0x1000000000000000ULL) ++failures;
    const Torus delta33 = my_ethmsb::delta_bits(33);
    const Torus guard33 = my_ethmsb::guard_value_for_parent_scale(33, 5);
    std::cout << "OFFSET_FORMULA k=33 delta_k=" << hex64(delta33)
              << " guard_value_delta_prime=" << hex64(guard33)
              << " pass="
              << (delta33 == 0x0000000080000000ULL &&
                  guard33 == 0x0400000000000000ULL)
              << "\n";
    if (delta33 != 0x0000000080000000ULL ||
        guard33 != 0x0400000000000000ULL)
        ++failures;
    return failures;
}

int trace_tests()
{
    int failures = 0;
    Context ctx = make_context();
    for (const std::string mode : {"trivial", "encrypted"}) {
        const auto ct = input_ct(mode, 0, 33, ctx);
        const auto out = trace_ethmsb(ctx, ct, 33, 5, 0, my_ethmsb::BOOL_ONE,
                                      0, mode, failures);
        const Torus out_phase = phase(out, ctx);
        const bool actual = decode_arith_phase(out_phase, my_ethmsb::BOOL_ONE);
        if (actual) {
            ++failures;
            std::cout << "OFFSET_FAILURE kind=final_decode mode=" << mode
                      << " expected=0 actual=1 output_phase="
                      << hex64(out_phase) << "\n";
        }
    }
    return failures;
}

}  // namespace

int main()
{
    const int failures = formula_tests() + trace_tests();
    std::cout << "MY_ETHMSB_OFFSET_SEMANTICS_RESULT failures=" << failures
              << "\n";
    return failures == 0 ? 0 : 1;
}
