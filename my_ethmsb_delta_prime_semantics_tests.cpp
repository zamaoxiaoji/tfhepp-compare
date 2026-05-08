#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "my_ethmsb_fixed.hpp"

namespace {

using Torus = std::uint64_t;

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

bool decode_arith_phase(const Torus phase, const Torus out_value)
{
    return my_ethmsb::torus_abs_centered(phase - out_value) <
           my_ethmsb::torus_abs_centered(phase);
}

int final_bit_decode_independence()
{
    int failures = 0;
    const Torus out_value = my_ethmsb::BOOL_ONE;
    const Torus delta_prime = 0x0400000000000000ULL;
    const Torus small_noise = 0x0000000000100000ULL;
    const std::vector<Torus> phases = {
        0,
        small_noise,
        out_value,
        out_value + small_noise,
        delta_prime,
        delta_prime / 2,
        delta_prime + small_noise,
    };
    for (const Torus phase : phases) {
        const Torus dist0 = my_ethmsb::torus_abs_centered(phase);
        const Torus dist1 = my_ethmsb::torus_abs_centered(phase - out_value);
        const bool actual = decode_arith_phase(phase, out_value);
        const bool expected = dist1 < dist0;
        const bool would_delta_prime_decode =
            (phase / delta_prime) != 0;  // deliberately not used.
        std::cout << "DELTA_PRIME_FINAL_DECODE phase=" << hex64(phase)
                  << " out_value=" << hex64(out_value)
                  << " delta_prime=" << hex64(delta_prime)
                  << " dist0=" << hex64(dist0)
                  << " dist1=" << hex64(dist1)
                  << " actual=" << actual << " expected=" << expected
                  << " coarse_delta_prime_division="
                  << would_delta_prime_decode
                  << " analysis_only=true pass=" << (actual == expected)
                  << "\n";
        if (actual != expected) ++failures;
        if (phase == delta_prime && actual) {
            ++failures;
            std::cout << "DELTA_PRIME_FAILURE phase_delta_prime_decoded_as_one"
                      << "\n";
        }
    }
    return failures;
}

int guard_decode_dynamic_out_value()
{
    int failures = 0;
    for (const int k : {8, 13, 18, 23, 28, 33}) {
        const Torus guard_value =
            my_ethmsb::guard_value_for_parent_scale(k, 5);
        const Torus phase0 = 0;
        const Torus phase1 = guard_value;
        const bool zero_actual = decode_arith_phase(phase0, guard_value);
        const bool one_actual = decode_arith_phase(phase1, guard_value);
        const bool one_with_wrong_bool =
            decode_arith_phase(phase1, my_ethmsb::BOOL_ONE);
        std::cout << "DELTA_PRIME_GUARD_DECODE k=" << k
                  << " guard_out_value=" << hex64(guard_value)
                  << " phase0_actual=" << zero_actual
                  << " phase1_actual=" << one_actual
                  << " phase1_wrong_BOOL_ONE_decode=" << one_with_wrong_bool
                  << " pass=" << (!zero_actual && one_actual) << "\n";
        if (zero_actual || !one_actual) ++failures;
    }
    return failures;
}

int margin_only_delta_prime()
{
    int failures = 0;
    for (const int k : {8, 13, 18, 23, 28, 33}) {
        const Torus delta_prime =
            my_ethmsb::guard_value_for_parent_scale(k, 5);
        const bool ok = delta_prime == 0x0400000000000000ULL;
        const Torus gap = my_ethmsb::gap_offset_for_current_layer(k, 5);
        std::cout << "DELTA_PRIME_MARGIN_ONLY k=" << k
                  << " delta_prime_hex=" << hex64(delta_prime)
                  << " local_gap_offset_hex=" << hex64(gap)
                  << " analysis_only=true correctness_decode_uses_delta_prime="
                  << "false pass=" << ok << "\n";
        if (!ok) ++failures;
    }
    return failures;
}

}  // namespace

int main()
{
    const int failures = final_bit_decode_independence() +
                         guard_decode_dynamic_out_value() +
                         margin_only_delta_prime();
    std::cout << "MY_ETHMSB_DELTA_PRIME_SEMANTICS_RESULT failures="
              << failures << "\n";
    return failures == 0 ? 0 : 1;
}
