#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "keyswitch.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"

namespace {

template <class T>
std::string hex_any(const T v)
{
    std::ostringstream os;
    os << "0x" << std::hex << +static_cast<std::uint64_t>(v);
    return os.str();
}

template <class P0>
typename P0::T expected_trunc(const std::uint64_t phase64)
{
    constexpr int bits = std::numeric_limits<typename P0::T>::digits;
    if constexpr (bits == 64) return static_cast<typename P0::T>(phase64);
    return static_cast<typename P0::T>(
        (phase64 + (std::uint64_t{1} << (64 - bits - 1))) >> (64 - bits));
}

void run_v10()
{
    using P2 = TFHEpp::lvl2param;
    using P0 = TFHEpp::lvl0param;
    using KS20 = TFHEpp::lvl20param;
    TFHEpp::SecretKey sk;
    auto iksk =
        std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    TFHEpp::ikskgen<KS20>(*iksk, sk.key.get<P2>(), sk.key.get<P0>());
    const std::vector<std::uint64_t> phases = {
        0,
        std::uint64_t{1} << 57,
        (std::uint64_t{1} << 57) + (std::uint64_t{1} << 46),
        std::uint64_t{1} << 58,
        std::uint64_t{1} << 61,
        std::uint64_t{1} << 63,
        ~std::uint64_t{0}};
    for (const auto phase64 : phases) {
        TFHEpp::TLWE<P2> in = {};
        in[P2::k * P2::n] = phase64;
        TFHEpp::TLWE<P0> out;
        TFHEpp::IdentityKeySwitch<KS20>(out, in, *iksk);
        const auto actual = TFHEpp::tlweSymPhase<P0>(out, sk.key.get<P0>());
        std::cout << "impl=v10_l20 phase64=" << hex_any(phase64)
                  << " expected_high_bits=" << hex_any(expected_trunc<P0>(phase64))
                  << " actual_lvl0_phase=" << hex_any(actual)
                  << " pass=" << (actual == expected_trunc<P0>(phase64))
                  << "\n";
    }
}

void run_h3compat()
{
    using namespace my_ethmsb_params;
    using P2 = my_h3_lvl2param;
    using P0 = my_h3_lvl0param;
    using KS20 = my_h3_lvl20param;
    std::mt19937_64 rng(0);
    auto keys = make_h3compat_keys(rng);
    auto iksk =
        std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    TFHEpp::ikskgen<KS20>(*iksk, keys.key2, keys.key0);
    const std::vector<std::uint64_t> phases = {
        0,
        std::uint64_t{1} << 57,
        (std::uint64_t{1} << 57) + (std::uint64_t{1} << 46),
        std::uint64_t{1} << 58,
        std::uint64_t{1} << 61,
        std::uint64_t{1} << 63,
        ~std::uint64_t{0}};
    for (const auto phase64 : phases) {
        TFHEpp::TLWE<P2> in = {};
        in[P2::k * P2::n] = phase64;
        TFHEpp::TLWE<P0> out;
        TFHEpp::IdentityKeySwitch<KS20>(out, in, *iksk);
        const auto actual = TFHEpp::tlweSymPhase<P0>(out, keys.key0);
        std::cout << "impl=h3compat_l20 phase64=" << hex_any(phase64)
                  << " expected_high_bits=" << hex_any(expected_trunc<P0>(phase64))
                  << " actual_lvl0_phase=" << hex_any(actual)
                  << " pass=" << (actual == expected_trunc<P0>(phase64))
                  << "\n";
    }
}

}  // namespace

int main()
{
    std::cout << "v10_lvl0_T_bits="
              << std::numeric_limits<typename TFHEpp::lvl0param::T>::digits
              << "\n";
    std::cout << "h3compat_lvl0_T_bits="
              << std::numeric_limits<
                     typename my_ethmsb_params::my_h3_lvl0param::T>::digits
              << "\n";
    run_v10();
    run_h3compat();
    return 0;
}
