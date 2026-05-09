#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

#include "HEDB/comparison/HomCompare.h"
#include "HEDB/utils/types.h"

namespace {

template <class P>
void RunGreaterThan(uint32_t plain_bits, int num_test)
{
    std::mt19937_64 engine(0x4845334442ULL + plain_bits);
    HEDB::TFHESecretKey sk;
    HEDB::TFHEEvalKey ek;
    ek.emplacebkfft<HEDB::Lvl01>(sk);
    if constexpr (std::is_same_v<P, HEDB::Lvl2>) ek.emplacebkfft<HEDB::Lvl02>(sk);
    ek.emplaceiksk<HEDB::Lvl10>(sk);
    if constexpr (std::is_same_v<P, HEDB::Lvl2>) {
        ek.emplaceiksk<HEDB::Lvl20>(sk);
        ek.emplaceiksk<HEDB::Lvl21>(sk);
    }

    const uint32_t scale_bits =
        std::numeric_limits<typename P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> message(
        0, (static_cast<typename P::T>(1) << (plain_bits - 1)) - 1);

    int errors = 0;
    double millis = 0;
    for (int t = 0; t < num_test; t++) {
        const typename P::T p0 = message(engine);
        const typename P::T p1 = message(engine);
        const bool expected = p0 > p1;
        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(
            p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(
            p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

        HEDB::TLWELvl1 cres;
        const auto start = std::chrono::high_resolution_clock::now();
        HEDB::greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        const auto end = std::chrono::high_resolution_clock::now();
        millis += std::chrono::duration<double, std::milli>(end - start).count();

        const bool got = TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(cres, sk.key.lvl1);
        if (got != expected) {
            errors++;
            if (errors <= 3) {
                std::cout << "  fail p0=" << p0 << " p1=" << p1
                          << " expected=" << expected << " got=" << got
                          << '\n';
            }
        }
    }

    std::cout << "plain_bits=" << plain_bits << " errors=" << errors << '/'
              << num_test << " avg_ms=" << millis / num_test << '\n';
}

}  // namespace

int main()
{
    constexpr int num_test = 20;
    std::cout << "HE3DB original comparison sources on local TFHEpp\n";
    RunGreaterThan<HEDB::Lvl1>(4, num_test);
    RunGreaterThan<HEDB::Lvl1>(8, num_test);
    RunGreaterThan<HEDB::Lvl2>(16, num_test);
    RunGreaterThan<HEDB::Lvl2>(32, num_test);
}
