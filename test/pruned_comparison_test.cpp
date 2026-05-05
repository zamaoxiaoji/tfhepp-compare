#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include "../src/HEDB/comparison/HomCompare.h"

using namespace HEDB;

void test_lvl1(uint32_t plain_bits, int num_test)
{
    const uint32_t p_eff = plain_bits + 1;
    std::cout << "\n====== lvl1 | bits=" << plain_bits
              << " | p_eff=" << p_eff << " | N=" << num_test << " ======\n";

    std::random_device rd;
    std::default_random_engine eng(rd());
    using P = Lvl1;

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> msg(0, (1 << (plain_bits-1))-1);

    int std_err = 0, eth_err = 0;
    double std_t = 0, eth_t = 0;

    for (int t = 0; t < num_test; t++) {
        auto p0 = msg(eng), p1 = msg(eng);
        bool expect = p0 > p1;
        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 cres;

        auto t0 = std::chrono::high_resolution_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        std_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
        if (TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1) != expect) std_err++;

        t0 = std::chrono::high_resolution_clock::now();
        ethmsb_greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        t1 = std::chrono::high_resolution_clock::now();
        eth_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
        if (TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1) != expect) eth_err++;
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Standard: err=" << std_err << "/" << num_test << "  avg=" << std_t/num_test << "ms\n";
    std::cout << "  ETHMSB:   err=" << eth_err << "/" << num_test << "  avg=" << eth_t/num_test << "ms\n";
}

void test_lvl2(uint32_t plain_bits, int num_test)
{
    const uint32_t p_eff = plain_bits + 1;
    std::cout << "\n====== lvl2 | bits=" << plain_bits
              << " | p_eff=" << p_eff << " | N=" << num_test << " ======\n";

    std::random_device rd;
    std::default_random_engine eng(rd());
    using P = Lvl2;

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> msg(0, (1ULL << (plain_bits-1))-1);

    int std_err = 0, eth_err = 0;
    double std_t = 0, eth_t = 0;

    for (int t = 0; t < num_test; t++) {
        auto p0 = msg(eng), p1 = msg(eng);
        bool expect = p0 > p1;
        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 cres;

        auto t0 = std::chrono::high_resolution_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        std_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
        if (TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1) != expect) std_err++;

        t0 = std::chrono::high_resolution_clock::now();
        ethmsb_greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        t1 = std::chrono::high_resolution_clock::now();
        eth_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
        if (TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1) != expect) eth_err++;
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Standard: err=" << std_err << "/" << num_test << "  avg=" << std_t/num_test << "ms\n";
    std::cout << "  ETHMSB:   err=" << eth_err << "/" << num_test << "  avg=" << eth_t/num_test << "ms\n";
}

int main()
{
    int N = 50;
    std::cout << "=== ETHMSB Comparison Test ===\n";
    test_lvl1(4, N);
    test_lvl1(8, N);
    test_lvl2(16, N);
    test_lvl2(32, N);
    std::cout << "\n=== Done ===\n";
    return 0;
}
