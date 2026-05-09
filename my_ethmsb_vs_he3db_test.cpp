/**
 * @file my_ethmsb_vs_he3db_test.cpp
 * @brief Benchmarks ETHMSB vs HE3DB-style HomMSB (both on same TFHEpp params)
 */
#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>

#include "ethmsb_compare.h"
#include "hommsb_he3db_style.h"

using ETHMSB_NS::tlweSymInt32Encrypt;

// ─────────────────────────────────────────────────────────
//  Lvl1 test (plain_bits ≤ 9)
// ─────────────────────────────────────────────────────────
void test_lvl1(uint32_t plain_bits, int num_test)
{
    std::cout << "\n=== Lvl1 Test  plain_bits=" << plain_bits
              << "  num_test=" << num_test << " ===" << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = TFHEpp::lvl1param;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> message(
        0, (1U << (plain_bits - 1)) - 1);

    uint32_t he3db_err = 0, ethmsb_err = 0;
    double he3db_time = 0, ethmsb_time = 0;

    for (int t = 0; t < num_test; t++) {
        typename P::T p0 = message(engine);
        typename P::T p1 = message(engine);
        uint32_t expected_gt = (p0 > p1) ? 1 : 0;

        auto c0 = tlweSymInt32Encrypt<P>(
            p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = tlweSymInt32Encrypt<P>(
            p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

        // --- HE3DB-style HomMSB ---
        TFHEpp::TLWE<P> res_he3db;
        auto c0h = c0, c1h = c1;
        auto t0 = std::chrono::high_resolution_clock::now();
        HomMSB_NS::greater_than<P>(c0h, c1h, res_he3db, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        he3db_time += std::chrono::duration<double, std::milli>(t1 - t0).count();
        auto d_he3db = TFHEpp::tlweSymDecrypt<P>(res_he3db, sk.key.lvl1);
        if (d_he3db != expected_gt) he3db_err++;

        // --- ETHMSB ---
        TFHEpp::TLWE<P> res_ethmsb;
        auto c0e = c0, c1e = c1;
        t0 = std::chrono::high_resolution_clock::now();
        ETHMSB_NS::greater_than<P>(c0e, c1e, res_ethmsb, plain_bits, ek, LOGIC);
        t1 = std::chrono::high_resolution_clock::now();
        ethmsb_time += std::chrono::duration<double, std::milli>(t1 - t0).count();
        auto d_ethmsb = TFHEpp::tlweSymDecrypt<P>(res_ethmsb, sk.key.lvl1);
        if (d_ethmsb != expected_gt) ethmsb_err++;

        if (t < 3) {
            std::cout << "  [" << t << "] p0=" << p0 << " p1=" << p1
                      << " exp=" << expected_gt
                      << " he3db=" << d_he3db
                      << " ethmsb=" << d_ethmsb << std::endl;
        }
    }

    double speedup = he3db_time / std::max(ethmsb_time, 0.001);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "| Method | Avg ms | Errors | Speedup |" << std::endl;
    std::cout << "|--------|--------|--------|---------|" << std::endl;
    std::cout << "| HE3DB  | " << he3db_time / num_test
              << " | " << he3db_err << "/" << num_test << " | 1.00x |" << std::endl;
    std::cout << "| ETHMSB | " << ethmsb_time / num_test
              << " | " << ethmsb_err << "/" << num_test
              << " | " << speedup << "x |" << std::endl;
}

// ─────────────────────────────────────────────────────────
//  Lvl2 test (plain_bits ≥ 10)
// ─────────────────────────────────────────────────────────
void test_lvl2(uint32_t plain_bits, int num_test)
{
    std::cout << "\n=== Lvl2 Test  plain_bits=" << plain_bits
              << "  num_test=" << num_test << " ===" << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = TFHEpp::lvl2param;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    ek.emplaceiksk<TFHEpp::lvl21param>(sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> message(
        0, (1ULL << (plain_bits - 1)) - 1);

    uint32_t he3db_err = 0, ethmsb_err = 0;
    double he3db_time = 0, ethmsb_time = 0;

    for (int t = 0; t < num_test; t++) {
        typename P::T p0 = message(engine);
        typename P::T p1 = message(engine);
        uint32_t expected_gt = (p0 > p1) ? 1 : 0;

        auto c0 = tlweSymInt32Encrypt<P>(
            p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = tlweSymInt32Encrypt<P>(
            p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

        // --- HE3DB-style HomMSB ---
        TFHEpp::TLWE<TFHEpp::lvl1param> res_he3db;
        auto c0h = c0, c1h = c1;
        auto t0 = std::chrono::high_resolution_clock::now();
        HomMSB_NS::greater_than<P>(c0h, c1h, res_he3db, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        he3db_time += std::chrono::duration<double, std::milli>(t1 - t0).count();
        auto d_he3db = TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(
            res_he3db, sk.key.lvl1);
        if (d_he3db != expected_gt) he3db_err++;

        // --- ETHMSB ---
        TFHEpp::TLWE<TFHEpp::lvl1param> res_ethmsb;
        auto c0e = c0, c1e = c1;
        t0 = std::chrono::high_resolution_clock::now();
        ETHMSB_NS::greater_than<P>(c0e, c1e, res_ethmsb, plain_bits, ek, LOGIC);
        t1 = std::chrono::high_resolution_clock::now();
        ethmsb_time += std::chrono::duration<double, std::milli>(t1 - t0).count();
        auto d_ethmsb = TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(
            res_ethmsb, sk.key.lvl1);
        if (d_ethmsb != expected_gt) ethmsb_err++;

        if (t < 3) {
            std::cout << "  [" << t << "] p0=" << p0 << " p1=" << p1
                      << " exp=" << expected_gt
                      << " he3db=" << d_he3db
                      << " ethmsb=" << d_ethmsb << std::endl;
        }
    }

    double speedup = he3db_time / std::max(ethmsb_time, 0.001);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "| Method | Avg ms | Errors | Speedup |" << std::endl;
    std::cout << "|--------|--------|--------|---------|" << std::endl;
    std::cout << "| HE3DB  | " << he3db_time / num_test
              << " | " << he3db_err << "/" << num_test << " | 1.00x |" << std::endl;
    std::cout << "| ETHMSB | " << ethmsb_time / num_test
              << " | " << ethmsb_err << "/" << num_test
              << " | " << speedup << "x |" << std::endl;
}

int main()
{
    int num_test = 100;

    std::cout << "========================================" << std::endl;
    std::cout << "  ETHMSB vs HE3DB Comparison Benchmark  " << std::endl;
    std::cout << "  (Same TFHEpp params for both)         " << std::endl;
    std::cout << "========================================" << std::endl;

    test_lvl1(4, num_test);
    test_lvl1(8, num_test);
    test_lvl2(16, num_test);
    test_lvl2(32, num_test);

    return 0;
}
