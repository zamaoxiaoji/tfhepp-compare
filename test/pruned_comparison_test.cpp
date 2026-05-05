#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/utils/utils.h"

using namespace HEDB;

// ============================================================
//  lvl1 comparison test (1–9 bits)
// ============================================================
void test_lvl1(uint32_t plain_bits, int num_test, uint32_t guard_k)
{
    const uint32_t p_eff = plain_bits + 1;
    // Validate guard_k range
    if (guard_k >= p_eff - 1) {
        std::cout << "\n[SKIP] lvl1 bits=" << plain_bits
                  << " guard_k=" << guard_k << " (invalid: guard_k >= p_eff-1)\n";
        return;
    }

    std::cout << "\n====== lvl1 Test | bits=" << plain_bits
              << " | p_eff=" << p_eff
              << " | guard_k=" << guard_k
              << " | w_k=" << (1ULL << (p_eff - 1 - guard_k))
              << " | N=" << num_test << " ======\n";

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);

    uint32_t scale_bits =
        std::numeric_limits<P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> message(
        0, (1 << (plain_bits - 1)) - 1);

    int std_errors = 0, gap_errors = 0;
    double std_time = 0, gap_time = 0;

    for (int t = 0; t < num_test; t++) {
        auto p0 = message(engine);
        auto p1 = message(engine);
        bool expected_gt = (p0 > p1) ? true : false;

        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(
            p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(
            p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 cres;

        // --- Standard comparison ---
        auto t0 = std::chrono::high_resolution_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        std_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        t1 - t0).count() / 1000.0;
        bool dec_std = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
        if (dec_std != expected_gt) std_errors++;

        // --- Gap MSB comparison ---
        t0 = std::chrono::high_resolution_clock::now();
        gap_greater_than<P>(c0, c1, cres, plain_bits, ek,
                            LOGIC, guard_k);
        t1 = std::chrono::high_resolution_clock::now();
        gap_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        t1 - t0).count() / 1000.0;
        bool dec_gap = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
        if (dec_gap != expected_gt) gap_errors++;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Standard:  errors=" << std_errors << "/" << num_test
              << "  avg=" << std_time / num_test << " ms\n";
    std::cout << "  GapMSB:    errors=" << gap_errors << "/" << num_test
              << "  avg=" << gap_time / num_test << " ms\n";
}

// ============================================================
//  lvl2 comparison test (10–32 bits)
// ============================================================
void test_lvl2(uint32_t plain_bits, int num_test, uint32_t guard_k)
{
    const uint32_t p_eff = plain_bits + 1;
    if (guard_k >= p_eff - 1) {
        std::cout << "\n[SKIP] lvl2 bits=" << plain_bits
                  << " guard_k=" << guard_k << " (invalid)\n";
        return;
    }

    std::cout << "\n====== lvl2 Test | bits=" << plain_bits
              << " | p_eff=" << p_eff
              << " | guard_k=" << guard_k
              << " | w_k=" << (1ULL << (p_eff - 1 - guard_k))
              << " | N=" << num_test << " ======\n";

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl2;

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    uint32_t scale_bits =
        std::numeric_limits<P::T>::digits - plain_bits - 1;
    std::uniform_int_distribution<typename P::T> message(
        0, (1ULL << (plain_bits - 1)) - 1);

    int std_errors = 0, gap_errors = 0;
    double std_time = 0, gap_time = 0;

    for (int t = 0; t < num_test; t++) {
        auto p0 = message(engine);
        auto p1 = message(engine);
        bool expected_gt = (p0 > p1) ? true : false;

        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(
            p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(
            p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 cres;

        // --- Standard comparison ---
        auto t0 = std::chrono::high_resolution_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        std_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        t1 - t0).count() / 1000.0;
        bool dec_std = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (dec_std != expected_gt) std_errors++;

        // --- Gap MSB comparison ---
        t0 = std::chrono::high_resolution_clock::now();
        gap_greater_than<P>(c0, c1, cres, plain_bits, ek,
                            LOGIC, guard_k);
        t1 = std::chrono::high_resolution_clock::now();
        gap_time += std::chrono::duration_cast<std::chrono::microseconds>(
                        t1 - t0).count() / 1000.0;
        bool dec_gap = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (dec_gap != expected_gt) gap_errors++;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Standard:  errors=" << std_errors << "/" << num_test
              << "  avg=" << std_time / num_test << " ms\n";
    std::cout << "  GapMSB:    errors=" << gap_errors << "/" << num_test
              << "  avg=" << gap_time / num_test << " ms\n";
}

int main()
{
    int N = 50;
    std::cout << "=== GapMSB Comparison Test ===\n";
    std::cout << "guard_k = p_eff - 5 (paper default, κ=5)\n";

    // 4-bit: p_eff=5, guard_k = max(1, p_eff-5) = 1
    // w_k = 2^(5-1-1) = 8, gap = 9*Δ, Δ = Q/32
    // offset = 4.5 * Q/32 ≈ Q/7, too large for 4-bit
    // For 4-bit, use guard_k = p_eff-2 = 3 (more conservative)
    test_lvl1(4, N, 3);

    // 8-bit: p_eff=9, guard_k = p_eff-5 = 4
    test_lvl1(8, N, 4);

    // 16-bit: p_eff=17, guard_k = p_eff-5 = 12
    test_lvl2(16, N, 12);

    // 32-bit: p_eff=33, guard_k = p_eff-5 = 28
    test_lvl2(32, N, 28);

    std::cout << "\n=== All tests complete ===\n";
    return 0;
}
