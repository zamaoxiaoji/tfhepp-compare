#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <vector>
#include "../src/HEDB/comparison/HomCompare.h"

using namespace HEDB;

std::vector<std::pair<uint64_t, uint64_t>> make_cases(uint32_t plain_bits,
                                                      int random_per_seed)
{
    const uint64_t max_msg = (1ULL << (plain_bits - 1)) - 1;
    const uint64_t half = max_msg / 2;
    const uint64_t margin = 16;
    std::vector<std::pair<uint64_t, uint64_t>> cases = {
        {margin, 0}, {0, margin}, {max_msg, max_msg - margin},
        {max_msg - margin, max_msg}, {max_msg, 0}, {0, max_msg},
        {half + margin, half}, {half, half + margin}, {half, 0}, {0, half},
    };

    const std::vector<uint64_t> seeds = {
        0x13579BDFULL, 0x2468ACE0ULL, 0xC0FFEE1234ULL,
    };
    std::uniform_int_distribution<uint64_t> msg(0, max_msg);
    for (uint64_t seed : seeds) {
        std::mt19937_64 eng(seed);
        for (int i = 0; i < random_per_seed; i++) {
            uint64_t p0 = msg(eng);
            uint64_t p1 = msg(eng);
            while (p0 > p1 ? p0 - p1 < margin : p1 - p0 < margin) {
                p0 = msg(eng);
                p1 = msg(eng);
            }
            cases.emplace_back(p0, p1);
        }
    }
    return cases;
}

bool test_lvl1(uint32_t plain_bits, int random_per_seed)
{
    const uint32_t p_eff = plain_bits + 1;
    const auto cases = make_cases(plain_bits, random_per_seed);
    std::cout << "\n====== lvl1 | bits=" << plain_bits
              << " | p_eff=" << p_eff
              << " | cases=" << cases.size() << " ======\n";

    using P = Lvl1;

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;

    int std_err = 0, eth_err = 0;
    double std_t = 0, eth_t = 0;

    size_t case_id = 0;
    for (auto [raw0, raw1] : cases) {
        auto p0 = static_cast<typename P::T>(raw0);
        auto p1 = static_cast<typename P::T>(raw1);
        bool expect = p0 > p1;
        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 cres;

        auto t0 = std::chrono::high_resolution_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        std_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
        const bool std_dec = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
        if (std_dec != expect) {
            std_err++;
            std::cout << "  Standard fail case#" << case_id
                      << " p0=" << raw0 << " p1=" << raw1
                      << " got=" << std_dec << " expect=" << expect << "\n";
        }

        t0 = std::chrono::high_resolution_clock::now();
        try {
            ethmsb_greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
            t1 = std::chrono::high_resolution_clock::now();
            eth_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
            const bool eth_dec = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
            if (eth_dec != expect) {
                eth_err++;
                std::cout << "  ETHMSB fail case#" << case_id
                          << " p0=" << raw0 << " p1=" << raw1
                          << " got=" << eth_dec << " expect=" << expect
                          << "\n";
            }
        } catch (const std::invalid_argument &e) {
            std::cout << "  ETHMSB skipped: " << e.what() << "\n";
            return false;
        }
        case_id++;
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Standard: err=" << std_err << "/" << cases.size()
              << "  avg=" << std_t/cases.size() << "ms\n";
    std::cout << "  ETHMSB:   err=" << eth_err << "/" << cases.size()
              << "  avg=" << eth_t/cases.size() << "ms\n";
    return std_err == 0 && eth_err == 0;
}

bool test_lvl2(uint32_t plain_bits, int random_per_seed)
{
    const uint32_t p_eff = plain_bits + 1;
    const auto cases = make_cases(plain_bits, random_per_seed);
    std::cout << "\n====== lvl2 | bits=" << plain_bits
              << " | p_eff=" << p_eff
              << " | cases=" << cases.size() << " ======\n";

    using P = Lvl2;

    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;

    int std_err = 0, eth_err = 0;
    double std_t = 0, eth_t = 0;

    size_t case_id = 0;
    for (auto [raw0, raw1] : cases) {
        auto p0 = static_cast<typename P::T>(raw0);
        auto p1 = static_cast<typename P::T>(raw1);
        bool expect = p0 > p1;
        auto c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        auto c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 cres;

        auto t0 = std::chrono::high_resolution_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        auto t1 = std::chrono::high_resolution_clock::now();
        std_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
        const bool std_dec = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (std_dec != expect) {
            std_err++;
            std::cout << "  Standard fail case#" << case_id
                      << " p0=" << raw0 << " p1=" << raw1
                      << " got=" << std_dec << " expect=" << expect << "\n";
        }

        t0 = std::chrono::high_resolution_clock::now();
        try {
            ethmsb_greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
            t1 = std::chrono::high_resolution_clock::now();
            eth_t += std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.0;
            const bool eth_dec = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
            if (eth_dec != expect) {
                eth_err++;
                std::cout << "  ETHMSB fail case#" << case_id
                          << " p0=" << raw0 << " p1=" << raw1
                          << " got=" << eth_dec << " expect=" << expect
                          << "\n";
            }
        } catch (const std::invalid_argument &e) {
            std::cout << "  ETHMSB skipped: " << e.what() << "\n";
            return false;
        }
        case_id++;
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Standard: err=" << std_err << "/" << cases.size()
              << "  avg=" << std_t/cases.size() << "ms\n";
    std::cout << "  ETHMSB:   err=" << eth_err << "/" << cases.size()
              << "  avg=" << eth_t/cases.size() << "ms\n";
    return std_err == 0 && eth_err == 0;
}

int main()
{
    int random_per_seed = 8;
    bool ok = true;
    std::cout << "=== ETHMSB Offset-only Comparison Test ===\n";
    ok &= test_lvl1(4, random_per_seed);
    ok &= test_lvl1(8, random_per_seed);
    ok &= test_lvl2(16, random_per_seed);
    ok &= test_lvl2(32, random_per_seed);
    std::cout << "\n=== Done ===\n";
    return ok ? 0 : 1;
}
