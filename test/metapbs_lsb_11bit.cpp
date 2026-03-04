// test/metapbs_lsb_11bit.cpp
// Test: MetaPBS with t=2N, period=2, LSB extraction LUT (0,q/2,0,q/2,...)
// over the full 11-bit input range [0, 2047].
//
// Build: make -j metapbs_lsb_11bit
// Run:   ./test/metapbs_lsb_11bit [--trials N]

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>

#include <metapbs.hpp>
#include <tfhe++.hpp>

int main(int argc, char **argv)
{
    int trials = 2048;  // default: full domain
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--trials" && i + 1 < argc)
            trials = std::stoi(argv[++i]);
    }

    using bkP = TFHEpp::lvlh1param;
    using TargetP = typename bkP::targetP;  // lvl1
    using DomainP = typename bkP::domainP;  // lvlhalf

    constexpr int N = static_cast<int>(TargetP::n);     // 1024
    constexpr std::int64_t t = 2 * N;                   // 2048 = 2N
    constexpr int nbit = static_cast<int>(TargetP::nbit); // 10

    std::cout << "[metapbs_lsb_11bit] N=" << N << " t=2N=" << t
              << " nbit=" << nbit << " trials=" << trials << "\n";

    // ── Key generation ──
    std::cout << "Generating keys..." << std::flush;
    TFHEpp::SecretKey sk;

    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<bkP>>();
    TFHEpp::bkfftgen<bkP>(*bkfft, sk);

    auto ahk = std::make_unique<TFHEpp::AnnihilateKey<TargetP>>();
    TFHEpp::annihilatekeygen<TargetP>(*ahk, sk);

    std::cout << " done.\n";

    // ── Build LSB test vector: tv[i] = (i & 1) ? q/2 : 0 ──
    TFHEpp::Polynomial<TargetP> tv{};
    tv.fill(0);
    constexpr auto half_q =
        static_cast<typename TargetP::T>(1)
        << (std::numeric_limits<typename TargetP::T>::digits - 1);
    for (int i = 0; i < N; i++)
        tv[i] = (i & 1) ? half_q : 0;

    // ── Iteration schedule (K=2) ──
    const std::vector<int> betas = {8, 8};
    const std::vector<int> Ts = {63, 63};

    // ── Run test ──
    std::mt19937 rng(0x4C534231U);  // "LSB1"
    std::uniform_int_distribution<std::uint32_t> dist(0, static_cast<std::uint32_t>(t - 1));

    int failures = 0;
    int total = 0;

    for (int trial = 0; trial < trials; trial++) {
        const std::uint32_t m = (trials <= static_cast<int>(t))
            ? static_cast<std::uint32_t>(trial)  // full domain sweep
            : dist(rng);                          // random if > 2N trials

        TFHEpp::TLWE<DomainP> cin{};
        TFHEpp::tlweSymIntEncrypt<DomainP, static_cast<std::uint32_t>(t)>(
            cin, m, sk);

        TFHEpp::TLWE<TargetP> cout{};

        // Call public MetaPBS with period=2 (LSB LUT is 2-periodic)
        TFHEpp::metapbs::MetaPBS<bkP>(
            cout, cin, *bkfft, *ahk, tv, betas, Ts, t,
            /*period=*/2);

        // Decode: closer to 0 or q/2?
        const auto phase_u = TFHEpp::tlweSymPhase<TargetP>(
            cout, sk.key.get<TargetP>());
        constexpr std::uint32_t q_over_4 = 1U << 30;
        constexpr std::uint32_t three_q_over_4 = 3U << 30;
        const bool got_one =
            (phase_u >= q_over_4) && (phase_u < three_q_over_4);
        const bool exp_one = (m & 1) != 0;

        if (got_one != exp_one) {
            failures++;
            if (failures <= 10) {
                std::cerr << "  MISMATCH: m=" << m
                          << " exp=" << (exp_one ? "odd" : "even")
                          << " got=" << (got_one ? "odd" : "even")
                          << " phase=0x" << std::hex << phase_u << std::dec
                          << "\n";
            }
        }
        total++;
    }

    std::cout << "\n[metapbs_lsb_11bit] Result: "
              << failures << "/" << total << " failures";
    if (total > 0) {
        double rate = 100.0 * failures / total;
        std::cout << " (" << std::fixed << std::setprecision(2) << rate << "%)";
    }
    std::cout << "\n";

    if (failures == 0) {
        std::cout << "PASS: All " << total
                  << " trials matched (full 11-bit range, t=2N, period=2).\n";
        return 0;
    } else {
        std::cout << "FAIL: " << failures << " mismatches detected.\n";
        return 1;
    }
}
