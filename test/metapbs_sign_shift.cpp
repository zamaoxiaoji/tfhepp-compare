// TFHEpp/test/metapbs_sign_shift.cpp
// Debugging test for MetaPBS sign detection.
//
// Tests MetaPBS with various TV constructions to find the correct
// encoding for sign detection over Z_{2N}.

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <metapbs.hpp>
#include <tfhe++.hpp>

int main(int argc, char **argv)
{
    using bkP     = TFHEpp::lvlh1param;
    using TargetP = typename bkP::targetP;   // lvl1
    using DomainP = typename bkP::domainP;   // lvlhalf

    constexpr int N    = static_cast<int>(TargetP::n);
    constexpr int nbit = static_cast<int>(TargetP::nbit);
    constexpr std::int64_t twoN = 2 * N;

    int trials = 100;
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--trials" && i + 1 < argc)
            trials = std::stoi(argv[++i]);
    }

    std::cout << "[metapbs_sign_shift] N=" << N << " nbit=" << nbit
              << " μ=0x" << std::hex << TargetP::μ << std::dec
              << " trials=" << trials << "\n";

    // Key generation
    std::cout << "Generating keys..." << std::flush;
    TFHEpp::SecretKey sk;
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<bkP>>();
    TFHEpp::bkfftgen<bkP>(*bkfft, sk);
    auto ahk = std::make_unique<TFHEpp::AnnihilateKey<TargetP>>();
    TFHEpp::annihilatekeygen<TargetP>(*ahk, sk);
    std::cout << " done.\n";

    const std::vector<int> betas = {8, 8};
    const std::vector<int> Ts = {(N / 8 - 1) / 2, (N / 8 - 1) / 2};

    // ─── Approach A: Use MetaPBSExtractBit2N directly (bit = nbit-1 = MSB)
    //     This is known to work — it extracts bit nbit-1 from Z_{2N}.
    //     Output: {0, q/2}
    {
        std::cout << "\n  ── A: MetaPBSExtractBit2N (bit=" << (nbit-1)
                  << ", MSB) ──\n";
        int failures = 0;

        for (int trial = 0; trial < trials; trial++) {
            const uint32_t m = static_cast<uint32_t>(trial % static_cast<int>(twoN));

            TFHEpp::TLWE<DomainP> cin{};
            TFHEpp::tlweSymIntEncrypt<DomainP, static_cast<uint32_t>(twoN)>(
                cin, static_cast<typename DomainP::T>(m), sk);

            TFHEpp::TLWE<TargetP> result{};
            TFHEpp::metapbs::detail::MetaPBSExtractBit2N<bkP>(
                result, cin, *bkfft, *ahk, nbit - 1);

            // Decode {0, q/2}: check if in [q/4, 3q/4)
            auto phase = TFHEpp::tlweSymPhase<TargetP>(
                result, sk.key.get<TargetP>());
            constexpr uint32_t q4 = 1U << 30;
            constexpr uint32_t q34 = 3U << 30;
            bool got_msb_set = (phase >= q4) && (phase < q34);
            bool exp_msb_set = (m >= static_cast<uint32_t>(N));

            if (got_msb_set != exp_msb_set) {
                failures++;
                if (failures <= 5)
                    std::cerr << "    MISMATCH: m=" << m
                              << " exp=" << exp_msb_set
                              << " got=" << got_msb_set
                              << " phase=0x" << std::hex << phase << std::dec
                              << "\n";
            }
        }
        std::cout << "  Result: " << failures << "/" << trials << " failures\n"
                  << (failures == 0 ? "  PASS ✅\n" : "  FAIL ❌\n");
    }

    // ─── Approach B: MetaPBS with t=2N, TV = BitExtractTV for MSB
    //     Same TV as approach A, but going through the public MetaPBS API.
    //     period = 2^{nbit} = N.
    {
        // Build the same TV as BitExtractTestVector2N(nbit-1)
        constexpr auto half_q = static_cast<typename TargetP::T>(1)
            << (std::numeric_limits<typename TargetP::T>::digits - 1);
        TFHEpp::Polynomial<TargetP> tv{};
        tv.fill(0);
        for (int i = 0; i < N; i++)
            tv[i] = ((i >> (nbit - 1)) & 1) ? half_q : 0;

        const int period = 1 << nbit;  // = N

        std::cout << "\n  ── B: MetaPBS (t=2N), BitExtract TV for MSB"
                  << " (period=" << period << ") ──\n";
        int failures = 0;

        for (int trial = 0; trial < trials; trial++) {
            const uint32_t m = static_cast<uint32_t>(trial % static_cast<int>(twoN));

            TFHEpp::TLWE<DomainP> cin{};
            TFHEpp::tlweSymIntEncrypt<DomainP, static_cast<uint32_t>(twoN)>(
                cin, static_cast<typename DomainP::T>(m), sk);

            TFHEpp::TLWE<TargetP> result{};
            TFHEpp::metapbs::MetaPBS<bkP>(
                result, cin, *bkfft, *ahk, tv, betas, Ts,
                /*t=*/twoN, /*period=*/period);

            auto phase = TFHEpp::tlweSymPhase<TargetP>(
                result, sk.key.get<TargetP>());
            constexpr uint32_t q4 = 1U << 30;
            constexpr uint32_t q34 = 3U << 30;
            bool got_msb_set = (phase >= q4) && (phase < q34);
            bool exp_msb_set = (m >= static_cast<uint32_t>(N));

            if (got_msb_set != exp_msb_set) {
                failures++;
                if (failures <= 5)
                    std::cerr << "    MISMATCH: m=" << m
                              << " exp=" << exp_msb_set
                              << " got=" << got_msb_set
                              << " phase=0x" << std::hex << phase << std::dec
                              << "\n";
            }
        }
        std::cout << "  Result: " << failures << "/" << trials << " failures\n"
                  << (failures == 0 ? "  PASS ✅\n" : "  FAIL ❌\n");
    }

    // ─── Approach C: Same as B, but test as comparison (a > b ?) ────────
    //     For 4-bit and 8-bit values.
    {
        const std::vector<int> bit_cases = {4, 8, 10};

        constexpr auto half_q = static_cast<typename TargetP::T>(1)
            << (std::numeric_limits<typename TargetP::T>::digits - 1);
        TFHEpp::Polynomial<TargetP> tv{};
        tv.fill(0);
        for (int i = 0; i < N; i++)
            tv[i] = ((i >> (nbit - 1)) & 1) ? half_q : 0;
        const int period = 1 << nbit;

        for (const int bits : bit_cases) {
            if (bits + 1 > nbit + 1) {
                std::cout << "\n  " << bits << "-bit: skipped (out of window)\n";
                continue;
            }

            std::cout << "\n  ── C: " << bits
                      << "-bit comparison via MetaPBS Sign ──\n";

            const uint32_t val_max = (1u << (bits - 1)) - 1;
            std::mt19937 rng(0x434D5033U + static_cast<uint32_t>(bits));
            std::uniform_int_distribution<uint32_t> dist(0, val_max);

            int failures = 0;
            double total_ms = 0.0;

            for (int trial = 0; trial < trials; trial++) {
                const uint32_t a = dist(rng);
                const uint32_t b = dist(rng);
                if (a == b) continue;
                const bool exp_gt = (a > b);

                // diff = b - a, wrap to [0, 2N)
                int32_t d = static_cast<int32_t>(b) - static_cast<int32_t>(a);
                int32_t d_mod = d % static_cast<int32_t>(twoN);
                if (d_mod < 0) d_mod += static_cast<int32_t>(twoN);

                TFHEpp::TLWE<DomainP> cin{};
                TFHEpp::tlweSymIntEncrypt<DomainP, static_cast<uint32_t>(twoN)>(
                    cin, static_cast<typename DomainP::T>(d_mod), sk);

                TFHEpp::TLWE<TargetP> result{};
                auto t0 = std::chrono::steady_clock::now();
                TFHEpp::metapbs::MetaPBS<bkP>(
                    result, cin, *bkfft, *ahk, tv, betas, Ts,
                    /*t=*/twoN, /*period=*/period);
                auto t1 = std::chrono::steady_clock::now();
                total_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

                // Decode {0, q/2}
                auto phase = TFHEpp::tlweSymPhase<TargetP>(
                    result, sk.key.get<TargetP>());
                constexpr uint32_t q4 = 1U << 30;
                constexpr uint32_t q34 = 3U << 30;
                bool got_msb_set = (phase >= q4) && (phase < q34);
                // a > b ⟺ d = b-a < 0 ⟺ d_mod ≥ N ⟺ MSB set
                bool got_gt = got_msb_set;

                if (got_gt != exp_gt) {
                    failures++;
                    if (failures <= 5) {
                        std::cerr << "    MISMATCH: a=" << a << " b=" << b
                                  << " d=" << d << " d_mod=" << d_mod
                                  << " got_msb=" << got_msb_set
                                  << " phase=0x" << std::hex << phase
                                  << std::dec << "\n";
                    }
                }
            }

            double avg = (trials > 0) ? total_ms / trials : 0.0;
            std::cout << "  Result: " << failures << "/" << trials
                      << " failures  avg=" << std::fixed
                      << std::setprecision(1) << avg << "ms\n"
                      << (failures == 0 ? "  PASS ✅\n" : "  FAIL ❌\n");
        }
    }

    std::cout << "\n[metapbs_sign_shift] done.\n";
    return 0;
}
