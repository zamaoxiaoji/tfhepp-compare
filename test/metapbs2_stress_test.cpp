// =============================================================
// Stress Test: Meta-PBS bit extraction
// Tests exhaustive combinations and random samples of bit extraction
// =============================================================
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <random>

#include "cloudkey.hpp"
#include "metapbs2/bit_extraction.hpp"
#include "metapbs2/paper_params.hpp"

using namespace MetaPBS2;

static double elapsed_ms(std::chrono::steady_clock::time_point start,
                         std::chrono::steady_clock::time_point stop) {
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

static int kth_bit_period(int k) {
    return 1 << (k + 1);
}

// ---------------------------------------------------------------
// Stress test 1: Exhaustive plaintext data for the supported paper row.
// This covers every message and every supported k without running tens of
// thousands of encrypted bootstraps.
// ---------------------------------------------------------------
void stress_test_exhaustive_plaintext_data() {
    printf("[Stress Test] Exhaustive plaintext data: t=4096, N=2048, k=0..10\n");
    using tgtP = TFHEpp::lvl2param;
    constexpr int N = tgtP::n;
    constexpr int t = 2 * N;
    const int max_k = MaxExtractableBit(N);

    int total_pass = 0;
    int total_fail = 0;
    for (int k = 0; k <= max_k; k++) {
        auto tv = BuildKthBitTestVector<tgtP>(k, t);
        int pass = 0;
        int fail = 0;
        for (int m = 0; m < t; m++) {
            int idx = m % N;
            uint64_t got = tv[idx];
            uint64_t want = ((m >> k) & 1) ? BinaryScale : 0;
            if (got == want) {
                pass++;
            } else {
                if (fail == 0)
                    printf("\n    FAIL: k=%d m=%d idx=%d got=%lu want=%lu",
                           k, m, idx, got, want);
                fail++;
            }
        }
        printf("  bit %d period=%d: %s (%d/%d)\n",
               k, kth_bit_period(k), fail == 0 ? "OK" : "FAILED", pass, t);
        total_pass += pass;
        total_fail += fail;
    }
    printf("  Subtotal: %d PASSED, %d FAILED\n", total_pass, total_fail);
    if (total_fail > 0) exit(1);
}

// Define a custom bootstrapping parameter for lvl2 -> lvl2
struct br_lvl22param {
    using domainP = TFHEpp::lvl2param;
    using targetP = TFHEpp::lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

// ---------------------------------------------------------------
// Stress test 2: Sampled test on Paper params (lvl22)
// Testing t=4096 with representative edge cases and random sampling.
// ---------------------------------------------------------------
void stress_test_paper_params() {
    printf("[Stress Test] Paper params: t=4096, N=2048, K=2 (lvl22)\n");
    using tgtP = TFHEpp::lvl2param;
    using domP = TFHEpp::lvl2param;
    using brP = br_lvl22param;

    auto cfg = PaperRowT2NConfig();

    // Private-key holder / setup side: excluded from server-side extraction time.
    TFHEpp::SecretKey sk;
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    TFHEpp::bkfftgen<brP>(*bkfft, sk);

    const int B1 = cfg.rounds[0].beta;
    const int B2 = cfg.rounds[1].beta;

    auto trkey1 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B1);
    auto trkey2 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B2);
    std::vector<TruncRepeatKey<tgtP>> trkeys = {std::move(trkey1), std::move(trkey2)};

    const int max_k = MaxExtractableBit(tgtP::n);
    std::mt19937 rng(1337);

    int total_pass = 0;
    int total_fail = 0;
    double total_baseline_extract_ms = 0.0;
    double total_pruned_extract_ms = 0.0;
    BlindRotatePruneStats all_stats{};

    for (int k = 0; k <= max_k; k++) {
        std::vector<int> msgs_to_test = {0, 1, cfg.t - 1, (1 << k), (1 << k) - 1};
        if (k < 11) msgs_to_test.push_back((1 << (k + 1)) - 1);
        for(int i=0; i<10; i++) {
            msgs_to_test.push_back(rng() % cfg.t);
        }

        struct EncryptedSample {
            int message;
            int expected_bit;
            TFHEpp::TLWE<domP> ct;
        };

        // Input ciphertext preparation is client/private-key-side test setup.
        // It is intentionally outside the timed server extraction window.
        std::vector<EncryptedSample> samples;
        samples.reserve(msgs_to_test.size());
        for (int m : msgs_to_test) {
            EncryptedSample sample{};
            sample.message = m;
            sample.expected_bit = (m >> k) & 1;
            TFHEpp::tlweSymEncrypt<domP>(
                sample.ct,
                static_cast<typename domP::T>(m) *
                    (((~typename domP::T(0)) / typename domP::T(cfg.t)) + 1),
                domP::α, sk.key.get<domP>());
            samples.push_back(sample);
        }

        int pass = 0, fail = 0;
        double baseline_extract_ms = 0.0;
        double pruned_extract_ms = 0.0;
        BlindRotatePruneStats bit_stats{};
        printf("  Testing bit %d period=%d (%zu samples)... ",
               k, kth_bit_period(k), samples.size());
        fflush(stdout);

        for (const auto& sample : samples) {
            // Server-side timed window: from input ciphertext to output bit ciphertext.
            // It includes bootstrap, TruncRepeat correction rounds, and optional pruning.
            auto t0 = std::chrono::steady_clock::now();
            auto cout_baseline = ExtractKthBit<brP>(
                sample.ct, k, cfg.t, *bkfft, trkeys, cfg);
            auto t1 = std::chrono::steady_clock::now();

            BlindRotatePruneStats stats;
            auto t2 = std::chrono::steady_clock::now();
            auto cout_pruned = ExtractKthBit<brP>(
                sample.ct, k, cfg.t, *bkfft, trkeys, cfg,
                kth_bit_period(k), &stats);
            auto t3 = std::chrono::steady_clock::now();

            baseline_extract_ms += elapsed_ms(t0, t1);
            pruned_extract_ms += elapsed_ms(t2, t3);
            bit_stats.total += stats.total;
            bit_stats.skipped += stats.skipped;

            // Verification/decryption is test-side only and is outside timing.
            typename tgtP::T baseline_phase =
                TFHEpp::tlweSymPhase<tgtP>(cout_baseline, sk.key.get<tgtP>());
            typename tgtP::T pruned_phase =
                TFHEpp::tlweSymPhase<tgtP>(cout_pruned, sk.key.get<tgtP>());
            int baseline_bit = DecodeBinaryPhase(baseline_phase);
            int pruned_bit = DecodeBinaryPhase(pruned_phase);

            if (baseline_bit == sample.expected_bit && pruned_bit == sample.expected_bit) {
                pass++;
            } else {
                printf("\n    FAIL: m=%d, expected=%d, baseline=%d pruned=%d "
                       "(baseline_phase=%lu pruned_phase=%lu Q/2=%lu)",
                       sample.message, sample.expected_bit, baseline_bit, pruned_bit,
                       baseline_phase, pruned_phase, BinaryScaleT<typename tgtP::T>);
                fail++;
            }
        }
        
        double speedup = pruned_extract_ms > 0.0
                             ? baseline_extract_ms / pruned_extract_ms
                             : 0.0;
        if (fail == 0)
            printf("OK prune=%.2f%% skipped=%lu/%lu "
                   "server_extract_baseline=%.2fms server_extract_pruned=%.2fms "
                   "speedup=%.2fx\n",
                   100.0 * bit_stats.prune_rate(),
                   bit_stats.skipped, bit_stats.total,
                   baseline_extract_ms, pruned_extract_ms, speedup);
        else
            printf("FAILED\n");

        total_pass += pass;
        total_fail += fail;
        total_baseline_extract_ms += baseline_extract_ms;
        total_pruned_extract_ms += pruned_extract_ms;
        all_stats.total += bit_stats.total;
        all_stats.skipped += bit_stats.skipped;
    }
    
    printf("  Subtotal: %d PASSED, %d FAILED\n", total_pass, total_fail);
    printf("  Total prune=%.2f%% skipped=%lu/%lu "
           "server_extract_baseline=%.2fms server_extract_pruned=%.2fms "
           "speedup=%.2fx\n",
           100.0 * all_stats.prune_rate(), all_stats.skipped, all_stats.total,
           total_baseline_extract_ms, total_pruned_extract_ms,
           total_pruned_extract_ms > 0.0
               ? total_baseline_extract_ms / total_pruned_extract_ms
               : 0.0);
    if (total_fail > 0) exit(1);
}

int main() {
    printf("=== MetaPBS2 Stress Test ===\n");
    stress_test_exhaustive_plaintext_data();
    stress_test_paper_params();
    printf("=== ALL STRESS TESTS PASSED ===\n");
    return 0;
}
