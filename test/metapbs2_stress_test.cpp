// =============================================================
// Stress Test: Meta-PBS bit extraction
// Tests exhaustive combinations and random samples of bit extraction
// =============================================================
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <stdexcept>
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

static int chapter_bit_from_lsb_bit(int p, int lsb_k) {
    return p - 1 - lsb_k;
}

static uint64_t abs_diff_u64(uint64_t a, uint64_t b) {
    return a >= b ? a - b : b - a;
}

static uint64_t binary_decision_margin(uint64_t phase) {
    const uint64_t q4 = BinaryScale / 2;
    const uint64_t three_q4 = BinaryScale + q4;
    uint64_t d1 = abs_diff_u64(phase, q4);
    uint64_t d2 = abs_diff_u64(phase, three_q4);
    return d1 < d2 ? d1 : d2;
}

static double normalized_binary_margin(uint64_t phase) {
    return static_cast<double>(binary_decision_margin(phase)) /
           static_cast<double>(BinaryScale / 2);
}

template <typename Fn>
static void expect_invalid_argument(const char* name, Fn&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        printf("  %s: rejected\n", name);
        return;
    }
    printf("  FAIL: %s was accepted\n", name);
    exit(1);
}

static void test_bitextract_parameter_robustness() {
    printf("[Robustness] BitExtract parameter checks\n");
    if (MessagePrecisionFromPowerOfTwoModulus(4096) != 12) {
        printf("  FAIL: p(log2 4096) should be 12\n");
        exit(1);
    }
    if (LSBIndexFromChapterBit(12, 11) != 0 ||
        BitExtractPeriodFromChapterBit(12, 11) != 2) {
        printf("  FAIL: Chapter bit k=11 should map to LSB and period=2\n");
        exit(1);
    }
    if (LSBIndexFromChapterBit(12, 10) != 1 ||
        BitExtractPeriodFromChapterBit(12, 10) != 4) {
        printf("  FAIL: Chapter bit k=10 should map to LSB-index 1 and period=4\n");
        exit(1);
    }

    expect_invalid_argument("non-power-of-two t", [] {
        (void)MessagePrecisionFromPowerOfTwoModulus(4095);
    });
    expect_invalid_argument("negative k", [] {
        (void)LSBIndexFromChapterBit(12, -1);
    });
    expect_invalid_argument("k >= p", [] {
        (void)LSBIndexFromChapterBit(12, 12);
    });
    expect_invalid_argument("unsupported top bit", [] {
        (void)BuildKthBitTestVector<TFHEpp::lvl2param>(11, 2 * TFHEpp::lvl2param::n);
    });
    printf("  PASSED\n");
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
        const int exact_period = ExactNegacyclicPeriod<tgtP>(tv);
        if (exact_period != kth_bit_period(k)) {
            printf("  FAIL: bit %d exact slot period=%d message-period=%d\n",
                   k, exact_period, kth_bit_period(k));
            exit(1);
        }
        for (int shift = 0; shift < 2 * N; shift += exact_period) {
            if (!NegacyclicRotationInvariant<tgtP>(tv, shift)) {
                printf("  FAIL: bit %d shift %d should be invariant\n", k, shift);
                exit(1);
            }
        }
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
        printf("  bit %d exact_period=%d: %s (%d/%d)\n",
               k, exact_period, fail == 0 ? "OK" : "FAILED", pass, t);
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

static void test_bitextract_wrapper_rejections() {
    printf("[Robustness] BitExtract wrapper rejection checks\n");
    using brP = br_lvl22param;
    using domP = TFHEpp::lvl2param;
    using tgtP = TFHEpp::lvl2param;

    auto cfg = PaperRowT2NConfig();
    const int p = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    TFHEpp::TLWE<domP> ct{};
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    std::vector<TruncRepeatKey<tgtP>> trkeys;

    expect_invalid_argument("unsupported Chapter/MSB top bit", [&] {
        (void)BitExtractBoolPruned<brP>(
            ct, *bkfft, trkeys, cfg,
            BitExtractOptions{.p = p, .k = 0});
    });
    expect_invalid_argument("period override not matching concrete TV", [&] {
        (void)BitExtractBoolPruned<brP>(
            ct, *bkfft, trkeys, cfg,
            BitExtractOptions{.p = p, .k = p - 1, .period = 3});
    });
    printf("  PASSED\n");
}

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
    const int p = MessagePrecisionFromPowerOfTwoModulus(cfg.t);

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
    double total_baseline_margin = 0.0;
    double total_pruned_margin = 0.0;
    double min_baseline_margin = 1.0;
    double min_pruned_margin = 1.0;
    int margin_count = 0;

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
        double bit_baseline_margin = 0.0;
        double bit_pruned_margin = 0.0;
        double bit_min_baseline_margin = 1.0;
        double bit_min_pruned_margin = 1.0;
        printf("  Testing bit %d period=%d (%zu samples)... ",
               k, kth_bit_period(k), samples.size());
        fflush(stdout);

        for (const auto& sample : samples) {
            // Server-side timed window: from input ciphertext to output bit ciphertext.
            // It includes bootstrap, TruncRepeat correction rounds, and optional pruning.
            auto t0 = std::chrono::steady_clock::now();
            auto cout_baseline = BitExtractBoolPruned<brP>(
                sample.ct, *bkfft, trkeys, cfg,
                BitExtractOptions{
                    .p = p,
                    .k = chapter_bit_from_lsb_bit(p, k),
                    .enable_periodic_pruning = false,
                });
            auto t1 = std::chrono::steady_clock::now();

            BlindRotatePruneStats stats;
            auto t2 = std::chrono::steady_clock::now();
            auto cout_pruned = BitExtractBoolPruned<brP>(
                sample.ct, *bkfft, trkeys, cfg,
                BitExtractOptions{
                    .p = p,
                    .k = chapter_bit_from_lsb_bit(p, k),
                    .enable_periodic_pruning = true,
                },
                &stats);
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
            double baseline_margin = normalized_binary_margin(baseline_phase);
            double pruned_margin = normalized_binary_margin(pruned_phase);
            bit_baseline_margin += baseline_margin;
            bit_pruned_margin += pruned_margin;
            bit_min_baseline_margin =
                baseline_margin < bit_min_baseline_margin ? baseline_margin : bit_min_baseline_margin;
            bit_min_pruned_margin =
                pruned_margin < bit_min_pruned_margin ? pruned_margin : bit_min_pruned_margin;

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
                   "speedup=%.2fx margin_min(base/pruned)=%.3f/%.3f "
                   "margin_avg(base/pruned)=%.3f/%.3f\n",
                   100.0 * bit_stats.prune_rate(),
                   bit_stats.skipped, bit_stats.total,
                   baseline_extract_ms, pruned_extract_ms, speedup,
                   bit_min_baseline_margin, bit_min_pruned_margin,
                   bit_baseline_margin / pass, bit_pruned_margin / pass);
        else
            printf("FAILED\n");

        total_pass += pass;
        total_fail += fail;
        total_baseline_margin += bit_baseline_margin;
        total_pruned_margin += bit_pruned_margin;
        min_baseline_margin = bit_min_baseline_margin < min_baseline_margin
                                  ? bit_min_baseline_margin
                                  : min_baseline_margin;
        min_pruned_margin = bit_min_pruned_margin < min_pruned_margin
                                ? bit_min_pruned_margin
                                : min_pruned_margin;
        margin_count += pass;
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
    printf("  Accuracy: %.6f%% (%d/%d)\n",
           100.0 * static_cast<double>(total_pass) /
               static_cast<double>(total_pass + total_fail),
           total_pass, total_pass + total_fail);
    printf("  Robustness margin: min(base/pruned)=%.3f/%.3f avg(base/pruned)=%.3f/%.3f\n",
           min_baseline_margin, min_pruned_margin,
           total_baseline_margin / margin_count, total_pruned_margin / margin_count);
    if (min_baseline_margin <= 0.05 || min_pruned_margin <= 0.05) {
        printf("  FAIL: binary decision margin below 5%% of half-cell radius\n");
        exit(1);
    }
    if (total_fail > 0) exit(1);
}

void robustness_test_multikey_samples() {
    printf("[Robustness] Multi-key sampled correctness (pruned BitExtract)\n");
    using tgtP = TFHEpp::lvl2param;
    using domP = TFHEpp::lvl2param;
    using brP = br_lvl22param;

    auto cfg = PaperRowT2NConfig();
    const int p = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    const int B1 = cfg.rounds[0].beta;
    const int B2 = cfg.rounds[1].beta;
    const std::vector<int> lsb_bits = {0, 1, 2, 5, 10};
    const std::vector<int> base_messages = {0, 1, cfg.t - 1, cfg.t / 2, cfg.t / 2 - 1};

    int total_pass = 0;
    int total_fail = 0;
    double min_margin = 1.0;
    double sum_margin = 0.0;

    for (int key_round = 0; key_round < 3; key_round++) {
        TFHEpp::SecretKey sk;
        auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
        TFHEpp::bkfftgen<brP>(*bkfft, sk);
        auto trkey1 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B1);
        auto trkey2 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B2);
        std::vector<TruncRepeatKey<tgtP>> trkeys = {std::move(trkey1), std::move(trkey2)};

        std::mt19937 rng(9000 + key_round);
        int round_pass = 0;
        int round_fail = 0;
        for (int lsb_k : lsb_bits) {
            std::vector<int> messages = base_messages;
            messages.push_back(1 << lsb_k);
            messages.push_back((1 << (lsb_k + 1)) - 1);
            messages.push_back(rng() % cfg.t);
            messages.push_back(rng() % cfg.t);

            for (int m : messages) {
                TFHEpp::TLWE<domP> ct;
                TFHEpp::tlweSymEncrypt<domP>(
                    ct,
                    static_cast<typename domP::T>(m) *
                        (((~typename domP::T(0)) / typename domP::T(cfg.t)) + 1),
                    domP::α, sk.key.get<domP>());

                BlindRotatePruneStats stats;
                auto cout = BitExtractBoolPruned<brP>(
                    ct, *bkfft, trkeys, cfg,
                    BitExtractOptions{
                        .p = p,
                        .k = chapter_bit_from_lsb_bit(p, lsb_k),
                        .enable_periodic_pruning = true,
                    },
                    &stats);
                auto phase = TFHEpp::tlweSymPhase<tgtP>(cout, sk.key.get<tgtP>());
                int got = DecodeBinaryPhase(phase);
                int want = (m >> lsb_k) & 1;
                double margin = normalized_binary_margin(phase);
                min_margin = margin < min_margin ? margin : min_margin;
                sum_margin += margin;
                if (got == want) {
                    round_pass++;
                } else {
                    printf("  FAIL key_round=%d bit=%d m=%d want=%d got=%d margin=%.3f\n",
                           key_round, lsb_k, m, want, got, margin);
                    round_fail++;
                }
            }
        }
        total_pass += round_pass;
        total_fail += round_fail;
        printf("  key_round=%d accuracy=%.6f%% (%d/%d)\n",
               key_round,
               100.0 * static_cast<double>(round_pass) /
                   static_cast<double>(round_pass + round_fail),
               round_pass, round_pass + round_fail);
    }

    printf("  Total accuracy=%.6f%% (%d/%d), margin_min=%.3f margin_avg=%.3f\n",
           100.0 * static_cast<double>(total_pass) /
               static_cast<double>(total_pass + total_fail),
           total_pass, total_pass + total_fail,
           min_margin,
           sum_margin / static_cast<double>(total_pass + total_fail));
    const double accuracy =
        static_cast<double>(total_pass) /
        static_cast<double>(total_pass + total_fail);
    if (accuracy < 0.99) {
        printf("  FAIL: multi-key accuracy below 99%% threshold\n");
        exit(1);
    }
    if (min_margin <= 0.05) {
        printf("  FAIL: multi-key margin below 5%% of half-cell radius\n");
        exit(1);
    }
}

int main() {
    printf("=== MetaPBS2 Stress Test ===\n");
    test_bitextract_parameter_robustness();
    test_bitextract_wrapper_rejections();
    stress_test_exhaustive_plaintext_data();
    stress_test_paper_params();
    robustness_test_multikey_samples();
    printf("=== ALL STRESS TESTS PASSED ===\n");
    return 0;
}
