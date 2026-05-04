// =============================================================
// Stress/robustness test for MetaPBS2 GapMSB.
//
// Covers the Chapter_3_revised (1).tex alg:et-hmsb-expanded flow:
//   BitExtract(0/Q/2) -> logical-to-arithmetic weight conversion ->
//   clear bit -> offset -> sign PBS.
// =============================================================
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "cloudkey.hpp"
#include "metapbs2/gap_msb.hpp"
#include "metapbs2/paper_params.hpp"

using namespace MetaPBS2;

struct br_lvl22param {
    using domainP = TFHEpp::lvl2param;
    using targetP = TFHEpp::lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

using brP = br_lvl22param;
using P = TFHEpp::lvl2param;

static double elapsed_ms(std::chrono::steady_clock::time_point start,
                         std::chrono::steady_clock::time_point stop) {
    return std::chrono::duration<double, std::milli>(stop - start).count();
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

template <class Param>
static typename Param::T encode_message(int m, int p) {
    const typename Param::T delta =
        TorusScaleForPrecision<typename Param::T>(p);
    return static_cast<typename Param::T>(m) * delta;
}

template <class Param>
static int decode_arithmetic_message(typename Param::T phase, int p) {
    const int shift = std::numeric_limits<typename Param::T>::digits - p;
    const typename Param::T half_delta =
        typename Param::T(1) << (shift - 1);
    const typename Param::T mask = (typename Param::T(1) << p) - 1;
    return static_cast<int>(((phase + half_delta) >> shift) & mask);
}

template <class Param>
static double arithmetic_cell_margin(
    typename Param::T phase, int expected_message, int p) {
    using T = typename Param::T;
    using SignedT = std::make_signed_t<T>;
    const int shift = std::numeric_limits<T>::digits - p;
    const T delta = T(1) << shift;
    const T target = static_cast<T>(expected_message) * delta;
    const SignedT diff = static_cast<SignedT>(phase - target);
    const auto abs_diff = diff < 0
                              ? static_cast<unsigned long long>(-diff)
                              : static_cast<unsigned long long>(diff);
    const double radius = static_cast<double>(delta >> 1);
    return 1.0 - static_cast<double>(abs_diff) / radius;
}

template <class Param>
static double sign_margin(typename Param::T phase) {
    using SignedT = std::make_signed_t<typename Param::T>;
    const SignedT signed_phase = static_cast<SignedT>(phase);
    const auto abs_phase = signed_phase < 0
                               ? static_cast<unsigned long long>(-signed_phase)
                               : static_cast<unsigned long long>(signed_phase);
    return static_cast<double>(abs_phase) /
           static_cast<double>(Param::μ);
}

static void test_parameter_helpers() {
    printf("[GapMSB helpers] parameter mapping\n");
    const int p = 12;
    if (ArithmeticWeightForChapterBit<uint64_t>(p, 7) !=
        uint64_t(16) * TorusScaleForPrecision<uint64_t>(p)) {
        printf("  FAIL: k=7 weight should be 16*Delta\n");
        exit(1);
    }
    if (HalfGapOffsetForChapterBit<uint64_t>(p, 7) !=
        uint64_t(17) * (TorusScaleForPrecision<uint64_t>(p) >> 1)) {
        printf("  FAIL: k=7 offset should be 17*Delta/2\n");
        exit(1);
    }
    expect_invalid_argument("GapMSB top bit k=0", [] {
        (void)ArithmeticWeightForChapterBit<uint64_t>(12, 0);
    });
    printf("  PASSED\n");
}

static void test_direct_logical_to_arithmetic(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft) {
    printf("[LOG_to_ARI equivalent] 0/Q/2 bit -> 0/w\n");
    const int p = 12;
    const int chapter_k = 7;
    const int lsb_k = LSBIndexFromChapterBit(p, chapter_k);
    const int weight = 1 << lsb_k;
    const auto weight_torus = ArithmeticWeightForChapterBit<typename P::T>(p, chapter_k);
    double min_margin = 1.0;

    for (int bit = 0; bit <= 1; bit++) {
        TFHEpp::TLWE<P> bit_ct{};
        TFHEpp::tlweSymEncrypt<P>(
            bit_ct,
            bit ? BinaryScaleT<typename P::T> : typename P::T(0),
            P::α, sk.key.get<P>());

        auto weighted_ct =
            LogicalBitToArithmeticWeight<brP>(bit_ct, weight_torus, bkfft);
        auto phase = TFHEpp::tlweSymPhase<P>(weighted_ct, sk.key.get<P>());
        const int got = decode_arithmetic_message<P>(phase, p);
        const int want = bit ? weight : 0;
        const double margin = arithmetic_cell_margin<P>(phase, want, p);
        min_margin = margin < min_margin ? margin : min_margin;
        if (got != want) {
            printf("  FAIL bit=%d got=%d want=%d margin=%.3f\n",
                   bit, got, want, margin);
            exit(1);
        }
        printf("  bit=%d -> %d OK margin=%.3f\n", bit, got, margin);
    }
    if (min_margin <= 0.05) {
        printf("  FAIL: arithmetic conversion margin below 5%%\n");
        exit(1);
    }
    printf("  PASSED\n");
}

static void test_gapmsb_parameter_rejections(
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<P>>& trkeys,
    const Algorithm1Config& cfg) {
    printf("[GapMSB robustness] parameter rejections\n");
    TFHEpp::TLWE<P> ct{};
    expect_invalid_argument("k=0 unsupported", [&] {
        (void)GapMSB<brP>(ct, bkfft, trkeys, cfg,
                          GapMSBOptions{.p = 12, .k = 0});
    });
    expect_invalid_argument("k>=p unsupported", [&] {
        (void)GapMSB<brP>(ct, bkfft, trkeys, cfg,
                          GapMSBOptions{.p = 12, .k = 12});
    });
    expect_invalid_argument("cfg.t != 2^p", [&] {
        (void)GapMSB<brP>(ct, bkfft, trkeys, cfg,
                          GapMSBOptions{.p = 11, .k = 7});
    });
    expect_invalid_argument("bad pruning period", [&] {
        (void)GapMSB<brP>(
            ct, bkfft, trkeys, cfg,
            GapMSBOptions{.p = 12, .k = 7, .period = 3});
    });
    printf("  PASSED\n");
}

static void test_weighted_extract_and_clear(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<P>>& trkeys,
    const Algorithm1Config& cfg) {
    printf("[GapMSB clear] weighted BitExtract and bit clearing\n");
    const int p = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    const int chapter_k = 7;  // w_k=16 for p=12, matching the tex example.
    const int lsb_k = LSBIndexFromChapterBit(p, chapter_k);
    const int weight = 1 << lsb_k;
    const std::vector<int> messages = {
        0, 1, weight - 1, weight, 2 * weight - 1,
        cfg.t / 2 - 1 - weight, cfg.t / 2 - 1, cfg.t / 2,
        cfg.t / 2 + weight, 3261, cfg.t - 1,
    };
    double min_weight_margin = 1.0;
    double min_clear_margin = 1.0;

    for (int m : messages) {
        TFHEpp::TLWE<P> ct{};
        TFHEpp::tlweSymEncrypt<P>(
            ct, encode_message<P>(m, p), P::α, sk.key.get<P>());

        BlindRotatePruneStats stats;
        auto weighted_ct = ExtractWeightedChapterBit<brP>(
            ct, bkfft, trkeys, cfg,
            GapMSBOptions{.p = p, .k = chapter_k},
            &stats);
        auto weighted_phase = TFHEpp::tlweSymPhase<P>(weighted_ct, sk.key.get<P>());
        const int want_weight = ((m >> lsb_k) & 1) ? weight : 0;
        const int got_weight = decode_arithmetic_message<P>(weighted_phase, p);
        const double weight_margin =
            arithmetic_cell_margin<P>(weighted_phase, want_weight, p);
        min_weight_margin =
            weight_margin < min_weight_margin ? weight_margin : min_weight_margin;
        if (got_weight != want_weight) {
            printf("  FAIL weight m=%d got=%d want=%d margin=%.3f\n",
                   m, got_weight, want_weight, weight_margin);
            exit(1);
        }

        TFHEpp::TLWE<P> cleared{};
        ClearChapterBitAssign<brP>(cleared, ct, weighted_ct);
        auto clear_phase = TFHEpp::tlweSymPhase<P>(cleared, sk.key.get<P>());
        const int want_clear = (m - want_weight) & (cfg.t - 1);
        const int got_clear = decode_arithmetic_message<P>(clear_phase, p);
        const double clear_margin =
            arithmetic_cell_margin<P>(clear_phase, want_clear, p);
        min_clear_margin =
            clear_margin < min_clear_margin ? clear_margin : min_clear_margin;
        if (got_clear != want_clear) {
            printf("  FAIL clear m=%d got=%d want=%d margin=%.3f\n",
                   m, got_clear, want_clear, clear_margin);
            exit(1);
        }
    }

    printf("  samples=%zu min_margin(weight/clear)=%.3f/%.3f\n",
           messages.size(), min_weight_margin, min_clear_margin);
    if (min_weight_margin <= 0.05 || min_clear_margin <= 0.05) {
        printf("  FAIL: weighted extraction/clear margin below 5%%\n");
        exit(1);
    }
    printf("  PASSED\n");
}

static void stress_gapmsb_correctness(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<P>>& trkeys,
    const Algorithm1Config& cfg) {
    printf("[GapMSB stress] boundary and random samples\n");
    const int p = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    const std::vector<int> chapter_bits = {1, 4, 7, 10, 11};
    std::mt19937 rng(424242);
    int total_pass = 0;
    int total_fail = 0;
    double min_sign_margin = 10.0;
    double total_ms = 0.0;
    BlindRotatePruneStats all_stats{};

    for (int chapter_k : chapter_bits) {
        const int lsb_k = LSBIndexFromChapterBit(p, chapter_k);
        const int weight = 1 << lsb_k;
        std::vector<int> messages = {
            0, 1, weight - 1, weight,
            cfg.t / 2 - 1 - weight, cfg.t / 2 - 1,
            cfg.t / 2, cfg.t / 2 + weight, cfg.t - 1,
        };
        for (int i = 0; i < 4; i++) messages.push_back(rng() % cfg.t);

        int pass = 0;
        int fail = 0;
        double bit_ms = 0.0;
        BlindRotatePruneStats bit_stats{};
        printf("  k=%d w=%d samples=%zu ... ", chapter_k, weight, messages.size());
        fflush(stdout);
        for (int m : messages) {
            TFHEpp::TLWE<P> ct{};
            TFHEpp::tlweSymEncrypt<P>(
                ct, encode_message<P>(m, p), P::α, sk.key.get<P>());

            BlindRotatePruneStats stats;
            auto t0 = std::chrono::steady_clock::now();
            auto out = GapMSB<brP>(
                ct, bkfft, trkeys, cfg,
                GapMSBOptions{.p = p, .k = chapter_k},
                &stats);
            auto t1 = std::chrono::steady_clock::now();
            bit_ms += elapsed_ms(t0, t1);
            bit_stats.total += stats.total;
            bit_stats.skipped += stats.skipped;

            auto phase = TFHEpp::tlweSymPhase<P>(out, sk.key.get<P>());
            const int got = DecodeSignPhase<P>(phase);
            const int want = m >= cfg.t / 2 ? 1 : 0;
            const double margin = sign_margin<P>(phase);
            min_sign_margin = margin < min_sign_margin ? margin : min_sign_margin;
            if (got == want) {
                pass++;
            } else {
                printf("\n    FAIL k=%d m=%d got=%d want=%d margin=%.3f",
                       chapter_k, m, got, want, margin);
                fail++;
            }
        }
        printf("%s prune=%.2f%% skipped=%lu/%lu time=%.2fms\n",
               fail == 0 ? "OK" : "FAILED",
               100.0 * bit_stats.prune_rate(),
               bit_stats.skipped, bit_stats.total, bit_ms);
        total_pass += pass;
        total_fail += fail;
        total_ms += bit_ms;
        all_stats.total += bit_stats.total;
        all_stats.skipped += bit_stats.skipped;
    }

    printf("  Total accuracy=%.6f%% (%d/%d) prune=%.2f%% time=%.2fms "
           "sign_margin_min=%.3f\n",
           100.0 * static_cast<double>(total_pass) /
               static_cast<double>(total_pass + total_fail),
           total_pass, total_pass + total_fail,
           100.0 * all_stats.prune_rate(), total_ms, min_sign_margin);
    if (total_fail > 0) exit(1);
    if (min_sign_margin <= 0.05) {
        printf("  FAIL: sign margin below 5%%\n");
        exit(1);
    }
}

static void robustness_multikey_recursive_gapmsb(
    const TFHEpp::BootstrappingKeyFFT<brP>&,
    const std::vector<TruncRepeatKey<P>>&,
    const Algorithm1Config&) {
    printf("[GapMSB recursion] multi-key sampled correctness\n");
    const auto cfg = PaperRowT2NConfig();
    const int p = MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    const int chapter_k = 7;
    const std::vector<int> messages = {
        0, 1, cfg.t / 2 - 17, cfg.t / 2 - 1,
        cfg.t / 2, cfg.t / 2 + 17, cfg.t - 1,
    };
    int total_pass = 0;
    int total_fail = 0;

    for (int key_round = 0; key_round < 2; key_round++) {
        TFHEpp::SecretKey sk;
        auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
        TFHEpp::bkfftgen<brP>(*bkfft, sk);
        auto trkey1 = GenerateTruncRepeatKey<P>(sk.key.get<P>(), cfg.rounds[0].beta);
        auto trkey2 = GenerateTruncRepeatKey<P>(sk.key.get<P>(), cfg.rounds[1].beta);
        std::vector<TruncRepeatKey<P>> trkeys = {std::move(trkey1), std::move(trkey2)};

        int pass = 0;
        int fail = 0;
        for (int m : messages) {
            TFHEpp::TLWE<P> ct{};
            TFHEpp::tlweSymEncrypt<P>(
                ct, encode_message<P>(m, p), P::α, sk.key.get<P>());
            auto out = RecursiveGapMSB<brP>(
                ct, *bkfft, trkeys, cfg,
                GapMSBOptions{.p = p, .k = chapter_k},
                /*rounds=*/2);
            auto phase = TFHEpp::tlweSymPhase<P>(out, sk.key.get<P>());
            const int got = DecodeSignPhase<P>(phase);
            const int want = m >= cfg.t / 2 ? 1 : 0;
            if (got == want) {
                pass++;
            } else {
                printf("  FAIL key_round=%d m=%d got=%d want=%d\n",
                       key_round, m, got, want);
                fail++;
            }
        }
        total_pass += pass;
        total_fail += fail;
        printf("  key_round=%d accuracy=%.6f%% (%d/%d)\n",
               key_round,
               100.0 * static_cast<double>(pass) /
                   static_cast<double>(pass + fail),
               pass, pass + fail);
    }
    printf("  Total recursive accuracy=%.6f%% (%d/%d)\n",
           100.0 * static_cast<double>(total_pass) /
               static_cast<double>(total_pass + total_fail),
           total_pass, total_pass + total_fail);
    if (total_fail > 0) exit(1);
}

int main() {
    printf("=== MetaPBS2 GapMSB Test ===\n");
    test_parameter_helpers();

    auto cfg = PaperRowT2NConfig();
    TFHEpp::SecretKey sk;
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    TFHEpp::bkfftgen<brP>(*bkfft, sk);
    auto trkey1 = GenerateTruncRepeatKey<P>(sk.key.get<P>(), cfg.rounds[0].beta);
    auto trkey2 = GenerateTruncRepeatKey<P>(sk.key.get<P>(), cfg.rounds[1].beta);
    std::vector<TruncRepeatKey<P>> trkeys = {std::move(trkey1), std::move(trkey2)};

    test_direct_logical_to_arithmetic(sk, *bkfft);
    test_gapmsb_parameter_rejections(*bkfft, trkeys, cfg);
    test_weighted_extract_and_clear(sk, *bkfft, trkeys, cfg);
    stress_gapmsb_correctness(sk, *bkfft, trkeys, cfg);
    robustness_multikey_recursive_gapmsb(*bkfft, trkeys, cfg);

    printf("=== ALL GAPMSB TESTS PASSED ===\n");
    return 0;
}
