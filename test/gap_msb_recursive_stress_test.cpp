// =============================================================
// gap_msb_recursive_stress_test.cpp — Optimized comparison test
//
// Optimizations:
//   1. LOG_to_ARI: IKS lvl2→lvl0 + GateBS lvl0→lvl2 (636 CMUXes)
//   2. Base case: IKS lvl2→lvl0 + GateBS lvl0→lvl1 (N=1024)
// =============================================================
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "keyswitch.hpp"
#include "metapbs2/hom_compare.hpp"
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

using brP_meta   = br_lvl22param;
using brP_logari = TFHEpp::lvl02param;  // lvl0→lvl2 (636 CMUXes)
using brP_base   = TFHEpp::lvl01param;  // lvl0→lvl1 (N=1024)
using iksP_t     = TFHEpp::lvl20param;  // lvl2→lvl0
using P_in       = TFHEpp::lvl2param;
using P_out      = TFHEpp::lvl1param;

static constexpr int KAPPA = 5;
static constexpr int TESTS_PER_BIT_WIDTH = 16;

template <class P>
static typename P::T encode_message(uint64_t m, int p) {
    return static_cast<typename P::T>(m) << (std::numeric_limits<typename P::T>::digits - p);
}

template <class P>
static int decode_sign(typename P::T phase) {
    return static_cast<std::make_signed_t<typename P::T>>(phase) < 0 ? 1 : 0;
}

static double elapsed_ms(std::chrono::steady_clock::time_point a,
                         std::chrono::steady_clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

static void stress_test(
    int plain_bits, int num_test,
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<brP_meta>& bk_meta,
    const std::vector<TruncRepeatKey<P_in>>& trkeys,
    const Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bk_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bk_base,
    const TFHEpp::KeySwitchingKey<iksP_t>& iksk) {

    int p_eff = plain_bits + 1;
    int rounds = 0, p_tmp = p_eff;
    while (p_tmp > KAPPA) { rounds++; p_tmp -= (KAPPA - 1); }

    printf("\n========================================\n");
    printf("  %d-bit comparison (p_eff=%d, %d rounds)\n", plain_bits, p_eff, rounds);
    printf("  LOG_to_ARI: lvl0→lvl2 (636 CMUXes)\n");
    printf("  Base case: lvl0→lvl1 (N=1024)\n");
    printf("========================================\n");

    std::random_device rd;
    std::default_random_engine eng(rd());
    uint64_t max_val = (uint64_t(1) << (plain_bits - 1)) - 1;
    std::uniform_int_distribution<uint64_t> dist(0, max_val);

    int pass = 0, fail = 0;
    double total_ms = 0;
    BlindRotatePruneStats all_stats{};
    int encode_p = plain_bits + 1;
    HomMSBOptions options;
    options.kappa = KAPPA;

    printf("  Running %d tests...\n", num_test);
    fflush(stdout);

    for (int t = 0; t < num_test; t++) {
        uint64_t a = dist(eng), b = dist(eng);
        int want = (a > b) ? 1 : 0;

        TFHEpp::TLWE<P_in> c0{}, c1{};
        TFHEpp::tlweSymEncrypt<P_in>(c0, encode_message<P_in>(a, encode_p), P_in::α, sk.key.get<P_in>());
        TFHEpp::tlweSymEncrypt<P_in>(c1, encode_message<P_in>(b, encode_p), P_in::α, sk.key.get<P_in>());

        BlindRotatePruneStats st{};
        auto t0 = std::chrono::steady_clock::now();
        auto res = HomGreaterThan<brP_meta, brP_logari, brP_base, iksP_t>(
            c0, c1, plain_bits, bk_meta, trkeys, cfg,
            bk_logari, bk_base, iksk, options, &st);
        auto t1 = std::chrono::steady_clock::now();

        double ms = elapsed_ms(t0, t1);
        total_ms += ms;
        all_stats.total += st.total;
        all_stats.skipped += st.skipped;

        auto phase = TFHEpp::tlweSymPhase<P_out>(res, sk.key.get<P_out>());
        int got = decode_sign<P_out>(phase);

        if (got == want) pass++;
        else {
            fail++;
            if (fail <= 5)
                printf("    FAIL #%d: a=%lu b=%lu want=%d got=%d (%.1fms)\n", t, a, b, want, got, ms);
        }
        if ((t+1) % 20 == 0 || t == num_test-1) {
            printf("    [%d/%d] pass=%d fail=%d (%.1fms avg)\r", t+1, num_test, pass, fail, total_ms/(t+1));
            fflush(stdout);
        }
    }

    printf("\n\n  Results for %d-bit:\n", plain_bits);
    printf("    Accuracy: %.4f%% (%d/%d)\n", 100.0*pass/(pass+fail), pass, pass+fail);
    printf("    Avg time: %.2f ms\n", total_ms/num_test);
    printf("    Total: %.1f s\n", total_ms/1000.0);
    if (all_stats.total > 0)
        printf("    Prune: %.2f%% (%lu/%lu)\n", 100.0*all_stats.prune_rate(), all_stats.skipped, all_stats.total);
    printf("    %s\n", fail == 0 ? "PASSED" : "FAILED");
}

int main() {
    printf("=== HomGapMSB Stress Test (v2: lightweight LOG_to_ARI + base) ===\n\n");

    auto cfg = PaperRowT2NConfig();
    TFHEpp::SecretKey sk;

    printf("Generating keys...\n");

    printf("  [1/5] BK lvl2→lvl2 (MetaPBS)...\n"); fflush(stdout);
    auto bk_meta = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_meta>>();
    TFHEpp::bkfftgen<brP_meta>(*bk_meta, sk);

    printf("  [2/5] BK lvl0→lvl2 (LOG_to_ARI)...\n"); fflush(stdout);
    auto bk_logari = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_logari>>();
    TFHEpp::bkfftgen<brP_logari>(*bk_logari, sk);

    printf("  [3/5] BK lvl0→lvl1 (base case)...\n"); fflush(stdout);
    auto bk_base = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_base>>();
    TFHEpp::bkfftgen<brP_base>(*bk_base, sk);

    printf("  [4/5] IKS lvl2→lvl0...\n"); fflush(stdout);
    auto iksk = std::make_unique<TFHEpp::KeySwitchingKey<iksP_t>>();
    TFHEpp::ikskgen<iksP_t>(*iksk, sk);

    printf("  [5/5] TruncRepeat keys...\n"); fflush(stdout);
    auto tr1 = GenerateTruncRepeatKey<P_in>(sk.key.get<P_in>(), cfg.rounds[0].beta);
    auto tr2 = GenerateTruncRepeatKey<P_in>(sk.key.get<P_in>(), cfg.rounds[1].beta);
    std::vector<TruncRepeatKey<P_in>> trkeys = {std::move(tr1), std::move(tr2)};

    printf("Keys generated.\n");

    for (int bits : {4, 8, 16, 24, 32})
        stress_test(bits, TESTS_PER_BIT_WIDTH, sk, *bk_meta, trkeys, cfg, *bk_logari, *bk_base, *iksk);

    printf("\n=== ALL TESTS COMPLETE ===\n");
    return 0;
}
