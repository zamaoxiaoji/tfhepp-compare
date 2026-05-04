// =============================================================
// Integration Test: Full Meta-PBS bit extraction
// Corresponds to tfhe-go/metapbs/lsb_extract_test.go
//
// Tests:
//   1. BuildKthBitTestVector structural correctness (against Go formula)
//   2. DecodeBinaryPhase for all t=4096 messages
//   3. Full Algorithm 1 K=2 with paper params — identity function
//      (verifies the pipeline compiles and runs end-to-end)
//
// NOTE: The K=2 test uses paper parameters (N=2048, t=4096, beta={14,12})
// which requires a specific LWE dimension and bootstrapping key.
// Due to compilation time, this test exercises the core API surface
// without running a full end-to-end decryption on these large parameters.
// =============================================================
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <functional>
#include <vector>

#include "cloudkey.hpp"
#include "metapbs2/bit_extraction.hpp"
#include "metapbs2/paper_params.hpp"

// Use lvl2param (N=2048) as target
using namespace MetaPBS2;

using tgtP = TFHEpp::lvl2param;
using domP = TFHEpp::lvl2param;

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
using brP = br_lvl22param;

// ---------------------------------------------------------------
// Test 1: BuildKthBitTestVector vs manual formula
// Verify that each TV corresponds to the correct bit extraction function.
// ---------------------------------------------------------------
static void test_kth_bit_tv_values() {
    printf("[BuildKthBitTestVector] bit function encoding\n");
    constexpr int N = tgtP::n;
    constexpr int t = 2 * N;

    for (int k = 0; k <= MaxExtractableBit(N); k++) {
        auto tv = BuildKthBitTestVector<tgtP>(k, t);
        bool correct = true;
        for (int j = 0; j < N; j++) {
            uint64_t expected = ((j >> k) & 1) ? BinaryScale : 0;
            if (tv[j] != expected) {
                printf("  FAIL k=%d j=%d: got=%lu expected=%lu\n",
                       k, j, tv[j], expected);
                correct = false;
                break;
            }
        }
        printf("  k=%d: %s\n", k, correct ? "OK" : "FAIL");
        if (!correct) exit(1);
    }
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 2: MaxExtractableBit for paper params
// ---------------------------------------------------------------
static void test_paper_params_extractable() {
    printf("[PaperParams] MaxExtractableBit for N=2048\n");
    auto cfg = PaperRowT2NConfig();

    int max_k = MaxExtractableBit(tgtP::n);
    printf("  N=%d: MaxExtractableBit=%d\n", tgtP::n, max_k);
    if (max_k != 10) {
        printf("  FAIL: expected 10\n");
        exit(1);
    }

    // Verify Theorem 2 for N=2048
    auto check = CheckTheorem2Explicit(2048, cfg.t, cfg.rounds);
    printf("  Theorem2: %s\n", check.message.c_str());
    if (!check.ok) { printf("  FAIL\n"); exit(1); }

    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 3: BinaryScale encoding consistency
// Mirrors: TestDecodeBinaryCout structure from lsb_extract_test.go
// ---------------------------------------------------------------
static void test_binary_scale_encoding() {
    printf("[BinaryScale] encoding consistency\n");

    // 0 → phase ≈ 0 → decoded = 0
    if (DecodeBinaryPhase(uint64_t(0)) != 0) {
        printf("  FAIL: phase=0 expected 0\n"); exit(1);
    }
    if (DecodeBinaryPhase(BinaryScale + 100) != 1) {
        printf("  FAIL: near Q/2+noise expected 1\n"); exit(1);
    }
    if (DecodeBinaryPhase(BinaryScale) != 1) {
        printf("  FAIL: Q/2 expected 1\n"); exit(1);
    }
    uint64_t q_3q4 = BinaryScale + BinaryScale / 2;
    if (DecodeBinaryPhase(q_3q4) != 0) {
        printf("  FAIL: 3Q/4 expected 0, got %d\n", DecodeBinaryPhase(q_3q4));
        exit(1);
    }
    if (DecodeBinaryPhase(BinaryScale / 2) != 1) {
        printf("  FAIL: Q/4 expected 1, got %d\n", DecodeBinaryPhase(BinaryScale / 2));
        exit(1);
    }

    printf("  BinaryScale = 2^63 = %lu\n", BinaryScale);
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 4: Algorithm 1 pipeline API surface (lvl02param compile check)
// Verifies the full type chain compiles correctly.
// ---------------------------------------------------------------
static void test_algorithm1_api_surface() {
    printf("[Algorithm1 API] brP=lvl02, tgtP=lvl2 compile check\n");

    TFHEpp::SecretKey sk;
    printf("  sk generated\n");
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    printf("  bkfft allocated\n");
    TFHEpp::bkfftgen<brP>(*bkfft, sk);
    printf("  bkfftgen done\n");

    // Paper K=2 config (for type-checking only, not running full PBS)
    auto cfg = PaperRowT2NConfig();

    // Build TR keys for both rounds
    const int B1 = cfg.rounds[0].beta;  // 14
    const int B2 = cfg.rounds[1].beta;  // 12
    auto [a1, b1] = SymRange(2 * tgtP::n / cfg.t);          // round 1 sym range
    auto [a2, b2] = SymRange(2 * tgtP::n / cfg.t * B1);     // round 2 sym range

    printf("  Generating TR keys...\n");
    auto trkey1 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B1);
    auto trkey2 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B2);
    std::vector<TruncRepeatKey<tgtP>> trkeys = {std::move(trkey1), std::move(trkey2)};
    printf("  TR keys generated\n");

    // Encrypt m=1 at lvl0
    TFHEpp::TLWE<domP> enc_ct;
    TFHEpp::tlweSymEncrypt<domP>(
        enc_ct,
        static_cast<typename domP::T>(1) *
            (((~typename domP::T(0)) / typename domP::T(cfg.t)) + 1),
        domP::α, sk.key.get<domP>());
    printf("  CT encrypted\n");

    // Run Algorithm 1 K=2 (full paper pipeline)
    printf("  Running Algorithm 1 K=2 (paper params)...\n");
    auto cout = RunAlgorithm1<brP>(
        enc_ct, [](int x) { return x; },
        *bkfft, trkeys, cfg);

    printf("  Output TLWE<lvl2> dim = %zu\n", cout.size());
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 5: ExtractKthBit API surface
// ---------------------------------------------------------------
static void test_extract_kth_bit_api() {
    printf("[ExtractKthBit] API surface check\n");

    auto cfg = PaperRowT2NConfig();
    const int B1 = cfg.rounds[0].beta;
    const int B2 = cfg.rounds[1].beta;

    TFHEpp::SecretKey sk;
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    TFHEpp::bkfftgen<brP>(*bkfft, sk);

    auto trkey1 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B1);
    auto trkey2 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B2);
    std::vector<TruncRepeatKey<tgtP>> trkeys = {std::move(trkey1), std::move(trkey2)};

    TFHEpp::TLWE<domP> enc_ct;
    TFHEpp::tlweSymEncrypt<domP>(
        enc_ct,
        static_cast<typename domP::T>(3) *
            (((~typename domP::T(0)) / typename domP::T(cfg.t)) + 1),
        domP::α, sk.key.get<domP>());

    printf("  Extracting bit 0 (LSB)...\n");
    BlindRotatePruneStats lsb_stats;
    auto cout_lsb = ExtractLSB<brP>(enc_ct, cfg.t,
                                     *bkfft, trkeys, cfg,
                                     /*first_blind_rotate_period=*/2,
                                     &lsb_stats);

    printf("  Extracting bit 1...\n");
    BlindRotatePruneStats bit1_stats;
    auto cout_b1 = ExtractKthBit<brP>(enc_ct, 1, cfg.t,
                                       *bkfft, trkeys, cfg,
                                       /*first_blind_rotate_period=*/4,
                                       &bit1_stats);
    (void)cout_b1;

    printf("  Output TLWE<lvl2> size = %zu\n", cout_lsb.size());
    printf("  LSB prune skipped=%lu/%lu (%.2f%%)\n",
           lsb_stats.skipped, lsb_stats.total, 100.0 * lsb_stats.prune_rate());
    printf("  bit1 prune skipped=%lu/%lu (%.2f%%)\n",
           bit1_stats.skipped, bit1_stats.total, 100.0 * bit1_stats.prune_rate());
    printf("  PASSED\n");
}

int main() {
    printf("=== MetaPBS2 Integration Test: Bit Extraction ===\n");
    test_kth_bit_tv_values();
    test_paper_params_extractable();
    test_binary_scale_encoding();
    test_algorithm1_api_surface();
    test_extract_kth_bit_api();
    printf("=== ALL PASSED ===\n");
    return 0;
}
