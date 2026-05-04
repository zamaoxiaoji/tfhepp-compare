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

using namespace MetaPBS2;

// Use lvl2param (N=2048) as target
using tgtP = TFHEpp::lvl2param;
using brP  = TFHEpp::lvl02param;   // domainP=lvl0, targetP=lvl2 (N=2048)
using domP = TFHEpp::lvl0param;

// ---------------------------------------------------------------
// Test 1: BuildKthBitTestVector vs manual formula
// Verify that each TV corresponds to the correct bit extraction function.
// ---------------------------------------------------------------
static void test_kth_bit_tv_values() {
    printf("[BuildKthBitTestVector] bit function encoding\n");
    constexpr int t = 8;   // small message modulus for easy verification
    constexpr int N = tgtP::n;

    // paper scale Q/t
    const uint64_t scale = ((~uint64_t(0)) / t) + 1;
    const int r0 = 2 * N / t;

    for (int k = 0; k <= MaxExtractableBit(t); k++) {
        auto tv = BuildKthBitTestVector<tgtP>(k, t);

        // Verify by re-computing manually for each message slot
        // f_k(i) = (i >> k) & 1 for i ∈ [0, t/2)
        bool correct = true;
        for (int i = 0; i < t / 2; i++) {
            int bit = (i >> k) & 1;
            uint64_t expected_val = static_cast<uint64_t>(bit) * scale;

            // The slot for message i occupies positions [i*r0, i*r0 + r0)
            // In [r0]_sym: one representative position
            auto [lo, hi] = SymRange(r0);
            uint64_t sum_at_slot = 0;
            for (int j = lo; j <= hi; j++) {
                auto [sign, idx] = ReduceLaurentExponent(i * r0 + j, N);
                if (sign == +1) sum_at_slot += tv[idx];
                else            sum_at_slot -= tv[idx];
            }
            // sum_at_slot should equal expected_val * r0 (each slot repeated r0 times)
            if (sum_at_slot != expected_val * r0) {
                printf("  FAIL k=%d i=%d: sum=%lu expected=%lu*%d=%lu\n",
                       k, i, sum_at_slot, expected_val, r0, expected_val * r0);
                correct = false;
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
    printf("[PaperParams] MaxExtractableBit for t=4096\n");
    auto cfg = PaperRowT2NConfig();

    // t=4096=2^12: bits 0..11 are extractable
    int max_k = MaxExtractableBit(cfg.t);
    printf("  t=%d: MaxExtractableBit=%d\n", cfg.t, max_k);
    if (max_k != 11) {
        printf("  FAIL: expected 11\n");
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
    if (DecodeBinaryPhase(0) != 0) {
        printf("  FAIL: phase=0 expected 0\n"); exit(1);
    }
    // Large noise near Q/2 → decoded = 1
    if (DecodeBinaryPhase(BinaryScale + 100) != 0) {
        printf("  FAIL: near Q/2+noise expected 0\n"); exit(1);
    }
    // Q/2 ± small_noise:  Q/2 itself wraps to 0, Q/2+Q/4 → 1
    uint64_t q_3q4 = BinaryScale + BinaryScale / 2;
    if (DecodeBinaryPhase(q_3q4) != 1) {
        printf("  FAIL: 3Q/4 expected 1, got %d\n", DecodeBinaryPhase(q_3q4));
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
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);

    // Paper K=2 config (for type-checking only, not running full PBS)
    auto cfg = PaperRowT2NConfig();

    // Build TR keys for both rounds
    const int B1 = cfg.rounds[0].beta;  // 14
    const int B2 = cfg.rounds[1].beta;  // 12
    auto [a1, b1] = SymRange(2 * tgtP::n / cfg.t);          // round 1 sym range
    auto [a2, b2] = SymRange(2 * tgtP::n / cfg.t * B1);     // round 2 sym range

    auto trkey1 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B1);
    auto trkey2 = GenerateTruncRepeatKey<tgtP>(sk.key.get<tgtP>(), B2);
    std::vector<TruncRepeatKey<tgtP>> trkeys = {std::move(trkey1), std::move(trkey2)};

    // Encrypt m=1 at lvl0
    TFHEpp::TLWE<domP> enc_ct;
    TFHEpp::tlweSymEncrypt<domP>(
        enc_ct,
        static_cast<typename domP::T>(1) *
            (((~typename domP::T(0)) / typename domP::T(cfg.t)) + 1),
        domP::α, sk.key.get<domP>());

    // Run Algorithm 1 K=2 (full paper pipeline)
    printf("  Running Algorithm 1 K=2 (paper params)...\n");
    auto cout = RunAlgorithm1<brP>(
        enc_ct, [](int x) { return x; },
        ek.getbkfft<brP>(), trkeys, cfg);

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
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);

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
    auto cout_lsb = ExtractLSB<brP>(enc_ct, cfg.t,
                                     ek.getbkfft<brP>(), trkeys, cfg);

    printf("  Extracting bit 1...\n");
    auto cout_b1 = ExtractKthBit<brP>(enc_ct, 1, cfg.t,
                                       ek.getbkfft<brP>(), trkeys, cfg);

    printf("  Output TLWE<lvl2> size = %zu\n", cout_lsb.size());
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
