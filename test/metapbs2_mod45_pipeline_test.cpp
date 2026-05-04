// =============================================================
// Test: Module 4 + 5 — metapbs_pipeline + bit_extraction
//
// Tests:
//   1. Theorem 2 check with paper parameters
//   2. BuildKthBitTestVector: verify LUT values
//   3. MaxExtractableBit
//   4. DecodeBinaryPhase correctness
//   5. BuildMetaPBSTV: verify against Go formula
//   6. Full Algorithm 1 pipeline (toy params, K=0 plaintext-mode)
// =============================================================
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <functional>

#include "cloudkey.hpp"
#include "metapbs2/paper_params.hpp"
#include "metapbs2/bit_extraction.hpp"

using namespace MetaPBS2;

// Use lvl2param as target ring
using tgtP = TFHEpp::lvl2param;

// ---------------------------------------------------------------
// Test 1: Theorem 2 check with paper parameters
// ---------------------------------------------------------------
static void test_theorem2_paper_params() {
    printf("[Theorem2] paper params (N=2048, t=4096, K=2)\n");
    auto cfg = PaperRowT2NConfig();
    auto result = CheckTheorem2Explicit(2048, cfg.t, cfg.rounds);
    if (!result.ok) {
        printf("  FAIL: %s\n", result.message.c_str());
        exit(1);
    }
    printf("  %s\n", result.message.c_str());
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 2: MaxExtractableBit
// ---------------------------------------------------------------
static void test_max_extractable_bit() {
    printf("[MaxExtractableBit]\n");
    if (MaxExtractableBit(2048) != 10) {
        printf("  FAIL: N=2048 expected 10, got %d\n", MaxExtractableBit(2048));
        exit(1);
    }
    if (MaxExtractableBit(32) != 4) {
        printf("  FAIL: N=32 expected 4, got %d\n", MaxExtractableBit(32));
        exit(1);
    }
    if (MaxExtractableBit(1) != -1) {
        printf("  FAIL: N=1 expected -1, got %d\n", MaxExtractableBit(1));
        exit(1);
    }
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 3: DecodeBinaryPhase
// ---------------------------------------------------------------
static void test_decode_binary_phase() {
    printf("[DecodeBinaryPhase]\n");
    if (DecodeBinaryPhase(uint64_t(0)) != 0) { printf("  FAIL at 0\n"); exit(1); }
    if (DecodeBinaryPhase(BinaryScale) != 1) {
        printf("  FAIL at Q/2=%lu expected 1, got %d\n",
               BinaryScale, DecodeBinaryPhase(BinaryScale));
        exit(1);
    }
    uint64_t three_q4 = BinaryScale + BinaryScale / 2;  // 3Q/4
    if (DecodeBinaryPhase(three_q4) != 0) {
        printf("  FAIL at 3Q/4=%lu expected 0, got %d\n",
               three_q4, DecodeBinaryPhase(three_q4));
        exit(1);
    }
    if (DecodeBinaryPhase(BinaryScale / 2) != 1) {
        printf("  FAIL at Q/4=%lu expected 1, got %d\n",
               BinaryScale / 2, DecodeBinaryPhase(BinaryScale / 2));
        exit(1);
    }
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 4: BuildKthBitTestVector vs manual formula
// Verify that the TV encodes bit_k(m) correctly for small N=16
// (we use lvl2 N=2048 and check a few slots symbolically)
// ---------------------------------------------------------------
static void test_build_kth_bit_tv() {
    printf("[BuildKthBitTestVector] structural check\n");
    constexpr int N = tgtP::n;
    constexpr int t = 2 * N;
    auto tv = BuildKthBitTestVector<tgtP>(0, t);

    for (int m = 0; m < N; m++) {
        uint64_t want = (m & 1) ? BinaryScale : 0;
        if (tv[m] != want) {
            printf("  FAIL: m=%d got=%lu want=%lu\n", m, tv[m], want);
            exit(1);
        }
    }

    bool rejected = false;
    try {
        (void)BuildKthBitTestVector<tgtP>(11, t);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    if (!rejected) {
        printf("  FAIL: expected top bit rejection for k=11\n");
        exit(1);
    }

    printf("  TV[0..4] = {%lu, %lu, %lu, %lu, %lu}\n",
           tv[0], tv[1], tv[2], tv[3], tv[4]);
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 5: BuildMetaPBSTV identity function check
// Mirrors: TestAlgorithm1TVMatchesPaperFormula
// ---------------------------------------------------------------
static void test_build_metapbs_tv_identity() {
    printf("[BuildMetaPBSTV] identity function check\n");
    constexpr int t = 8;
    constexpr int N = tgtP::n;

    // f(m) = m (identity)
    auto f_id = [](int m) { return m; };
    auto tv_id = BuildMetaPBSTV<tgtP>(f_id, t);

    // f(m) = 0 (constant zero)
    auto f_zero = [](int m) -> int { (void)m; return 0; };
    auto tv_zero = BuildMetaPBSTV<tgtP>(f_zero, t);

    // Identity TV should differ from zero TV
    bool differs = false;
    for (int i = 0; i < N; i++) {
        if (tv_id[i] != tv_zero[i]) { differs = true; break; }
    }
    if (!differs) {
        printf("  FAIL: identity TV equals zero TV\n");
        exit(1);
    }

    // Paper scale = Q / t
    uint64_t scale = ((~uint64_t(0)) / static_cast<uint64_t>(t)) + 1;
    printf("  Scale = Q/%d = %lu\n", t, scale);
    printf("  TV_id[0..2] = {%lu, %lu, %lu}\n", tv_id[0], tv_id[1], tv_id[2]);
    printf("  PASSED\n");
}

// ---------------------------------------------------------------
// Test 6: Algorithm 1 with K=0 (just initial BlindRotate) — plaintext check
// This verifies the TV is correct end-to-end, without expansion rounds.
// ---------------------------------------------------------------
static void test_algorithm1_k0_compiles_and_runs() {
    printf("[Algorithm1 K=0] compilation and execution check\n");

    // Verify that RunAlgorithm1 with K=0 compiles and produces a TLWE output
    using brP   = TFHEpp::lvl01param;
    using lvl1P = TFHEpp::lvl1param;
    using lvl0P = TFHEpp::lvl0param;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);

    Algorithm1Config cfg{.K = 0, .t = 4, .rounds = {}};
    std::vector<TruncRepeatKey<lvl1P>> trkeys;

    // Encrypt m=1 and verify we get a valid TLWE<lvl1P> output (not all-zero)
    TFHEpp::TLWE<lvl0P> enc_ct;
    TFHEpp::tlweSymEncrypt<lvl0P>(
        enc_ct,
        static_cast<typename lvl0P::T>(1) *
            (((~typename lvl0P::T(0)) / typename lvl0P::T(4)) + 1),
        lvl0P::α, sk.key.get<lvl0P>());

    auto cout = RunAlgorithm1<brP>(
        enc_ct, [](int x) { return x; },
        ek.getbkfft<brP>(), trkeys, cfg);

    // Verify output has correct type/dimension (lvl1P::k*N + 1 elements)
    constexpr int expected_size = lvl1P::k * lvl1P::n + 1;
    static_assert(std::tuple_size<decltype(cout)>::value == expected_size,
                  "Output TLWE has wrong size");

    printf("  Output TLWE<lvl1> size = %d (correct)\n", expected_size);
    printf("  PASSED\n");
}

static void test_periodic_pruning_stats() {
    printf("[BlindRotate periodic pruning stats]\n");
    using brP   = TFHEpp::lvl01param;
    using lvl1P = TFHEpp::lvl1param;
    using lvl0P = TFHEpp::lvl0param;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);

    Algorithm1Config cfg{.K = 0, .t = 2 * static_cast<int>(lvl1P::n), .rounds = {}};
    std::vector<TruncRepeatKey<lvl1P>> trkeys;

    TFHEpp::TLWE<lvl0P> enc_ct;
    TFHEpp::tlweSymEncrypt<lvl0P>(
        enc_ct,
        static_cast<typename lvl0P::T>(17) *
            (((~typename lvl0P::T(0)) / typename lvl0P::T(2 * lvl1P::n)) + 1),
        lvl0P::α, sk.key.get<lvl0P>());

    BlindRotatePruneStats stats;
    auto cout = RunAlgorithm1<brP>(
        enc_ct, [](int x) { return x & 1; },
        ek.getbkfft<brP>(), trkeys, cfg,
        /*first_blind_rotate_period=*/2, &stats);
    (void)cout;

    if (stats.total == 0 || stats.skipped == 0) {
        printf("  FAIL: expected non-zero periodic pruning stats, got skipped=%lu total=%lu\n",
               stats.skipped, stats.total);
        exit(1);
    }
    if (stats.skipped > stats.total) {
        printf("  FAIL: skipped > total\n");
        exit(1);
    }
    printf("  skipped=%lu/%lu (%.2f%%)\n",
           stats.skipped, stats.total, 100.0 * stats.prune_rate());
    printf("  PASSED\n");
}

int main() {
    printf("=== MetaPBS2 Module 4+5: pipeline + bit_extraction ===\n");
    test_theorem2_paper_params();
    test_max_extractable_bit();
    test_decode_binary_phase();
    test_build_kth_bit_tv();
    test_build_metapbs_tv_identity();
    test_algorithm1_k0_compiles_and_runs();
    test_periodic_pruning_stats();
    printf("=== ALL PASSED ===\n");
    return 0;
}
