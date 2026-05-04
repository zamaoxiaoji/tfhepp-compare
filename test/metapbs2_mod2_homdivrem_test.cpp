// =============================================================
// Test: Module 2 — homdivrem.hpp
// Mirrors tfhe-go/metapbs/metapbs_d1.go tests.
//
// Test HomDivRemLWE with exact multiples and sub-delta values.
// =============================================================
#include <cassert>
#include <cstdint>
#include <cstdio>

#include "metapbs2/homdivrem.hpp"
#include "params.hpp"

using namespace MetaPBS2;

// Use lvl0param for a simple LWE domain
using domP = TFHEpp::lvl0param;
using T    = typename domP::T;  // uint32_t
constexpr int N = 1;  // just test scalar logic via b-term

// For a scalar TLWE (only the b-term), we test HomDivRemLWE directly
// using a helper that creates a zero-mask TLWE.
static TFHEpp::TLWE<domP> make_tlwe(T b) {
    TFHEpp::TLWE<domP> ct = {};
    ct[domP::k * domP::n] = b;
    return ct;
}

static void test_homdivrem_exact_multiples() {
    printf("[HomDivRemLWE] exact multiples\n");
    // Use lvl2param-like scale: q' = 2^(nbit+1) = 2^12 for lvl0 nbit=10
    // But since we want a simple test, use q_prime = 4 (so scale = 2^30 for uint32)
    constexpr T q_prime = 4;
    constexpr T scale = ((~T(0)) / q_prime) + T(1);  // = 2^32 / 4 = 2^30 = 0x40000000

    // Test: b = m * scale  → quo = m, rem = 0
    for (int m = 0; m < 8; m++) {
        T b = static_cast<T>(m) * scale;
        auto ct = make_tlwe(b);
        TFHEpp::TLWE<domP> cq = {}, cr = {};
        HomDivRemLWE<domP>(cq, cr, ct, q_prime);
        assert(cq[domP::k * domP::n] == static_cast<T>(m));
        assert(cr[domP::k * domP::n] == 0);
    }
    printf("  PASSED\n");
}

static void test_homdivrem_sub_delta() {
    printf("[HomDivRemLWE] sub-delta values\n");
    // q' = 8, scale = 2^32/8 = 2^29 = 0x20000000
    constexpr T q_prime = 8;
    constexpr T scale = ((~T(0)) / q_prime) + T(1);
    constexpr T sub = scale / 4;  // sub-unit = scale/4

    // b = r * sub  (r=0..3) → quo = 0, rem = r
    for (int r = 0; r < 4; r++) {
        T b = static_cast<T>(r) * sub;
        auto ct = make_tlwe(b);
        TFHEpp::TLWE<domP> cq = {}, cr = {};
        HomDivRemLWE<domP>(cq, cr, ct, q_prime);
        assert(cq[domP::k * domP::n] == 0);
        // rem = CenteredRem(r*sub, scale) = r*sub if r*sub <= scale/2
        // scale = 2^29, scale/2 = 2^28
        // For r <= 2 (since sub = scale/4): r*sub = r*2^27 <= 2^28 → rem = r*sub
        if (r <= 2) {
            assert(cr[domP::k * domP::n] == static_cast<T>(r) * sub);
        }
    }
    printf("  PASSED\n");
}

static void test_homdivrem_raw_phase_identity() {
    printf("[HomDivRemLWE] raw-phase identity: c = scale*quo + rem\n");
    constexpr T q_prime = 4;
    constexpr T scale = ((~T(0)) / q_prime) + T(1);

    // Test several b values
    for (int trial = 0; trial < 20; trial++) {
        T b = static_cast<T>(trial * 7 + 3) * (scale / 8);
        auto ct = make_tlwe(b);
        TFHEpp::TLWE<domP> cq = {}, cr = {};
        HomDivRemLWE<domP>(cq, cr, ct, q_prime);
        // Identity: b == scale * quo + rem
        T reconstructed = scale * cq[domP::k * domP::n] + cr[domP::k * domP::n];
        assert(reconstructed == b);
    }
    printf("  PASSED\n");
}

static void test_homdivrem_at_scale() {
    printf("[HomDivRemAtScale] cascade correctness\n");
    // Simulate two-level decomposition as in Algorithm 1
    // currentMod = 2*N = let's say 8, beta = 4 → q_prime = 32
    constexpr int current_mod = 8;
    constexpr int beta = 4;
    constexpr T q_prime_eff = static_cast<T>(current_mod * beta);  // 32
    constexpr T scale = ((~T(0)) / q_prime_eff) + T(1);

    for (int m = 0; m < 32; m++) {
        T b = static_cast<T>(m) * scale;
        auto ct = make_tlwe(b);
        TFHEpp::TLWE<domP> cq1 = {}, cr1 = {};
        HomDivRemAtScale<domP>(cq1, cr1, ct, current_mod, beta);
        // Identity must hold
        T reconstructed = scale * cq1[domP::k * domP::n] + cr1[domP::k * domP::n];
        assert(reconstructed == b);
    }
    printf("  PASSED\n");
}

int main() {
    printf("=== MetaPBS2 Module 2: homdivrem ===\n");
    test_homdivrem_exact_multiples();
    test_homdivrem_sub_delta();
    test_homdivrem_raw_phase_identity();
    test_homdivrem_at_scale();
    printf("=== ALL PASSED ===\n");
    return 0;
}
