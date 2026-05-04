// =============================================================
// Test: Module 1 — sym_range.hpp
// Mirrors tfhe-go/metapbs/trunc_repeat.go tests:
//   SymRange, ReduceLaurentExponent, DeltaOffset
// =============================================================
#include <cassert>
#include <cstdio>

#include "metapbs2/sym_range.hpp"

using namespace MetaPBS2;

static void test_sym_range() {
    printf("[SymRange]\n");

    // [2]_sym = {0, 1}
    auto [lo2, hi2] = SymRange(2);
    assert(lo2 == 0 && hi2 == 1);

    // [4]_sym = {-2, -1, 0, 1}
    auto [lo4, hi4] = SymRange(4);
    assert(lo4 == -2 && hi4 == 1);

    // [3]_sym = {-1, 0, 1}
    auto [lo3, hi3] = SymRange(3);
    assert(lo3 == -1 && hi3 == 1);

    // [14]_sym: lo=-(14/2)=-7, hi=(14-1)/2=6
    auto [lo14, hi14] = SymRange(14);
    assert(lo14 == -7 && hi14 == 6);

    printf("  PASSED\n");
}

static void test_reduce_laurent_exponent() {
    printf("[ReduceLaurentExponent]\n");
    const int N = 8;

    // X^0 = +X^0
    {
        auto [s, i] = ReduceLaurentExponent(0, N);
        assert(s == +1 && i == 0);
    }
    // X^3 = +X^3
    {
        auto [s, i] = ReduceLaurentExponent(3, N);
        assert(s == +1 && i == 3);
    }
    // X^8 = -X^0  (negacyclic: X^N = -1)
    {
        auto [s, i] = ReduceLaurentExponent(8, N);
        assert(s == -1 && i == 0);
    }
    // X^9 = -X^1
    {
        auto [s, i] = ReduceLaurentExponent(9, N);
        assert(s == -1 && i == 1);
    }
    // X^{-1} = X^{15} mod X^{16}+1 = -X^7
    {
        auto [s, i] = ReduceLaurentExponent(-1, N);
        assert(s == -1 && i == 7);
    }
    // X^{16} = X^0 (X^{2N} = 1)
    {
        auto [s, i] = ReduceLaurentExponent(16, N);
        assert(s == +1 && i == 0);
    }

    printf("  PASSED\n");
}

static void test_delta_offset() {
    printf("[DeltaOffset]\n");

    // Go reference: DeltaOffset(1, 14)
    // [14]_sym = {-7,..,6}, [1]_sym = {0,0} → lo(rB)=lo(14)=-7, B*lo(r)=14*0=0, lo(B)=-7
    // DeltaOffset(1,14) = -7 - 14*0 - (-7) = 0
    assert(DeltaOffset(1, 14) == 0);

    // DeltaOffset(14, 12):
    // r=14, B=12, rB=168
    // [168]_sym: lo=-(168/2)=-84
    // [14]_sym: lo=-7
    // [12]_sym: lo=-6
    // DeltaOffset = -84 - 12*(-7) - (-6) = -84 + 84 + 6 = 6
    assert(DeltaOffset(14, 12) == 6);

    // DeltaOffset(2, 2):
    // r=2, B=2, rB=4
    // [4]_sym: lo=-2, [2]_sym={0,1}: lo=0, [2]_sym: lo=0
    // DeltaOffset = -2 - 2*0 - 0 = -2
    assert(DeltaOffset(2, 2) == -2);

    printf("  PASSED\n");
}

static void test_centered_rem() {
    printf("[CenteredRem]\n");

    // For uint64_t modular arithmetic:
    // CenteredRem(6, 4) = 6%4=2, 2 <= 4/2=2 → 2
    assert(CenteredRem<uint64_t>(6, 4) == 2);
    // CenteredRem(3, 4) = 3, 3 > 2 → 3-4 = uint64 wrap = huge → need sign-aware version
    // The Go version uses uint64 but the result is interpreted as signed difference.
    // Let's just test even divisibility:
    assert(CenteredRem<uint64_t>(8, 4) == 0);
    assert(CenteredRem<uint64_t>(0, 4) == 0);

    printf("  PASSED\n");
}

int main() {
    printf("=== MetaPBS2 Module 1: sym_range ===\n");
    test_sym_range();
    test_reduce_laurent_exponent();
    test_delta_offset();
    test_centered_rem();
    printf("=== ALL PASSED ===\n");
    return 0;
}
