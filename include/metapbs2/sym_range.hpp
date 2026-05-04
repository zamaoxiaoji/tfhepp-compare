#pragma once
// =============================================================
// sym_range.hpp — Laurent ring helper functions
//
// Mirrors tfhe-go/metapbs/trunc_repeat.go:
//   SymRange(), ReduceLaurentExponent(), DeltaOffset()
// and tfhe-go/metapbs/hom_divrem.go:
//   CenteredRem()
// =============================================================

#include <cstdint>
#include <limits>
#include <utility>

namespace MetaPBS2 {

// ---------------------------------------------------------------
// SymRange: [q]_sym = Z ∩ [-floor(q/2), ceil(q/2)-1]
// Exception: [2]_sym = {0, 1} (NOT {-1, 0}).
//
// Mirrors: func SymRange(q int) (lo, hi int)
// ---------------------------------------------------------------
inline std::pair<int, int> SymRange(int q) {
    if (q == 2) return {0, 1};
    return {-(q / 2), (q - 1) / 2};
}

// ---------------------------------------------------------------
// ReduceLaurentExponent: reduce X^exp mod (X^N + 1)
// Returns (sign, idx) such that X^exp ≡ sign * X^idx
// with 0 <= idx < N, sign ∈ {+1, -1}.
//
// Mirrors: func ReduceLaurentExponent(exp, N int) (sign int, idx int)
// ---------------------------------------------------------------
inline std::pair<int, int> ReduceLaurentExponent(int exp, int N) {
    int twoN = 2 * N;
    // Normalize to [0, 2N)
    int r = ((exp % twoN) + twoN) % twoN;
    if (r < N) return {+1, r};
    return {-1, r - N};
}

// ---------------------------------------------------------------
// DeltaOffset: Δ_{r,B} = min([rB]_sym) - B*min([r]_sym) - min([B]_sym)
//
// Mirrors: func DeltaOffset(r, B int) int
// ---------------------------------------------------------------
inline int DeltaOffset(int r, int B) {
    auto [lorB, _1] = SymRange(r * B);
    auto [lor,  _2] = SymRange(r);
    auto [loB,  _3] = SymRange(B);
    (void)_1; (void)_2; (void)_3;
    return lorB - B * lor - loB;
}

// ---------------------------------------------------------------
// CenteredRem: centered remainder of x mod d.
// Result is in [-(d-1)/2, d/2].
//
// Mirrors: func CenteredRem[T TorusInt](x T, d T) T
// Works for unsigned types: interpret result unsigned.
// ---------------------------------------------------------------
template <typename T>
inline T CenteredRem(T x, T d) {
    T r = x % d;
    if (r > d / 2) r -= d;
    return r;
}

}  // namespace MetaPBS2
