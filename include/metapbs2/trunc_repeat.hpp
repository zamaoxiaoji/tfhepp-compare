#pragma once
// =============================================================
// trunc_repeat.hpp — Plaintext and Homomorphic TruncRepeat
//
// Mirrors tfhe-go/metapbs/trunc_repeat.go:
//   PlainTruncPadByLemma5(), PlainTruncRepeatByLemma5()
//   BuildSumPoly(), SymRange(), ReduceLaurentExponent()
// and tfhe-go/metapbs/hom_trunc_pad.go:
//   TruncPadKey, GenTruncPadKey, TruncPadGLWE
// and tfhe-go/metapbs/hom_trunc_repeat.go:
//   TruncRepeatKey, GenTruncRepeatKey, TruncRepeatGLWE
//
// The homomorphic TruncRepeat uses the Column Method (paper §4.3):
//   TruncRepeat(C,[a,b],β) = Σ_{i∈[N]} <G⁻¹(truncPad(A·Xⁱ)), RLev(s[i]·sumPoly)>
//                            + truncRepeat(B, [a,b], β)
//
// Key: RLev uses HalfTRGSWFFT in TFHEpp (= FourierGLevCiphertext in tfhe-go).
// =============================================================

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "mulfft.hpp"
#include "params.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"
#include "metapbs2/sym_range.hpp"

namespace MetaPBS2 {

// =============================================================
// Plaintext Operations
// =============================================================

// ---------------------------------------------------------------
// PolyMulByXk: Negacyclic multiplication by X^k (plaintext)
//   result[i] = poly[i-k] with negacyclic sign flips
//
// Mirrors: polyEval.MonomialMulPoly(poly, k) from tfhe-go
// ---------------------------------------------------------------
template <class P>
TFHEpp::Polynomial<P> PolyMulByXk(const TFHEpp::Polynomial<P>& poly, int k) {
    constexpr int N = P::n;
    k = ((k % (2 * N)) + 2 * N) % (2 * N);
    TFHEpp::Polynomial<P> res = {};
    for (int i = 0; i < N; i++) {
        int src = i - k;
        if (src >= 0 && src < N) {
            res[i] = poly[src];
        } else if (src < 0) {
            src += 2 * N;
            if (src >= N)
                res[i] = -poly[src - N];
            else
                res[i] = poly[src];
        } else {
            // src >= N
            res[i] = -poly[src - N];
        }
    }
    return res;
}

// ---------------------------------------------------------------
// PlainTruncPadByLemma5: truncPad(M, [a,b], B) via Lemma 5.
//
// Formula (Lemma 5):
//   j >= 0: coeff = M[j],   exp = j*B
//   j <  0: coeff = M[N+j], exp = N+j*B   (NO negation!)
// Reduce each X^exp mod (X^N+1) via ReduceLaurentExponent.
//
// Mirrors: func PlainTruncPadByLemma5[T](M, a, b, B, N)
// ---------------------------------------------------------------
template <class P>
TFHEpp::Polynomial<P> PlainTruncPadByLemma5(const TFHEpp::Polynomial<P>& M,
                                              int a, int b, int B) {
    constexpr int N = P::n;
    TFHEpp::Polynomial<P> result = {};
    for (int j = a; j <= b; j++) {
        typename P::T coeff;
        int exp;
        if (j >= 0) {
            coeff = M[j];
            exp = j * B;
        } else {
            coeff = M[N + j];  // positive, per Lemma 5 (no negation)
            exp = N + j * B;
        }
        auto [sign, idx] = ReduceLaurentExponent(exp, N);
        if (sign == +1)
            result[idx] += coeff;
        else
            result[idx] -= coeff;
    }
    return result;
}

// ---------------------------------------------------------------
// BuildSumPoly: sumPoly = Σ_{k ∈ [B]_sym} X^k  mod (X^N+1)
//
// Mirrors: func BuildSumPoly[T](B, N int) poly.Poly[T]
// ---------------------------------------------------------------
template <class P>
TFHEpp::Polynomial<P> BuildSumPoly(int B) {
    constexpr int N = P::n;
    TFHEpp::Polynomial<P> result = {};
    auto [lo, hi] = SymRange(B);
    for (int k = lo; k <= hi; k++) {
        auto [sign, idx] = ReduceLaurentExponent(k, N);
        if (sign == +1)
            result[idx] += 1;
        else
            result[idx] -= 1;
    }
    return result;
}

// ---------------------------------------------------------------
// NegacyclicPolyMul: schoolbook negacyclic convolution (plaintext)
// ---------------------------------------------------------------
template <class P>
TFHEpp::Polynomial<P> NegacyclicPolyMul(const TFHEpp::Polynomial<P>& a,
                                          const TFHEpp::Polynomial<P>& b) {
    constexpr int N = P::n;
    TFHEpp::Polynomial<P> res = {};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i + j;
            if (idx < N)
                res[idx] += a[i] * b[j];
            else
                res[idx - N] -= a[i] * b[j];
        }
    }
    return res;
}

// ---------------------------------------------------------------
// PlainTruncRepeatByLemma5: truncRepeat(M,[a,b],B) via Lemma 5.
//
// truncRepeat = Σ_{j∈[a,b], j>=0} Σ_{k∈[B]_sym} M[j]   * X^{j*B+k}
//            + Σ_{j∈[a,b], j<0}  Σ_{k∈[B]_sym} M[N+j] * X^{N+j*B+k}
//
// Mirrors: func PlainTruncRepeatByLemma5[T](M, a, b, B, N)
// ---------------------------------------------------------------
template <class P>
TFHEpp::Polynomial<P> PlainTruncRepeatByLemma5(const TFHEpp::Polynomial<P>& M,
                                                 int a, int b, int B) {
    constexpr int N = P::n;
    TFHEpp::Polynomial<P> result = {};
    auto [lo, hi] = SymRange(B);

    for (int j = a; j <= b; j++) {
        typename P::T coeff;
        int baseExp;
        if (j >= 0) {
            coeff = M[j];
            baseExp = j * B;
        } else {
            coeff = M[N + j];
            baseExp = N + j * B;
        }
        for (int k = lo; k <= hi; k++) {
            auto [sign, idx] = ReduceLaurentExponent(baseExp + k, N);
            if (sign == +1)
                result[idx] += coeff;
            else
                result[idx] -= coeff;
        }
    }
    return result;
}

// =============================================================
// Homomorphic TruncPad Key (mirrors TruncPadKey + GenTruncPadKey)
//
// Keys[i] = HalfTRGSWFFT(s[i])   for i = 0..N-1
// (RLev in paper = HalfTRGSW in TFHEpp = FourierGLev in tfhe-go)
// =============================================================
template <class P>
using TruncPadKey = std::unordered_map<uint32_t, TFHEpp::HalfTRGSWFFT<P>>;

// ---------------------------------------------------------------
// GenerateTruncPadKey: encrypts HalfTRGSWFFT(s[i]) for each i.
//
// Mirrors: func GenTruncPadKey[T](enc, gadgetParams) TruncPadKey[T]
// ---------------------------------------------------------------
template <class P>
TruncPadKey<P> GenerateTruncPadKey(const TFHEpp::Key<P>& key) {
    TruncPadKey<P> trkey;
    for (uint32_t i = 0; i < P::k * P::n; i++) {
        // Encrypt s[i] as a constant polynomial (just s[i] at index 0)
        TFHEpp::Polynomial<P> si_poly = {};
        si_poly[0] = static_cast<typename P::T>(key[i]);
        TFHEpp::HalfTRGSW<P> halftrgsw;
        TFHEpp::halftrgswSymEncrypt<P>(halftrgsw, si_poly, P::α, key);
        trkey[i] = TFHEpp::ApplyFFT2halftrgsw<P>(halftrgsw);
    }
    return trkey;
}

// =============================================================
// Homomorphic TruncRepeat Key (mirrors TruncRepeatKey + GenTruncRepeatKey)
//
// Keys[i] = HalfTRGSWFFT(s[i] * sumPoly)  for i = 0..N-1
// =============================================================
template <class P>
using TruncRepeatKey = std::unordered_map<uint32_t, TFHEpp::HalfTRGSWFFT<P>>;

// ---------------------------------------------------------------
// GenerateTruncRepeatKey: encrypts HalfTRGSWFFT(s[i] * sumPoly).
//
// Mirrors: func GenTruncRepeatKey[T](enc, gadgetParams, B) TruncRepeatKey[T]
// ---------------------------------------------------------------
template <class P>
TruncRepeatKey<P> GenerateTruncRepeatKey(const TFHEpp::Key<P>& key, int B) {
    TruncRepeatKey<P> trkey;
    auto sumPoly = BuildSumPoly<P>(B);

    for (uint32_t i = 0; i < P::k * P::n; i++) {
        auto si = static_cast<typename P::T>(key[i]);
        if (si == 0) continue;  // Skip zero key coefficients (sparse optimization)

        // s[i] * sumPoly: scalar multiply sumPoly by s[i]
        TFHEpp::Polynomial<P> si_sum = {};
        for (int j = 0; j < (int)P::n; j++)
            si_sum[j] = si * sumPoly[j];

        TFHEpp::HalfTRGSW<P> halftrgsw;
        TFHEpp::halftrgswSymEncrypt<P>(halftrgsw, si_sum, P::α, key);
        trkey[i] = TFHEpp::ApplyFFT2halftrgsw<P>(halftrgsw);
    }
    return trkey;
}

// =============================================================
// Homomorphic TruncPad (column method)
//
// TruncPad(C,[a,b],B) = Σ_i <G⁻¹(truncPad(A·Xⁱ)), RLev(s[i])> + truncPad(B)
//
// Mirrors: TruncPadEvaluator.TruncPadGLWEAssign in tfhe-go
// =============================================================
template <class P>
void HomTruncPad(TFHEpp::TRLWE<P>& result,
                 const TFHEpp::TRLWE<P>& input,
                 int a, int b, int B,
                 const TruncPadKey<P>& padkey) {
    constexpr int N = P::n;

    // Initialize result to zero
    result = {};

    // Body term: truncPad(B_body, [a,b], B)  (B_body = input[P::k] = body poly)
    result[P::k] = PlainTruncPadByLemma5<P>(input[P::k], a, b, B);

    // Mask columns: for each mask index ki, for each position i
    for (uint32_t ki = 0; ki < P::k; ki++) {
        for (int i = 0; i < N; i++) {
            uint32_t key_idx = i + ki * N;
            auto it = padkey.find(key_idx);
            if (it == padkey.end()) continue;  // s[i]=0, skip

            // Rotate A[ki] by X^i
            auto rotA = PolyMulByXk<P>(input[ki], i);

            // truncPad on rotated polynomial
            auto tp = PlainTruncPadByLemma5<P>(rotA, a, b, B);

            // Check all-zero
            bool all_zero = true;
            for (int j = 0; j < N; j++) {
                if (tp[j] != 0) { all_zero = false; break; }
            }
            if (all_zero) continue;

            // External product: <G⁻¹(tp), RLev(s[i])>
            TFHEpp::TRLWE<P> contribution;
            TFHEpp::ExternalProduct<P>(contribution, tp, it->second);

            // Accumulate (subtract: TFHEpp uses B = A·s + M + e convention)
            for (int comp = 0; comp <= (int)P::k; comp++)
                for (int j = 0; j < N; j++)
                    result[comp][j] -= contribution[comp][j];
        }
    }
}

// =============================================================
// Homomorphic TruncRepeat (column method) — MAIN implementation
//
// TruncRepeat(C,[a,b],B) = Σ_i <G⁻¹(truncPad(A·Xⁱ)), RLev(s[i]·sumPoly)>
//                          + truncRepeat(B_body, [a,b], B)
//
// This is more efficient than HomTruncPad + multiply by sumPoly because
// the sumPoly is baked into the key material.
//
// Mirrors: TruncRepeatEvaluator.TruncRepeatGLWEAssign in tfhe-go
// =============================================================
template <class P>
void HomTruncRepeat(TFHEpp::TRLWE<P>& result,
                    const TFHEpp::TRLWE<P>& input,
                    int a, int b, int B,
                    const TruncRepeatKey<P>& repkey) {
    constexpr int N = P::n;

    // Initialize result to zero
    result = {};

    // Body term: truncRepeat(B_body, [a,b], B)
    result[P::k] = PlainTruncRepeatByLemma5<P>(input[P::k], a, b, B);

    // Mask columns
    for (uint32_t ki = 0; ki < P::k; ki++) {
        for (int i = 0; i < N; i++) {
            uint32_t key_idx = i + ki * N;
            auto it = repkey.find(key_idx);
            if (it == repkey.end()) continue;  // s[i]=0, skip

            // Rotate A[ki] by X^i  (monomial rotation)
            auto rotA = PolyMulByXk<P>(input[ki], i);

            // truncPad (NOT truncRepeat) on rotated polynomial
            auto tp = PlainTruncPadByLemma5<P>(rotA, a, b, B);

            // Zero-skip optimization
            bool all_zero = true;
            for (int j = 0; j < N; j++) {
                if (tp[j] != 0) { all_zero = false; break; }
            }
            if (all_zero) continue;

            // External product: <G⁻¹(tp), RLev(s[i]·sumPoly)>
            TFHEpp::TRLWE<P> contribution;
            TFHEpp::ExternalProduct<P>(contribution, tp, it->second);

            // Accumulate (subtract)
            for (int comp = 0; comp <= (int)P::k; comp++)
                for (int j = 0; j < N; j++)
                    result[comp][j] -= contribution[comp][j];
        }
    }
}

// =============================================================
// HomTruncRepeatShifted: TruncRepeat + X^delta_shift rotation
//
// C' = TruncRepeat(C, [a,b], B) · X^{delta_shift}
//
// Mirrors: TruncRepeatEvaluator.TruncRepeatGLWE + applyShift in tfhe-go
// =============================================================
template <class P>
void HomTruncRepeatShifted(TFHEpp::TRLWE<P>& result,
                            const TFHEpp::TRLWE<P>& input,
                            int a, int b, int B, int delta_shift,
                            const TruncRepeatKey<P>& repkey) {
    HomTruncRepeat<P>(result, input, a, b, B, repkey);

    // Apply X^{delta_shift} to each polynomial component
    for (int comp = 0; comp <= (int)P::k; comp++) {
        auto rotated = PolyMulByXk<P>(result[comp], delta_shift);
        result[comp] = rotated;
    }
}

}  // namespace MetaPBS2
