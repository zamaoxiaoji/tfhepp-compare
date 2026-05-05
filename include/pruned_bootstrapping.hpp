#pragma once

#include <cmath>
#include <limits>

#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"
#include "gatebootstrapping.hpp"

namespace TFHEpp {

// ============================================================
//  Periodic Test Vector construction for bit extraction
// ============================================================
//
// For extracting bit k (MSB = bit 0) from a p-bit message encoded
// at scale Δ = Q / 2^p, we need a LUT with period M_k = 2^(p-k)
// in the extended 2N representation.
//
// Using negacyclic-safe {0, Q/2} boolean encoding:
//   bit_k = 0  →  0
//   bit_k = 1  →  Q/2
// Since -(Q/2) ≡ Q/2 (mod Q), the encoding is self-negating
// and naturally compatible with the negacyclic ring X^N + 1.
//
// The polynomial stored has N coefficients; the negacyclic
// extension to 2N slots is implicit:  slot[j+N] = -slot[j].
//
// Each period M_k consists of w_k = M_k/2 zeros followed by
// w_k copies of Q/2.  Since -(Q/2) = Q/2, the negacyclic
// wraparound preserves this pattern, giving a true period
// of M_k over the full 2N extended representation.
//
// period_M must divide 2N.

template <class P>
Polynomial<P> PeriodicTestVector(uint32_t period_M)
{
    Polynomial<P> poly;
    const uint32_t w_k = period_M / 2;
    // Q/2 in the type's modular arithmetic
    const typename P::T half_Q =
        static_cast<typename P::T>(1)
            << (std::numeric_limits<typename P::T>::digits - 1);

    for (uint32_t i = 0; i < P::n; i++) {
        // Position within one period
        uint32_t pos_in_period = i % period_M;
        // First half of the period → 0 ; second half → Q/2
        poly[i] = (pos_in_period < w_k) ? 0 : half_Q;
    }
    return poly;
}

// ============================================================
//  Pruned Blind Rotation
// ============================================================
//
// Identical to the standard BlindRotate, except that the CMUX
// for index i is skipped when the mod-switched mask ā_i satisfies
//   ā_i ≡ 0  (mod period_M)
// because a rotation by a full period of the test vector leaves
// the accumulator unchanged.
//
// When period_M == 0, no periodic pruning is performed (only the
// standard ā == 0 skip is used, same as vanilla BlindRotate).

template <class P, uint32_t num_out = 1>
void PrunedBlindRotate(TRLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const BootstrappingKeyFFT<P> &bkfft,
                       const Polynomial<typename P::targetP> &testvector,
                       uint32_t period_M = 0)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                            << bitwidth);
    res = {};
    PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], testvector, b̄);

    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;

        // Standard skip: ā == 0
        if (ā == 0) continue;

        // Periodic pruning: ā is a multiple of the LUT period
        if (period_M > 0 && (ā % period_M) == 0) continue;

        CMUXFFTwithPolynomialMulByXaiMinusOne<P>(res, bkfft[i], ā);
    }
}

// ============================================================
//  Pruned Gate Bootstrapping  (IKS → PrunedBR → SampleExtract)
// ============================================================

template <class P>
void PrunedGateBootstrappingTLWE2TLWEFFT(
    TLWE<typename P::targetP> &res,
    const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector,
    uint32_t period_M = 0)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    PrunedBlindRotate<P>(acc, tlwe, bkfft, testvector, period_M);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

}  // namespace TFHEpp
