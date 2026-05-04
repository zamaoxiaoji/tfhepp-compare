#pragma once
// =============================================================
// homdivrem.hpp — HomDivRem (Definition 1) and HomDivRemAtScale
//
// Mirrors tfhe-go/metapbs/metapbs_d1.go:
//   HomDivRemLWE[T](ct, qPrime) -> (ctQuo, ctRem)
// and tfhe-go/metapbs/divrem_scale_fix.go:
//   HomDivRemAtScale[T](ct, currentMod, beta) -> (ctQuo, ctRem)
//
// Definition 1 (paper):
//   Each coordinate c[i] is split as:
//     c_rem[i] = CenteredRem(c[i], Q/q')
//     c_quo[i] = (c[i] - c_rem[i]) / (Q/q')
//   Raw-phase identity: ψ*(c) = (Q/q') · ψ*(c_quo) + ψ*(c_rem)
//   c_quo is NOT reduced mod q'.
//
// For Q = 2^64 and power-of-two q':
//   scale = Q / q' = 2^(64 - log2(q'))
//   CenteredRem(c, scale): remainder in [-scale/2, scale/2]
// =============================================================

#include <cstdint>
#include <limits>

#include "metapbs2/sym_range.hpp"
#include "params.hpp"
#include "tlwe.hpp"

namespace MetaPBS2 {

// ---------------------------------------------------------------
// HomDivRemLWE: Definition 1 applied to LWE ciphertext.
//
// domP = LWE parameter type (determines T and dimension)
// q_prime = target modulus (power of two)
//
// Mirrors: func HomDivRemLWE[T TorusInt](ct LWECiphertext[T], qPrime T)
//          -> (ctQuo, ctRem LWECiphertext[T])
// ---------------------------------------------------------------
template <class domP>
void HomDivRemLWE(TFHEpp::TLWE<domP>& c_quo, TFHEpp::TLWE<domP>& c_rem,
                  const TFHEpp::TLWE<domP>& c, typename domP::T q_prime) {
    using T = typename domP::T;
    constexpr uint32_t bits = std::numeric_limits<T>::digits;
    constexpr uint32_t len = domP::k * domP::n + 1;

    // scale = Q / q' = 2^bits / q_prime
    // For T=uint64, Q=2^64: scale = (2^64 - 1) / q_prime + 1 when q_prime | 2^64
    T scale = ((~T(0)) / q_prime) + T(1);

    for (uint32_t i = 0; i < len; i++) {
        T ci = c[i];
        T rem = CenteredRem<T>(ci, scale);
        c_rem[i] = rem;
        c_quo[i] = (ci - rem) / scale;
    }
}

// ---------------------------------------------------------------
// HomDivRemAtScale: split remainder at cumulative Meta-PBS scale.
//
// Algorithm 1 recursively applies HomDivRem to the previous remainder.
// Round k must split at Q/(currentMod * beta), where
//   currentMod = 2N * Π_{j<k} β_j
//
// Mirrors: func HomDivRemAtScale[T TorusInt](ct, currentMod, beta) -> (ctQuo, ctRem)
//          which calls HomDivRemLWE(ct, T(currentMod*beta))
// ---------------------------------------------------------------
template <class domP>
void HomDivRemAtScale(TFHEpp::TLWE<domP>& c_quo, TFHEpp::TLWE<domP>& c_rem,
                      const TFHEpp::TLWE<domP>& c,
                      int current_mod, int beta) {
    using T = typename domP::T;
    HomDivRemLWE<domP>(c_quo, c_rem, c, static_cast<T>(current_mod * beta));
}

}  // namespace MetaPBS2
