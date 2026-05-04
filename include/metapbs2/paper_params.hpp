#pragma once
// =============================================================
// paper_params.hpp — Paper parameter configurations and Theorem 2 checks
//
// Mirrors tfhe-go/metapbs/paper_row_config.go:
//   PaperRowT2NConfig() -> (ParametersLiteral, Algorithm1ExplicitConfig)
// and tfhe-go/metapbs/algorithm1.go / algorithm1_explicit.go:
//   CheckTheorem2Explicit(N, t, rounds)
//   SymRange, DeltaOffset (already in sym_range.hpp)
//
// Paper (LNCS 2025) Table 3 row: N=2048, t=4096, K=2
//   Round 1: beta=14, T=66,  DeltaBound=66
//   Round 2: beta=12, T=73,  DeltaBound=66
//   TR gadget: base=2^11, level=3
// =============================================================

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "metapbs2/sym_range.hpp"
#include "metapbs2/metapbs_pipeline.hpp"

namespace MetaPBS2 {

// =============================================================
// Theorem 2 condition check
//
// For Algorithm 1 with K rounds, the correctness condition requires:
//   For each round k: the [r_k]_sym range does not overlap the "error zone"
//   after the quotient extraction.
//
// Simplified check (as in Go CheckTheorem2Explicit):
//   For round k:
//     redundancy r_k = r0 * Π_{j<k} β_j
//     sym range = [lo, hi] with |lo|, |hi| ≤ r_k/2
//     DeltaBound must be ≥ DeltaOffset(r_k, β_k)
//
// Returns true if all conditions are satisfied.
// =============================================================
struct Theorem2Result {
    bool ok;
    std::string message;
};

inline Theorem2Result CheckTheorem2Explicit(
    int N, int t, const std::vector<RoundConfig>& rounds) {
    if (t <= 0 || (2 * N) % t != 0)
        return {false, "t must be positive and divide 2N"};

    int r = 2 * N / t;  // r0 = initial redundancy
    for (int k = 0; k < (int)rounds.size(); k++) {
        const auto& rnd = rounds[k];
        if (rnd.beta <= 1)
            return {false, "beta must be > 1 for round " + std::to_string(k)};

        int delta = DeltaOffset(r, rnd.beta);

        // Delta bound check: |delta| must not exceed delta_bound
        if (rnd.delta_bound > 0 && std::abs(delta) > rnd.delta_bound) {
            return {false, "Round " + std::to_string(k) +
                               ": |delta|=" + std::to_string(std::abs(delta)) +
                               " exceeds delta_bound=" + std::to_string(rnd.delta_bound)};
        }

        r *= rnd.beta;
    }
    return {true, "Theorem 2 satisfied"};
}

// =============================================================
// Paper Table 3 row: N=2048, t=4096, K=2
//
// TruncRepeat parameters (for key generation):
//   base = 2^11, level = 3
//   (These are baked into the TruncRepeatKey generation call)
//
// Mirrors: func PaperRowT2NConfig() (ParametersLiteral, Algorithm1ExplicitConfig)
// =============================================================
inline Algorithm1Config PaperRowT2NConfig() {
    return Algorithm1Config{
        .K = 2,
        .t = 4096,
        .rounds = {
            // Round 1: beta=14, T=66, delta_bound=66
            RoundConfig{.beta = 14, .T = 66, .delta_bound = 66},
            // Round 2: beta=12, T=73, delta_bound=66
            RoundConfig{.beta = 12, .T = 73, .delta_bound = 66},
        },
    };
}

}  // namespace MetaPBS2
