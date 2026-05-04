#pragma once
// =============================================================
// bit_extraction.hpp — LSB/K-th bit extraction via Meta-PBS
//
// Mirrors tfhe-go/metapbs/lsb_extract.go:
//   BinaryScale = 1 << 63     — binary output encoding (0 or Q/2)
//   MaxExtractableBit(N)      — max k such that 2^{k+1} | 2N
//   BuildKthBitTestVector(N, k) -> poly
//   BuildLSBTestVector(N) -> poly
//   ExtractKthBit(ct, brEval, trEvals, cfg, k) -> LWECiphertext
//   ExtractLSB(...) -> LWECiphertext
//   DecodeBinaryCout / DecodeBinaryGLWE
//
// Binary encoding: m_bit ∈ {0,1} is encoded as 0 or Q/2 = BinaryScale.
// This allows direct use as a boolean ciphertext in subsequent computation.
//
// Algorithm:
//   1. Extract bit k from m ∈ Z_t: bit_k(m) = floor(m / 2^k) mod 2
//   2. The function f(m) = bit_k(m) ∈ {0, 1}
//   3. Use Meta-PBS Algorithm 1 to evaluate f
//   4. Output is encoded in binary scale (0 or Q/2)
// =============================================================

#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>

#include "params.hpp"
#include "tlwe.hpp"
#include "trlwe.hpp"

#include "metapbs2/sym_range.hpp"
#include "metapbs2/metapbs_pipeline.hpp"

namespace MetaPBS2 {

// =============================================================
// BinaryScale = Q/2 = 2^63  (for uint64_t torus)
// Binary encoding: 0 → 0,  1 → BinaryScale = Q/2
// =============================================================
constexpr uint64_t BinaryScale = uint64_t(1) << 63;

// =============================================================
// MaxExtractableBit: largest k such that 2^k divides t
// Equivalently: the number of trailing zeros in t (in binary).
//
// For t = 4096 = 2^12: MaxExtractableBit = 11 (bits 0..11).
//
// Mirrors: func MaxExtractableBit(N int) int
// =============================================================
inline int MaxExtractableBit(int t) {
    if (t <= 0) return -1;
    int count = 0;
    while ((t & 1) == 0) {
        count++;
        t >>= 1;
    }
    return count - 1;  // bits 0..count-1 are extractable
}

// =============================================================
// BuildKthBitTestVector: LUT for extracting the k-th bit of m ∈ Z_t.
//
// f(m) = floor(m / 2^k) mod 2
//
// The output is encoded in BinaryScale (0 or Q/2) for binary representation.
// The LUT is built with paper scale Q/t to ensure correct window coverage.
//
// Mirrors: func BuildKthBitTestVector(N, k int) Polynomial[uint64]
// =============================================================
template <class targetP>
TFHEpp::Polynomial<targetP> BuildKthBitTestVector(int k, int t) {
    static_assert(std::is_same_v<typename targetP::T, uint64_t>,
                  "BuildKthBitTestVector requires uint64_t torus type");

    if (k < 0 || k > MaxExtractableBit(t))
        throw std::invalid_argument("k out of range for given t");

    // f(m) = (m >> k) & 1, encoded as 0 or Q/2 = BinaryScale
    auto f = [k, t](int m) -> int {
        m = ((m % t) + t) % t;
        return (m >> k) & 1;
    };

    // BuildMetaPBSTV uses paper scale Q/t
    return BuildMetaPBSTV<targetP>(f, t);
}

// =============================================================
// BuildLSBTestVector: LUT for extracting bit 0 (LSB).
//
// Mirrors: func BuildLSBTestVector(N int) Polynomial[uint64]
// =============================================================
template <class targetP>
TFHEpp::Polynomial<targetP> BuildLSBTestVector(int t) {
    return BuildKthBitTestVector<targetP>(0, t);
}

// =============================================================
// DecodeBinaryOutput: decode a uint64_t phase to 0 or 1.
// Phase ≈ 0 → 0,  Phase ≈ BinaryScale = Q/2 → 1.
// =============================================================
inline int DecodeBinaryPhase(uint64_t phase) {
    // Q/4 = BinaryScale/2: midpoint threshold
    uint64_t threshold = BinaryScale / 2;
    // Wrap phase to [0, Q/2) range
    if (phase >= BinaryScale)
        phase -= BinaryScale;
    return (phase >= threshold) ? 1 : 0;
}

// =============================================================
// DecodeBinaryCout: decrypt and decode a GLWE-key LWE ciphertext
// (sample-extracted from TRLWE) to a binary value {0, 1}.
//
// Mirrors: func DecodeBinaryCout(enc, cout, params) int
// =============================================================
template <class P>
int DecodeBinaryCout(
    const TFHEpp::TLWE<P>& cout,
    const TFHEpp::Key<P>& key) {
    constexpr int N = P::n;
    constexpr int k_rank = P::k;

    uint64_t phase = cout[k_rank * N];  // b term
    for (int ki = 0; ki < k_rank; ki++)
        for (int i = 0; i < N; i++)
            phase += cout[ki * N + i] *
                     static_cast<uint64_t>(key[ki * N + i]);

    return DecodeBinaryPhase(phase);
}

// =============================================================
// DecodeBinaryGLWE: decrypt and decode a TRLWE ciphertext's
// constant coefficient to a binary value {0, 1}.
//
// Mirrors: func DecodeBinaryGLWE(enc, ct) int
// =============================================================
template <class P>
int DecodeBinaryGLWE(
    const TFHEpp::TRLWE<P>& ct,
    const TFHEpp::Key<P>& key) {
    constexpr int N = P::n;

    // Phase = b[0] - sum_i A_i[0] * s_i[0]
    uint64_t phase = ct[P::k][0];  // body[0]
    for (uint32_t ki = 0; ki < P::k; ki++)
        phase -= ct[ki][0] * static_cast<uint64_t>(key[ki * N + 0]);

    return DecodeBinaryPhase(phase);
}

// =============================================================
// ExtractKthBit: extract the k-th bit of the message using Meta-PBS.
//
// Input:  LWE ciphertext encrypting m ∈ Z_t (encoded at paper scale Q/t)
// Output: LWE ciphertext encrypting bit_k(m) ∈ {0,1} (encoded at Q/2)
//
// Mirrors:
//   func ExtractKthBit(ct, brEval, trEvals, cfg, k) LWECiphertext
//
// Note: The output is a TLWE<tgtP> (GLWE-key LWE, i.e., sample-extracted).
//       Call SampleExtract from tgtP to convert to standard LWE if needed.
// =============================================================
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
ExtractKthBit(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    int k, int t,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg) {
    using tgtP = typename brP::targetP;

    // Build the k-th bit test vector
    auto tv = BuildKthBitTestVector<tgtP>(k, t);

    // Run Algorithm 1 with this TV
    // We use the internal pipeline directly with the custom TV:
    using domP = typename brP::domainP;
    constexpr int N = tgtP::n;

    // Compute round parameters
    std::vector<int> redundancies(cfg.K + 1);
    std::vector<int> deltas(cfg.K);
    redundancies[0] = 2 * N / cfg.t;
    for (int rnd = 0; rnd < cfg.K; rnd++) {
        deltas[rnd] = DeltaOffset(redundancies[rnd], cfg.rounds[rnd].beta);
        redundancies[rnd + 1] = redundancies[rnd] * cfg.rounds[rnd].beta;
    }

    // Initial BlindRotate with k-th bit TV
    TFHEpp::TRLWE<tgtP> C0;
    TFHEpp::BlindRotate<brP>(C0, ct, bkfft, tv);

    // Initial HomDivRem
    TFHEpp::TLWE<domP> cquo0, crem0;
    HomDivRemLWE<domP>(cquo0, crem0, ct, static_cast<typename domP::T>(2 * N));

    int current_mod = 2 * N;
    TFHEpp::TRLWE<tgtP> prev_GLWE = C0;
    TFHEpp::TLWE<domP> prev_crem = crem0;

    for (int rnd = 0; rnd < cfg.K; rnd++) {
        const auto& round = cfg.rounds[rnd];
        auto [a_k, b_k] = SymRange(redundancies[rnd]);
        int delta_k = deltas[rnd];

        TFHEpp::TRLWE<tgtP> Ck_prime;
        HomTruncRepeatShifted<tgtP>(Ck_prime, prev_GLWE,
                                     a_k, b_k, round.beta, delta_k,
                                     trkeys[rnd]);

        TFHEpp::TLWE<domP> cquo_k, crem_k;
        HomDivRemAtScale<domP>(cquo_k, crem_k, prev_crem, current_mod, round.beta);

        TFHEpp::TRLWE<tgtP> Ck;
        BlindRotateGLWEFromQuotient<brP>(Ck, cquo_k, Ck_prime, bkfft, N);

        current_mod *= round.beta;
        prev_GLWE = Ck;
        prev_crem = crem_k;
    }

    TFHEpp::TLWE<tgtP> cout;
    SampleExtractIndex0<tgtP>(cout, prev_GLWE);
    return cout;
}

// =============================================================
// ExtractLSB: convenience wrapper for k=0 (LSB extraction)
// =============================================================
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
ExtractLSB(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    int t,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg) {
    return ExtractKthBit<brP>(ct, 0, t, bkfft, trkeys, cfg);
}

}  // namespace MetaPBS2
