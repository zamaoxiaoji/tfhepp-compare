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
#include <type_traits>

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
template <typename T>
constexpr T BinaryScaleT = T(1) << (std::numeric_limits<T>::digits - 1);

constexpr uint64_t BinaryScale = BinaryScaleT<uint64_t>;

// =============================================================
// MaxExtractableBit: largest supported k for t=2N binary 0/Q/2 output.
//
// The top bit of Z_{2N} is not compatible with negacyclic binary 0/Q/2
// encoding because bit_k(m+N) must equal bit_k(m). For N=2048, this returns
// 10, matching tfhe-go's MaxExtractableBit(N).
//
// Mirrors: func MaxExtractableBit(N int) int
// =============================================================
inline int MaxExtractableBit(int N) {
    if (N <= 0) return -1;
    int width = 0;
    unsigned int x = static_cast<unsigned int>(N);
    while (x != 0) {
        width++;
        x >>= 1;
    }
    return width - 2;
}

// =============================================================
// BuildKthBitTestVector: LUT for extracting the k-th bit of m ∈ Z_t.
//
// f(m) = floor(m / 2^k) mod 2
//
// The output is encoded in BinaryScale (0 or Q/2). This C++ overload keeps the
// existing t argument, but it intentionally supports only the Go-validated
// t=2N paper-row semantics.
//
// Mirrors: func BuildKthBitTestVector(N, k int) Polynomial[uint64]
// =============================================================
template <class targetP>
TFHEpp::Polynomial<targetP> BuildKthBitTestVector(int k, int t) {
    constexpr int N = targetP::n;
    if (t != 2 * N)
        throw std::invalid_argument("BuildKthBitTestVector currently requires t=2N");
    if (k < 0 || k > MaxExtractableBit(N))
        throw std::invalid_argument("k out of range for given t");

    TFHEpp::Polynomial<targetP> tv = {};
    for (int j = 0; j < N; j++)
        if (((j >> k) & 1) == 1)
            tv[j] = BinaryScaleT<typename targetP::T>;
    return tv;
}

template <class targetP>
TFHEpp::Polynomial<targetP> BuildKthBitTestVector(int k) {
    return BuildKthBitTestVector<targetP>(k, 2 * targetP::n);
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
// DecodeBinaryOutput: decode a T phase to 0 or 1.
// Phase ≈ 0 → 0,  Phase ≈ BinaryScale = Q/2 → 1.
// =============================================================
template <typename T>
inline int DecodeBinaryPhase(T phase) {
    static_assert(std::is_unsigned_v<T>, "DecodeBinaryPhase expects unsigned torus type");
    T binScale = BinaryScaleT<T>;
    return static_cast<int>(((phase + binScale / 2) / binScale) & T(1));
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

    typename P::T phase = cout[k_rank * N];  // b term
    for (int ki = 0; ki < k_rank; ki++)
        for (int i = 0; i < N; i++)
            phase -= cout[ki * N + i] *
                     static_cast<typename P::T>(key[ki * N + i]);

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
    auto phase = TFHEpp::trlwePhase<P>(ct, key);
    return DecodeBinaryPhase(phase[0]);
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
    const Algorithm1Config& cfg,
    int first_blind_rotate_period = 0,
    BlindRotatePruneStats* prune_stats = nullptr) {
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

    // Initial HomDivRem
    TFHEpp::TLWE<domP> cquo0, crem0;
    HomDivRemLWE<domP>(cquo0, crem0, ct, static_cast<typename domP::T>(2 * N));

    // Initial BlindRotate with k-th bit TV from quotient coordinates.
    TFHEpp::TRLWE<tgtP> C0;
    BlindRotateGLWEFromQuotient<brP>(
        C0, cquo0, MakeAccumulator<tgtP>(tv), bkfft, N,
        first_blind_rotate_period, prune_stats);

    int current_mod = 2 * N;
    TFHEpp::TRLWE<tgtP> prev_GLWE = C0;
    TFHEpp::TLWE<domP> prev_crem = crem0;

    for (int rnd = 0; rnd < cfg.K; rnd++) {
        const auto& round = cfg.rounds[rnd];
        int a_k = -round.T;
        int b_k = round.T;
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
    const Algorithm1Config& cfg,
    int first_blind_rotate_period = 0,
    BlindRotatePruneStats* prune_stats = nullptr) {
    return ExtractKthBit<brP>(
        ct, 0, t, bkfft, trkeys, cfg,
        first_blind_rotate_period, prune_stats);
}

}  // namespace MetaPBS2
