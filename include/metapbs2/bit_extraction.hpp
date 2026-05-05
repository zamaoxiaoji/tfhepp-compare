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

struct BitExtractOptions {
    int p;       // Message precision. For the current exact path, cfg.t = 2^p = 2N.
    int k;       // Chapter-3 bit index: k=0 is the top/sign bit, k=p-1 is LSB.
    bool enable_periodic_pruning = true;
    int period = 0;  // 0 means derive M_k from p,k; otherwise use this pruning period.
};

inline int CheckedPowerOfTwo(int exponent) {
    if (exponent < 0 || exponent >= static_cast<int>(std::numeric_limits<int>::digits))
        throw std::invalid_argument("power-of-two exponent out of int range");
    return int(1) << exponent;
}

inline int LSBIndexFromChapterBit(int p, int k) {
    if (p <= 0)
        throw std::invalid_argument("BitExtract requires positive message precision p");
    if (k < 0 || k >= p)
        throw std::invalid_argument("BitExtract bit index k out of range for p");
    return p - 1 - k;
}

inline int BitExtractPeriodFromChapterBit(int p, int k) {
    return CheckedPowerOfTwo(LSBIndexFromChapterBit(p, k) + 1);
}

inline int MessagePrecisionFromPowerOfTwoModulus(int t) {
    if (t <= 0 || (t & (t - 1)) != 0)
        throw std::invalid_argument("message modulus must be a positive power of two");
    int p = 0;
    while ((int(1) << p) < t) p++;
    return p;
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
    auto tv = BuildKthBitTestVector<tgtP>(k, t);
    return RunAlgorithm1WithTV<brP>(
        ct, tv, bkfft, trkeys, cfg, first_blind_rotate_period, prune_stats);
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

// =============================================================
// BitExtract: Chapter-3 periodic-pruned iterative bit extraction.
//
// This is the paper-facing wrapper for algorithm
// alg:full_periodic_pruned_iterative_pbs:
//   1. Interpret k in the Chapter-3/MSB-side convention.
//   2. Compute w_k=2^{p-1-k} and M_k=2w_k automatically.
//   3. Run the Meta-PBS exact extraction path.
//   4. Apply M_k-periodic pruning only in the first blind rotation.
//
// The current exact implementation supports the t=2N paper row. Thus cfg.t
// must equal 2^p and the converted LSB-side bit must be in the supported range.
// =============================================================
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
BitExtract(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const BitExtractOptions& options,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    constexpr int N = tgtP::n;

    if (cfg.t != CheckedPowerOfTwo(options.p))
        throw std::invalid_argument("BitExtract requires cfg.t == 2^p");
    if (cfg.t != 2 * N)
        throw std::invalid_argument("BitExtract currently supports only cfg.t=2N");

    int lsb_k = LSBIndexFromChapterBit(options.p, options.k);
    if (lsb_k > MaxExtractableBit(N))
        throw std::invalid_argument("BitExtract requested bit is not supported by 0/Q/2 negacyclic encoding");

    int first_period = 0;
    if (options.enable_periodic_pruning)
        first_period = options.period > 0
                           ? options.period
                           : BitExtractPeriodFromChapterBit(options.p, options.k);

    return ExtractKthBit<brP>(
        ct, lsb_k, cfg.t, bkfft, trkeys, cfg, first_period, prune_stats);
}

template <class brP>
TFHEpp::TLWE<typename brP::targetP>
BitExtract(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    int p, int k,
    BlindRotatePruneStats* prune_stats = nullptr) {
    return BitExtract<brP>(
        ct, bkfft, trkeys, cfg,
        BitExtractOptions{.p = p, .k = k},
        prune_stats);
}

template <typename T>
inline bool IsTorusSelfNegating(T x) {
    static_assert(std::is_unsigned_v<T>, "IsTorusSelfNegating expects unsigned torus");
    return x == T(0) || x == BinaryScaleT<T>;
}

inline bool KthBitInvariantUnderHalfTurn(int k, int N) {
    if (k < 0) return false;
    if (k + 1 >= static_cast<int>(std::numeric_limits<unsigned int>::digits))
        return false;
    return (N % (1u << (k + 1))) == 0;
}

template <typename T>
inline bool WeightedBitNegacyclicCompatible(int k, int N, T weight_torus) {
    return KthBitInvariantUnderHalfTurn(k, N) && IsTorusSelfNegating(weight_torus);
}

// =============================================================
// BuildWeightedBitTestVector: guarded direct arithmetic-weight LUT.
//
// This path is only valid when the resulting TV satisfies the negacyclic
// constraint required by blind rotation. For the usual GapMSB clear-bit
// weights, that condition is false unless weight_torus is self-negating
// (0 or Q/2), so callers should be prepared to fall back to 0/Q/2 extraction
// followed by LOG_to_ARI.
// =============================================================
template <class targetP>
TFHEpp::Polynomial<targetP> BuildWeightedBitTestVector(
    int k, int t, typename targetP::T weight_torus) {
    constexpr int N = targetP::n;
    if (t != 2 * N)
        throw std::invalid_argument("BuildWeightedBitTestVector currently requires t=2N");
    if (k < 0 || k > MaxExtractableBit(N))
        throw std::invalid_argument("k out of range for given t");
    if (!WeightedBitNegacyclicCompatible(k, N, weight_torus))
        throw std::invalid_argument("direct weighted bit LUT violates negacyclic encoding");

    TFHEpp::Polynomial<targetP> tv = {};
    for (int j = 0; j < N; j++)
        if (((j >> k) & 1) == 1)
            tv[j] = weight_torus;
    return tv;
}

// =============================================================
// ExtractKthBitWeighted: extract bit k with direct arithmetic weight output.
//
// Same Meta-PBS pipeline as ExtractKthBit, but the LUT outputs
// 0/weight_torus instead of 0/BinaryScale.
//
// Output: TLWE encrypting 0 or weight_torus (arithmetic encoding).
// This can be directly subtracted from the original ciphertext
// to clear the bit, without a separate LOG_to_ARI bootstrap.
// =============================================================
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
ExtractKthBitWeighted(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    int k, int t,
    typename brP::targetP::T weight_torus,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    int first_blind_rotate_period = 0,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    auto tv = BuildWeightedBitTestVector<tgtP>(k, t, weight_torus);
    return RunAlgorithm1WithTV<brP>(
        ct, tv, bkfft, trkeys, cfg, first_blind_rotate_period, prune_stats);
}

// =============================================================
// BitExtractWeighted: Chapter-3 wrapper for weighted bit extraction.
//
// Same as BitExtract but outputs 0/weight_torus directly.
// Eliminates the need for LogicalBitToArithmeticWeight PBS.
// =============================================================
template <class brP>
TFHEpp::TLWE<typename brP::targetP>
BitExtractWeighted(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<TruncRepeatKey<typename brP::targetP>>& trkeys,
    const Algorithm1Config& cfg,
    const BitExtractOptions& options,
    typename brP::targetP::T weight_torus,
    BlindRotatePruneStats* prune_stats = nullptr) {
    using tgtP = typename brP::targetP;
    constexpr int N = tgtP::n;

    if (cfg.t != CheckedPowerOfTwo(options.p))
        throw std::invalid_argument("BitExtractWeighted requires cfg.t == 2^p");
    if (cfg.t != 2 * N)
        throw std::invalid_argument("BitExtractWeighted currently supports only cfg.t=2N");

    int lsb_k = LSBIndexFromChapterBit(options.p, options.k);
    if (lsb_k > MaxExtractableBit(N))
        throw std::invalid_argument("BitExtractWeighted requested bit is not supported by negacyclic encoding");

    int first_period = 0;
    if (options.enable_periodic_pruning)
        first_period = options.period > 0
                           ? options.period
                           : BitExtractPeriodFromChapterBit(options.p, options.k);

    return ExtractKthBitWeighted<brP>(
        ct, lsb_k, cfg.t, weight_torus, bkfft, trkeys, cfg, first_period, prune_stats);
}

}  // namespace MetaPBS2
