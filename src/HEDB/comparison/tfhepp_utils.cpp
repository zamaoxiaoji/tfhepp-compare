#include "tfhepp_utils.h"
#include "HEDB/utils/types.h"

namespace TFHEpp
{

// ================================================================
//  LOG ↔ ARI conversions  (identical to HE3DB)
// ================================================================

void ARI_to_LOG(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                const EvalKey &ek)
{
    lvl1param::T μ = lvl1param::μ;
    TLWE<lvl1param> tlweoffset = tlwe;
    tlweoffset[lvl1param::k * lvl1param::n] += μ;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
}

void LOG_to_ARI(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                const EvalKey &ek)
{
    lvl1param::T μ = lvl1param::μ;
    μ = μ << 1;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(-μ));
    res[lvl1param::k * lvl1param::n] += lvl1param::μ;
}

void log_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                 uint32_t scale_bits, const EvalKey &ek)
{
    lvl1param::T μ = 1ULL << (scale_bits - 1);
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(-μ));
    res[lvl1param::k * lvl1param::n] += μ;
}

void ari_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                 uint32_t scale_bits, const EvalKey &ek)
{
    lvl1param::T μ = 1ULL << (scale_bits - 1);
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
    TLWE<lvl1param> tlweoffset = tlwe;
    tlweoffset[lvl1param::k * lvl1param::n] += offset;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
    res[lvl1param::k * lvl1param::n] += μ;
}

// ================================================================
//  Standard MSB Gate Bootstrapping  (HE3DB original)
// ================================================================

void MSBGateBootstrapping(TLWE<lvl1param> &res,
                          const TLWE<lvl1param> &tlwe,
                          const EvalKey &ek, bool result_type)
{
    uint32_t μ = 1U << 29;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
    TLWE<lvl1param> tlweoffset = tlwe;
    tlweoffset[lvl1param::k * lvl1param::n] += offset;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
    if (IS_ARITHMETIC(result_type))
        res[lvl1param::k * lvl1param::n] += μ;
}

void MSBGateBootstrapping(TLWE<lvl1param> &res,
                          const TLWE<lvl2param> &tlwe,
                          const EvalKey &ek, bool result_type)
{
    uint32_t μ = 1U << 29;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl2param::T>::digits - 6);
    TLWE<lvl2param> tlweoffset = tlwe;
    tlweoffset[lvl2param::k * lvl2param::n] += offset;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
    if (IS_ARITHMETIC(result_type))
        res[lvl1param::k * lvl1param::n] += μ;
}

void MSBGateBootstrapping(TLWE<lvl2param> &res,
                          const TLWE<lvl2param> &tlwe,
                          const EvalKey &ek, bool result_type)
{
    uint64_t μ = 1ULL << 61;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl2param::T>::digits - 7);
    TLWE<lvl2param> tlweoffset = tlwe;
    tlweoffset[lvl2param::k * lvl2param::n] += offset;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
    GateBootstrappingTLWE2TLWEFFT<lvl02param>(
        res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<lvl2param>(μ));
    if (IS_ARITHMETIC(result_type))
        res[lvl2param::k * lvl2param::n] += μ;
}

// ================================================================
//  Identity Gate Bootstrapping  (rescale)
// ================================================================

void IdeGateBootstrapping(TLWE<lvl1param> &res,
                          const TLWE<lvl1param> &tlwe,
                          uint32_t scale_bits, const EvalKey &ek)
{
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
    TLWE<lvl1param> tlweoffset = tlwe;
    tlweoffset[lvl1param::k * lvl1param::n] += offset;
    constexpr uint32_t plain_bits = 4;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01,
        gpolygen<lvl1param>(plain_bits, scale_bits));
}

void IdeGateBootstrapping(TLWE<lvl1param> &res,
                          const TLWE<lvl2param> &tlwe,
                          uint32_t scale_bits, const EvalKey &ek)
{
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl2param::T>::digits - 6);
    TLWE<lvl2param> tlweoffset = tlwe;
    tlweoffset[lvl2param::k * lvl2param::n] += offset;
    constexpr uint32_t plain_bits = 4;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01,
        gpolygen<lvl1param>(plain_bits, scale_bits - 32));
}

void IdeGateBootstrapping(TLWE<lvl2param> &res,
                          const TLWE<lvl2param> &tlwe,
                          uint32_t scale_bits, const EvalKey &ek)
{
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl2param::T>::digits - 7);
    TLWE<lvl2param> tlweoffset = tlwe;
    tlweoffset[lvl2param::k * lvl2param::n] += offset;
    constexpr uint32_t plain_bits = 5;
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
    GateBootstrappingTLWE2TLWEFFT<lvl02param>(
        res, tlwelvl0, *ek.bkfftlvl02,
        gpolygen<lvl2param>(plain_bits, scale_bits));
}

// ================================================================
//  Pruned MSB Gate Bootstrapping
// ================================================================
//
//  Uses a periodic test vector (period = period_M) and skips
//  CMUX when ā_i ≡ 0 (mod period_M).
//
//  The test vector encodes the standard MSB function using
//  {0, Q/2} boolean encoding when period_M divides 2N.
//  For period_M == 0, falls back to the standard MSB LUT.

void PrunedMSBGateBootstrapping(TLWE<lvl1param> &res,
                                const TLWE<lvl1param> &tlwe,
                                const EvalKey &ek,
                                bool result_type,
                                uint32_t period_M)
{
    // μ: output encoding value
    uint32_t μ = 1U << 29;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;

    // Standard offset (same as HE3DB)
    constexpr uint32_t plain_bits_eff = 5;
    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl1param::T>::digits -
                 plain_bits_eff - 1);

    TLWE<lvl1param> tlweoffset = tlwe;
    tlweoffset[lvl1param::k * lvl1param::n] += offset;

    // Key switch lvl1 → lvl0
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);

    if (period_M > 0) {
        // Use periodic test vector + pruned blind rotation
        auto tv = PeriodicTestVector<lvl1param>(period_M);
        PrunedGateBootstrappingTLWE2TLWEFFT<lvl01param>(
            res, tlwelvl0, *ek.bkfftlvl01, tv, period_M);
    } else {
        // Standard MSB LUT
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
    }

    if (IS_ARITHMETIC(result_type))
        res[lvl1param::k * lvl1param::n] += μ;
}

void PrunedMSBGateBootstrapping(TLWE<lvl2param> &res,
                                const TLWE<lvl2param> &tlwe,
                                const EvalKey &ek,
                                bool result_type,
                                uint32_t period_M)
{
    uint64_t μ = 1ULL << 61;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;

    constexpr uint64_t offset =
        1ULL << (std::numeric_limits<lvl2param::T>::digits - 7);
    TLWE<lvl2param> tlweoffset = tlwe;
    tlweoffset[lvl2param::k * lvl2param::n] += offset;

    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);

    if (period_M > 0) {
        auto tv = PeriodicTestVector<lvl2param>(period_M);
        PrunedGateBootstrappingTLWE2TLWEFFT<lvl02param>(
            res, tlwelvl0, *ek.bkfftlvl02, tv, period_M);
    } else {
        GateBootstrappingTLWE2TLWEFFT<lvl02param>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<lvl2param>(μ));
    }

    if (IS_ARITHMETIC(result_type))
        res[lvl2param::k * lvl2param::n] += μ;
}

// ================================================================
//  Gap-MSB Gate Bootstrapping
// ================================================================
//
//  After the caller has cleared bit guard_k from the ciphertext
//  (by subtracting w_k · ct_k), this function performs the final
//  MSB extraction with the gap-aligned offset:
//
//      offset' = (w_k + 1) / 2 * Δ
//
//  where w_k = 2^(plain_bits - 1 - guard_k) and Δ = Q / 2^plain_bits.
//
//  The client must decrypt using the shifted boundary Δ'.

void GapMSBGateBootstrapping(TLWE<lvl1param> &res,
                             const TLWE<lvl1param> &tlwe,
                             const EvalKey &ek,
                             bool result_type,
                             uint32_t guard_k,
                             uint32_t plain_bits)
{
    uint32_t μ = 1U << 29;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;

    // Compute gap offset:
    //   w_k = 2^(plain_bits - 1 - guard_k)
    //   Δ   = Q / 2^plain_bits = 2^(q - plain_bits)  where q = 32
    //   offset' = (w_k + 1) / 2 * Δ
    //           = (2^(p-1-k) + 1) * 2^(q - p - 1)
    constexpr uint32_t q = std::numeric_limits<lvl1param::T>::digits;
    const uint64_t w_k = 1ULL << (plain_bits - 1 - guard_k);
    const uint64_t delta = 1ULL << (q - plain_bits);
    const uint64_t gap_offset = ((w_k + 1) * delta) / 2;

    TLWE<lvl1param> tlweoffset = tlwe;
    tlweoffset[lvl1param::k * lvl1param::n] +=
        static_cast<lvl1param::T>(gap_offset);

    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));

    if (IS_ARITHMETIC(result_type))
        res[lvl1param::k * lvl1param::n] += μ;
}

void GapMSBGateBootstrapping(TLWE<lvl2param> &res,
                             const TLWE<lvl2param> &tlwe,
                             const EvalKey &ek,
                             bool result_type,
                             uint32_t guard_k,
                             uint32_t plain_bits)
{
    uint64_t μ = 1ULL << 61;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;

    constexpr uint32_t q = std::numeric_limits<lvl2param::T>::digits;
    const uint64_t w_k = 1ULL << (plain_bits - 1 - guard_k);
    const uint64_t delta = 1ULL << (q - plain_bits);
    const uint64_t gap_offset = ((w_k + 1) * delta) / 2;

    TLWE<lvl2param> tlweoffset = tlwe;
    tlweoffset[lvl2param::k * lvl2param::n] +=
        static_cast<lvl2param::T>(gap_offset);

    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
    GateBootstrappingTLWE2TLWEFFT<lvl02param>(
        res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<lvl2param>(μ));

    if (IS_ARITHMETIC(result_type))
        res[lvl2param::k * lvl2param::n] += μ;
}

} // namespace TFHEpp
