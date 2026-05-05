#include "ethmsb.h"
using namespace TFHEpp;

namespace HEDB
{

static constexpr uint32_t KAPPA_LVL1 = 5;
static constexpr uint32_t KAPPA_LVL2 = 6;

// ================================================================
//  ETHMSB lvl1 (Algorithm 2, κ=5)
// ================================================================

void ETHMSB_lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                 uint32_t k, const TFHEEvalKey &ek,
                 bool result_type)
{
    constexpr uint32_t kappa = KAPPA_LVL1;
    constexpr uint32_t q = std::numeric_limits<Lvl1::T>::digits;

    if (k <= kappa) {
        // Base case: direct PBS-MSB with offset = Q/2^(κ+1)
        Lvl1::T mu_out = 1U << 29;
        if (IS_ARITHMETIC(result_type)) mu_out = mu_out << 1;
        constexpr Lvl1::T offset = 1ULL << (q - kappa - 1);
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;
        TLWELvl0 tlwelvl0;
        IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01,
            μ_polygen<Lvl1>(mu_out));
        if (IS_ARITHMETIC(result_type))
            res[Lvl1::k * Lvl1::n] += mu_out;
        return;
    }

    // Recursive case: k > κ
    // The ciphertext encrypts m at scale Δ = Q/2^k.
    // Bit b_κ has weight 2^(k-1-κ) in message space,
    // which corresponds to phase 2^(k-1-κ) · Δ = Q/2^(κ+1).
    const uint32_t shift = k - kappa;

    // Step 1: left-shift (k-κ) bits to align b_κ to MSB of κ-bit window
    TLWELvl1 shift_tlwe;
    for (size_t i = 0; i <= Lvl1::n; i++)
        shift_tlwe[i] = tlwe[i] << shift;
    // Shifted phase: 2^shift · Δ · m = (Q/2^κ) · m
    // The shifted κ-bit window has scale Q/2^κ.

    // Step 2: LUT to extract b_κ at its original weight
    // Paper: coefficients ±2^(k-κ-1).
    // But these are in MESSAGE-space units of the original encoding.
    // In Q-space (phase): 2^(k-κ-1) · Δ = 2^(k-κ-1) · Q/2^k = Q/2^(κ+1)
    // So μ for the LUT = Q/2^(κ+1):
    const Lvl1::T mu_guard = 1U << (q - kappa - 1); // Q/2^(κ+1)
    // μ_polygen(mu_guard): all coefficients = -mu_guard
    // PBS output: MSB=0 (positive phase) → -mu_guard
    //             MSB=1 (negative phase) → +mu_guard
    // Add mu_guard: MSB=0 → 0, MSB=1 → 2·mu_guard = Q/2^κ
    // But we want Q/2^(κ+1) for MSB=1...
    // So use half: mu_half = Q/2^(κ+2)
    const Lvl1::T mu_half = 1U << (q - kappa - 2); // Q/2^(κ+2)
    // PBS output + mu_half: MSB=0 → 0, MSB=1 → 2·mu_half = Q/2^(κ+1) ✓

    // Step 3: offset₂ = (2^(κ-1)+1)/2 · Δ
    // Δ = Q/2^k, offset₂ = (2^(κ-1)+1) · Q/2^(k+1)
    // This is added to the SHIFTED ciphertext body.
    const uint64_t delta = 1ULL << (q - k);
    const uint64_t w_kappa = 1ULL << (kappa - 1);
    const Lvl1::T offset2 =
        static_cast<Lvl1::T>(((w_kappa + 1) * delta) / 2);
    shift_tlwe[Lvl1::k * Lvl1::n] += offset2;

    // PBS
    TLWELvl0 tlwelvl0;
    IdentityKeySwitch<Lvl10>(tlwelvl0, shift_tlwe, *ek.iksklvl10);
    TLWELvl1 guard_ct;
    GateBootstrappingTLWE2TLWEFFT<Lvl01>(
        guard_ct, tlwelvl0, *ek.bkfftlvl01,
        μ_polygen<Lvl1>(mu_half));
    guard_ct[Lvl1::k * Lvl1::n] += mu_half;
    // guard_ct now encrypts b_κ · Q/2^(κ+1) = b_κ · 2^(k-κ-1) · Δ

    // Step 4: subtract → m' = m - b_κ · 2^(k-κ-1)
    TLWELvl1 guarded;
    for (size_t i = 0; i <= Lvl1::n; i++)
        guarded[i] = tlwe[i] - guard_ct[i];

    // Step 5: recurse with k-κ
    ETHMSB_lvl1(res, guarded, k - kappa, ek, result_type);
}

// ================================================================
//  ETHMSB lvl2
// ================================================================

void ETHMSB_lvl2(TLWELvl1 &res, const TLWELvl2 &tlwe,
                 uint32_t k, const TFHEEvalKey &ek,
                 bool result_type)
{
    constexpr uint32_t kappa = KAPPA_LVL2;
    constexpr uint32_t q = std::numeric_limits<Lvl2::T>::digits;

    if (k <= KAPPA_LVL1) {
        IdentityKeySwitch<lvl21param>(res, tlwe, *ek.iksklvl21);
        ETHMSB_lvl1(res, res, k, ek, result_type);
        return;
    }

    if (k <= kappa) {
        Lvl1::T mu_out = 1U << 29;
        if (IS_ARITHMETIC(result_type)) mu_out = mu_out << 1;
        constexpr Lvl2::T offset = 1ULL << (q - kappa - 1);
        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;
        TLWELvl0 tlwelvl0;
        IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01,
            μ_polygen<Lvl1>(mu_out));
        if (IS_ARITHMETIC(result_type))
            res[Lvl1::k * Lvl1::n] += mu_out;
        return;
    }

    const uint32_t shift = k - kappa;

    TLWELvl2 shift_tlwe;
    for (size_t i = 0; i <= Lvl2::n; i++)
        shift_tlwe[i] = tlwe[i] << shift;

    // Guard bit weight in Q-space: Q/2^(κ+1)
    const Lvl2::T mu_half = 1ULL << (q - kappa - 2);

    const uint64_t delta = 1ULL << (q - k);
    const uint64_t w_kappa = 1ULL << (kappa - 1);
    const Lvl2::T offset2 =
        static_cast<Lvl2::T>(((w_kappa + 1) * delta) / 2);
    shift_tlwe[Lvl2::k * Lvl2::n] += offset2;

    TLWELvl0 tlwelvl0;
    IdentityKeySwitch<Lvl20>(tlwelvl0, shift_tlwe, *ek.iksklvl20);
    TLWELvl2 guard_ct;
    GateBootstrappingTLWE2TLWEFFT<Lvl02>(
        guard_ct, tlwelvl0, *ek.bkfftlvl02,
        μ_polygen<Lvl2>(mu_half));
    guard_ct[Lvl2::k * Lvl2::n] += mu_half;

    TLWELvl2 guarded;
    for (size_t i = 0; i <= Lvl2::n; i++)
        guarded[i] = tlwe[i] - guard_ct[i];

    ETHMSB_lvl2(res, guarded, k - kappa, ek, result_type);
}

// Dispatchers
void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
               uint32_t k, const TFHEEvalKey &ek, bool result_type)
{ ETHMSB_lvl1(res, tlwe, k, ek, result_type); }

void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
               uint32_t k, const TFHEEvalKey &ek, bool result_type)
{ ETHMSB_lvl2(res, tlwe, k, ek, result_type); }

} // namespace HEDB
