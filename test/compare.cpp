// TFHEpp/test/compare.cpp
// Comparison test: HE3DB vs Meta-PBS
//
// Following HE3DB comparison_test.cpp, we test all 5 operators:
//   >, >=, <, <=, ==
//
// Two paths:
//   HE3DB:    recursive MSB clearing (reimplemented in TFHEpp APIs)
//   MetaPBS:  in-window  → MetaPBSExtractBit2N MSB extraction (periodic pruning)
//             out-of-window → HE3DB reduction + extract LSB (period=2) + clear +
//                             offset + normal PBS
//
// Build & Run:
//   cmake .. -DENABLE_TEST=ON && make -j compare
//   ./test/compare --trials 100

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <metapbs.hpp>
#include <tfhe++.hpp>

namespace {

// ─── helpers ────────────────────────────────────────────────────────────────

/// Torus left-shift: multiply every coefficient by 2^shift.
template <class P>
void TLWELeftShift(TFHEpp::TLWE<P> &ct, const int shift)
{
    if (shift <= 0) return;
    for (auto &x : ct)
        x = static_cast<typename P::T>(
            static_cast<typename P::T>(x) << shift);
}

/// ct_out = ct_a − ct_b, element-wise.
template <class P>
void TLWESub(TFHEpp::TLWE<P> &out,
             const TFHEpp::TLWE<P> &a,
             const TFHEpp::TLWE<P> &b)
{
    for (size_t i = 0; i < out.size(); i++)
        out[i] = a[i] - b[i];
}

// ─── HE3DB-style primitives (reimplemented with TFHEpp APIs) ─────────────

/// MSBGateBootstrapping: extract MSB as ±μ (ARITHMETIC) or ±μ/4 (LOGIC).
template <class iksP, class brP>
void HE3DB_MSBGateBootstrapping(
    TFHEpp::TLWE<typename brP::targetP> &res,
    const TFHEpp::TLWE<typename iksP::domainP> &tlwe,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    bool arithmetic)
{
    using TargetP = typename brP::targetP;
    using DomainP = typename iksP::domainP;
    using DomBRP  = typename brP::domainP;

    constexpr int digits = std::numeric_limits<typename DomainP::T>::digits;
    const auto offset =
        static_cast<typename DomainP::T>(
            static_cast<typename DomainP::T>(1) << (digits - 6));

    TFHEpp::TLWE<DomainP> tlwe_off = tlwe;
    tlwe_off[DomainP::k * DomainP::n] += offset;

    TFHEpp::TLWE<DomBRP> dom{};
    TFHEpp::IdentityKeySwitch<iksP>(dom, tlwe_off, ksk);

    typename TargetP::T mu = TargetP::μ;
    if (arithmetic) mu = static_cast<typename TargetP::T>(mu << 1);

    TFHEpp::Polynomial<TargetP> tv{};
    for (auto &p : tv) p = static_cast<typename TargetP::T>(-mu);

    TFHEpp::GateBootstrappingTLWE2TLWE<brP>(res, dom, bk, tv);

    if (arithmetic)
        res[TargetP::k * TargetP::n] += mu;
}

/// IdeGateBootstrapping: identity PBS that maps m → m (rescaled).
template <class iksP, class brP>
void HE3DB_IdeGateBootstrapping(
    TFHEpp::TLWE<typename brP::targetP> &res,
    const TFHEpp::TLWE<typename iksP::domainP> &tlwe,
    uint32_t scale_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk)
{
    using TargetP = typename brP::targetP;
    using DomainP = typename iksP::domainP;
    using DomBRP  = typename brP::domainP;

    constexpr int digits = std::numeric_limits<typename DomainP::T>::digits;
    const auto offset =
        static_cast<typename DomainP::T>(
            static_cast<typename DomainP::T>(1) << (digits - 6));

    TFHEpp::TLWE<DomainP> tlwe_off = tlwe;
    tlwe_off[DomainP::k * DomainP::n] += offset;

    TFHEpp::TLWE<DomBRP> dom{};
    TFHEpp::IdentityKeySwitch<iksP>(dom, tlwe_off, ksk);

    constexpr uint32_t plain_bits_lut = 4;
    constexpr uint32_t padding_bits = TargetP::nbit - plain_bits_lut;
    TFHEpp::Polynomial<TargetP> tv{};
    for (uint32_t i = 0; i < TargetP::n; i++)
        tv[i] = static_cast<typename TargetP::T>(
            static_cast<typename TargetP::T>(1) << scale_bits) *
            static_cast<typename TargetP::T>(i >> padding_bits);

    TFHEpp::GateBootstrappingTLWE2TLWE<brP>(res, dom, bk, tv);
}

// ─── HE3DB recursive HomMSB ─────────────────────────────────────────────

template <class iksP, class brP>
void HE3DB_HomMSB_Lvl1(
    TFHEpp::TLWE<typename brP::targetP> &res,
    const TFHEpp::TLWE<typename brP::targetP> &tlwe,
    uint32_t plain_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    bool result_type)
{
    using P = typename brP::targetP;

    if (plain_bits <= 5) {
        HE3DB_MSBGateBootstrapping<iksP, brP>(res, tlwe, ksk, bk, result_type);
        return;
    }

    constexpr int digits = std::numeric_limits<typename P::T>::digits;
    uint32_t scale_bits = digits - plain_bits;

    TFHEpp::TLWE<P> shift_tlwe = tlwe;
    TLWELeftShift<P>(shift_tlwe, static_cast<int>(plain_bits - 5));

    TFHEpp::TLWE<P> sign_tlwe{};
    HE3DB_MSBGateBootstrapping<iksP, brP>(
        sign_tlwe, shift_tlwe, ksk, bk, /*arithmetic=*/true);

    for (size_t i = 0; i <= P::k * P::n; i++)
        shift_tlwe[i] -= sign_tlwe[i];

    TFHEpp::TLWE<P> cleaned{};
    HE3DB_IdeGateBootstrapping<iksP, brP>(
        cleaned, shift_tlwe, scale_bits, ksk, bk);

    for (size_t i = 0; i <= P::k * P::n; i++)
        res[i] = tlwe[i] - cleaned[i];

    uint32_t new_bits = (plain_bits <= 9) ? 5 : (plain_bits - 4);
    HE3DB_HomMSB_Lvl1<iksP, brP>(res, res, new_bits, ksk, bk, result_type);
}

// ─── HE3DB comparison operators ─────────────────────────────────────────

template <class P, class iksP, class brP>
bool he3db_greater_than(
    const TFHEpp::TLWE<P> &a_ct, const TFHEpp::TLWE<P> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    const TFHEpp::SecretKey &sk)
{
    // a > b ⟺ MSB(b − a) = 1
    TFHEpp::TLWE<P> diff{};
    TLWESub<P>(diff, b_ct, a_ct);
    TFHEpp::TLWE<typename brP::targetP> msb{};
    HE3DB_HomMSB_Lvl1<iksP, brP>(msb, diff, plain_bits + 1, ksk, bk, false);
    return TFHEpp::tlweSymDecrypt<typename brP::targetP>(msb, sk);
}

template <class P, class iksP, class brP>
bool he3db_greater_equal(
    const TFHEpp::TLWE<P> &a_ct, const TFHEpp::TLWE<P> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    const TFHEpp::SecretKey &sk)
{
    // a >= b ⟺ NOT MSB(a − b)
    TFHEpp::TLWE<P> diff{};
    TLWESub<P>(diff, a_ct, b_ct);
    TFHEpp::TLWE<typename brP::targetP> msb{};
    HE3DB_HomMSB_Lvl1<iksP, brP>(msb, diff, plain_bits + 1, ksk, bk, false);
    return !TFHEpp::tlweSymDecrypt<typename brP::targetP>(msb, sk);
}

template <class P, class iksP, class brP>
bool he3db_less_than(
    const TFHEpp::TLWE<P> &a_ct, const TFHEpp::TLWE<P> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    const TFHEpp::SecretKey &sk)
{
    // a < b ⟺ MSB(a − b) = 1
    TFHEpp::TLWE<P> diff{};
    TLWESub<P>(diff, a_ct, b_ct);
    TFHEpp::TLWE<typename brP::targetP> msb{};
    HE3DB_HomMSB_Lvl1<iksP, brP>(msb, diff, plain_bits + 1, ksk, bk, false);
    return TFHEpp::tlweSymDecrypt<typename brP::targetP>(msb, sk);
}

template <class P, class iksP, class brP>
bool he3db_less_equal(
    const TFHEpp::TLWE<P> &a_ct, const TFHEpp::TLWE<P> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    const TFHEpp::SecretKey &sk)
{
    // a <= b ⟺ NOT MSB(b − a)
    TFHEpp::TLWE<P> diff{};
    TLWESub<P>(diff, b_ct, a_ct);
    TFHEpp::TLWE<typename brP::targetP> msb{};
    HE3DB_HomMSB_Lvl1<iksP, brP>(msb, diff, plain_bits + 1, ksk, bk, false);
    return !TFHEpp::tlweSymDecrypt<typename brP::targetP>(msb, sk);
}

template <class P, class iksP, class brP>
bool he3db_equal(
    const TFHEpp::TLWE<P> &a_ct, const TFHEpp::TLWE<P> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk,
    const TFHEpp::SecretKey &sk)
{
    // a == b ⟺ (a >= b) AND (a <= b)
    return he3db_greater_equal<P, iksP, brP>(a_ct, b_ct, plain_bits, ksk, bk, sk) &&
           he3db_less_equal<P, iksP, brP>(a_ct, b_ct, plain_bits, ksk, bk, sk);
}

// ─── Meta-PBS HomMSB ─────────────────────────────────────────────────────
//
// Core MSB extraction using MetaPBS with periodic pruning.
//
// Architecture:
//   In-window  (p+1 ≤ nbit+1): single MetaPBSExtractBit2N at bit=nbit-1
//   Out-of-window (p+1 > nbit+1): HE3DB reduction → extract LSB (period=2)
//                                  → clear → offset → normal sign PBS

template <class TargetP, class DomainP,
          class brP,       // DomainP → TargetP (extraction)
          class ksToDomP,  // TargetP → DomainP (key switch to BS domain)
          class weightBkP, // DomainP → lvl2 (weight bootstrap)
          class ksDownP,   // lvl2 → TargetP (key switch down)
          class iksP_he3db, class brP_he3db>  // for HE3DB MSB reduction
void HomMSB_MetaPBS(
    TFHEpp::TLWE<TargetP> &res,
    const TFHEpp::TLWE<TargetP> &tlwe,
    uint32_t plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP> &bkfft_extract,
    const TFHEpp::AnnihilateKey<TargetP> &ahk,
    const TFHEpp::KeySwitchingKey<ksToDomP> &ksk_to_dom,
    const TFHEpp::BootstrappingKeyFFT<weightBkP> &bkfft_weight,
    const TFHEpp::KeySwitchingKey<ksDownP> &ksk_down,
    const TFHEpp::KeySwitchingKey<iksP_he3db> &ksk_he3db,
    const TFHEpp::BootstrappingKeyFFT<brP_he3db> &bk_he3db)
{
    constexpr int nbit = static_cast<int>(TargetP::nbit);
    constexpr int stable_bits = nbit + 1;

    TFHEpp::TLWE<TargetP> diff = tlwe;
    uint32_t current_bits = plain_bits;

    if (current_bits <= static_cast<uint32_t>(stable_bits)) {
        // ── In-window: direct MSB extraction via MetaPBS ──
        //
        // MetaPBSExtractBit2N at bit=nbit-1 (LSB-0 indexed, i.e., the MSB of
        // the Z_N portion) with period = 2^{nbit}.
        // Output: {0, q/2} → decrypt via tlweSymDecrypt (+μ convention).
        //
        // The sign LUT with TV = [-μ, ...] and standard GateBootstrapping:
        //   m ∈ [0, N)  → -μ → false
        //   m ∈ [N, 2N) → +μ → true
        TFHEpp::TLWE<DomainP> dom{};
        TFHEpp::IdentityKeySwitch<ksToDomP>(dom, diff, ksk_to_dom);

        TFHEpp::Polynomial<TargetP> tv_sign{};
        tv_sign.fill(static_cast<typename TargetP::T>(-TargetP::μ));

        TFHEpp::GateBootstrappingTLWE2TLWE<brP>(res, dom, bkfft_extract, tv_sign);
        return;
    }

    // ── Out-of-window: HE3DB reduction → extract LSB → clear → offset → PBS ──

    // Phase 1: HE3DB-style MSB reduction
    while (current_bits > static_cast<uint32_t>(stable_bits)) {
        constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
        uint32_t scale_bits = digits - current_bits;

        TFHEpp::TLWE<TargetP> shift_tlwe = diff;
        TLWELeftShift<TargetP>(shift_tlwe,
                               static_cast<int>(current_bits - 5));

        TFHEpp::TLWE<TargetP> sign_ct{};
        HE3DB_MSBGateBootstrapping<iksP_he3db, brP_he3db>(
            sign_ct, shift_tlwe, ksk_he3db, bk_he3db, true);

        for (size_t i = 0; i <= TargetP::k * TargetP::n; i++)
            shift_tlwe[i] -= sign_ct[i];

        TFHEpp::TLWE<TargetP> cleaned{};
        HE3DB_IdeGateBootstrapping<iksP_he3db, brP_he3db>(
            cleaned, shift_tlwe, scale_bits, ksk_he3db, bk_he3db);

        for (size_t i = 0; i <= TargetP::k * TargetP::n; i++)
            diff[i] -= cleaned[i];

        current_bits = (current_bits <= 9) ? 5 : (current_bits - 4);
    }

    // Phase 2: Extract LSB (period=2) and clear it
    //
    // After reduction, data fills the full window. The data's LSB in Z_{2N}
    // is at bit position 0 (LSB-0 indexed). MetaPBS extraction at bit=0
    // has period = 2, giving maximum CMUX pruning.
    constexpr int bit_msb_extract = nbit - 1;  // MSB-0 for bit_lsb=0

    TFHEpp::TLWE<DomainP> diff_dom{};
    TFHEpp::IdentityKeySwitch<ksToDomP>(diff_dom, diff, ksk_to_dom);

    TFHEpp::TLWE<TargetP> extracted{};
    TFHEpp::metapbs::ExtractBitInPlaceViaLvl2<
        brP, weightBkP, ksToDomP, ksDownP>(
        extracted, diff_dom, bkfft_extract, ahk,
        ksk_to_dom, bkfft_weight, ksk_down, bit_msb_extract);

    // Clear the LSB: diff -= extracted
    for (size_t i = 0; i <= TargetP::k * TargetP::n; i++)
        diff[i] -= extracted[i];

    // Phase 3: Offset + sign PBS
    //
    // After clearing bit 0, odd values near the sign boundary N become even,
    // creating a gap. The gap-center offset improves noise tolerance.
    // Offset = Δ/2 (half a Z_{2N} unit) to center the gap.
    constexpr auto half_delta = static_cast<typename TargetP::T>(
        static_cast<typename TargetP::T>(1)
        << (std::numeric_limits<typename TargetP::T>::digits - nbit - 2));
    diff[TargetP::k * TargetP::n] += half_delta;

    TFHEpp::TLWE<DomainP> dom{};
    TFHEpp::IdentityKeySwitch<ksToDomP>(dom, diff, ksk_to_dom);

    TFHEpp::Polynomial<TargetP> tv_sign{};
    tv_sign.fill(static_cast<typename TargetP::T>(-TargetP::μ));

    TFHEpp::GateBootstrappingTLWE2TLWE<brP>(res, dom, bkfft_extract, tv_sign);
}

// ─── Meta-PBS comparison operators ──────────────────────────────────────

template <class TargetP, class DomainP,
          class brP, class ksToDomP, class weightBkP, class ksDownP,
          class iksP_he3db, class brP_he3db>
bool metapbs_greater_than(
    const TFHEpp::TLWE<TargetP> &a_ct, const TFHEpp::TLWE<TargetP> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP> &bkfft,
    const TFHEpp::AnnihilateKey<TargetP> &ahk,
    const TFHEpp::KeySwitchingKey<ksToDomP> &ksk_to_dom,
    const TFHEpp::BootstrappingKeyFFT<weightBkP> &bkfft_weight,
    const TFHEpp::KeySwitchingKey<ksDownP> &ksk_down,
    const TFHEpp::KeySwitchingKey<iksP_he3db> &ksk_he3db,
    const TFHEpp::BootstrappingKeyFFT<brP_he3db> &bk_he3db,
    const TFHEpp::SecretKey &sk)
{
    // a > b ⟺ MSB(b − a) = 1
    TFHEpp::TLWE<TargetP> diff{};
    TLWESub<TargetP>(diff, b_ct, a_ct);
    TFHEpp::TLWE<TargetP> msb{};
    HomMSB_MetaPBS<TargetP, DomainP, brP, ksToDomP, weightBkP, ksDownP,
                   iksP_he3db, brP_he3db>(
        msb, diff, plain_bits + 1,
        bkfft, ahk, ksk_to_dom, bkfft_weight, ksk_down, ksk_he3db, bk_he3db);
    return TFHEpp::tlweSymDecrypt<TargetP>(msb, sk);
}

template <class TargetP, class DomainP,
          class brP, class ksToDomP, class weightBkP, class ksDownP,
          class iksP_he3db, class brP_he3db>
bool metapbs_greater_equal(
    const TFHEpp::TLWE<TargetP> &a_ct, const TFHEpp::TLWE<TargetP> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP> &bkfft,
    const TFHEpp::AnnihilateKey<TargetP> &ahk,
    const TFHEpp::KeySwitchingKey<ksToDomP> &ksk_to_dom,
    const TFHEpp::BootstrappingKeyFFT<weightBkP> &bkfft_weight,
    const TFHEpp::KeySwitchingKey<ksDownP> &ksk_down,
    const TFHEpp::KeySwitchingKey<iksP_he3db> &ksk_he3db,
    const TFHEpp::BootstrappingKeyFFT<brP_he3db> &bk_he3db,
    const TFHEpp::SecretKey &sk)
{
    // a >= b ⟺ NOT MSB(a − b)
    TFHEpp::TLWE<TargetP> diff{};
    TLWESub<TargetP>(diff, a_ct, b_ct);
    TFHEpp::TLWE<TargetP> msb{};
    HomMSB_MetaPBS<TargetP, DomainP, brP, ksToDomP, weightBkP, ksDownP,
                   iksP_he3db, brP_he3db>(
        msb, diff, plain_bits + 1,
        bkfft, ahk, ksk_to_dom, bkfft_weight, ksk_down, ksk_he3db, bk_he3db);
    return !TFHEpp::tlweSymDecrypt<TargetP>(msb, sk);
}

template <class TargetP, class DomainP,
          class brP, class ksToDomP, class weightBkP, class ksDownP,
          class iksP_he3db, class brP_he3db>
bool metapbs_less_than(
    const TFHEpp::TLWE<TargetP> &a_ct, const TFHEpp::TLWE<TargetP> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP> &bkfft,
    const TFHEpp::AnnihilateKey<TargetP> &ahk,
    const TFHEpp::KeySwitchingKey<ksToDomP> &ksk_to_dom,
    const TFHEpp::BootstrappingKeyFFT<weightBkP> &bkfft_weight,
    const TFHEpp::KeySwitchingKey<ksDownP> &ksk_down,
    const TFHEpp::KeySwitchingKey<iksP_he3db> &ksk_he3db,
    const TFHEpp::BootstrappingKeyFFT<brP_he3db> &bk_he3db,
    const TFHEpp::SecretKey &sk)
{
    // a < b ⟺ MSB(a − b) = 1
    TFHEpp::TLWE<TargetP> diff{};
    TLWESub<TargetP>(diff, a_ct, b_ct);
    TFHEpp::TLWE<TargetP> msb{};
    HomMSB_MetaPBS<TargetP, DomainP, brP, ksToDomP, weightBkP, ksDownP,
                   iksP_he3db, brP_he3db>(
        msb, diff, plain_bits + 1,
        bkfft, ahk, ksk_to_dom, bkfft_weight, ksk_down, ksk_he3db, bk_he3db);
    return TFHEpp::tlweSymDecrypt<TargetP>(msb, sk);
}

template <class TargetP, class DomainP,
          class brP, class ksToDomP, class weightBkP, class ksDownP,
          class iksP_he3db, class brP_he3db>
bool metapbs_less_equal(
    const TFHEpp::TLWE<TargetP> &a_ct, const TFHEpp::TLWE<TargetP> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP> &bkfft,
    const TFHEpp::AnnihilateKey<TargetP> &ahk,
    const TFHEpp::KeySwitchingKey<ksToDomP> &ksk_to_dom,
    const TFHEpp::BootstrappingKeyFFT<weightBkP> &bkfft_weight,
    const TFHEpp::KeySwitchingKey<ksDownP> &ksk_down,
    const TFHEpp::KeySwitchingKey<iksP_he3db> &ksk_he3db,
    const TFHEpp::BootstrappingKeyFFT<brP_he3db> &bk_he3db,
    const TFHEpp::SecretKey &sk)
{
    // a <= b ⟺ NOT MSB(b − a)
    TFHEpp::TLWE<TargetP> diff{};
    TLWESub<TargetP>(diff, b_ct, a_ct);
    TFHEpp::TLWE<TargetP> msb{};
    HomMSB_MetaPBS<TargetP, DomainP, brP, ksToDomP, weightBkP, ksDownP,
                   iksP_he3db, brP_he3db>(
        msb, diff, plain_bits + 1,
        bkfft, ahk, ksk_to_dom, bkfft_weight, ksk_down, ksk_he3db, bk_he3db);
    return !TFHEpp::tlweSymDecrypt<TargetP>(msb, sk);
}

template <class TargetP, class DomainP,
          class brP, class ksToDomP, class weightBkP, class ksDownP,
          class iksP_he3db, class brP_he3db>
bool metapbs_equal(
    const TFHEpp::TLWE<TargetP> &a_ct, const TFHEpp::TLWE<TargetP> &b_ct,
    uint32_t plain_bits,
    const TFHEpp::BootstrappingKeyFFT<brP> &bkfft,
    const TFHEpp::AnnihilateKey<TargetP> &ahk,
    const TFHEpp::KeySwitchingKey<ksToDomP> &ksk_to_dom,
    const TFHEpp::BootstrappingKeyFFT<weightBkP> &bkfft_weight,
    const TFHEpp::KeySwitchingKey<ksDownP> &ksk_down,
    const TFHEpp::KeySwitchingKey<iksP_he3db> &ksk_he3db,
    const TFHEpp::BootstrappingKeyFFT<brP_he3db> &bk_he3db,
    const TFHEpp::SecretKey &sk)
{
    // a == b ⟺ (a >= b) AND (a <= b)
    return metapbs_greater_equal<TargetP, DomainP, brP, ksToDomP, weightBkP,
                                 ksDownP, iksP_he3db, brP_he3db>(
               a_ct, b_ct, plain_bits, bkfft, ahk, ksk_to_dom,
               bkfft_weight, ksk_down, ksk_he3db, bk_he3db, sk) &&
           metapbs_less_equal<TargetP, DomainP, brP, ksToDomP, weightBkP,
                              ksDownP, iksP_he3db, brP_he3db>(
               a_ct, b_ct, plain_bits, bkfft, ahk, ksk_to_dom,
               bkfft_weight, ksk_down, ksk_he3db, bk_he3db, sk);
}

}  // namespace

int main(int argc, char **argv)
{
    // ─── type aliases ───────────────────────────────────────────────────────
    using brP = TFHEpp::lvlh1param;
    using ksToDomP = TFHEpp::lvl1hparam;
    using weightBkP = TFHEpp::lvlh2param;
    using ksDownP = TFHEpp::lvl21param;

    using iksP_he3db = TFHEpp::lvl10param;
    using brP_he3db  = TFHEpp::lvl01param;

    using DomainP = typename brP::domainP;  // lvlhalf
    using TargetP = typename brP::targetP;  // lvl1

    // ─── CLI ────────────────────────────────────────────────────────────────
    uint32_t trials = 100;
    bool verify = true;

    for (int i = 1; i < argc; i++) {
        const std::string arg(argv[i]);
        if (arg == "--trials" && i + 1 < argc)
            trials = static_cast<uint32_t>(std::stoul(argv[++i]));
        else if (arg == "--no-check")
            verify = false;
    }

    // ─── key generation ─────────────────────────────────────────────────────
    std::cout << "[compare] generating keys..." << std::flush;

    TFHEpp::SecretKey sk;

    auto bkfft_extract = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    TFHEpp::bkfftgen<brP>(*bkfft_extract, sk);

    auto ksk_to_dom = std::make_unique<TFHEpp::KeySwitchingKey<ksToDomP>>();
    TFHEpp::ikskgen<ksToDomP>(*ksk_to_dom, sk);

    auto ahk = std::make_unique<TFHEpp::AnnihilateKey<TargetP>>();
    TFHEpp::annihilatekeygen<TargetP>(*ahk, sk);

    auto bkfft_weight =
        std::make_unique<TFHEpp::BootstrappingKeyFFT<weightBkP>>();
    TFHEpp::bkfftgen<weightBkP>(*bkfft_weight, sk);

    auto ksk_down = std::make_unique<TFHEpp::KeySwitchingKey<ksDownP>>();
    TFHEpp::ikskgen<ksDownP>(*ksk_down, sk);

    auto ksk_he3db = std::make_unique<TFHEpp::KeySwitchingKey<iksP_he3db>>();
    TFHEpp::ikskgen<iksP_he3db>(*ksk_he3db, sk);

    auto bk_he3db = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_he3db>>();
    TFHEpp::bkfftgen<brP_he3db>(*bk_he3db, sk);

    std::cout << " done.\n";

    // ─── experiment ─────────────────────────────────────────────────────────
    std::mt19937 rng(0x434D5031U);
    const std::vector<int> bit_cases = {4, 8, 16};
    const char* op_names[5] = {"  >  ", " >= ", "  <  ", " <= ", " == "};

    for (const int bits : bit_cases) {
        constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
        if (bits + 1 > digits - 1) {
            std::cout << "  " << bits << "-bit: skipped (needs Lvl2 torus)\n";
            continue;
        }

        const uint32_t scale_bits_val = digits - (bits + 1);
        const typename TargetP::T scale =
            static_cast<typename TargetP::T>(1) << scale_bits_val;

        const uint32_t vmax = (1u << (bits - 1)) - 1u;
        std::uniform_int_distribution<uint32_t> dist(0, vmax);

        // [0]=GT, [1]=GE, [2]=LT, [3]=LE, [4]=EQ
        uint32_t he3db_err[5] = {};
        uint32_t mpbs_err[5] = {};
        double he3db_ms[5] = {};
        double mpbs_ms[5] = {};

        std::cout << "\n  ───── " << bits << "-bit comparison (trials=" << trials
                  << ", range=[0," << vmax << "]) ─────\n";
        std::cout << "    Op  | HE3DB err | MetaPBS err | HE3DB ms | MetaPBS ms\n"
                  << "  ------+-----------+-------------+----------+-----------\n";

        for (uint32_t rep = 0; rep < trials; rep++) {
            const uint32_t a = dist(rng);
            const uint32_t b = dist(rng);

            const bool exp_gt = (a > b);
            const bool exp_ge = (a >= b);
            const bool exp_lt = (a < b);
            const bool exp_le = (a <= b);
            const bool exp_eq = (a == b);
            const bool expected[5] = {exp_gt, exp_ge, exp_lt, exp_le, exp_eq};

            TFHEpp::TLWE<TargetP> ca{}, cb{};
            {
                auto pa = static_cast<typename TargetP::T>(a * scale);
                TFHEpp::tlweSymEncrypt<TargetP>(ca, pa, sk);
            }
            {
                auto pb = static_cast<typename TargetP::T>(b * scale);
                TFHEpp::tlweSymEncrypt<TargetP>(cb, pb, sk);
            }

            // ── HE3DB path ──
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = he3db_greater_than<TargetP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                he3db_ms[0] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[0]) he3db_err[0]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = he3db_greater_equal<TargetP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                he3db_ms[1] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[1]) he3db_err[1]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = he3db_less_than<TargetP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                he3db_ms[2] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[2]) he3db_err[2]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = he3db_less_equal<TargetP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                he3db_ms[3] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[3]) he3db_err[3]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = he3db_equal<TargetP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                he3db_ms[4] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[4]) he3db_err[4]++;
            }

            // ── Meta-PBS path ──
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = metapbs_greater_than<TargetP, DomainP, brP, ksToDomP,
                    weightBkP, ksDownP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *bkfft_extract, *ahk, *ksk_to_dom,
                    *bkfft_weight, *ksk_down, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                mpbs_ms[0] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[0]) mpbs_err[0]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = metapbs_greater_equal<TargetP, DomainP, brP, ksToDomP,
                    weightBkP, ksDownP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *bkfft_extract, *ahk, *ksk_to_dom,
                    *bkfft_weight, *ksk_down, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                mpbs_ms[1] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[1]) mpbs_err[1]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = metapbs_less_than<TargetP, DomainP, brP, ksToDomP,
                    weightBkP, ksDownP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *bkfft_extract, *ahk, *ksk_to_dom,
                    *bkfft_weight, *ksk_down, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                mpbs_ms[2] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[2]) mpbs_err[2]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = metapbs_less_equal<TargetP, DomainP, brP, ksToDomP,
                    weightBkP, ksDownP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *bkfft_extract, *ahk, *ksk_to_dom,
                    *bkfft_weight, *ksk_down, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                mpbs_ms[3] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[3]) mpbs_err[3]++;
            }
            {
                auto t0 = std::chrono::steady_clock::now();
                bool got = metapbs_equal<TargetP, DomainP, brP, ksToDomP,
                    weightBkP, ksDownP, iksP_he3db, brP_he3db>(
                    ca, cb, bits, *bkfft_extract, *ahk, *ksk_to_dom,
                    *bkfft_weight, *ksk_down, *ksk_he3db, *bk_he3db, sk);
                auto t1 = std::chrono::steady_clock::now();
                mpbs_ms[4] += std::chrono::duration<double, std::milli>(t1-t0).count();
                if (verify && got != expected[4]) mpbs_err[4]++;
            }
        }

        // Print results
        for (int op = 0; op < 5; op++) {
            std::cout << "  " << op_names[op] << " | ";
            if (verify)
                std::cout << std::setw(4) << he3db_err[op] << "/"
                          << std::setw(3) << trials;
            else
                std::cout << "   n/a   ";
            std::cout << "  |  ";
            if (verify)
                std::cout << std::setw(4) << mpbs_err[op] << "/"
                          << std::setw(3) << trials;
            else
                std::cout << "    n/a    ";
            std::cout << "  | " << std::fixed << std::setprecision(1)
                      << std::setw(6) << he3db_ms[op] / trials
                      << "ms | " << std::setw(6)
                      << mpbs_ms[op] / trials << "ms\n";
        }
    }

    std::cout << "\n[compare] done.\n";
    return 0;
}
