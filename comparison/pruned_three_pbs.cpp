#include "pruned_three_pbs.h"

namespace tfhepp_compare::three_pbs
{
    // ══════════════════════════════════════════════════════════════
    //  Pruned blind-rotate primitives
    //
    //  Per Chapter 3 §subsec:periodic_pruned_pbs, the BitExtract step uses a
    //  Q/2 LUT with period M_κ = 2·w_κ. CMUX rounds with ā_i ≡ 0 (mod M_κ)
    //  leave the accumulator unchanged and can therefore be skipped.
    // ══════════════════════════════════════════════════════════════

    template <class P>
    static void PrunedBlindRotate(
        TFHEpp::TRLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t period)
    {
        TFHEpp::ModswitchTLWE<typename P::domainP> moded;
        TFHEpp::BRModSwitch<P, 1>(moded, tlwe);
        res = {};
        TFHEpp::PolynomialMulByXai<typename P::targetP>(
            res[P::targetP::k], testvector,
            moded[P::domainP::k * P::domainP::n]);
        for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
            const uint32_t a_bar = moded[i];
            if (a_bar == 0) continue;
            if (period > 1 && (a_bar % period == 0)) continue;
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(
                res, bkfft[i], a_bar);
        }
    }

    template <class P>
    static void PrunedGateBootstrappingTLWE2TLWEFFT(
        TFHEpp::TLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t period)
    {
        alignas(64) TFHEpp::TRLWE<typename P::targetP> acc;
        PrunedBlindRotate<P>(acc, tlwe, bkfft, testvector, period);
        TFHEpp::SampleExtractIndex<typename P::targetP>(res, acc, 0);
    }

    // Q/2 step LUT (boolean encoding from Eq. (eq:boolean_encoding)).
    // Coefficient j is Q/2 iff (j mod period) ≥ period/2.
    template <class P>
    static TFHEpp::Polynomial<P> QHalfBitPolygen(uint32_t period)
    {
        TFHEpp::Polynomial<P> poly = {};
        const uint32_t        half_period = period / 2;
        const typename P::T   q_half =
            typename P::T(1) << (std::numeric_limits<typename P::T>::digits - 1);
        for (uint32_t i = 0; i < P::n; i++)
            poly[i] = ((i % period) >= half_period) ? q_half : 0;
        return poly;
    }

    // Pick a power-of-two period inside [2, 2N] roughly aligned with M_κ.
    template <class P>
    static uint32_t NormalizePeriod(uint32_t period)
    {
        if (period < 2) period = 32;
        // floor to power of two
        uint32_t pow2 = 1;
        while ((pow2 << 1) != 0 && (pow2 << 1) <= period) pow2 <<= 1;
        period = pow2;
        const uint32_t max_period = 2 * P::n;
        if (period > max_period) {
            uint32_t cap = 1;
            while ((cap << 1) != 0 && (cap << 1) <= max_period) cap <<= 1;
            period = cap;
        }
        if (period < 2) period = 2;
        return period;
    }

    // ── Step 1: BitExtract — pruned PBS that emits boolean Enc(b_κ ∈ {0,1})
    //   in {0, Q/2} encoding.                                              ──
    static void BitExtract_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                const TFHEEvalKey &ek, uint32_t period)
    {
        constexpr Lvl1::T offset =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, QHalfBitPolygen<Lvl1>(period),
            period);
    }

    static void BitExtract_Lvl2(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                const TFHEEvalKey &ek, uint32_t period)
    {
        constexpr Lvl2::T offset =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - 7);
        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, QHalfBitPolygen<Lvl2>(period),
            period);
    }

    // ── Step 2: B2A weight PBS — convert boolean Enc(b_κ) → arithmetic
    //   Enc(b_κ · weight) using the standard ETHMSB ±weight LUT.        ──
    static void Boolean2Weight_Lvl1(TLWELvl1 &res, const TLWELvl1 &boolean_bit,
                                    Lvl1::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl1::T offset =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = boolean_bit;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(weight));
        res[Lvl1::k * Lvl1::n] += weight;
    }

    static void Boolean2Weight_Lvl2(TLWELvl2 &res, const TLWELvl2 &boolean_bit,
                                    Lvl2::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl2::T offset =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - 7);
        TLWELvl2 tlweoffset = boolean_bit;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(weight));
        res[Lvl2::k * Lvl2::n] += weight;
    }

    // Default LUT period: matches the κ=5 reliable window.
    static constexpr uint32_t kDefaultPeriod = 4;

    // ══════════════════════════════════════════════════════════════
    //  One pruned-3-PBS recursion level
    // ══════════════════════════════════════════════════════════════

    static void three_pbs_lvl1_recursive(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                         uint32_t plain_bits,
                                         const TFHEEvalKey &ek,
                                         bool result_type)
    {
        constexpr uint32_t kappa = 5;
        if (plain_bits <= kappa) {
            MSBGateBootstrapping(res, tlwe, ek, result_type);
            return;
        }

        // Step 1: align b_κ to MSB of κ-bit window.
        TLWELvl1 shift_tlwe;
        for (size_t i = 0; i <= Lvl1::n; i++) shift_tlwe[i] = tlwe[i] << kappa;

        // Step 2 (PBS #1): pruned BitExtract → boolean Enc(b_κ).
        const uint32_t period = NormalizePeriod<Lvl1>(kDefaultPeriod);
        TLWELvl1       boolean_bit;
        BitExtract_Lvl1(boolean_bit, shift_tlwe, ek, period);

        // Step 3 (PBS #2): convert boolean → weighted ciphertext.
        constexpr Lvl1::T weight =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - kappa - 2);
        TLWELvl1 weight_bit;
        Boolean2Weight_Lvl1(weight_bit, boolean_bit, weight, ek);

        // Step 4: subtract from original to clear b_κ.
        TLWELvl1 guarded;
        for (size_t i = 0; i <= Lvl1::n; i++)
            guarded[i] = tlwe[i] - weight_bit[i];

        three_pbs_lvl1_recursive(res, guarded, plain_bits - kappa, ek,
                                 result_type);
    }

    static void three_pbs_lvl2_one_level(TLWELvl2 &out, const TLWELvl2 &tlwe,
                                         const TFHEEvalKey &ek)
    {
        constexpr uint32_t kappa = 5;
        TLWELvl2 shift_tlwe;
        for (size_t i = 0; i <= Lvl2::n; i++) shift_tlwe[i] = tlwe[i] << kappa;

        const uint32_t period = NormalizePeriod<Lvl2>(kDefaultPeriod);
        TLWELvl2       boolean_bit;
        BitExtract_Lvl2(boolean_bit, shift_tlwe, ek, period);

        constexpr Lvl2::T weight =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - kappa - 2);
        TLWELvl2 weight_bit;
        Boolean2Weight_Lvl2(weight_bit, boolean_bit, weight, ek);

        for (size_t i = 0; i <= Lvl2::n; i++)
            out[i] = tlwe[i] - weight_bit[i];
    }

    static void three_pbs_lvl2_dispatch(TLWELvl1 &res, const TLWELvl2 &tlwe,
                                        uint32_t plain_bits,
                                        const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 5) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            MSBGateBootstrapping(res, res, ek, result_type);
            return;
        }
        if (plain_bits <= 9) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            three_pbs_lvl1_recursive(res, res, plain_bits, ek, result_type);
            return;
        }
        TLWELvl2 guarded;
        three_pbs_lvl2_one_level(guarded, tlwe, ek);
        three_pbs_lvl2_dispatch(res, guarded, plain_bits - 5, ek, result_type);
    }

    // ══════════════════════════════════════════════════════════════
    //  Public API
    // ══════════════════════════════════════════════════════════════

    void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     const TFHEEvalKey &ek, bool result_type)
    {
        MSBGateBootstrapping(res, tlwe, ek, result_type);
    }

    void ExtractMSB10(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                      const TFHEEvalKey &ek, bool result_type)
    {
        three_pbs_lvl1_recursive(res, tlwe, plain_bits, ek, result_type);
    }

    void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type)
    {
        (void) plain_bits;
        three_pbs_lvl2_dispatch(res, tlwe, 5, ek, result_type);
    }

    void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type);
    }

    void ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type);
    }

    void ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type);
    }

    void ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type);
    }

    void ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type);
    }

    void ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, result_type);
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 5) ExtractMSB5(res, tlwe, ek, result_type);
        else if (plain_bits <= 10)
            ExtractMSB10(res, tlwe, plain_bits, ek, result_type);
        else
            throw std::invalid_argument(
                "Pruned3PBS: Lvl1 plain_bits out of range (max 10).");
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 5)
            ImExtractMSB5(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 9)
            ImExtractMSB9(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 14)
            ImExtractMSB14(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 19)
            ImExtractMSB19(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 24)
            ImExtractMSB24(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 29)
            ImExtractMSB29(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 33)
            ImExtractMSB33(res, tlwe, plain_bits, ek, result_type);
        else
            throw std::invalid_argument(
                "Pruned3PBS: Lvl2 plain_bits out of range (max 33).");
    }

} // namespace tfhepp_compare::three_pbs
