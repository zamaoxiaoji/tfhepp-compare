#include "pruned_three_pbs.h"

#include <algorithm>

#include "micro_pbs.h"

namespace tfhepp_compare::three_pbs
{
    FastB2AEvalKeyPack GenerateFastB2AEvalKeyPack(const TFHESecretKey &sk,
                                                  bool with_lvl2)
    {
        FastB2AEvalKeyPack pack;
        pack.lvl1 =
            micro_pbs::GenerateMicroEvalKeyPack<FastB2ALvl1Candidate>(sk,
                                                                      false);
        if (with_lvl2) {
            pack.lvl2 = std::make_unique<FastB2ALvl2EvalKeyPack>(
                micro_pbs::GenerateMicroEvalKeyPack<FastB2ALvl2Candidate>(
                    sk, false));
        }
        return pack;
    }

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
        res = {};
        const uint32_t br_index =
            2 * P::targetP::n -
            (tlwe[P::domainP::k * P::domainP::n] >>
             (std::numeric_limits<typename P::domainP::T>::digits - 1 -
              P::targetP::nbit));
        TFHEpp::PolynomialMulByXai<typename P::targetP>(
            res[P::targetP::k], testvector, br_index);
        for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
            constexpr typename P::domainP::T roundoffset =
                typename P::domainP::T(1)
                << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                    P::targetP::nbit);
            const uint32_t a_bar =
                (tlwe[i] + roundoffset) >>
                (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                 P::targetP::nbit);
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

    template <class P>
    static typename P::T PhaseDelta(uint32_t plain_bits)
    {
        constexpr uint32_t digits = std::numeric_limits<typename P::T>::digits;
        return plain_bits >= digits
                   ? typename P::T(1)
                   : (typename P::T(1) << (digits - plain_bits));
    }

    template <class P>
    static typename P::T CenterOffset(uint32_t plain_bits)
    {
        return PhaseDelta<P>(plain_bits) >> 1;
    }

    template <class P>
    static uint32_t WindowLocalPeriod(uint32_t window_local_k)
    {
        return uint32_t(1) << (P::nbit + 1 - window_local_k);
    }

    template <class P>
    static uint32_t BitExtractPreScaleShift(uint32_t plain_bits,
                                            uint32_t window_local_k)
    {
        const uint32_t desired =
            plain_bits > window_local_k + 1
                ? plain_bits - (window_local_k + 1)
                : 0;
        const uint32_t base_period = WindowLocalPeriod<P>(window_local_k);
        // A {0,Q/2} LUT is invariant under negation, so a period-2N step
        // cannot be represented by the negacyclic polynomial: the implicit
        // [N,2N) half is just the negation of [0,N). Keep the expanded
        // BitExtract plateau within N slots.
        const uint32_t max_period = P::n;
        uint32_t max_shift = 0;
        while (max_shift < desired &&
               base_period <= (max_period >> (max_shift + 1)))
            max_shift++;
        return std::min(desired, max_shift);
    }

    template <class P>
    static uint32_t WindowLocalPeriod(uint32_t window_local_k,
                                      uint32_t prescale_shift)
    {
        return WindowLocalPeriod<P>(window_local_k) << prescale_shift;
    }

    static Lvl1::T GapOffsetLvl1(uint32_t guarded_plain_bits)
    {
        constexpr uint32_t kappa = 5;
        const uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        if (guarded_plain_bits <= kappa + 1 || guarded_plain_bits >= digits)
            return Lvl1::T(1) << (digits - 6);
        return (Lvl1::T(1) << (digits - kappa - 2)) +
               (Lvl1::T(1) << (digits - guarded_plain_bits - 1));
    }

    static void GapMSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                        uint32_t guarded_plain_bits,
                                        const TFHEEvalKey &ek,
                                        bool result_type)
    {
        Lvl1::T mu = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) mu = mu << 1;
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] +=
            GapOffsetLvl1(guarded_plain_bits);

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(mu));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += mu;
    }

    // ── Step 1: BitExtract — pruned PBS that emits boolean Enc(b_κ ∈ {0,1})
    //   in {0, Q/2} encoding.                                              ──
    static void BitExtract_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                const TFHEEvalKey &ek, uint32_t plain_bits,
                                uint32_t window_local_k)
    {
        const uint32_t prescale_shift =
            BitExtractPreScaleShift<Lvl1>(plain_bits, window_local_k);
        const uint32_t effective_plain_bits = plain_bits - prescale_shift;
        const uint32_t period =
            WindowLocalPeriod<Lvl1>(window_local_k, prescale_shift);
        TLWELvl1 tlweoffset = tlwe;
        if (prescale_shift != 0) {
            for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
                tlweoffset[i] <<= prescale_shift;
        }
        tlweoffset[Lvl1::k * Lvl1::n] +=
            CenterOffset<Lvl1>(effective_plain_bits);

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, QHalfBitPolygen<Lvl1>(period),
            period);
    }

    static void BitExtract_Lvl2(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                const TFHEEvalKey &ek, uint32_t plain_bits,
                                uint32_t window_local_k)
    {
        const uint32_t prescale_shift =
            BitExtractPreScaleShift<Lvl2>(plain_bits, window_local_k);
        const uint32_t effective_plain_bits = plain_bits - prescale_shift;
        const uint32_t period =
            WindowLocalPeriod<Lvl2>(window_local_k, prescale_shift);
        TLWELvl2 tlweoffset = tlwe;
        if (prescale_shift != 0) {
            for (size_t i = 0; i <= Lvl2::k * Lvl2::n; i++)
                tlweoffset[i] <<= prescale_shift;
        }
        tlweoffset[Lvl2::k * Lvl2::n] +=
            CenterOffset<Lvl2>(effective_plain_bits);

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

    static void Boolean2Weight_Lvl1_Fast(
        TLWELvl1 &res, const TLWELvl1 &boolean_bit, Lvl1::T A,
        const TFHEEvalKey &ek, const FastB2AEvalKeyPack *micro_pack)
    {
        const Lvl1::T guard_value = A << 1;
        if (micro_pack != nullptr) {
            const micro_pbs::MicroUnsafeBestPlan plan =
                micro_pbs::SelectMicroUnsafeBestPlan(5);
            if (plan.supported) {
                micro_pbs::MicroUnsafeKsPreQHalfToGuardValue_Lvl1<
                    FastB2ALvl1Candidate>(res, boolean_bit, guard_value,
                                          plan.offset_value, micro_pack->lvl1);
                return;
            }
        }
        // Compatibility fallback keeps the old API usable when the caller has
        // not generated the explicit unsafe micro key pack.
        Boolean2Weight_Lvl1(res, boolean_bit, A, ek);
    }

    static void Boolean2Weight_Lvl2_Fast(
        TLWELvl2 &res, const TLWELvl2 &boolean_bit, Lvl2::T A,
        const TFHEEvalKey &ek, const FastB2AEvalKeyPack *micro_pack)
    {
        if (micro_pack != nullptr && micro_pack->lvl2 != nullptr) {
            constexpr uint32_t digits = std::numeric_limits<Lvl2::T>::digits;
            constexpr Lvl2::T offset = Lvl2::T(1) << (digits - 5);
            const Lvl2::T guard_value = A << 1;
            micro_pbs::MicroUnsafeKsPreQHalfToGuardValue_Lvl2<
                FastB2ALvl2Candidate>(res, boolean_bit, guard_value, offset,
                                      *micro_pack->lvl2);
            return;
        }
        Boolean2Weight_Lvl2(res, boolean_bit, A, ek);
    }

    // ══════════════════════════════════════════════════════════════
    //  One pruned-3-PBS recursion level
    // ══════════════════════════════════════════════════════════════

    static void three_pbs_lvl1_recursive(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                         uint32_t plain_bits,
                                         const TFHEEvalKey &ek,
                                         const FastB2AEvalKeyPack *micro_pack,
                                         bool result_type,
                                         uint32_t gap_parent_bits)
    {
        constexpr uint32_t kappa = 5;
        const bool gap_allows_early_final =
            gap_parent_bits != 0 && plain_bits <= kappa + 3;
        if (plain_bits <= kappa || gap_allows_early_final) {
            if (gap_parent_bits == 0)
                MSBGateBootstrapping(res, tlwe, plain_bits, ek, result_type);
            else
                GapMSBGateBootstrapping(res, tlwe, gap_parent_bits, ek,
                                        result_type);
            return;
        }

        // Step 1 (PBS #1): window-local BitExtract → boolean Enc(b_κ).
        // kappa is a window-local MSB-first bit index here, not a global
        // plaintext bit index. The Q/2 LUT period and the BR representative
        // offset are selected from the current recursion window.
        TLWELvl1 boolean_bit;
        BitExtract_Lvl1(boolean_bit, tlwe, ek, plain_bits, kappa);

        // Step 2 (PBS #2): convert boolean → weighted ciphertext.
        constexpr Lvl1::T A =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - kappa - 2);
        TLWELvl1 weight_bit;
        if (micro_pack != nullptr)
            Boolean2Weight_Lvl1_Fast(weight_bit, boolean_bit, A, ek,
                                     micro_pack);
        else
            Boolean2Weight_Lvl1(weight_bit, boolean_bit, A, ek);

        // Step 3: subtract from original to clear b_κ.
        TLWELvl1 guarded;
        for (size_t i = 0; i <= Lvl1::n; i++)
            guarded[i] = tlwe[i] - weight_bit[i];

        three_pbs_lvl1_recursive(res, guarded, plain_bits - kappa, ek,
                                 micro_pack, result_type, plain_bits);
    }

    static void three_pbs_lvl2_one_level(TLWELvl2 &out, const TLWELvl2 &tlwe,
                                         uint32_t plain_bits,
                                         const TFHEEvalKey &ek,
                                         const FastB2AEvalKeyPack *micro_pack)
    {
        constexpr uint32_t kappa = 5;
        TLWELvl2 boolean_bit;
        BitExtract_Lvl2(boolean_bit, tlwe, ek, plain_bits, kappa);

        constexpr Lvl2::T A =
            Lvl2::T(1) << (std::numeric_limits<Lvl2::T>::digits - kappa - 2);
        TLWELvl2 weight_bit;
        Boolean2Weight_Lvl2_Fast(weight_bit, boolean_bit, A, ek, micro_pack);

        for (size_t i = 0; i <= Lvl2::n; i++)
            out[i] = tlwe[i] - weight_bit[i];
    }

    static void three_pbs_lvl2_dispatch(TLWELvl1 &res, const TLWELvl2 &tlwe,
                                        uint32_t plain_bits,
                                        const TFHEEvalKey &ek,
                                        const FastB2AEvalKeyPack *micro_pack,
                                        bool result_type,
                                        uint32_t gap_parent_bits)
    {
        constexpr uint32_t kappa = 5;
        if (gap_parent_bits != 0 && plain_bits <= kappa + 9) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            three_pbs_lvl1_recursive(res, res, plain_bits, ek, micro_pack,
                                     result_type, gap_parent_bits);
            return;
        }
        if (plain_bits <= 6) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            if (gap_parent_bits == 0)
                MSBGateBootstrapping(res, res, plain_bits, ek, result_type);
            else
                GapMSBGateBootstrapping(res, res, gap_parent_bits, ek,
                                        result_type);
            return;
        }
        if (plain_bits <= 9) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            three_pbs_lvl1_recursive(res, res, plain_bits, ek, micro_pack,
                                     result_type, gap_parent_bits);
            return;
        }
        TLWELvl2 guarded;
        three_pbs_lvl2_one_level(guarded, tlwe, plain_bits, ek, micro_pack);
        three_pbs_lvl2_dispatch(res, guarded, plain_bits - 5, ek, micro_pack,
                                result_type, plain_bits);
    }

    // ══════════════════════════════════════════════════════════════
    //  Public API
    // ══════════════════════════════════════════════════════════════

    void ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                     const TFHEEvalKey &ek, bool result_type)
    {
        MSBGateBootstrapping(res, tlwe, 5, ek, result_type);
    }

    void ExtractMSB10(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                      const TFHEEvalKey &ek, bool result_type)
    {
        three_pbs_lvl1_recursive(res, tlwe, plain_bits, ek, nullptr,
                                 result_type, 0);
    }

    void ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type)
    {
        (void) plain_bits;
        three_pbs_lvl2_dispatch(res, tlwe, 5, ek, nullptr, result_type, 0);
    }

    void ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                       const TFHEEvalKey &ek, bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, nullptr, result_type,
                                0);
    }

    void ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, nullptr, result_type,
                                0);
    }

    void ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, nullptr, result_type,
                                0);
    }

    void ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, nullptr, result_type,
                                0);
    }

    void ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, nullptr, result_type,
                                0);
    }

    void ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        bool result_type)
    {
        three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, nullptr, result_type,
                                0);
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 6)
            MSBGateBootstrapping(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 10)
            ExtractMSB10(res, tlwe, plain_bits, ek, result_type);
        else
            throw std::invalid_argument(
                "Pruned3PBS: Lvl1 plain_bits out of range (max 10).");
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, bool result_type)
    {
        if (plain_bits <= 6)
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

    void HomMSB(TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, const FastB2AEvalKeyPack &micro_pack,
                bool result_type)
    {
        if (plain_bits <= 6)
            MSBGateBootstrapping(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 10)
            three_pbs_lvl1_recursive(res, tlwe, plain_bits, ek, &micro_pack,
                                     result_type, 0);
        else
            throw std::invalid_argument(
                "Pruned3PBS+MicroPBS: Lvl1 plain_bits out of range (max 10).");
    }

    void HomMSB(TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
                const TFHEEvalKey &ek, const FastB2AEvalKeyPack &micro_pack,
                bool result_type)
    {
        if (plain_bits <= 6) {
            TFHEpp::IdentityKeySwitch<Lvl21>(res, tlwe, *ek.iksklvl21);
            MSBGateBootstrapping(res, res, plain_bits, ek, result_type);
        }
        else if (plain_bits <= 33) {
            three_pbs_lvl2_dispatch(res, tlwe, plain_bits, ek, &micro_pack,
                                    result_type, 0);
        }
        else {
            throw std::invalid_argument(
                "Pruned3PBS+MicroPBS: Lvl2 plain_bits out of range (max 33).");
        }
    }

} // namespace tfhepp_compare::three_pbs
