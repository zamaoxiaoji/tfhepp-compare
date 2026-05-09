#include "ethmsb.h"

namespace ETHMSB_NS
{
    // ══════════════════════════════════════════════════════════════
    //  Low-level PBS wrappers (identical to HE3DB's MSBGateBootstrapping)
    // ══════════════════════════════════════════════════════════════

    void MSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                              const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T μ = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        constexpr uint64_t offset =
            1ULL << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(μ));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += μ;
    }

    void MSBGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                              const TFHEEvalKey &ek, bool result_type)
    {
        Lvl2::T μ = Lvl2::μ;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        constexpr uint64_t offset =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - 7);
        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(μ));
        if (IS_ARITHMETIC(result_type)) res[Lvl2::k * Lvl2::n] += μ;
    }

    static Lvl1::T GapOffsetLvl1(uint32_t guarded_plain_bits)
    {
        constexpr uint32_t kappa = 5;
        const uint32_t digits = std::numeric_limits<Lvl1::T>::digits;
        if (guarded_plain_bits <= kappa + 1 || guarded_plain_bits >= digits)
            return 1U << (digits - 6);
        return (1U << (digits - kappa - 2)) +
               (1U << (digits - guarded_plain_bits - 1));
    }

    static void GapMSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                        uint32_t guarded_plain_bits,
                                        const TFHEEvalKey &ek,
                                        bool result_type)
    {
        Lvl1::T μ = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;

        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += GapOffsetLvl1(guarded_plain_bits);

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(μ));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += μ;
    }

    // ══════════════════════════════════════════════════════════════
    //  Guard-bit extraction PBS (ETHMSB's key innovation)
    //
    //  Instead of MSB-BS + Identity-BS (2 PBS), we do a single PBS
    //  with LUT coefficients ±weight. This PBS uses the same small-window
    //  offset as HE3DB's reliable MSB extraction; the gap-specific offset is
    //  applied only at the final GapMSB decision.
    // ══════════════════════════════════════════════════════════════

    void GuardBitExtractBS_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                Lvl1::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl1::T offset =
            1U << (std::numeric_limits<Lvl1::T>::digits - 6);

        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);

        // LUT with ±weight (ARITHMETIC mode always for guard-bit extraction)
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(weight));
        res[Lvl1::k * Lvl1::n] += weight;
    }

    void GuardBitExtractBS_Lvl2(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                Lvl2::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl2::T offset =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - 7);

        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);

        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(weight));
        res[Lvl2::k * Lvl2::n] += weight;
    }

    static uint32_t NormalizeSkipPeriod(uint32_t period, uint32_t fallback)
    {
        if (period < 2) return fallback;
        return period;
    }

    static uint32_t FloorPowerOfTwo(uint32_t value)
    {
        uint32_t pow2 = 1;
        while ((pow2 << 1) != 0 && (pow2 << 1) <= value) pow2 <<= 1;
        return pow2;
    }

    template <class P>
    static uint32_t NormalizePaperPeriod(uint32_t period)
    {
        if (period < 2) period = 32;
        period = FloorPowerOfTwo(period);
        const uint32_t max_period = 2 * P::n;
        if (period > max_period) period = FloorPowerOfTwo(max_period);
        if (period < 2) period = 2;
        return period;
    }

    static bool PeriodicSkip(uint32_t a_bar, uint32_t period)
    {
        return period > 1 && (a_bar % period == 0);
    }

    static void RecordPruneTerm(PruneStats *stats, uint32_t a_bar,
                                uint32_t period, bool &skip)
    {
        skip = false;
        if (stats != nullptr) stats->total_terms++;
        if (a_bar == 0) {
            if (stats != nullptr) stats->zero_skipped++;
            skip = true;
            return;
        }
        if (PeriodicSkip(a_bar, period)) {
            if (stats != nullptr) stats->periodic_skipped++;
            skip = true;
            return;
        }
        if (stats != nullptr) stats->executed_terms++;
    }

    template <class P, uint32_t num_out = 1>
    static void PrunedBlindRotate(
        TFHEpp::TRLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t period, PruneStats *stats)
    {
#ifdef USE_HE3DB_COMPAT
        constexpr uint32_t bitwidth = TFHEpp::bits_needed<num_out - 1>();
        const uint32_t b_bar =
            2 * P::targetP::n -
            ((tlwe[P::domainP::k * P::domainP::n] >>
              (std::numeric_limits<typename P::domainP::T>::digits - 1 -
               P::targetP::nbit + bitwidth))
             << bitwidth);
        res = {};
        TFHEpp::PolynomialMulByXai<typename P::targetP>(
            res[P::targetP::k], testvector, b_bar);
        for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
            constexpr typename P::domainP::T roundoffset =
                1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                         2 - P::targetP::nbit + bitwidth);
            const uint32_t a_bar =
                (tlwe[i] + roundoffset) >>
                (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                 P::targetP::nbit + bitwidth)
                    << bitwidth;
            bool skip = false;
            RecordPruneTerm(stats, a_bar, period, skip);
            if (skip) continue;
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(
                res, bkfft[i], a_bar);
        }
#else
        TFHEpp::ModswitchTLWE<typename P::domainP> moded;
        TFHEpp::BRModSwitch<P, num_out>(moded, tlwe);
        res = {};
        TFHEpp::PolynomialMulByXai<typename P::targetP>(
            res[P::targetP::k], testvector,
            moded[P::domainP::k * P::domainP::n]);
        for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
            const uint32_t a_bar = moded[i];
            bool skip = false;
            RecordPruneTerm(stats, a_bar, period, skip);
            if (skip) continue;
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(
                res, bkfft[i], a_bar);
        }
#endif
    }

    template <class P>
    static void PrunedGateBootstrappingTLWE2TLWEFFT(
        TFHEpp::TLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t period, PruneStats *stats)
    {
        alignas(64) TFHEpp::TRLWE<typename P::targetP> acc;
        PrunedBlindRotate<P>(acc, tlwe, bkfft, testvector, period, stats);
        TFHEpp::SampleExtractIndex<typename P::targetP>(res, acc, 0);
    }

    template <class P>
    static TFHEpp::Polynomial<P> QHalfBitPolygen(uint32_t period)
    {
        TFHEpp::Polynomial<P> poly = {};
        period = NormalizePaperPeriod<P>(period);
        const uint32_t half_period = period / 2;
        const typename P::T q_half =
            typename P::T(1) << (std::numeric_limits<typename P::T>::digits - 1);
        for (uint32_t i = 0; i < P::n; i++)
            poly[i] = ((i % period) >= half_period) ? q_half : 0;
        return poly;
    }

    template <class P>
    static typename P::T PaperWeightCoeff(uint32_t period)
    {
        period = NormalizePaperPeriod<P>(period);
        uint32_t log_period = 0;
        while ((1U << log_period) < period) log_period++;
        const uint32_t digits = std::numeric_limits<typename P::T>::digits;
        const uint32_t shift = (log_period + 2 >= digits)
                                   ? 0
                                   : digits - log_period - 2;
        return typename P::T(1) << shift;
    }

    static void PrunedWeightedGuardBitExtractBS_Lvl1(
        TLWELvl1 &res, const TLWELvl1 &tlwe, Lvl1::T weight,
        const TFHEEvalKey &ek, uint32_t period, PruneStats *stats)
    {
        constexpr Lvl1::T offset =
            1U << (std::numeric_limits<Lvl1::T>::digits - 6);

        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(weight),
            NormalizeSkipPeriod(period, 4), stats);
        res[Lvl1::k * Lvl1::n] += weight;
    }

    static void PrunedWeightedGuardBitExtractBS_Lvl2(
        TLWELvl2 &res, const TLWELvl2 &tlwe, Lvl2::T weight,
        const TFHEEvalKey &ek, uint32_t period, PruneStats *stats)
    {
        constexpr Lvl2::T offset =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - 7);

        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(weight),
            NormalizeSkipPeriod(period, 4), stats);
        res[Lvl2::k * Lvl2::n] += weight;
    }

    static void PrunedQHalfBitExtractBS_Lvl1(
        TLWELvl1 &res, const TLWELvl1 &tlwe, const TFHEEvalKey &ek,
        uint32_t period, PruneStats *stats)
    {
        constexpr Lvl1::T offset =
            1U << (std::numeric_limits<Lvl1::T>::digits - 6);

        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        const uint32_t paper_period = NormalizePaperPeriod<Lvl1>(period);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01,
            QHalfBitPolygen<Lvl1>(paper_period), paper_period, stats);
    }

    static void PrunedQHalfBitExtractBS_Lvl2(
        TLWELvl2 &res, const TLWELvl2 &tlwe, const TFHEEvalKey &ek,
        uint32_t period, PruneStats *stats)
    {
        constexpr Lvl2::T offset =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - 7);

        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        const uint32_t paper_period = NormalizePaperPeriod<Lvl2>(period);
        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02,
            QHalfBitPolygen<Lvl2>(paper_period), paper_period, stats);
    }

    static void QHalfBitToWeightPBS_Lvl1(TLWELvl1 &res,
                                         const TLWELvl1 &qhalf_bit,
                                         Lvl1::T weight,
                                         const TFHEEvalKey &ek)
    {
        GuardBitExtractBS_Lvl1(res, qhalf_bit, weight, ek);
    }

    static void QHalfBitToWeightPBS_Lvl2(TLWELvl2 &res,
                                         const TLWELvl2 &qhalf_bit,
                                         Lvl2::T weight,
                                         const TFHEEvalKey &ek)
    {
        GuardBitExtractBS_Lvl2(res, qhalf_bit, weight, ek);
    }

    // ══════════════════════════════════════════════════════════════
    //  ETHMSB: Lvl1 path
    // ══════════════════════════════════════════════════════════════

    // Base case: direct PBS-MSB for ≤5 bits (same as HE3DB's ExtractMSB5)
    void ETHMSB_ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                            const TFHEEvalKey &ek, bool result_type)
    {
        MSBGateBootstrapping(res, tlwe, ek, result_type);
    }

    static void ethmsb_extract_msb_lvl1_impl(TLWELvl1 &res,
                                             const TLWELvl1 &tlwe,
                                             uint32_t plain_bits,
                                             const TFHEEvalKey &ek,
                                             bool result_type,
                                             uint32_t gap_parent_bits)
    {
        if (plain_bits <= 5) {
            if (gap_parent_bits == 0)
                ETHMSB_ExtractMSB5(res, tlwe, ek, result_type);
            else
                GapMSBGateBootstrapping(res, tlwe, gap_parent_bits, ek,
                                        result_type);
            return;
        }

        constexpr uint32_t kappa = 5;
        uint32_t shift = kappa;

        // Step 1: Left-shift to align b_κ to the MSB of a κ-bit window
        TLWELvl1 shift_tlwe;
        for (size_t i = 0; i <= Lvl1::n; i++)
            shift_tlwe[i] = tlwe[i] << shift;

        // Step 2: guard-bit extraction PBS.
        // The arithmetic PBS pattern with TV ±w and final +w encrypts 0 or 2w.
        // Algorithm 2 uses plaintext weight 2^{k-κ-1}; under the k-bit scale
        // this is always Q/64 for κ=5, so w=Q/128.
        constexpr Lvl1::T weight =
            1U << (std::numeric_limits<Lvl1::T>::digits - kappa - 2);
        TLWELvl1 guard_bit;
        GuardBitExtractBS_Lvl1(guard_bit, shift_tlwe, weight, ek);

        // Step 3: Subtract to zero the guard bit
        TLWELvl1 guarded;
        for (size_t i = 0; i <= Lvl1::n; i++)
            guarded[i] = tlwe[i] - guard_bit[i];

        // Step 4: Recurse on the guarded ciphertext
        uint32_t remaining = plain_bits - kappa;
        ethmsb_extract_msb_lvl1_impl(res, guarded, remaining, ek, result_type,
                                     plain_bits);
    }

    // Recursive: 5 < plain_bits ≤ 10
    void ETHMSB_ExtractMSB_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                uint32_t plain_bits, const TFHEEvalKey &ek,
                                bool result_type)
    {
        ethmsb_extract_msb_lvl1_impl(res, tlwe, plain_bits, ek, result_type,
                                     0);
    }

    // ══════════════════════════════════════════════════════════════
    //  ETHMSB: Lvl2 → Lvl1 path (mirrors HE3DB's ImExtractMSB*)
    // ══════════════════════════════════════════════════════════════

    static void ethmsb_im_extract_msb5_impl(TLWELvl1 &res,
                                            const TLWELvl2 &tlwe,
                                            uint32_t plain_bits,
                                            const TFHEEvalKey &ek,
                                            bool result_type,
                                            uint32_t gap_parent_bits)
    {
        (void) plain_bits;
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        if (gap_parent_bits == 0)
            ETHMSB_ExtractMSB5(res, res, ek, result_type);
        else
            GapMSBGateBootstrapping(res, res, gap_parent_bits, ek, result_type);
    }

    void ETHMSB_ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type)
    {
        ethmsb_im_extract_msb5_impl(res, tlwe, plain_bits, ek, result_type,
                                    0);
    }

    static void ethmsb_im_extract_msb9_impl(TLWELvl1 &res,
                                            const TLWELvl2 &tlwe,
                                            uint32_t plain_bits,
                                            const TFHEEvalKey &ek,
                                            bool result_type,
                                            uint32_t gap_parent_bits)
    {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        ethmsb_extract_msb_lvl1_impl(res, res, plain_bits, ek, result_type,
                                     gap_parent_bits);
    }

    void ETHMSB_ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type)
    {
        ethmsb_im_extract_msb9_impl(res, tlwe, plain_bits, ek, result_type,
                                    0);
    }

    // Helper: one ETHMSB guard-bit level on Lvl2.
    static void ethmsb_lvl2_one_level(TLWELvl2 &res_lvl2, const TLWELvl2 &tlwe,
                                      uint32_t plain_bits, const TFHEEvalKey &ek)
    {
        (void) plain_bits;
        constexpr uint32_t kappa = 5;
        uint32_t shift = kappa;

        // Step 1: Left-shift
        TLWELvl2 shift_tlwe;
        for (size_t i = 0; i <= Lvl2::n; i++)
            shift_tlwe[i] = tlwe[i] << shift;

        // Step 2: guard-bit extraction PBS. The fixed guard position kappa=5
        // has phase weight Q/64, so the arithmetic PBS coefficient is Q/128.
        constexpr Lvl2::T weight =
            1ULL << (std::numeric_limits<Lvl2::T>::digits - kappa - 2);
        TLWELvl2 guard_bit;
        GuardBitExtractBS_Lvl2(guard_bit, shift_tlwe, weight, ek);

        // Step 3: Subtract to zero guard bit
        for (size_t i = 0; i <= Lvl2::n; i++)
            res_lvl2[i] = tlwe[i] - guard_bit[i];
    }

    static void ethmsb_im_extract_dispatch(TLWELvl1 &res,
                                           const TLWELvl2 &tlwe,
                                           uint32_t plain_bits,
                                           const TFHEEvalKey &ek,
                                           bool result_type,
                                           uint32_t gap_parent_bits)
    {
        if (plain_bits <= 5)
            ethmsb_im_extract_msb5_impl(res, tlwe, plain_bits, ek, result_type,
                                        gap_parent_bits);
        else if (plain_bits <= 9)
            ethmsb_im_extract_msb9_impl(res, tlwe, plain_bits, ek, result_type,
                                        gap_parent_bits);
        else {
            TLWELvl2 guarded;
            ethmsb_lvl2_one_level(guarded, tlwe, plain_bits, ek);
            ethmsb_im_extract_dispatch(res, guarded, plain_bits - 5, ek,
                                       result_type, plain_bits);
        }
    }

    void ETHMSB_ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        ethmsb_im_extract_dispatch(res, tlwe, plain_bits, ek, result_type,
                                   0);
    }

    void ETHMSB_ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        ethmsb_im_extract_dispatch(res, tlwe, plain_bits, ek, result_type,
                                   0);
    }

    void ETHMSB_ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        ethmsb_im_extract_dispatch(res, tlwe, plain_bits, ek, result_type,
                                   0);
    }

    void ETHMSB_ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        ethmsb_im_extract_dispatch(res, tlwe, plain_bits, ek, result_type,
                                   0);
    }

    void ETHMSB_ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        ethmsb_im_extract_dispatch(res, tlwe, plain_bits, ek, result_type,
                                   0);
    }

    // ══════════════════════════════════════════════════════════════
    //  Dispatch (mirrors HE3DB's HomMSB)
    // ══════════════════════════════════════════════════════════════

    void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type)
    {
        if (plain_bits <= 5) ETHMSB_ExtractMSB5(res, tlwe, ek, result_type);
        else if (plain_bits <= 10)
            ETHMSB_ExtractMSB_Lvl1(res, tlwe, plain_bits, ek, result_type);
        else throw std::invalid_argument("Lvl1 plain_bits out of range.");
    }

    void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type)
    {
        if (plain_bits <= 5)
            ETHMSB_ImExtractMSB5(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 9)
            ETHMSB_ImExtractMSB9(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 14)
            ETHMSB_ImExtractMSB14(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 19)
            ETHMSB_ImExtractMSB19(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 24)
            ETHMSB_ImExtractMSB24(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 29)
            ETHMSB_ImExtractMSB29(res, tlwe, plain_bits, ek, result_type);
        else if (plain_bits <= 33)
            ETHMSB_ImExtractMSB33(res, tlwe, plain_bits, ek, result_type);
        else throw std::invalid_argument("Plain bits out of range.");
    }

    // ══════════════════════════════════════════════════════════════
    //  Explicit fast/pruned path
    // ══════════════════════════════════════════════════════════════

    static uint32_t FastWindowBits(const PrunedETHMSBOptions &options)
    {
        (void) options;
        return 5;
    }

    template <class P>
    static typename P::T FixedGuardCoeff(uint32_t kappa)
    {
        const uint32_t digits = std::numeric_limits<typename P::T>::digits;
        return static_cast<typename P::T>(1) << (digits - kappa - 2);
    }

    static void ethmsb_pruned_lvl1_one_level(
        TLWELvl1 &guarded, const TLWELvl1 &tlwe, const TFHEEvalKey &ek,
        const PrunedETHMSBOptions &options, PruneStats *stats)
    {
        const uint32_t kappa = FastWindowBits(options);
        TLWELvl1 shift_tlwe;
        for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
            shift_tlwe[i] = tlwe[i] << kappa;

        TLWELvl1 guard_bit;
        if (options.mode == PrunedETHMSBMode::WeightedApprox) {
            PrunedWeightedGuardBitExtractBS_Lvl1(
                guard_bit, shift_tlwe, FixedGuardCoeff<Lvl1>(kappa), ek,
                options.period_idx, stats);
        }
        else {
            const uint32_t period =
                NormalizePaperPeriod<Lvl1>(options.period_idx);
            TLWELvl1 qhalf_bit;
            PrunedQHalfBitExtractBS_Lvl1(qhalf_bit, shift_tlwe, ek, period,
                                         stats);
            QHalfBitToWeightPBS_Lvl1(guard_bit, qhalf_bit,
                                     PaperWeightCoeff<Lvl1>(period), ek);
        }

        for (size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
            guarded[i] = tlwe[i] - guard_bit[i];
    }

    static void ethmsb_pruned_extract_msb_lvl1_impl(
        TLWELvl1 &res, const TLWELvl1 &tlwe, uint32_t plain_bits,
        const TFHEEvalKey &ek, bool result_type,
        const PrunedETHMSBOptions &options, PruneStats *stats,
        uint32_t gap_parent_bits)
    {
        if (plain_bits <= 5) {
            if (gap_parent_bits == 0)
                ETHMSB_ExtractMSB5(res, tlwe, ek, result_type);
            else
                GapMSBGateBootstrapping(res, tlwe, gap_parent_bits, ek,
                                        result_type);
            return;
        }

        TLWELvl1 guarded;
        ethmsb_pruned_lvl1_one_level(guarded, tlwe, ek, options, stats);
        ethmsb_pruned_extract_msb_lvl1_impl(
            res, guarded, plain_bits - FastWindowBits(options), ek,
            result_type, options, stats, plain_bits);
    }

    static void ethmsb_pruned_im_extract_msb5_impl(
        TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
        const TFHEEvalKey &ek, bool result_type,
        const PrunedETHMSBOptions &options, PruneStats *stats,
        uint32_t gap_parent_bits)
    {
        (void) plain_bits;
        (void) options;
        (void) stats;
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        if (gap_parent_bits == 0)
            ETHMSB_ExtractMSB5(res, res, ek, result_type);
        else
            GapMSBGateBootstrapping(res, res, gap_parent_bits, ek, result_type);
    }

    static void ethmsb_pruned_im_extract_msb9_impl(
        TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
        const TFHEEvalKey &ek, bool result_type,
        const PrunedETHMSBOptions &options, PruneStats *stats,
        uint32_t gap_parent_bits)
    {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        ethmsb_pruned_extract_msb_lvl1_impl(
            res, res, plain_bits, ek, result_type, options, stats,
            gap_parent_bits);
    }

    static void ethmsb_pruned_lvl2_one_level(
        TLWELvl2 &guarded, const TLWELvl2 &tlwe, const TFHEEvalKey &ek,
        const PrunedETHMSBOptions &options, PruneStats *stats)
    {
        const uint32_t kappa = FastWindowBits(options);
        TLWELvl2 shift_tlwe;
        for (size_t i = 0; i <= Lvl2::k * Lvl2::n; i++)
            shift_tlwe[i] = tlwe[i] << kappa;

        TLWELvl2 guard_bit;
        if (options.mode == PrunedETHMSBMode::WeightedApprox) {
            PrunedWeightedGuardBitExtractBS_Lvl2(
                guard_bit, shift_tlwe, FixedGuardCoeff<Lvl2>(kappa), ek,
                options.period_idx, stats);
        }
        else {
            const uint32_t period =
                NormalizePaperPeriod<Lvl2>(options.period_idx);
            TLWELvl2 qhalf_bit;
            PrunedQHalfBitExtractBS_Lvl2(qhalf_bit, shift_tlwe, ek, period,
                                         stats);
            QHalfBitToWeightPBS_Lvl2(guard_bit, qhalf_bit,
                                     PaperWeightCoeff<Lvl2>(period), ek);
        }

        for (size_t i = 0; i <= Lvl2::k * Lvl2::n; i++)
            guarded[i] = tlwe[i] - guard_bit[i];
    }

    static void ethmsb_pruned_im_extract_dispatch(
        TLWELvl1 &res, const TLWELvl2 &tlwe, uint32_t plain_bits,
        const TFHEEvalKey &ek, bool result_type,
        const PrunedETHMSBOptions &options, PruneStats *stats,
        uint32_t gap_parent_bits)
    {
        if (plain_bits <= 5)
            ethmsb_pruned_im_extract_msb5_impl(
                res, tlwe, plain_bits, ek, result_type, options, stats,
                gap_parent_bits);
        else if (plain_bits <= 9)
            ethmsb_pruned_im_extract_msb9_impl(
                res, tlwe, plain_bits, ek, result_type, options, stats,
                gap_parent_bits);
        else {
            TLWELvl2 guarded;
            ethmsb_pruned_lvl2_one_level(guarded, tlwe, ek, options, stats);
            ethmsb_pruned_im_extract_dispatch(
                res, guarded, plain_bits - FastWindowBits(options), ek,
                result_type, options, stats, plain_bits);
        }
    }

    void HomETHMSBPrunedFast(TLWELvl1 &res, const TLWELvl1 &tlwe,
                             uint32_t plain_bits, const TFHEEvalKey &ek,
                             bool result_type,
                             const PrunedETHMSBOptions &options,
                             PruneStats *stats)
    {
        if (plain_bits <= 5) ETHMSB_ExtractMSB5(res, tlwe, ek, result_type);
        else if (plain_bits <= 10)
            ethmsb_pruned_extract_msb_lvl1_impl(
                res, tlwe, plain_bits, ek, result_type, options, stats, 0);
        else throw std::invalid_argument("Lvl1 plain_bits out of range.");
    }

    void HomETHMSBPrunedFast(TLWELvl1 &res, const TLWELvl2 &tlwe,
                             uint32_t plain_bits, const TFHEEvalKey &ek,
                             bool result_type,
                             const PrunedETHMSBOptions &options,
                             PruneStats *stats)
    {
        if (plain_bits <= 5)
            ethmsb_pruned_im_extract_msb5_impl(
                res, tlwe, plain_bits, ek, result_type, options, stats, 0);
        else if (plain_bits <= 9)
            ethmsb_pruned_im_extract_msb9_impl(
                res, tlwe, plain_bits, ek, result_type, options, stats, 0);
        else if (plain_bits <= 33)
            ethmsb_pruned_im_extract_dispatch(
                res, tlwe, plain_bits, ek, result_type, options, stats, 0);
        else throw std::invalid_argument("Plain bits out of range.");
    }

    // ══════════════════════════════════════════════════════════════
    //  ARI ↔ LOG conversion & gate operators (same as HE3DB)
    // ══════════════════════════════════════════════════════════════

    void ARI_to_LOG(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    const TFHEEvalKey &ek)
    {
        Lvl1::T μ = Lvl1::μ;
        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += μ;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(μ));
    }

    void LOG_to_ARI(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    const TFHEEvalKey &ek)
    {
        Lvl1::T μ = Lvl1::μ;
        μ = μ << 1;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlwe, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-μ));
        res[Lvl1::k * Lvl1::n] += μ;
    }

    void HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
                const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T offset = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = ca[i] + cb[i];
        res[Lvl1::k * Lvl1::n] -= Lvl1::μ >> 1;
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-offset));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
    }

    void HomOR(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
               const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T offset = Lvl1::μ;
        if (IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = ca[i] + cb[i];
        res[Lvl1::k * Lvl1::n] += (Lvl1::μ >> 1);
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(-offset));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
    }

} // namespace ETHMSB_NS
