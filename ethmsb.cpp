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

    // ══════════════════════════════════════════════════════════════
    //  Guard-bit extraction PBS (ETHMSB's key innovation)
    //
    //  Instead of MSB-BS + Identity-BS (2 PBS), we do a single PBS
    //  with LUT coefficients ±weight and the ETHMSB offset
    //  offset = ((2^{κ-1} + 1) / 2) * Δ
    //
    //  For Lvl1: κ=5, offset = (17/2) * (Q/2^5) = 17 * Q/64
    //  For Lvl2: κ=6, offset = (33/2) * (Q/2^6) = 33 * Q/128
    // ══════════════════════════════════════════════════════════════

    void GuardBitExtractBS_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                Lvl1::T weight, const TFHEEvalKey &ek)
    {
        // offset = ((2^4 + 1) / 2) * Δ where Δ = Q/2^5
        // = 17 * Q / 64 = 17 * (1 << (32-6))
        // But more precisely: (2^{κ-1}+1) * Q / 2^{κ+1}
        // κ=5: (16+1) * 2^32 / 2^6 = 17 * 2^26
        constexpr Lvl1::T offset = 17U * (1U << (std::numeric_limits<Lvl1::T>::digits - 6));

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
        // κ=6 for Lvl2: offset = (2^5+1) * Q / 2^7 = 33 * 2^(64-7)
        // But we use κ=5 recursion step (peel 5 bits), matching HE3DB
        // HE3DB uses offset = 1 << (digits - 7) for Lvl2 MSB
        // ETHMSB offset for guard bit: (2^4+1) * (Q/2^6) = 17 * 2^(64-6) = 17 * 2^58
        // Actually for Lvl2, κ_effective for the window = 6 (since lvl2 has N=2048, nbit=11)
        // So offset = (2^5+1)/2 * (Q/2^6) = 33/2 * 2^58 ≈ 16.5 * 2^58
        // Since we're using integer arithmetic: (33 * (1ULL << (64-7)))
        constexpr Lvl2::T offset = 17ULL * (1ULL << (std::numeric_limits<Lvl2::T>::digits - 6));

        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);

        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(weight));
        res[Lvl2::k * Lvl2::n] += weight;
    }

    // ══════════════════════════════════════════════════════════════
    //  Pruned Blind Rotation (Chapter 3: periodic CMUX pruning)
    //
    //  Identical to TFHEpp::BlindRotate but adds one extra check:
    //  if (ā % Mk == 0) continue;
    //  This skips CMUX steps where the periodic LUT is invariant
    //  under the rotation, i.e., Rot_{ā}(v) = v.
    // ══════════════════════════════════════════════════════════════

    template <class P>
    static void PrunedBlindRotate(
        TFHEpp::TRLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t Mk)
    {
        constexpr uint32_t bitwidth = TFHEpp::bits_needed<0>();
        const uint32_t b_bar = 2 * P::targetP::n -
            ((tlwe[P::domainP::k * P::domainP::n] >>
              (std::numeric_limits<typename P::domainP::T>::digits -
               1 - P::targetP::nbit + bitwidth))
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
                (std::numeric_limits<typename P::domainP::T>::digits -
                 1 - P::targetP::nbit + bitwidth)
                    << bitwidth;
            if (a_bar == 0) continue;                        // standard: skip zero
            if (Mk > 1 && (a_bar % Mk == 0)) continue;      // periodic pruning
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(
                res, bkfft[i], a_bar);
        }
    }

    // Pruned GateBootstrapping: PrunedBlindRotate + SampleExtract
    template <class P>
    static void PrunedGateBootstrappingTLWE2TLWEFFT(
        TFHEpp::TLWE<typename P::targetP> &res,
        const TFHEpp::TLWE<typename P::domainP> &tlwe,
        const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
        const TFHEpp::Polynomial<typename P::targetP> &testvector,
        uint32_t Mk)
    {
        alignas(64) TFHEpp::TRLWE<typename P::targetP> acc;
        PrunedBlindRotate<P>(acc, tlwe, bkfft, testvector, Mk);
        TFHEpp::SampleExtractIndex<typename P::targetP>(res, acc, 0);
    }

    // ══════════════════════════════════════════════════════════════
    //  Pruned Guard-bit extraction PBS
    // ══════════════════════════════════════════════════════════════

    void PrunedGuardBitExtractBS_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                      Lvl1::T weight, uint32_t Mk,
                                      const TFHEEvalKey &ek)
    {
        constexpr Lvl1::T offset =
            17U * (1U << (std::numeric_limits<Lvl1::T>::digits - 6));

        TLWELvl1 tlweoffset = tlwe;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);

        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(weight), Mk);
        res[Lvl1::k * Lvl1::n] += weight;
    }

    void PrunedGuardBitExtractBS_Lvl2(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                      Lvl2::T weight, uint32_t Mk,
                                      const TFHEEvalKey &ek)
    {
        constexpr Lvl2::T offset =
            17ULL * (1ULL << (std::numeric_limits<Lvl2::T>::digits - 6));

        TLWELvl2 tlweoffset = tlwe;
        tlweoffset[Lvl2::k * Lvl2::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl20>(tlwelvl0, tlweoffset, *ek.iksklvl20);

        PrunedGateBootstrappingTLWE2TLWEFFT<Lvl02>(
            res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<Lvl2>(weight), Mk);
        res[Lvl2::k * Lvl2::n] += weight;
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

    // Recursive: 5 < plain_bits ≤ 10
    // Uses PRUNED guard-bit extraction (Chapter 3 optimization)
    void ETHMSB_ExtractMSB_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                uint32_t plain_bits, const TFHEEvalKey &ek,
                                bool result_type)
    {
        if (plain_bits <= 5) {
            ETHMSB_ExtractMSB5(res, tlwe, ek, result_type);
            return;
        }

        constexpr uint32_t kappa = 5;
        uint32_t shift = plain_bits - kappa;

        // Step 1: Left-shift to align b_κ to the MSB of a κ-bit window
        TLWELvl1 shift_tlwe;
        for (size_t i = 0; i <= Lvl1::n; i++)
            shift_tlwe[i] = tlwe[i] << shift;

        // Step 2: Pruned guard-bit extraction PBS
        // weight = 2^{shift-1}, Mk = 2^shift (LUT period after shift)
        Lvl1::T weight = 1U << (shift - 1);
        uint32_t Mk = 1U << shift;
        TLWELvl1 guard_bit;
        PrunedGuardBitExtractBS_Lvl1(guard_bit, shift_tlwe, weight, Mk, ek);

        // Step 3: Subtract to zero the guard bit
        TLWELvl1 guarded;
        for (size_t i = 0; i <= Lvl1::n; i++)
            guarded[i] = tlwe[i] - guard_bit[i];

        // Step 4: Recurse on the guarded ciphertext
        uint32_t remaining = plain_bits - kappa;
        ETHMSB_ExtractMSB_Lvl1(res, guarded, remaining, ek, result_type);
    }

    // ══════════════════════════════════════════════════════════════
    //  ETHMSB: Lvl2 → Lvl1 path (mirrors HE3DB's ImExtractMSB*)
    // ══════════════════════════════════════════════════════════════

    void ETHMSB_ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type)
    {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        ETHMSB_ExtractMSB5(res, res, ek, result_type);
    }

    void ETHMSB_ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type)
    {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res, tlwe, *ek.iksklvl21);
        ETHMSB_ExtractMSB_Lvl1(res, res, plain_bits, ek, result_type);
    }

    // Helper: one ETHMSB guard-bit level on Lvl2 with PRUNING
    static void ethmsb_lvl2_one_level(TLWELvl2 &res_lvl2, const TLWELvl2 &tlwe,
                                      uint32_t plain_bits, const TFHEEvalKey &ek)
    {
        constexpr uint32_t kappa = 5;
        uint32_t shift = plain_bits - kappa;

        // Step 1: Left-shift
        TLWELvl2 shift_tlwe;
        for (size_t i = 0; i <= Lvl2::n; i++)
            shift_tlwe[i] = tlwe[i] << shift;

        // Step 2: Pruned guard-bit extraction
        Lvl2::T weight = 1ULL << (shift - 1);
        uint32_t Mk = 1U << shift;  // LUT period
        TLWELvl2 guard_bit;
        PrunedGuardBitExtractBS_Lvl2(guard_bit, shift_tlwe, weight, Mk, ek);

        // Step 3: Subtract to zero guard bit
        for (size_t i = 0; i <= Lvl2::n; i++)
            res_lvl2[i] = tlwe[i] - guard_bit[i];
    }

    void ETHMSB_ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 guarded;
        ethmsb_lvl2_one_level(guarded, tlwe, plain_bits, ek);
        ETHMSB_ImExtractMSB9(res, guarded, plain_bits - 5, ek, result_type);
    }

    void ETHMSB_ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 guarded;
        ethmsb_lvl2_one_level(guarded, tlwe, plain_bits, ek);
        ETHMSB_ImExtractMSB14(res, guarded, plain_bits - 5, ek, result_type);
    }

    void ETHMSB_ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 guarded;
        ethmsb_lvl2_one_level(guarded, tlwe, plain_bits, ek);
        ETHMSB_ImExtractMSB19(res, guarded, plain_bits - 5, ek, result_type);
    }

    void ETHMSB_ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 guarded;
        ethmsb_lvl2_one_level(guarded, tlwe, plain_bits, ek);
        ETHMSB_ImExtractMSB24(res, guarded, plain_bits - 5, ek, result_type);
    }

    void ETHMSB_ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type)
    {
        TLWELvl2 guarded;
        ethmsb_lvl2_one_level(guarded, tlwe, plain_bits, ek);
        ETHMSB_ImExtractMSB29(res, guarded, plain_bits - 5, ek, result_type);
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
