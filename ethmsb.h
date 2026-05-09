#pragma once
/**
 * @file ethmsb.h
 * @brief Error-Truncating Homomorphic MSB Extraction (ETHMSB)
 *
 * Implements the ETHMSB algorithm following the recursive architecture of
 * HE3DB's HomMSB, but replacing the two-PBS-per-level pattern
 * (MSB-BS + Identity-BS) with a single guard-bit extraction PBS per level.
 *
 * Key idea: At each recursion level, left-shift to align the boundary bit,
 * extract it at its original weight via a coefficient-weighted PBS,
 * subtract to zero the guard bit, then recurse on the guarded ciphertext.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"
#include "gatebootstrapping.hpp"

// Re-use HE3DB type conventions
#define LOGIC true
#define ARITHMETIC false
#define IS_LOGIC(a) a
#define IS_ARITHMETIC(a) !a

namespace ETHMSB_NS
{
    // ── Type aliases (matching HE3DB conventions) ──
    using Lvl0  = TFHEpp::lvl0param;
    using Lvl1  = TFHEpp::lvl1param;
    using Lvl2  = TFHEpp::lvl2param;
    using Lvl01 = TFHEpp::lvl01param;
    using Lvl02 = TFHEpp::lvl02param;
    using Lvl10 = TFHEpp::lvl10param;
    using Lvl20 = TFHEpp::lvl20param;
    using Lvl21 = TFHEpp::lvl21param;

    using TLWELvl0 = TFHEpp::TLWE<Lvl0>;
    using TLWELvl1 = TFHEpp::TLWE<Lvl1>;
    using TLWELvl2 = TFHEpp::TLWE<Lvl2>;
    using TFHEEvalKey  = TFHEpp::EvalKey;
    using TFHESecretKey = TFHEpp::SecretKey;

    enum class PrunedETHMSBMode {
        WeightedApprox,
        PaperQHalf3PBS
    };

    struct PruneStats {
        uint64_t total_terms = 0;
        uint64_t zero_skipped = 0;
        uint64_t periodic_skipped = 0;
        uint64_t executed_terms = 0;

        void Add(const PruneStats &other)
        {
            total_terms += other.total_terms;
            zero_skipped += other.zero_skipped;
            periodic_skipped += other.periodic_skipped;
            executed_terms += other.executed_terms;
        }
    };

    struct PrunedETHMSBOptions {
        PrunedETHMSBMode mode = PrunedETHMSBMode::WeightedApprox;
        uint32_t period_idx = 4;
        uint32_t window_bits = 5;
    };

    // ── Utility: μ polynomial ──
    template <class P>
    constexpr TFHEpp::Polynomial<P> μ_polygen(typename P::T μ)
    {
        TFHEpp::Polynomial<P> poly;
        for (typename P::T &p : poly) p = -μ;
        return poly;
    }

    // ── Low-level PBS wrappers (same as HE3DB's MSBGateBootstrapping) ──
    void MSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                              const TFHEEvalKey &ek, bool result_type);

    void MSBGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                              const TFHEEvalKey &ek, bool result_type);

    // ── Guard-bit extraction PBS ──
    // Extracts the boundary bit at its original weight via a PBS
    // with LUT coefficients ±weight and ETHMSB offset.
    void GuardBitExtractBS_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                Lvl1::T weight, const TFHEEvalKey &ek);

    void GuardBitExtractBS_Lvl2(TLWELvl2 &res, const TLWELvl2 &tlwe,
                                Lvl2::T weight, const TFHEEvalKey &ek);

    // ── ETHMSB core: Lvl1 path ──
    // Base case: plain_bits <= 5
    void ETHMSB_ExtractMSB5(TLWELvl1 &res, const TLWELvl1 &tlwe,
                            const TFHEEvalKey &ek, bool result_type);

    // Recursive: 5 < plain_bits <= 10
    void ETHMSB_ExtractMSB_Lvl1(TLWELvl1 &res, const TLWELvl1 &tlwe,
                                uint32_t plain_bits, const TFHEEvalKey &ek,
                                bool result_type);

    // ── ETHMSB core: Lvl2 path (output to Lvl1) ──
    // Base case: plain_bits <= 5, Lvl2 → IKS21 → Lvl1 → MSB-BS
    void ETHMSB_ImExtractMSB5(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type);

    // Recursive: 5 < plain_bits <= 9, Lvl2 → IKS21 → Lvl1 → ETHMSB_Lvl1
    void ETHMSB_ImExtractMSB9(TLWELvl1 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type);

    // Recursive: 9 < plain_bits <= 14
    void ETHMSB_ImExtractMSB14(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type);

    // Recursive: 14 < plain_bits <= 19
    void ETHMSB_ImExtractMSB19(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type);

    // Recursive: 19 < plain_bits <= 24
    void ETHMSB_ImExtractMSB24(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type);

    // Recursive: 24 < plain_bits <= 29
    void ETHMSB_ImExtractMSB29(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type);

    // Recursive: 29 < plain_bits <= 33
    void ETHMSB_ImExtractMSB33(TLWELvl1 &res, const TLWELvl2 &tlwe,
                               uint32_t plain_bits, const TFHEEvalKey &ek,
                               bool result_type);

    // ── Dispatch (mirrors HE3DB's HomMSB) ──
    void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type);

    void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &tlwe,
                   uint32_t plain_bits, const TFHEEvalKey &ek,
                   bool result_type);

    // ── Explicit fast/pruned path ──
    // This is deliberately not wired into HomETHMSB. WeightedApprox is an
    // experimental speed-first route; PaperQHalf3PBS keeps the Chapter 3
    // 0/Q2 periodic LUT shape as an audit path.
    void HomETHMSBPrunedFast(TLWELvl1 &res, const TLWELvl1 &tlwe,
                             uint32_t plain_bits, const TFHEEvalKey &ek,
                             bool result_type,
                             const PrunedETHMSBOptions &options,
                             PruneStats *stats = nullptr);

    void HomETHMSBPrunedFast(TLWELvl1 &res, const TLWELvl2 &tlwe,
                             uint32_t plain_bits, const TFHEEvalKey &ek,
                             bool result_type,
                             const PrunedETHMSBOptions &options,
                             PruneStats *stats = nullptr);

    // ── ARI ↔ LOG conversion (re-implemented from HE3DB) ──
    void ARI_to_LOG(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    const TFHEEvalKey &ek);

    void LOG_to_ARI(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    const TFHEEvalKey &ek);

    // ── Gate operators ──
    void HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
                const TFHEEvalKey &ek, bool result_type);

    void HomOR(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
               const TFHEEvalKey &ek, bool result_type);

    template <typename P>
    inline void HomNOT(TFHEpp::TLWE<P> &res, const TFHEpp::TLWE<P> &tlwe)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = -tlwe[i];
    }

    // ── Encrypt / Decrypt helpers (same as HE3DB) ──
    template <class P>
    TFHEpp::TLWE<P> tlweSymInt32Encrypt(const typename P::T p, const double α,
                                        const double scale,
                                        const TFHEpp::Key<P> &key)
    {
        std::uniform_int_distribution<typename P::T> Torusdist(
            0, std::numeric_limits<typename P::T>::max());
        TFHEpp::TLWE<P> res = {};
        res[P::k * P::n] =
            TFHEpp::ModularGaussian<P>(static_cast<typename P::T>(p * scale), α);
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++) {
                res[k * P::n + i] = Torusdist(TFHEpp::generator);
                res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
            }
        return res;
    }

    template <class P>
    typename P::T tlweSymInt32Decrypt(const TFHEpp::TLWE<P> &c,
                                     const double scale,
                                     const TFHEpp::Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus =
            (1ULL << (std::numeric_limits<typename P::T>::digits - 1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        typename P::T res =
            static_cast<typename P::T>(std::round(phase / scale)) %
            plain_modulus;
        return res;
    }

} // namespace ETHMSB_NS
