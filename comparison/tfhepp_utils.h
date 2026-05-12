#pragma once
/**
 * @file tfhepp_utils.h
 * @brief Shared TFHEpp aliases, encoding helpers, gate ops, and rescale
 *        utilities used by both ETHMSB and 3-PBS comparison pipelines.
 *
 * Mirrors the layout of HE3DB's HEDB/comparison/tfhepp_utils.h so that this
 * directory can be used as a drop-in replacement.
 */

#include <cmath>
#include <cstdint>
#include <limits>

#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace tfhepp_compare
{
    // Boolean / arithmetic encoding selector (HE3DB convention).
    static constexpr bool LOGIC      = true;
    static constexpr bool ARITHMETIC = false;
    inline constexpr bool IS_LOGIC(bool a) { return a; }
    inline constexpr bool IS_ARITHMETIC(bool a) { return !a; }

    // Parameter aliases.
    using Lvl0  = TFHEpp::lvl0param;
    using Lvl1  = TFHEpp::lvl1param;
    using Lvl2  = TFHEpp::lvl2param;
    using Lvl01 = TFHEpp::lvl01param;
    using Lvl02 = TFHEpp::lvl02param;
    using Lvl10 = TFHEpp::lvl10param;
    using Lvl20 = TFHEpp::lvl20param;
    using Lvl21 = TFHEpp::lvl21param;

    using TLWELvl0    = TFHEpp::TLWE<Lvl0>;
    using TLWELvl1    = TFHEpp::TLWE<Lvl1>;
    using TLWELvl2    = TFHEpp::TLWE<Lvl2>;
    using TFHEEvalKey   = TFHEpp::EvalKey;
    using TFHESecretKey = TFHEpp::SecretKey;

    // μ test-vector polynomial (every coefficient = -μ; LUT shape used by
    // arithmetic-style PBSes).
    template <class P>
    constexpr TFHEpp::Polynomial<P> μ_polygen(typename P::T μ)
    {
        TFHEpp::Polynomial<P> poly;
        for (typename P::T &p : poly) p = -μ;
        return poly;
    }

    // Identity-LUT polynomial used by IdeGateBootstrapping.
    template <class P>
    TFHEpp::Polynomial<P> g_polygen(uint32_t plain_bits, uint32_t scale_bits)
    {
        TFHEpp::Polynomial<P> poly;
        const uint32_t padding_bits = P::nbit - plain_bits;
        for (int i = 0; i < (int)P::n; i++)
            poly[i] = (typename P::T(1) << scale_bits) * (i >> padding_bits);
        return poly;
    }

    // Encrypt a p-bit integer at the given Δ scale (HE3DB convention).
    template <class P>
    TFHEpp::TLWE<P> tlweSymInt32Encrypt(typename P::T p, double alpha,
                                        double scale,
                                        const TFHEpp::Key<P> &key)
    {
        std::uniform_int_distribution<typename P::T> Torusdist(
            0, std::numeric_limits<typename P::T>::max());
        TFHEpp::TLWE<P> res = {};
        res[P::k * P::n] =
            TFHEpp::ModularGaussian<P>(static_cast<typename P::T>(p * scale),
                                       alpha);
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++) {
                res[k * P::n + i] = Torusdist(TFHEpp::generator);
                res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
            }
        return res;
    }

    template <class P>
    typename P::T tlweSymInt32Decrypt(const TFHEpp::TLWE<P> &c, double scale,
                                     const TFHEpp::Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus =
            (typename P::T(1) << (std::numeric_limits<typename P::T>::digits - 1)) /
            static_cast<typename P::T>(scale);
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        return static_cast<typename P::T>(std::round(phase / scale)) %
               plain_modulus;
    }

    // ── Building-block PBSes shared by both algorithms ──

    // Direct PBS-MSB at level 1 / level 2.  The overloads that take
    // plain_bits use the cell-center offset for that signed-difference width.
    void MSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type);

    void MSBGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                              const TFHEEvalKey &ek, bool result_type);

    void MSBGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                              uint32_t plain_bits, const TFHEEvalKey &ek,
                              bool result_type);

    void MSBGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                              const TFHEEvalKey &ek, bool result_type);

    // Identity PBS used by HE3DB-style HomMSB and the 3-PBS pipeline.
    void IdeGateBootstrapping(TLWELvl1 &res, const TLWELvl1 &tlwe,
                              uint32_t scale_bits, const TFHEEvalKey &ek);

    void IdeGateBootstrapping(TLWELvl2 &res, const TLWELvl2 &tlwe,
                              uint32_t scale_bits, const TFHEEvalKey &ek);

    // ARI ↔ LOG result encoding conversion at level 1.
    void ARI_to_LOG(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    const TFHEEvalKey &ek);

    void LOG_to_ARI(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    const TFHEEvalKey &ek);

    // ── Boolean gate operators (mirrors HE3DB's operators.h subset) ──
    void HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
                const TFHEEvalKey &ek, bool result_type);

    void HomOR(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,
               const TFHEEvalKey &ek, bool result_type);

    template <typename P>
    inline void HomNOT(TFHEpp::TLWE<P> &res, const TFHEpp::TLWE<P> &tlwe)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = -tlwe[i];
    }

} // namespace tfhepp_compare
