#pragma once
#include <cmath>
#include <limits>
#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"
#include "gatebootstrapping.hpp"
#include "pruned_bootstrapping.hpp"

namespace TFHEpp
{
    // ============================================================
    //  Encrypt / Decrypt helpers (arithmetic integer encoding)
    // ============================================================

    // Encrypt integer p at given scale: phase = p * scale + noise
    template <class P>
    TLWE<P> tlweSymInt32Encrypt(const typename P::T p, const double α,
                                const double scale, const Key<P> &key)
    {
        std::uniform_int_distribution<typename P::T> Torusdist(
            0, std::numeric_limits<typename P::T>::max());
        TLWE<P> res = {};
        res[P::k * P::n] =
            ModularGaussian<P>(static_cast<typename P::T>(p * scale), α);
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++) {
                res[k * P::n + i] = Torusdist(generator);
                res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
            }
        return res;
    }

    // Decrypt to integer: round(phase / scale) mod plain_modulus
    template <class P>
    typename P::T tlweSymInt32Decrypt(const TLWE<P> &c, const double scale,
                                      const Key<P> &key)
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

    // ============================================================
    //  TRLWE encrypt / decrypt
    // ============================================================

    template <class P>
    TRLWE<P> trlweSymInt32Encrypt(const std::array<typename P::T, P::n> &p,
                                   const double α, const double scale,
                                   const Key<P> &key)
    {
        TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
        for (int i = 0; i < P::n; i++)
            c[P::k][i] += static_cast<typename P::T>(scale * p[i]);
        return c;
    }

    template <class P>
    Polynomial<P> trlweSymInt32Decrypt(const TRLWE<P> &c, double scale,
                                        const Key<P> &key)
    {
        Polynomial<P> phase = c[P::k];
        typename P::T plain_modulus =
            (1ULL << (std::numeric_limits<typename P::T>::digits - 1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++) {
            Polynomial<P> mulres;
            std::array<typename P::T, P::n> partkey;
            for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
            PolyMul<P>(mulres, c[k], partkey);
            for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
        }
        Polynomial<P> p;
        for (int i = 0; i < P::n; i++)
            p[i] = static_cast<typename P::T>(std::round(phase[i] / scale)) %
                   plain_modulus;
        return p;
    }

    // ============================================================
    //  Polynomial generators
    // ============================================================

    // All-constant polynomial with value -μ (used as standard MSB LUT)
    template <class P>
    constexpr Polynomial<P> μ_polygen(typename P::T μ)
    {
        Polynomial<P> poly;
        for (typename P::T &p : poly) p = -μ;
        return poly;
    }

    // Identity polynomial: f(m) = m * 2^scale_bits
    template <class P>
    Polynomial<P> gpolygen(uint32_t plain_bits, uint32_t scale_bits)
    {
        Polynomial<P> poly;
        uint32_t padding_bits = P::nbit - plain_bits;
        for (int i = 0; i < P::n; i++)
            poly[i] = (1ULL << scale_bits) * (i >> padding_bits);
        return poly;
    }

    // ============================================================
    //  TLWE / TRLWE arithmetic helpers
    // ============================================================

    template <class P>
    void TLWEAdd(TLWE<P> &ca, TLWE<P> &cb, TLWE<P> &res)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i] + cb[i];
    }

    template <class P>
    void TLWESub(TLWE<P> &ca, TLWE<P> &cb, TLWE<P> &res)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i] - cb[i];
    }

    // ============================================================
    //  LOG ↔ ARI conversions  (HE3DB style)
    // ============================================================

    // ARI {0, Q/2} → LOG {-Q/8, Q/8}
    void ARI_to_LOG(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                    const EvalKey &ek);

    // LOG {-Q/8, Q/8} → ARI {0, Q/2}
    void LOG_to_ARI(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                    const EvalKey &ek);

    void log_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                     uint32_t scale_bits, const EvalKey &ek);

    void ari_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                     uint32_t scale_bits, const EvalKey &ek);

    // ============================================================
    //  Standard MSB Gate Bootstrapping  (HE3DB original)
    // ============================================================

    // lvl1 → lvl1
    void MSBGateBootstrapping(TLWE<lvl1param> &res,
                              const TLWE<lvl1param> &tlwe,
                              const EvalKey &ek, bool result_type);

    // lvl2 → lvl1 (key-switch + bootstrap)
    void MSBGateBootstrapping(TLWE<lvl1param> &res,
                              const TLWE<lvl2param> &tlwe,
                              const EvalKey &ek, bool result_type);

    // lvl2 → lvl2
    void MSBGateBootstrapping(TLWE<lvl2param> &res,
                              const TLWE<lvl2param> &tlwe,
                              const EvalKey &ek, bool result_type);

    // ============================================================
    //  Identity Gate Bootstrapping  (rescale)
    // ============================================================

    void IdeGateBootstrapping(TLWE<lvl1param> &res,
                              const TLWE<lvl1param> &tlwe,
                              uint32_t scale_bits, const EvalKey &ek);

    void IdeGateBootstrapping(TLWE<lvl1param> &res,
                              const TLWE<lvl2param> &tlwe,
                              uint32_t scale_bits, const EvalKey &ek);

    void IdeGateBootstrapping(TLWE<lvl2param> &res,
                              const TLWE<lvl2param> &tlwe,
                              uint32_t scale_bits, const EvalKey &ek);

    // ============================================================
    //  Pruned MSB Gate Bootstrapping  (with periodic CMUX skip)
    // ============================================================
    //
    // Extracts the MSB using a periodic test vector and skips CMUX
    // steps where ā_i ≡ 0 (mod period_M).
    //
    // period_M:  LUT period in [0, 2N].  0 = no pruning.
    //
    // lvl1 → lvl1
    void PrunedMSBGateBootstrapping(TLWE<lvl1param> &res,
                                    const TLWE<lvl1param> &tlwe,
                                    const EvalKey &ek,
                                    bool result_type,
                                    uint32_t period_M = 0);

    // lvl2 → lvl2
    void PrunedMSBGateBootstrapping(TLWE<lvl2param> &res,
                                    const TLWE<lvl2param> &tlwe,
                                    const EvalKey &ek,
                                    bool result_type,
                                    uint32_t period_M = 0);

    // ============================================================
    //  Gap-MSB Gate Bootstrapping
    // ============================================================
    //
    // After clearing bit guard_k from the ciphertext (creating an
    // unreachable interval near the MSB threshold), this function
    // performs MSB extraction with an adjusted offset that aligns
    // the decision boundary to the midpoint of the gap.
    //
    // offset' = (w_k + 1) / 2 * Δ   where w_k = 2^(p-1-guard_k)
    //
    // Decryption must use the shifted Δ' boundary.
    //
    // lvl1 → lvl1
    void GapMSBGateBootstrapping(TLWE<lvl1param> &res,
                                 const TLWE<lvl1param> &tlwe,
                                 const EvalKey &ek,
                                 bool result_type,
                                 uint32_t guard_k,
                                 uint32_t plain_bits);

    // lvl2 → lvl2
    void GapMSBGateBootstrapping(TLWE<lvl2param> &res,
                                 const TLWE<lvl2param> &tlwe,
                                 const EvalKey &ek,
                                 bool result_type,
                                 uint32_t guard_k,
                                 uint32_t plain_bits);

} // namespace TFHEpp
