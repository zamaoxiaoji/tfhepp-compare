#pragma once

// Adapted from /home/wrn/AE_submit/SOTA/ArcEDB/src/ARCEDB/comparison/comparable.h
// and comparable.cpp. This header isolates ArcEDB's TFHEpp-native
// exponent/equality/ordering operators for the paper enumerated baselines.

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cloudkey.hpp"
#include "gate.hpp"
#include "keyswitch.hpp"
#include "tlwe.hpp"
#include "trlwe.hpp"
#include "trgsw.hpp"

namespace PaperReview::ArcEDBNative {

using Lvl0 = TFHEpp::lvl0param;
using Lvl1 = TFHEpp::lvl1param;
using Lvl01 = TFHEpp::lvl01param;
using Lvl10 = TFHEpp::lvl10param;
using TLWELvl0 = TFHEpp::TLWE<Lvl0>;
using TLWELvl1 = TFHEpp::TLWE<Lvl1>;
using TRLWELvl1 = TFHEpp::TRLWE<Lvl1>;
using TRGSWLvl1 = TFHEpp::TRGSWFFT<Lvl1>;
using TFHEEvalKey = TFHEpp::EvalKey;
using TFHESecretKey = TFHEpp::SecretKey;

inline std::size_t BitsCount(std::uint64_t v) {
    return v == 0 ? 1 : static_cast<std::size_t>(std::log2(v)) + 1;
}

template <class P>
inline TFHEpp::Polynomial<P> MuPoly(typename P::T mu) {
    TFHEpp::Polynomial<P> poly{};
    for (auto& p : poly) p = -mu;
    return poly;
}

template <class P>
inline void exponent_encrypt(std::uint64_t value, TFHEpp::TRLWE<P>& cipher,
                             TFHESecretKey& sk) {
    if (value >= P::n) throw std::invalid_argument("ArcEDB exponent value out of range");
    std::array<typename P::T, P::n> plain{};
    plain[value] = P::μ;
    TFHEpp::trlweSymEncrypt<P>(cipher, plain, P::α, sk.key.get<P>());
}

template <class P>
inline void exponent_encrypt(std::uint64_t value, std::uint32_t precision,
                             std::vector<TFHEpp::TRLWE<P>>& ciphers,
                             TFHESecretKey& sk) {
    if (precision != 64 && value >= (std::uint64_t{1} << precision))
        throw std::invalid_argument("ArcEDB exponent value exceeds precision");
    const auto digit_count = static_cast<std::size_t>(
        std::ceil(static_cast<double>(precision) / std::log2(P::n)));
    ciphers.resize(std::max<std::size_t>(1, digit_count));
    for (auto& cipher : ciphers) {
        const auto digit = value % P::n;
        exponent_encrypt<P>(digit, cipher, sk);
        value = (value - digit) / P::n;
    }
}

template <class P>
inline void exponent_encrypt_rgsw(std::uint64_t value,
                                  TFHEpp::TRGSWFFT<P>& cipher,
                                  TFHESecretKey& sk,
                                  bool is_negative) {
    if (value >= P::n) throw std::invalid_argument("ArcEDB RGSW exponent value out of range");
    constexpr std::array<typename P::T, P::l> h = TFHEpp::hgen<P, false>();
    TFHEpp::TRGSW<P> trgsw;
    for (auto& trlwe : trgsw)
        TFHEpp::trlweSymEncryptZero<P>(trlwe, P::α, sk.key.get<P>());
    if (value != 0) {
        for (int i = 0; i < P::l; i++) {
            for (int k = 0; k < P::k + 1; k++) {
                if (is_negative)
                    trgsw[i + k * P::l][k][P::n - value] -= h[i];
                else
                    trgsw[i + k * P::l][k][value] += h[i];
            }
        }
    } else {
        for (int i = 0; i < P::l; i++) {
            for (int k = 0; k < P::k + 1; k++)
                trgsw[i + k * P::l][k][value] += h[i];
        }
    }
    cipher = TFHEpp::ApplyFFT2trgsw<P>(trgsw);
}

template <class P>
inline void exponent_encrypt_rgsw(std::uint64_t value, std::uint32_t precision,
                                  std::vector<TFHEpp::TRGSWFFT<P>>& ciphers,
                                  TFHESecretKey& sk, bool is_negative) {
    if (precision != 64 && value >= (std::uint64_t{1} << precision))
        throw std::invalid_argument("ArcEDB RGSW exponent value exceeds precision");
    const auto digit_count = static_cast<std::size_t>(
        std::ceil(static_cast<double>(precision) / std::log2(P::n)));
    ciphers.resize(std::max<std::size_t>(1, digit_count));
    for (auto& cipher : ciphers) {
        const auto digit = value % P::n;
        exponent_encrypt_rgsw<P>(digit, cipher, sk, is_negative);
        value = (value - digit) / P::n;
    }
}

inline void greater_than_tfhepp(TRLWELvl1& cipher1, TRGSWLvl1& cipher2,
                                TLWELvl1& res, TFHESecretKey&) {
    TRLWELvl1 trlwe_mul, trlwe_product;
    TFHEpp::ExternalProduct<Lvl1>(trlwe_mul, cipher1, cipher2);
    TFHEpp::Polynomial<Lvl1> test_plaintext;
    for (auto& x : test_plaintext) x = 1;
    TFHEpp::PolyMul<Lvl1>(trlwe_product[0], trlwe_mul[0], test_plaintext);
    TFHEpp::PolyMul<Lvl1>(trlwe_product[1], trlwe_mul[1], test_plaintext);
    TFHEpp::SampleExtractIndex<Lvl1>(res, trlwe_product, 0);
    for (std::size_t i = 0; i <= Lvl1::n; i++) res[i] = -res[i];
}

inline void equality_tfhepp(TRLWELvl1& cipher1, TRGSWLvl1& cipher2,
                            TLWELvl1& res, TFHESecretKey&) {
    TRLWELvl1 trlwe_mul;
    TFHEpp::ExternalProduct<Lvl1>(trlwe_mul, cipher1, cipher2);
    TFHEpp::SampleExtractIndex<Lvl1>(res, trlwe_mul, 0);
    for (std::size_t i = 0; i <= Lvl1::n; i++) res[i] = 2 * res[i];
    res[Lvl1::n] -= Lvl1::μ;
}

inline void less_than_tfhepp(TRLWELvl1& cipher1, TRGSWLvl1& cipher2,
                             TLWELvl1& res, TFHESecretKey&) {
    TRLWELvl1 trlwe_mul, trlwe_product;
    TFHEpp::ExternalProduct<Lvl1>(trlwe_mul, cipher1, cipher2);
    TFHEpp::Polynomial<Lvl1> test_plaintext{};
    test_plaintext[0] = Lvl1::plain_modulus - 1;
    for (std::size_t i = 1; i < Lvl1::n; i++) test_plaintext[i] = 1;
    TFHEpp::PolyMul<Lvl1>(trlwe_product[0], trlwe_mul[0], test_plaintext);
    TFHEpp::PolyMul<Lvl1>(trlwe_product[1], trlwe_mul[1], test_plaintext);
    TFHEpp::SampleExtractIndex<Lvl1>(res, trlwe_product, 0);
}

inline void equality_tfhepp(std::vector<TRLWELvl1>& ciphers1,
                            std::vector<TRGSWLvl1>& ciphers2,
                            std::size_t cipher_size,
                            TLWELvl1& res,
                            TFHEEvalKey& ek,
                            TFHESecretKey& sk) {
    if (cipher_size == 0) throw std::invalid_argument("empty ArcEDB equality input");
    if (cipher_size == 1) {
        equality_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        return;
    }
    TLWELvl1 low, high;
    equality_tfhepp(ciphers1, ciphers2, cipher_size - 1, low, ek, sk);
    equality_tfhepp(ciphers1[cipher_size - 1], ciphers2[cipher_size - 1], high, sk);
    TFHEpp::HomAND(res, low, high, ek);
}

inline void less_than_tfhepp(std::vector<TRLWELvl1>& ciphers1,
                             std::vector<TRGSWLvl1>& ciphers2,
                             std::size_t cipher_size,
                             TLWELvl1& res,
                             TFHEEvalKey& ek,
                             TFHESecretKey& sk) {
    if (cipher_size == 0) throw std::invalid_argument("empty ArcEDB less-than input");
    if (cipher_size == 1) {
        less_than_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        return;
    }
    TLWELvl1 low, high, equal_res;
    TRLWELvl1 trlwe_mul;
    less_than_tfhepp(ciphers1, ciphers2, cipher_size - 1, low, ek, sk);
    TFHEpp::ExternalProduct<Lvl1>(
        trlwe_mul, ciphers1[cipher_size - 1], ciphers2[cipher_size - 1]);
    TFHEpp::SampleExtractIndex<Lvl1>(equal_res, trlwe_mul, 0);
    less_than_tfhepp(ciphers1[cipher_size - 1], ciphers2[cipher_size - 1], high, sk);
    for (std::size_t i = 0; i <= Lvl1::n; i++) high[i] += high[i];

    TLWELvl1 selector;
    const std::uint32_t offset = Lvl1::μ >> 1;
    for (std::size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
        selector[i] = equal_res[i] + high[i] + low[i];
    selector[Lvl1::n] += offset;
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, selector, ek.getiksk<Lvl10>());
    TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
        res, lvl0, ek.getbkfft<Lvl01>(),
        TFHEpp::μpolygen<Lvl1, Lvl1::μ>());
}

inline void greater_than_tfhepp(std::vector<TRLWELvl1>& ciphers1,
                                std::vector<TRGSWLvl1>& ciphers2,
                                std::size_t cipher_size,
                                TLWELvl1& res,
                                TFHEEvalKey& ek,
                                TFHESecretKey& sk) {
    if (cipher_size == 0) throw std::invalid_argument("empty ArcEDB greater-than input");
    if (cipher_size == 1) {
        greater_than_tfhepp(ciphers1[0], ciphers2[0], res, sk);
        return;
    }
    TLWELvl1 low, high, equal_res;
    TRLWELvl1 trlwe_mul;
    greater_than_tfhepp(ciphers1, ciphers2, cipher_size - 1, low, ek, sk);
    TFHEpp::ExternalProduct<Lvl1>(
        trlwe_mul, ciphers1[cipher_size - 1], ciphers2[cipher_size - 1]);
    TFHEpp::SampleExtractIndex<Lvl1>(equal_res, trlwe_mul, 0);
    greater_than_tfhepp(ciphers1[cipher_size - 1], ciphers2[cipher_size - 1], high, sk);
    for (std::size_t i = 0; i <= Lvl1::n; i++) high[i] += high[i];

    TLWELvl1 selector;
    const std::uint32_t offset = Lvl1::μ >> 1;
    for (std::size_t i = 0; i <= Lvl1::k * Lvl1::n; i++)
        selector[i] = equal_res[i] + high[i] + low[i];
    selector[Lvl1::n] += offset;
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, selector, ek.getiksk<Lvl10>());
    TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
        res, lvl0, ek.getbkfft<Lvl01>(),
        TFHEpp::μpolygen<Lvl1, Lvl1::μ>());
}

inline int DecryptLogic(const TLWELvl1& ct, const TFHESecretKey& sk) {
    return TFHEpp::tlweSymDecrypt<Lvl1>(ct, sk.key.get<Lvl1>()) ? 1 : 0;
}

}  // namespace PaperReview::ArcEDBNative
