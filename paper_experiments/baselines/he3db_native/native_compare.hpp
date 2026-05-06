#pragma once

// Adapted from /home/wrn/HE3DB/src/HEDB/comparison/HomCompare.h.
// This file keeps the HE3DB comparison formulas in the current TFHEpp project
// so the paper baselines can call native HE3DB-style comparison operators
// without adding HE3DB as an external build target.

#include <cstdint>
#include <limits>
#include <stdexcept>

#include "cloudkey.hpp"
#include "gate.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "tlwe.hpp"

namespace PaperReview::HE3DBNative {

using Lvl0 = TFHEpp::lvl0param;
using Lvl1 = TFHEpp::lvl1param;
using Lvl01 = TFHEpp::lvl01param;
using Lvl10 = TFHEpp::lvl10param;
using TLWELvl0 = TFHEpp::TLWE<Lvl0>;
using TLWELvl1 = TFHEpp::TLWE<Lvl1>;
using TFHEEvalKey = TFHEpp::EvalKey;
using TFHESecretKey = TFHEpp::SecretKey;

constexpr bool LOGIC = true;
constexpr bool ARITHMETIC = false;
inline bool IsArithmetic(bool result_type) { return !result_type; }

template <class P>
inline TFHEpp::Polynomial<P> MuPoly(typename P::T mu) {
    TFHEpp::Polynomial<P> poly{};
    for (auto& p : poly) p = -mu;
    return poly;
}

template <class P>
inline TFHEpp::Polynomial<P> GPoly(std::uint32_t plain_bits,
                                   std::uint32_t scale_bits) {
    TFHEpp::Polynomial<P> poly{};
    const std::uint32_t padding_bits = P::nbit - plain_bits;
    for (int i = 0; i < P::n; i++)
        poly[i] = (typename P::T{1} << scale_bits) * (i >> padding_bits);
    return poly;
}

inline void MsbGateBootstrapping(TLWELvl1& res, const TLWELvl1& tlwe,
                                 const TFHEEvalKey& ek, bool result_type) {
    Lvl1::T mu = Lvl1::μ;
    if (IsArithmetic(result_type)) mu <<= 1;
    constexpr auto offset =
        Lvl1::T{1} << (std::numeric_limits<Lvl1::T>::digits - 6);
    auto shifted = tlwe;
    shifted[Lvl1::k * Lvl1::n] += offset;
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, shifted, ek.getiksk<Lvl10>());
    TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
        res, lvl0, ek.getbkfft<Lvl01>(), MuPoly<Lvl1>(mu));
    if (IsArithmetic(result_type)) res[Lvl1::k * Lvl1::n] += mu;
}

inline void IdeGateBootstrapping(TLWELvl1& res, const TLWELvl1& tlwe,
                                 std::uint32_t scale_bits,
                                 const TFHEEvalKey& ek) {
    constexpr auto offset =
        Lvl1::T{1} << (std::numeric_limits<Lvl1::T>::digits - 6);
    auto shifted = tlwe;
    shifted[Lvl1::k * Lvl1::n] += offset;
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, shifted, ek.getiksk<Lvl10>());
    TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
        res, lvl0, ek.getbkfft<Lvl01>(), GPoly<Lvl1>(4, scale_bits));
}

inline void ExtractMSB5(TLWELvl1& res, const TLWELvl1& tlwe,
                        const TFHEEvalKey& ek, bool result_type) {
    MsbGateBootstrapping(res, tlwe, ek, result_type);
}

inline void ExtractMSB9(TLWELvl1& res, const TLWELvl1& tlwe,
                        std::uint32_t plain_bits, const TFHEEvalKey& ek,
                        bool result_type) {
    TLWELvl1 shifted, sign5;
    const std::uint32_t scale_bits =
        std::numeric_limits<Lvl1::T>::digits - plain_bits;
    for (std::size_t i = 0; i <= Lvl1::n; i++) shifted[i] = tlwe[i] << (plain_bits - 5);
    MsbGateBootstrapping(sign5, shifted, ek, ARITHMETIC);
    for (std::size_t i = 0; i <= Lvl1::n; i++) shifted[i] -= sign5[i];
    IdeGateBootstrapping(shifted, shifted, scale_bits, ek);
    for (std::size_t i = 0; i <= Lvl1::n; i++) res[i] = tlwe[i] - shifted[i];
    ExtractMSB5(res, res, ek, result_type);
}

inline void ExtractMSB10(TLWELvl1& res, const TLWELvl1& tlwe,
                         std::uint32_t plain_bits, const TFHEEvalKey& ek,
                         bool result_type) {
    TLWELvl1 shifted, sign5;
    const std::uint32_t scale_bits =
        std::numeric_limits<Lvl1::T>::digits - plain_bits;
    for (std::size_t i = 0; i <= Lvl1::n; i++) shifted[i] = tlwe[i] << (plain_bits - 5);
    MsbGateBootstrapping(sign5, shifted, ek, ARITHMETIC);
    for (std::size_t i = 0; i <= Lvl1::n; i++) shifted[i] -= sign5[i];
    IdeGateBootstrapping(shifted, shifted, scale_bits, ek);
    for (std::size_t i = 0; i <= Lvl1::n; i++) res[i] = tlwe[i] - shifted[i];
    ExtractMSB9(res, res, plain_bits - 4, ek, result_type);
}

inline void HomMSB(TLWELvl1& res, const TLWELvl1& tlwe,
                   std::uint32_t plain_bits, const TFHEEvalKey& ek,
                   bool result_type) {
    if (plain_bits <= 5) return ExtractMSB5(res, tlwe, ek, result_type);
    if (plain_bits <= 9) return ExtractMSB9(res, tlwe, plain_bits, ek, result_type);
    if (plain_bits == 10) return ExtractMSB10(res, tlwe, plain_bits, ek, result_type);
    throw std::invalid_argument("HE3DB lvl1 HomMSB supports at most 10 bits");
}

inline void HomAND(TLWELvl1& res, const TLWELvl1& ca, const TLWELvl1& cb,
                   const TFHEEvalKey& ek, bool result_type) {
    Lvl1::T offset = Lvl1::μ;
    if (IsArithmetic(result_type)) offset <<= 1;
    for (int i = 0; i <= Lvl1::k * Lvl1::n; i++) res[i] = ca[i] + cb[i];
    res[Lvl1::k * Lvl1::n] -= Lvl1::μ >> 1;
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, res, ek.getiksk<Lvl10>());
    TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
        res, lvl0, ek.getbkfft<Lvl01>(), MuPoly<Lvl1>(-offset));
    if (IsArithmetic(result_type)) res[Lvl1::k * Lvl1::n] += offset;
}

template <typename P>
inline void greater_than(TFHEpp::TLWE<P>& cipher1, TFHEpp::TLWE<P>& cipher2,
                         TLWELvl1& res, std::uint32_t plain_bits,
                         TFHEEvalKey& ek, bool result_type) {
    TFHEpp::TLWE<P> sub{};
    for (std::size_t i = 0; i <= P::k * P::n; i++) sub[i] = cipher2[i] - cipher1[i];
    HomMSB(res, sub, plain_bits + 1, ek, result_type);
}

template <typename P>
inline void less_than(TFHEpp::TLWE<P>& cipher1, TFHEpp::TLWE<P>& cipher2,
                      TLWELvl1& res, std::uint32_t plain_bits,
                      TFHEEvalKey& ek, bool result_type) {
    TFHEpp::TLWE<P> sub{};
    for (std::size_t i = 0; i <= P::k * P::n; i++) sub[i] = cipher1[i] - cipher2[i];
    HomMSB(res, sub, plain_bits + 1, ek, result_type);
}

template <typename P>
inline void greater_than_equal(TFHEpp::TLWE<P>& cipher1,
                               TFHEpp::TLWE<P>& cipher2,
                               TLWELvl1& res, std::uint32_t plain_bits,
                               TFHEEvalKey& ek, bool result_type) {
    TFHEpp::TLWE<P> sub{};
    for (std::size_t i = 0; i <= P::k * P::n; i++) sub[i] = cipher1[i] - cipher2[i];
    HomMSB(res, sub, plain_bits + 1, ek, LOGIC);
    TFHEpp::HomNOT<Lvl1>(res, res);
    if (IsArithmetic(result_type)) {
        TLWELvl0 lvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, res, ek.getiksk<Lvl10>());
        TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
            res, lvl0, ek.getbkfft<Lvl01>(), MuPoly<Lvl1>(-(Lvl1::μ << 1)));
        res[Lvl1::k * Lvl1::n] += Lvl1::μ;
    }
}

template <typename P>
inline void less_than_equal(TFHEpp::TLWE<P>& cipher1,
                            TFHEpp::TLWE<P>& cipher2,
                            TLWELvl1& res, std::uint32_t plain_bits,
                            TFHEEvalKey& ek, bool result_type) {
    TFHEpp::TLWE<P> sub{};
    for (std::size_t i = 0; i <= P::k * P::n; i++) sub[i] = cipher2[i] - cipher1[i];
    HomMSB(res, sub, plain_bits + 1, ek, LOGIC);
    TFHEpp::HomNOT<Lvl1>(res, res);
    if (IsArithmetic(result_type)) {
        TLWELvl0 lvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, res, ek.getiksk<Lvl10>());
        TFHEpp::GateBootstrappingTLWE2TLWE<Lvl01>(
            res, lvl0, ek.getbkfft<Lvl01>(), MuPoly<Lvl1>(-(Lvl1::μ << 1)));
        res[Lvl1::k * Lvl1::n] += Lvl1::μ;
    }
}

template <typename P>
inline void equal(TFHEpp::TLWE<P>& cipher1, TFHEpp::TLWE<P>& cipher2,
                  TLWELvl1& res, std::uint32_t plain_bits,
                  TFHEEvalKey& ek, bool result_type) {
    TLWELvl1 ge, le;
    greater_than_equal<P>(cipher1, cipher2, ge, plain_bits, ek, LOGIC);
    less_than_equal<P>(cipher1, cipher2, le, plain_bits, ek, LOGIC);
    HomAND(res, ge, le, ek, result_type);
}

template <class P>
inline TFHEpp::TLWE<P> EncryptInt(std::uint64_t value, std::uint32_t bits,
                                  const TFHESecretKey& sk) {
    if (bits == 0 || bits >= std::numeric_limits<typename P::T>::digits)
        throw std::invalid_argument("invalid HE3DB plaintext bit width");
    TFHEpp::TLWE<P> out;
    switch (bits) {
    case 1:
        TFHEpp::tlweSymIntEncrypt<P, 2>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 2:
        TFHEpp::tlweSymIntEncrypt<P, 4>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 3:
        TFHEpp::tlweSymIntEncrypt<P, 8>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 4:
        TFHEpp::tlweSymIntEncrypt<P, 16>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 5:
        TFHEpp::tlweSymIntEncrypt<P, 32>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 6:
        TFHEpp::tlweSymIntEncrypt<P, 64>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 7:
        TFHEpp::tlweSymIntEncrypt<P, 128>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 8:
        TFHEpp::tlweSymIntEncrypt<P, 256>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 9:
        TFHEpp::tlweSymIntEncrypt<P, 512>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    case 10:
        TFHEpp::tlweSymIntEncrypt<P, 1024>(out, static_cast<typename P::T>(value), P::α, sk.key.get<P>());
        break;
    default:
        throw std::invalid_argument("HE3DB lvl1 baseline supports up to 10-bit integer encryption");
    }
    return out;
}

inline int DecryptLogic(const TLWELvl1& ct, const TFHESecretKey& sk) {
    return TFHEpp::tlweSymDecrypt<Lvl1>(ct, sk.key.get<Lvl1>()) ? 1 : 0;
}

}  // namespace PaperReview::HE3DBNative
