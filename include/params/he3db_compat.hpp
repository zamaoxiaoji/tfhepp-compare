#pragma once
/**
 * @file he3db_compat.hpp
 * @brief TFHEpp parameters matching HE3DB's TFHEpp fork for cross-validation.
 *        Key differences from default 128bit.hpp:
 *        - lvl0param::T = uint32_t (not uint16_t)
 *        - lvl0param::n = 672 (not 636)
 *        - lvl1param::l = 3, Bgbit = 6 (not l=2, Bgbit=8)
 *        - lvl10param: t=2, basebit=10 (not t=4, basebit=3)
 */

#include <cmath>
#include <cstdint>
#include <limits>

// ═════════════ LWE / RLWE core params ═════════════

struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 672;  // HE3DB uses 672
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -16);  // HE3DB noise
    using T = uint32_t;  // 32-bit torus (critical for precision)
    static constexpr std::make_signed_t<T> μ =
        1U << (std::numeric_limits<T>::digits - 2);
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvlhalfparam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 760;
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -17);
    using T = uint32_t;
    static constexpr T μ = 1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 32;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 3;       // HE3DB uses 3
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 6;   // HE3DB uses 6
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -25);
    using T = uint32_t;
    static constexpr std::make_signed_t<T> μ = 1 << 29;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct AHlvl1param {
    using baseP = lvl1param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 5;
    static constexpr std::uint32_t Bgₐbit = 5;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;
    using T = typename baseP::T;
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
};

struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static const std::uint32_t nbit = 11;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;        // HE3DB uses 4
    static constexpr std::uint32_t Bgbit = 9;    // HE3DB uses 9
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -52);  // HE3DB uses -52
    using T = uint64_t;
    static constexpr std::make_signed_t<T> μ = 1LL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;  // HE3DB: Δ = μ
};

struct cblvl2param {
    using baseP = lvl2param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;
    static constexpr std::uint32_t k = baseP::k;
    static_assert(n*k==baseP::n*baseP::k);
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;
    using T = typename baseP::T;
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
};

struct AHlvl2param {
    using baseP = cblvl2param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 5;
    static constexpr std::uint32_t l = 5;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;
    using T = typename baseP::T;
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
};

struct lvl3param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static const std::uint32_t nbit = 13;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -51);
    using T = uint64_t;
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr uint64_t plain_modulus = 1ULL << plain_modulusbit;
    static constexpr double Δ = 1ULL << (64 - plain_modulusbit - 1);
};

// ═════════════ Key Switching params ═════════════

struct lvl10param {
    static constexpr std::uint32_t t = 2;       // HE3DB uses 2
    static constexpr std::uint32_t basebit = 10; // HE3DB uses 10
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;
    using domainP = lvl1param;
    using targetP = lvl0param;
};

struct lvl1hparam {
    static constexpr std::uint32_t t = 10;
    static constexpr std::uint32_t basebit = 3;
    static const inline double α = lvlhalfparam::α;
    using domainP = lvl1param;
    using targetP = lvlhalfparam;
};

struct lvl11param {
    static constexpr std::uint32_t t = 6;
    static constexpr std::uint32_t basebit = 4;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl1param::α;
    using domainP = lvl1param;
    using targetP = lvl1param;
};

struct lvl20param {
    static constexpr std::uint32_t t = 2;       // HE3DB uses 2
    static constexpr std::uint32_t basebit = 10; // HE3DB uses 10
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;
    using domainP = lvl2param;
    using targetP = lvl0param;
};

struct lvl2hparam {
    static constexpr std::uint32_t t = 7;
    static constexpr std::uint32_t basebit = 2;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvlhalfparam::α;
    using domainP = lvl2param;
    using targetP = lvlhalfparam;
};

struct lvl21param {
    static constexpr std::uint32_t t = 10;      // HE3DB uses 10
    static constexpr std::uint32_t basebit = 3;  // HE3DB uses 3
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl1param::α;
    using domainP = lvl2param;
    using targetP = lvl1param;
};

struct lvl22param {
    static constexpr std::uint32_t t = 8;       // HE3DB uses 8
    static constexpr std::uint32_t basebit = 4;  // HE3DB uses 4
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl2param::α;
    using domainP = lvl2param;
    using targetP = lvl2param;
};

struct lvl31param {
    static constexpr std::uint32_t t = 7;
    static constexpr std::uint32_t basebit = 2;
    static const inline double α = lvl1param::α;
    using domainP = lvl3param;
    using targetP = lvl1param;
};
