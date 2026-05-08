#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>

#include "evalkeygens.hpp"
#include "key.hpp"
#include "params.hpp"

namespace my_ethmsb_params {

struct my_h3_lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = 1;
    static constexpr std::uint32_t n = 672;
    static constexpr std::uint32_t k = 1;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double alpha = std::pow(2.0, -16);
    static const inline double α = alpha;
    using T = std::uint32_t;
    static constexpr T mu = T{1} << (std::numeric_limits<T>::digits - 2);
    static constexpr T μ = mu;
    static constexpr std::uint32_t plain_modulus = 2;
    static constexpr double Delta =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
    static constexpr double Δ = Delta;
};

struct my_h3_lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_diff = 2;
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t l_a = l;
    static constexpr std::uint32_t lₐ = l_a;
    static constexpr std::uint32_t Bgbit = 6;
    static constexpr std::uint32_t Bg_a_bit = Bgbit;
    static constexpr std::uint32_t Bgₐbit = Bg_a_bit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bg_a = Bg;
    static constexpr std::uint32_t Bgₐ = Bg_a;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double alpha = std::pow(2.0, -25);
    static const inline double α = alpha;
    using T = std::uint32_t;
    static constexpr T mu = T{1} << 29;
    static constexpr T μ = mu;
    static constexpr std::uint32_t plain_modulus = 2;
    static constexpr double Delta =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
    static constexpr double Δ = Delta;
};

struct my_h3_lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_diff = 2;
    static constexpr std::uint32_t nbit = 11;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t l_a = l;
    static constexpr std::uint32_t lₐ = l_a;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bg_a_bit = Bgbit;
    static constexpr std::uint32_t Bgₐbit = Bg_a_bit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bg_a = Bg;
    static constexpr std::uint32_t Bgₐ = Bg_a;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double alpha = std::pow(2.0, -52);
    static const inline double α = alpha;
    using T = std::uint64_t;
    static constexpr T mu = T{1} << 61;
    static constexpr T μ = mu;
    static constexpr std::uint32_t plain_modulus = 8;
    static constexpr double Delta = static_cast<double>(mu);
    static constexpr double Δ = Delta;
};

struct my_h3_lvl20param {
    static constexpr std::uint32_t t = 2;
    static constexpr std::uint32_t basebit = 10;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double alpha = my_h3_lvl0param::alpha;
    static const inline double α = alpha;
    using domainP = my_h3_lvl2param;
    using targetP = my_h3_lvl0param;
};

struct my_h3_lvl10param {
    static constexpr std::uint32_t t = 2;
    static constexpr std::uint32_t basebit = 10;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double alpha = my_h3_lvl0param::alpha;
    static const inline double α = alpha;
    using domainP = my_h3_lvl1param;
    using targetP = my_h3_lvl0param;
};

struct my_h3_lvl21param {
    static constexpr std::uint32_t t = 10;
    static constexpr std::uint32_t basebit = 3;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double alpha = my_h3_lvl1param::alpha;
    static const inline double α = alpha;
    using domainP = my_h3_lvl2param;
    using targetP = my_h3_lvl1param;
};

struct my_h3_lvl01param {
    using domainP = my_h3_lvl0param;
    using targetP = my_h3_lvl1param;
#ifdef USE_KEY_BUNDLE
    static constexpr std::uint32_t Addends = 2;
#else
    static constexpr std::uint32_t Addends = 1;
#endif
};

struct my_h3_lvl02param {
    using domainP = my_h3_lvl0param;
    using targetP = my_h3_lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr std::uint32_t Addends = 2;
#else
    static constexpr std::uint32_t Addends = 1;
#endif
};

struct my_v10_lvl2param : TFHEpp::lvl2param {
    static constexpr int32_t key_value_diff =
        TFHEpp::lvl2param::key_value_max - TFHEpp::lvl2param::key_value_min;
    static constexpr std::uint32_t l_a = TFHEpp::lvl2param::lₐ;
    static constexpr std::uint32_t Bg_a_bit = TFHEpp::lvl2param::Bgₐbit;
    static constexpr std::uint32_t Bg_a = TFHEpp::lvl2param::Bgₐ;
};

struct my_v10_lvl22param {
    using domainP = my_v10_lvl2param;
    using targetP = my_v10_lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr std::uint32_t Addends = 2;
#else
    static constexpr std::uint32_t Addends = 1;
#endif
};

template <class P, class URBG>
TFHEpp::Key<P> make_key(URBG& rng)
{
    TFHEpp::Key<P> key;
    std::uniform_int_distribution<int32_t> dist(P::key_value_min,
                                                P::key_value_max);
    for (auto& v : key) v = static_cast<typename P::T>(dist(rng));
    return key;
}

template <class P>
TFHEpp::TLWE<P> trivial_tlwe(const typename P::T body)
{
    TFHEpp::TLWE<P> ct = {};
    ct[P::k * P::n] = body;
    return ct;
}

template <class P>
typename P::T phase(const TFHEpp::TLWE<P>& ct, const TFHEpp::Key<P>& key)
{
    return TFHEpp::tlweSymPhase<P>(ct, key);
}

struct H3CompatKeys {
    TFHEpp::Key<my_h3_lvl0param> key0;
    TFHEpp::Key<my_h3_lvl1param> key1;
    TFHEpp::Key<my_h3_lvl2param> key2;
};

struct V10DirectKeys {
    TFHEpp::Key<my_v10_lvl2param> key2;
};

template <class URBG>
H3CompatKeys make_h3compat_keys(URBG& rng)
{
    return {make_key<my_h3_lvl0param>(rng), make_key<my_h3_lvl1param>(rng),
            make_key<my_h3_lvl2param>(rng)};
}

template <class URBG>
V10DirectKeys make_v10_direct_keys(URBG& rng)
{
    return {make_key<my_v10_lvl2param>(rng)};
}

}  // namespace my_ethmsb_params
