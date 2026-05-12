#pragma once
/**
 * @file micro_params.h
 * @brief Local UNSAFE_PERFORMANCE_ONLY parameter structs for isolated
 *        Boolean Q/2 -> guard_value conversion sweeps.
 *
 * These structs are deliberately kept outside TFHEpp vendor headers.  They are
 * used only through explicit keys passed to TFHEpp templates, so EvalKey /
 * SecretKey registries do not need to know about them.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "tfhepp_utils.h"

namespace tfhepp_compare::micro_pbs
{
    template <std::uint32_t N>
    struct MicroNBit;

    template <>
    struct MicroNBit<128> {
        static constexpr std::uint32_t value = 7;
    };

    template <>
    struct MicroNBit<256> {
        static constexpr std::uint32_t value = 8;
    };

    template <>
    struct MicroNBit<512> {
        static constexpr std::uint32_t value = 9;
    };

    template <std::uint32_t NIn, int AlphaLog2 = -45>
    struct MicroInParam {
        static constexpr int32_t key_value_max = 1;
        static constexpr int32_t key_value_min = 0;
        static constexpr int32_t key_value_diff = key_value_max - key_value_min;
        static constexpr std::uint32_t n = NIn;
        static constexpr std::uint32_t k = 1;
        static constexpr TFHEpp::ErrorDistribution errordist =
            TFHEpp::ErrorDistribution::ModularGaussian;
        static const inline double α = std::pow(2.0, AlphaLog2);
        using T = typename Lvl1::T;
        static constexpr std::make_signed_t<T> μ = 1 << 29;
        static constexpr std::uint32_t plain_modulus = 8;
        static constexpr double Δ =
            std::ldexp(1.0, std::numeric_limits<T>::digits) / plain_modulus;
    };

    template <std::uint32_t NIn, int AlphaLog2 = -60>
    struct MicroIn64Param {
        static constexpr int32_t key_value_max = 1;
        static constexpr int32_t key_value_min = 0;
        static constexpr int32_t key_value_diff = key_value_max - key_value_min;
        static constexpr std::uint32_t n = NIn;
        static constexpr std::uint32_t k = 1;
        static constexpr TFHEpp::ErrorDistribution errordist =
            TFHEpp::ErrorDistribution::ModularGaussian;
        static const inline double α = std::pow(2.0, AlphaLog2);
        using T = typename Lvl2::T;
        static constexpr std::make_signed_t<T> μ = Lvl2::μ;
        static constexpr std::uint32_t plain_modulus = Lvl2::plain_modulus;
        static constexpr double Δ =
            std::ldexp(1.0, std::numeric_limits<T>::digits) / plain_modulus;
    };

    template <std::uint32_t NOut, std::uint32_t Level,
              std::uint32_t Basebit, int AlphaLog2 = -45>
    struct MicroOutParam {
        static constexpr int32_t key_value_max = 1;
        static constexpr int32_t key_value_min = 0;
        static constexpr int32_t key_value_diff = key_value_max - key_value_min;
        static constexpr std::uint32_t nbit = MicroNBit<NOut>::value;
        static constexpr std::uint32_t n = NOut;
        static constexpr std::uint32_t k = 1;
        static constexpr std::uint32_t l = Level;
        static constexpr std::uint32_t lₐ = Level;
        static constexpr std::uint32_t Bgbit = Basebit;
        static constexpr std::uint32_t Bgₐbit = Basebit;
        static constexpr std::uint32_t Bg = 1U << Bgbit;
        static constexpr std::uint32_t Bgₐ = 1U << Bgₐbit;
        static constexpr TFHEpp::ErrorDistribution errordist =
            TFHEpp::ErrorDistribution::ModularGaussian;
        static const inline double α = std::pow(2.0, AlphaLog2);
        using T = typename Lvl1::T;
        static constexpr std::make_signed_t<T> μ = 1 << 29;
        static constexpr std::uint32_t plain_modulus = 8;
        static constexpr double Δ =
            static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
            plain_modulus;
    };

    template <class DomainP, class TargetP>
    struct MicroPBSParam {
        using domainP = DomainP;
        using targetP = TargetP;
        static constexpr std::uint32_t Addends = 1;
    };

    struct MicroLvl1DirectDomainParam {
        static constexpr int32_t key_value_max = Lvl1::key_value_max;
        static constexpr int32_t key_value_min = Lvl1::key_value_min;
        static constexpr int32_t key_value_diff = key_value_max - key_value_min;
        static constexpr std::uint32_t nbit = Lvl1::nbit;
        static constexpr std::uint32_t n = Lvl1::n;
        static constexpr std::uint32_t k = Lvl1::k;
        static constexpr std::uint32_t l = Lvl1::l;
        static constexpr std::uint32_t lₐ = Lvl1::lₐ;
        static constexpr std::uint32_t Bgbit = Lvl1::Bgbit;
        static constexpr std::uint32_t Bgₐbit = Lvl1::Bgₐbit;
        static constexpr std::uint32_t Bg = Lvl1::Bg;
        static constexpr std::uint32_t Bgₐ = Lvl1::Bgₐ;
        static constexpr TFHEpp::ErrorDistribution errordist = Lvl1::errordist;
        static const inline double α = Lvl1::α;
        using T = typename Lvl1::T;
        static constexpr std::make_signed_t<T> μ = Lvl1::μ;
        static constexpr std::uint32_t plain_modulus = Lvl1::plain_modulus;
        static constexpr double Δ = Lvl1::Δ;
    };

    template <class TargetP>
    struct MicroDirectPBSParam {
        using domainP = MicroLvl1DirectDomainParam;
        using targetP = TargetP;
        static constexpr std::uint32_t Addends = 1;
    };

    template <class DomainP, class TargetP, std::uint32_t Level,
              std::uint32_t Basebit, int AlphaLog2 = -45>
    struct MicroIKSParam {
        static constexpr std::uint32_t t = Level;
        static constexpr std::uint32_t basebit = Basebit;
        static constexpr TFHEpp::ErrorDistribution errordist =
            TFHEpp::ErrorDistribution::ModularGaussian;
        static const inline double α = std::pow(2.0, AlphaLog2);
        using domainP = DomainP;
        using targetP = TargetP;
    };

    template <class MicroInP, class MicroOutP, std::uint32_t IKSLevel = 7,
              std::uint32_t IKSBasebit = 4>
    struct MicroCandidateBase {
        using MicroIn = MicroInP;
        using MicroOut = MicroOutP;
        using PreIKS = MicroIKSParam<Lvl1, MicroIn, IKSLevel, IKSBasebit>;
        using PBS = MicroPBSParam<MicroIn, MicroOut>;
        using DirectPBS = MicroDirectPBSParam<MicroOut>;
        using PostIKS = MicroIKSParam<MicroOut, Lvl1, IKSLevel, IKSBasebit>;
        static constexpr std::uint32_t iks_level = IKSLevel;
        static constexpr std::uint32_t iks_basebit = IKSBasebit;
        static constexpr bool needs_post_iks =
            !(std::is_same_v<MicroOut, Lvl1> ||
              std::is_same_v<MicroOut, TFHEpp::AHlvl1param>);
    };

    template <class MicroInP, class MicroOutP, std::uint32_t IKSLevel = 10,
              std::uint32_t IKSBasebit = 3>
    struct MicroLvl2CandidateBase {
        using MicroIn = MicroInP;
        using MicroOut = MicroOutP;
        using PreIKS = MicroIKSParam<Lvl2, MicroIn, IKSLevel, IKSBasebit>;
        using PBS = MicroPBSParam<MicroIn, MicroOut>;
        using DirectPBS = MicroPBSParam<MicroIn, MicroOut>;
        using PostIKS = MicroIKSParam<MicroOut, Lvl2, IKSLevel, IKSBasebit>;
        static constexpr std::uint32_t iks_level = IKSLevel;
        static constexpr std::uint32_t iks_basebit = IKSBasebit;
        static constexpr bool needs_post_iks =
            !(std::is_same_v<MicroOut, Lvl2> ||
              std::is_same_v<MicroOut, TFHEpp::AHlvl2param>);
    };

    // Feasible fallback family for this checkout: TFHEpp's SPQLIOS FFT backend
    // only has uint32_t torus FFT specializations for lvl1param/AHlvl1param.
    // We still shrink the input LWE dimension aggressively, but the output
    // GLWE degree remains N=1024 unless vendor FFT registry code is extended.
    struct micro_n64_N1024_l2_b8
        : MicroCandidateBase<MicroInParam<64>, Lvl1> {
#ifdef USE_HE3DB_COMPAT
        static constexpr const char *name = "micro_n64_N1024_l3_b6";
#else
        static constexpr const char *name = "micro_n64_N1024_l2_b8";
#endif
    };

    struct micro_n96_N1024_l2_b8
        : MicroCandidateBase<MicroInParam<96>, Lvl1> {
#ifdef USE_HE3DB_COMPAT
        static constexpr const char *name = "micro_n96_N1024_l3_b6";
#else
        static constexpr const char *name = "micro_n96_N1024_l2_b8";
#endif
    };

    struct micro_n128_N1024_l2_b8
        : MicroCandidateBase<MicroInParam<128>, Lvl1> {
#ifdef USE_HE3DB_COMPAT
        static constexpr const char *name = "micro_n128_N1024_l3_b6";
#else
        static constexpr const char *name = "micro_n128_N1024_l2_b8";
#endif
    };

    struct micro_n192_N1024_l2_b8
        : MicroCandidateBase<MicroInParam<192>, Lvl1> {
#ifdef USE_HE3DB_COMPAT
        static constexpr const char *name = "micro_n192_N1024_l3_b6";
#else
        static constexpr const char *name = "micro_n192_N1024_l2_b8";
#endif
    };

    struct micro_n256_N1024_l2_b8
        : MicroCandidateBase<MicroInParam<256>, Lvl1> {
#ifdef USE_HE3DB_COMPAT
        static constexpr const char *name = "micro_n256_N1024_l3_b6";
#else
        static constexpr const char *name = "micro_n256_N1024_l2_b8";
#endif
    };

    struct micro_n64_N1024_ahlvl1_l4_b5
        : MicroCandidateBase<MicroInParam<64>, TFHEpp::AHlvl1param> {
        static constexpr const char *name = "micro_n64_N1024_ahlvl1_l4_b5";
    };

    struct micro2_n32_N2048_l4_b9
        : MicroLvl2CandidateBase<MicroIn64Param<32>, Lvl2> {
        static constexpr const char *name = "micro2_n32_N2048_l4_b9";
    };

    struct micro2_n64_N2048_l4_b9
        : MicroLvl2CandidateBase<MicroIn64Param<64>, Lvl2> {
        static constexpr const char *name = "micro2_n64_N2048_l4_b9";
    };

} // namespace tfhepp_compare::micro_pbs
