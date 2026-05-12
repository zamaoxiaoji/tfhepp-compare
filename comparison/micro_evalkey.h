#pragma once
/**
 * @file micro_evalkey.h
 * @brief Local key-pack and conversion helpers for unsafe micro-parameter PBS.
 */

#include <chrono>
#include <memory>
#include <random>

#include "micro_params.h"

namespace tfhepp_compare::micro_pbs
{
    struct MicroTiming {
        double t_pre_ks_ms = 0.0;
        double t_micro_pbs_ms = 0.0;
        double t_post_ks_ms = 0.0;

        [[nodiscard]] double total_ms() const
        {
            return t_pre_ks_ms + t_micro_pbs_ms + t_post_ks_ms;
        }
    };

    template <class Candidate>
    struct MicroEvalKeyPack {
        using MicroIn = typename Candidate::MicroIn;
        using MicroOut = typename Candidate::MicroOut;
        using PreIKS = typename Candidate::PreIKS;
        using PBS = typename Candidate::PBS;
        using DirectPBS = typename Candidate::DirectPBS;
        using PostIKS = typename Candidate::PostIKS;

        TFHEpp::Key<MicroIn> microin_key = {};
        TFHEpp::Key<MicroOut> microout_key = {};

        std::unique_ptr<TFHEpp::KeySwitchingKey<PreIKS>>
            iksk_lvl1_to_microin;
        std::unique_ptr<TFHEpp::BootstrappingKeyFFT<PBS>>
            bkfft_microin_to_microout;
        std::unique_ptr<TFHEpp::KeySwitchingKey<PostIKS>>
            iksk_microout_to_lvl1;

        std::unique_ptr<TFHEpp::BootstrappingKeyFFT<DirectPBS>>
            bkfft_lvl1_to_microout;
    };

    template <class P>
    TFHEpp::Key<P> GenerateBinaryKey()
    {
        TFHEpp::Key<P> key;
        std::uniform_int_distribution<int32_t> dist(P::key_value_min,
                                                    P::key_value_max);
        for (auto &x : key) x = static_cast<typename P::T>(dist(TFHEpp::generator));
        return key;
    }

    template <class P>
    TFHEpp::Polynomial<P> MinusAPolygen(typename P::T A)
    {
        TFHEpp::Polynomial<P> poly;
        for (typename P::T &p : poly) p = -A;
        return poly;
    }

    template <class P>
    constexpr bool SameKeyShapeAsLvl1()
    {
        return std::is_same_v<typename P::T, typename Lvl1::T> &&
               P::k == Lvl1::k && P::n == Lvl1::n;
    }

    template <class P>
    constexpr bool SameKeyShapeAsLvl2()
    {
        return std::is_same_v<typename P::T, typename Lvl2::T> &&
               P::k == Lvl2::k && P::n == Lvl2::n;
    }

    template <class Candidate>
    MicroEvalKeyPack<Candidate> GenerateMicroEvalKeyPack(const TFHESecretKey &sk,
                                                         bool with_direct)
    {
        MicroEvalKeyPack<Candidate> pack;
        pack.microin_key = GenerateBinaryKey<typename Candidate::MicroIn>();
        if constexpr (std::is_same_v<typename Candidate::MicroOut, Lvl1> ||
                      std::is_same_v<typename Candidate::MicroOut,
                                     TFHEpp::AHlvl1param>)
            pack.microout_key = sk.key.lvl1;
        else if constexpr (std::is_same_v<typename Candidate::MicroOut,
                                          Lvl2> ||
                           std::is_same_v<typename Candidate::MicroOut,
                                          TFHEpp::AHlvl2param>)
            pack.microout_key = sk.key.lvl2;
        else
            pack.microout_key =
                GenerateBinaryKey<typename Candidate::MicroOut>();

        pack.iksk_lvl1_to_microin =
            std::make_unique_for_overwrite<
                TFHEpp::KeySwitchingKey<typename Candidate::PreIKS>>();
        if constexpr (SameKeyShapeAsLvl1<typename Candidate::PreIKS::domainP>())
            TFHEpp::ikskgen<typename Candidate::PreIKS>(
                *pack.iksk_lvl1_to_microin, sk.key.lvl1, pack.microin_key);
        else if constexpr (SameKeyShapeAsLvl2<
                               typename Candidate::PreIKS::domainP>())
            TFHEpp::ikskgen<typename Candidate::PreIKS>(
                *pack.iksk_lvl1_to_microin, sk.key.lvl2, pack.microin_key);
        else
            static_assert(SameKeyShapeAsLvl1<
                              typename Candidate::PreIKS::domainP>() ||
                              SameKeyShapeAsLvl2<
                                  typename Candidate::PreIKS::domainP>(),
                          "unsupported micro pre-IKS domain key");

        pack.bkfft_microin_to_microout =
            std::make_unique_for_overwrite<
                TFHEpp::BootstrappingKeyFFT<typename Candidate::PBS>>();
        TFHEpp::bkfftgen<typename Candidate::PBS>(
            *pack.bkfft_microin_to_microout, pack.microin_key,
            pack.microout_key);

        if constexpr (Candidate::needs_post_iks) {
            pack.iksk_microout_to_lvl1 =
                std::make_unique_for_overwrite<
                    TFHEpp::KeySwitchingKey<typename Candidate::PostIKS>>();
            TFHEpp::ikskgen<typename Candidate::PostIKS>(
                *pack.iksk_microout_to_lvl1, pack.microout_key, sk.key.lvl1);
        }

        if (with_direct) {
            pack.bkfft_lvl1_to_microout =
                std::make_unique_for_overwrite<
                    TFHEpp::BootstrappingKeyFFT<typename Candidate::DirectPBS>>();
            if constexpr (SameKeyShapeAsLvl1<
                              typename Candidate::DirectPBS::domainP>())
                TFHEpp::bkfftgen<typename Candidate::DirectPBS>(
                    *pack.bkfft_lvl1_to_microout, sk.key.lvl1,
                    pack.microout_key);
            else if constexpr (SameKeyShapeAsLvl2<
                                   typename Candidate::DirectPBS::domainP>())
                TFHEpp::bkfftgen<typename Candidate::DirectPBS>(
                    *pack.bkfft_lvl1_to_microout, sk.key.lvl2,
                    pack.microout_key);
            else
                TFHEpp::bkfftgen<typename Candidate::DirectPBS>(
                    *pack.bkfft_lvl1_to_microout, pack.microin_key,
                    pack.microout_key);
        }
        return pack;
    }

    template <class Clock = std::chrono::steady_clock>
    inline double elapsed_ms(typename Clock::time_point a,
                             typename Clock::time_point b)
    {
        return std::chrono::duration<double, std::milli>(b - a).count();
    }

    template <class Candidate>
    void MicroUnsafeKsPreQHalfToGuardValue_Lvl1(
        TLWELvl1 &ct_guard, const TLWELvl1 &bit_qhalf, Lvl1::T guard_value,
        Lvl1::T offset, const MicroEvalKeyPack<Candidate> &pack,
        MicroTiming *timing = nullptr)
    {
        using MicroIn = typename Candidate::MicroIn;
        using MicroOut = typename Candidate::MicroOut;
        static_assert(std::is_same_v<typename MicroOut::T, Lvl1::T>,
                      "micro output torus must match Lvl1");
        const typename MicroOut::T A = guard_value >> 1;

        const auto t0 = std::chrono::steady_clock::now();
        TFHEpp::TLWE<MicroIn> ct_microin;
        TFHEpp::IdentityKeySwitch<typename Candidate::PreIKS>(
            ct_microin, bit_qhalf, *pack.iksk_lvl1_to_microin);
        ct_microin[MicroIn::k * MicroIn::n] += offset;
        const auto t1 = std::chrono::steady_clock::now();

        TFHEpp::TLWE<MicroOut> ct_microout;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<typename Candidate::PBS>(
            ct_microout, ct_microin, *pack.bkfft_microin_to_microout,
            MinusAPolygen<MicroOut>(A));
        ct_microout[MicroOut::k * MicroOut::n] += A;
        const auto t2 = std::chrono::steady_clock::now();

        if constexpr (Candidate::needs_post_iks)
            TFHEpp::IdentityKeySwitch<typename Candidate::PostIKS>(
                ct_guard, ct_microout, *pack.iksk_microout_to_lvl1);
        else
            ct_guard = ct_microout;
        const auto t3 = std::chrono::steady_clock::now();

        if (timing != nullptr) {
            timing->t_pre_ks_ms += elapsed_ms(t0, t1);
            timing->t_micro_pbs_ms += elapsed_ms(t1, t2);
            timing->t_post_ks_ms += elapsed_ms(t2, t3);
        }
    }

    template <class Candidate>
    void MicroUnsafeDirectQHalfToGuardValue_Lvl1(
        TLWELvl1 &ct_guard, const TLWELvl1 &bit_qhalf, Lvl1::T guard_value,
        Lvl1::T offset, const MicroEvalKeyPack<Candidate> &pack,
        MicroTiming *timing = nullptr)
    {
        using MicroOut = typename Candidate::MicroOut;
        static_assert(std::is_same_v<typename MicroOut::T, Lvl1::T>,
                      "micro output torus must match Lvl1");
        const typename MicroOut::T A = guard_value >> 1;

        const auto t0 = std::chrono::steady_clock::now();
        TLWELvl1 ct_offset = bit_qhalf;
        ct_offset[Lvl1::k * Lvl1::n] += offset;
        const auto t1 = std::chrono::steady_clock::now();

        TFHEpp::TLWE<MicroOut> ct_microout;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<typename Candidate::DirectPBS>(
            ct_microout, ct_offset, *pack.bkfft_lvl1_to_microout,
            MinusAPolygen<MicroOut>(A));
        ct_microout[MicroOut::k * MicroOut::n] += A;
        const auto t2 = std::chrono::steady_clock::now();

        if constexpr (Candidate::needs_post_iks)
            TFHEpp::IdentityKeySwitch<typename Candidate::PostIKS>(
                ct_guard, ct_microout, *pack.iksk_microout_to_lvl1);
        else
            ct_guard = ct_microout;
        const auto t3 = std::chrono::steady_clock::now();

        if (timing != nullptr) {
            timing->t_pre_ks_ms += elapsed_ms(t0, t1);
            timing->t_micro_pbs_ms += elapsed_ms(t1, t2);
            timing->t_post_ks_ms += elapsed_ms(t2, t3);
        }
    }

    template <class Candidate>
    void MicroUnsafeKsPreQHalfToGuardValue_Lvl2(
        TLWELvl2 &ct_guard, const TLWELvl2 &bit_qhalf, Lvl2::T guard_value,
        Lvl2::T offset, const MicroEvalKeyPack<Candidate> &pack,
        MicroTiming *timing = nullptr)
    {
        using MicroIn = typename Candidate::MicroIn;
        using MicroOut = typename Candidate::MicroOut;
        static_assert(std::is_same_v<typename MicroOut::T, Lvl2::T>,
                      "micro output torus must match Lvl2");
        const typename MicroOut::T A = guard_value >> 1;

        const auto t0 = std::chrono::steady_clock::now();
        TFHEpp::TLWE<MicroIn> ct_microin;
        TFHEpp::IdentityKeySwitch<typename Candidate::PreIKS>(
            ct_microin, bit_qhalf, *pack.iksk_lvl1_to_microin);
        ct_microin[MicroIn::k * MicroIn::n] += offset;
        const auto t1 = std::chrono::steady_clock::now();

        TFHEpp::TLWE<MicroOut> ct_microout;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<typename Candidate::PBS>(
            ct_microout, ct_microin, *pack.bkfft_microin_to_microout,
            MinusAPolygen<MicroOut>(A));
        ct_microout[MicroOut::k * MicroOut::n] += A;
        const auto t2 = std::chrono::steady_clock::now();

        if constexpr (Candidate::needs_post_iks)
            TFHEpp::IdentityKeySwitch<typename Candidate::PostIKS>(
                ct_guard, ct_microout, *pack.iksk_microout_to_lvl1);
        else
            ct_guard = ct_microout;
        const auto t3 = std::chrono::steady_clock::now();

        if (timing != nullptr) {
            timing->t_pre_ks_ms += elapsed_ms(t0, t1);
            timing->t_micro_pbs_ms += elapsed_ms(t1, t2);
            timing->t_post_ks_ms += elapsed_ms(t2, t3);
        }
    }

} // namespace tfhepp_compare::micro_pbs
