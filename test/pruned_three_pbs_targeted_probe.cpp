#include <cmath>
#include <iostream>
#include <limits>

#include "../comparison/comparison.h"

using namespace tfhepp_compare;

namespace
{
    template <class P>
    TFHEpp::TLWE<P> encrypt_message(typename P::T value, uint32_t plain_bits,
                                    const TFHESecretKey &sk)
    {
        const uint32_t scale_bits =
            std::numeric_limits<typename P::T>::digits - plain_bits - 1;
        return tlweSymInt32Encrypt<P>(value, P::α, std::pow(2., scale_bits),
                                      sk.key.get<P>());
    }

    template <class P>
    TFHEpp::TLWE<P> subtract(const TFHEpp::TLWE<P> &a,
                             const TFHEpp::TLWE<P> &b)
    {
        TFHEpp::TLWE<P> out;
        for (size_t i = 0; i <= P::k * P::n; i++) out[i] = a[i] - b[i];
        return out;
    }

    template <class P>
    TFHEpp::TLWE<P> trivial_difference(typename P::T p0, typename P::T p1,
                                       uint32_t plain_bits)
    {
        const uint32_t scale_bits =
            std::numeric_limits<typename P::T>::digits - plain_bits - 1;
        TFHEpp::TLWE<P> out = {};
        out[P::k * P::n] = (p0 - p1) * (typename P::T(1) << scale_bits);
        return out;
    }

    template <class P>
    bool check_decode(const char *label, const TFHEpp::TLWE<P> &sub,
                      uint32_t hom_bits, typename P::T expected,
                      TFHEEvalKey &ek, const TFHESecretKey &sk,
                      const three_pbs::FastB2AEvalKeyPack &micro_pack)
    {
        bool     ok = true;
        TLWELvl1 out;

        ethmsb::HomMSB(out, sub, hom_bits, ek, LOGIC);
        const auto ethmsb_dec =
            TFHEpp::tlweSymDecrypt<Lvl1>(out, sk.key.lvl1);
        ok = ok && ethmsb_dec == expected;

        three_pbs::HomMSB(out, sub, hom_bits, ek, LOGIC);
        const auto std_dec = TFHEpp::tlweSymDecrypt<Lvl1>(out, sk.key.lvl1);
        ok = ok && std_dec == expected;

        three_pbs::HomMSB(out, sub, hom_bits, ek, micro_pack, LOGIC);
        const auto opt_dec = TFHEpp::tlweSymDecrypt<Lvl1>(out, sk.key.lvl1);
        ok = ok && opt_dec == expected;

        std::cout << "    " << label << " ethmsb=" << ethmsb_dec
                  << " three_pbs_std=" << std_dec
                  << " three_pbs_opt=" << opt_dec << "\n";
        return ok;
    }

    template <class P>
    bool run_case(const char *name, typename P::T p0, typename P::T p1,
                  uint32_t plain_bits)
    {
        TFHESecretKey sk;
        TFHEEvalKey   ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplacebkfft<Lvl02>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        ek.emplaceiksk<Lvl20>(sk);
        ek.emplaceiksk<Lvl21>(sk);
        const auto micro_pack = three_pbs::GenerateFastB2AEvalKeyPack(sk, true);

        const uint32_t hom_bits = plain_bits + 1;
        const auto     expected_msb = typename P::T(p0 < p1);
        std::cout << name << " plain_bits=" << plain_bits << " p0=" << p0
                  << " p1=" << p1 << " expected_msb=" << expected_msb
                  << "\n";

        const auto trivial = trivial_difference<P>(p0, p1, plain_bits);
        bool ok = check_decode<P>("trivial ", trivial, hom_bits, expected_msb,
                                  ek, sk, micro_pack);

        const auto c0 = encrypt_message<P>(p0, plain_bits, sk);
        const auto c1 = encrypt_message<P>(p1, plain_bits, sk);
        const auto encrypted = subtract<P>(c0, c1);
        ok = check_decode<P>("encrypted", encrypted, hom_bits, expected_msb, ek,
                             sk, micro_pack) &&
             ok;
        return ok;
    }
} // namespace

int main()
{
    bool ok = true;
    ok = run_case<Lvl2>("case16", 10012, 10661, 16) && ok;
    ok = run_case<Lvl2>("case32", 1693249469ULL, 1706842735ULL, 32) && ok;
    std::cout << "EDGE_REGRESSION " << (ok ? "ok" : "failed") << "\n";
    return ok ? 0 : 1;
}
