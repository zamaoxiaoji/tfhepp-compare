// TFHEpp/test/metapbs_k2.cpp
// Minimal Unit H: K=2 Meta-PBS end-to-end (sign LUT) by calling library metapbs::MetaPBS.

#include <tfhe++.hpp>
#include <metapbs.hpp>

#include <iostream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <limits>
#include <type_traits>

// sign LUT expectation
static inline int f_sign(int m, int t) { return (m < (t / 2)) ? +1 : -1; }

template <class T>
static inline int sign_of_torus(T x)
{
    using ST = std::make_signed_t<T>;
    return (static_cast<ST>(x) > 0) ? +1 : -1;
}

int main()
{
    using DomainP = TFHEpp::lvl0param;
    using TargetP = TFHEpp::lvl1param;
    using brP = TFHEpp::lvl01param;

    constexpr int64_t t = 8;

    // K=2 parameters
    std::vector<int> betas = {4, 4};
    std::vector<int> Ts    = {2, 2};

    std::cout << "[metapbs_k2] building keys...\n";
    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);

    TFHEpp::AnnihilateKey<TargetP> ahk;
    TFHEpp::annihilatekeygen<TargetP>(ahk, sk);

    // Build test vector using library helper to avoid any mismatch:
    // It expects f_half of size t/2 giving f(0..t/2-1).
    std::vector<int64_t> f_half((size_t)(t/2), +1); // sign LUT: +1 on first half, implicit -1 on second half
    const auto TV =
        TFHEpp::metapbs::GenerateTestVectorNegacyclic<TargetP, int64_t>(f_half,
                                                                       t);

    // input encode (q/t)*m in torus word
    constexpr int digits0 = std::numeric_limits<typename DomainP::T>::digits;
    const int logt = __builtin_ctz((unsigned)t);
    const typename DomainP::T scale_in = (typename DomainP::T)1 << (digits0 - logt);

    bool ok = true;

    for (int m = 0; m < (int)t; m++) {
        TFHEpp::TLWE<DomainP> cin{};
        TFHEpp::tlweSymEncrypt<DomainP>(cin, (typename DomainP::T)(m * scale_in), sk);

        TFHEpp::TLWE<TargetP> cout{};

        // Call library MetaPBS (Algorithm 1)
        TFHEpp::metapbs::MetaPBS<brP>(cout, cin, ek.getbkfft<brP>(), ahk, TV,
                                      betas, Ts, t);

        auto phase = TFHEpp::tlweSymPhase<TargetP>(cout, sk.key.get<TargetP>());
        int got = sign_of_torus<typename TargetP::T>(phase);
        int expect = f_sign(m, (int)t);

        if (got != expect) {
            ok = false;
            std::cerr << "FAIL m=" << m << " expect=" << expect << " got=" << got << "\n";
            break;
        }
    }

    if (ok) {
        std::cout << "PASS: K=2 metapbs::MetaPBS end-to-end (sign LUT) works.\n";
        return 0;
    } else {
        std::cout << "FAIL\n";
        return 1;
    }
}
