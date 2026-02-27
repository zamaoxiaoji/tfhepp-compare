// TFHEpp/test/metapbs_k0.cpp
// Minimal unit: K=0 PBS path in TFHEpp using Meta-PBS test vector (Algorithm 1 line 1)
// - Build TV polynomial (LUT) in R_q = Torus[X]/(X^N+1)
// - Encrypt message m in Z_t into TLWE(lvl0)
// - Run GateBootstrappingTLWE2TLWE<lvl01param>(...) with our TV
//   (internally: BlindRotate + SampleExtractIndex)
// - Decrypt TLWE(lvl1) and check sign matches a negacyclic sign LUT
//
// Build&Run: (from TFHEpp/build)
//   cmake .. -DENABLE_TEST=ON
//   make -j --target metapbs_k0
//   ./test/metapbs_k0

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>

#include <tfhe++.hpp>

// ---------- Helpers for Torus polynomial in negacyclic ring ----------
template <class P>
static inline void add_monomial_TV(TFHEpp::Polynomial<P> &poly,
                                  typename P::T coeff,
                                  int64_t exp)
{
    // Work in Z_{2N} exponents then map into ring degree [0..N-1] with sign flip
    // for exp>=N (X^{N+k} == -X^k).
    constexpr int64_t N = P::n;
    int64_t e = exp % (2 * N);
    if (e < 0) e += 2 * N;

    if (e >= N) {
        poly[e - N] -= coeff;
    }
    else {
        poly[e] += coeff;
    }
}

// f(m): +1 for m in [0, t/2), -1 otherwise. (negacyclic)
static inline int f_sign(int m, int t) { return (m < (t / 2)) ? +1 : -1; }

// Build TV exactly like Algorithm 1 line 1 (Meta-PBS paper):
// TV = X^{-N/t} * sum_{i=0}^{t/2-1} ( (q/t) f(i) * X^{2N i/t} * sum_{j=0}^{2N/t-1} X^j )
// Here q is implicit Torus modulus 2^{digits}. We realize (q/t) as 2^{digits}/t.
template <class TargetP>
static TFHEpp::Polynomial<TargetP> build_test_vector_tv(int t)
{
    constexpr int64_t N = TargetP::n;
    assert((2 * N) % t == 0);
    assert((N % t) == 0);  // so -N/t is integer

    TFHEpp::Polynomial<TargetP> tv{};
    tv.fill(0);

    constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
    // scale = q/t = 2^digits / t (require t is power of 2 for exact shift)
    assert((t & (t - 1)) == 0);
    const int logt = __builtin_ctz((unsigned)t);
    const typename TargetP::T scale =
        (typename TargetP::T)1 << (digits - logt);  // 2^{digits-log2(t)}

    const int64_t r0 = (2 * N) / t;   // 2N/t
    const int64_t shift0 = -(N / t);  // -N/t

    for (int i = 0; i < t / 2; i++) {
        const int s = f_sign(i, t);
        const typename TargetP::T vm =
            (s > 0) ? scale : (typename TargetP::T)(-scale);

        const int64_t base_exp = (2 * N * (int64_t)i) / t;
        for (int64_t j = 0; j < r0; j++) {
            const int64_t exp = shift0 + base_exp + j;
            add_monomial_TV<TargetP>(tv, vm, exp);
        }
    }
    return tv;
}

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
    using BRP = TFHEpp::lvl01param;  // domain=lvl0, target=lvl1

    constexpr int t = 8;  // small, power-of-two (required by build_test_vector_tv)
    std::cout << "[metapbs_k0] building keys...\n";

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<BRP>(sk);

    // Build TV polynomial in lvl1 ring.
    const auto tv = build_test_vector_tv<TargetP>(t);

    // Torus value for m/t when t=2^k is: m * 2^{digits-k}.
    constexpr int digits = std::numeric_limits<typename DomainP::T>::digits;
    const int logt = __builtin_ctz((unsigned)t);
    const typename DomainP::T scale_in =
        (typename DomainP::T)1 << (digits - logt);

    for (int m = 0; m < t; m++) {
        TFHEpp::TLWE<DomainP> cin{};
        const typename DomainP::T p = (typename DomainP::T)(m * scale_in);
        TFHEpp::tlweSymEncrypt<DomainP>(cin, p, sk);

        TFHEpp::TLWE<TargetP> cout{};
        TFHEpp::GateBootstrappingTLWE2TLWE<BRP>(cout, cin, ek.getbkfft<BRP>(),
                                                tv);

        const auto phase =
            TFHEpp::tlweSymPhase<TargetP>(cout, sk.key.get<TargetP>());
        const int got = sign_of_torus<typename TargetP::T>(phase);
        const int expect = f_sign(m, t);
        if (got != expect) {
            std::cerr << "FAIL at m=" << m << " expect_sign=" << expect
                      << " got_sign=" << got << "\n";
            return 1;
        }
    }

    std::cout << "PASS: K=0 PBS path with Meta-PBS TV (sign LUT) works.\n";
    return 0;
}

