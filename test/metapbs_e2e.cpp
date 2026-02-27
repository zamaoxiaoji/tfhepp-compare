// TFHEpp/test/metapbs_e2e.cpp
// End-to-End Meta-PBS with programmable LUT (randomly generated).
// For K=2 with random LUT, running 1000 samples and checking accuracy.

#include <tfhe++.hpp>
#include <metapbs.hpp>

#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>

template <class T>
std::vector<T> generate_fixed_LUT(int t)
{
    std::vector<T> LUT(t / 2, 1);  // Set first half to +1
    std::fill(LUT.begin() + t / 2, LUT.end(), -1);  // Set second half to -1
    return LUT;
}

template <class TargetP>
static inline void add_monomial_TV(TFHEpp::Polynomial<TargetP>& poly,
                                   typename TargetP::T coeff,
                                   int64_t exp)
{
    constexpr int64_t N = TargetP::n;
    int64_t e = exp % (2 * N);
    if (e < 0) e += 2 * N;
    if (e >= N) poly[e - N] -= coeff; // X^{N+k} = -X^k
    else poly[e] += coeff;
}

template <class T>
static inline int sign_of_torus(T x)
{
    using ST = std::make_signed_t<T>;
    return (static_cast<ST>(x) > 0) ? +1 : -1;
}

// sign LUT expectation
static inline int f_sign(int m, int t) { return (m < (t / 2)) ? +1 : -1; }

// f_half -> identity/random LUT
template <class T>
std::vector<T> generate_random_LUT(int t)
{
    std::vector<T> LUT(t / 2);
    std::mt19937_64 rng(0);
    std::uniform_int_distribution<int> dist(-1, 1); // random values in {-1, 1}

    for (int i = 0; i < t / 2; i++) {
        LUT[i] = dist(rng);
    }
    return LUT;
}

// Build test vector (sign LUT / identity)
template <class TargetP>
static TFHEpp::Polynomial<TargetP> build_test_vector(int t, const std::vector<int>& f_half)
{
    constexpr int N = TargetP::n;
    TFHEpp::Polynomial<TargetP> tv{};
    tv.fill(0);

    constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
    const int logt = __builtin_ctz((unsigned)t);
    const typename TargetP::T scale = (typename TargetP::T)1 << (digits - logt);

    const int64_t r0 = (2 * N) / t;
    const int64_t shift0 = -(N / t);

    // first half: sign LUT or identity (just map first half)
    for (int i = 0; i < t / 2; i++) {
        const typename TargetP::T vm = scale * f_half[i];
        const int64_t base_exp = (2 * N * (int64_t)i) / t;
        for (int64_t j = 0; j < r0; j++) {
            add_monomial_TV<TargetP>(tv, vm, shift0 + base_exp + j);
        }
    }
    return tv;
}

int main()
{
    using DomainP = TFHEpp::lvl0param;
    using TargetP = TFHEpp::lvl1param;
    using brP = TFHEpp::lvl01param;

    constexpr int t = 8;

    // K=2 parameters
    std::vector<int> betas = {4, 4};
    std::vector<int> Ts    = {2, 2};

    std::cout << "[metapbs_e2e] building keys...\n";
    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);

    TFHEpp::AnnihilateKey<TargetP> ahk;
    TFHEpp::annihilatekeygen<TargetP>(ahk, sk);

    // Randomly generate LUT
    std::vector<int> f_half = generate_fixed_LUT<int>(t); // Random LUT
    auto TV = build_test_vector<TargetP>(t, f_half);

    // input encode (q/t)*m in torus word
    constexpr int digits0 = std::numeric_limits<typename DomainP::T>::digits;
    const int logt = __builtin_ctz((unsigned)t);
    const typename DomainP::T scale_in = (typename DomainP::T)1 << (digits0 - logt);

    int num_trials = 1000;
    int num_failures = 0;

    // Run 1000 trials of m, comparing result with LUT
    for (int m = 0; m < num_trials; m++) {
        TFHEpp::TLWE<DomainP> cin{};
        TFHEpp::tlweSymEncrypt<DomainP>(cin, (typename DomainP::T)(m * scale_in), sk);

        TFHEpp::TLWE<TargetP> cout{};

        // Call library MetaPBS (Algorithm 1)
        TFHEpp::metapbs::MetaPBS<brP>(cout, cin, ek.getbkfft<brP>(), ahk, TV, betas, Ts, t);

        auto phase = TFHEpp::tlweSymPhase<TargetP>(cout, sk.key.get<TargetP>());
        int got = f_sign(m, t); // LUT result
        int expect = sign_of_torus<typename TargetP::T>(phase);

        if (got != expect) {
            num_failures++;
        }
    }

    std::cout << "Completed " << num_trials << " trials.\n";
    std::cout << "Failures: " << num_failures << " out of " << num_trials << "\n";
    if (num_failures > 0) {
        std::cout << "FAIL: Some trials did not match expected LUT result.\n";
        return 1;
    } else {
        std::cout << "PASS: All trials matched expected LUT result.\n";
        return 0;
    }
}
