// TFHEpp/test/metapbs_truncrepeat_cipher_coeffs.cpp
// Minimal Unit F2 (fixed): cipher-domain truncRepeat via TLWE coefficient extraction + linear combination,
// verified by comparing decrypted PHASE high-bits to expected plaintext torus coefficient.
// (No TLWE->TRLWE packing yet.)

#include <tfhe++.hpp>

#include <iostream>
#include <vector>
#include <random>
#include <cstdint>
#include <cassert>
#include <limits>
#include <algorithm>

// ------------------ TruncRepeatPlan (same as F1) ------------------

static std::vector<int> B_sym(int B)
{
    int lo = -(B / 2);
    int hi = (B + 1) / 2 - 1;
    std::vector<int> v;
    for (int i = lo; i <= hi; i++) v.push_back(i);
    return v;
}

struct TruncRepeatPlan {
    int N, a, b, B;
    std::vector<std::vector<std::pair<int,int>>> out_terms; // (in_idx, sign)

    static int wrap_exp_to_index(int N, int64_t exp, int& sign)
    {
        int64_t e = exp % (2LL * N);
        if (e < 0) e += 2LL * N;
        if (e >= N) { sign = -1; return (int)(e - N); }
        sign = +1; return (int)e;
    }

    TruncRepeatPlan(int N_, int a_, int b_, int B_)
        : N(N_), a(a_), b(b_), B(B_), out_terms(N_)
    {
        auto sym = B_sym(B);

        for (int j = a; j <= b; j++) {
            int in_idx = (j >= 0) ? j : (N + j);
            int64_t base_exp = (j >= 0) ? (int64_t)j * B : (int64_t)N + (int64_t)j * B;

            for (int i : sym) {
                int sgn = +1;
                int out_idx = wrap_exp_to_index(N, base_exp + i, sgn);
                out_terms[out_idx].push_back({in_idx, sgn});
            }
        }
    }
};

// ------------------ TLWE helpers (generic add/sub) ------------------

template<class P>
static void tlwe_set_zero(TFHEpp::TLWE<P>& c)
{
    for (auto& x : c) x = 0;
}

template<class P>
static void tlwe_add_inplace(TFHEpp::TLWE<P>& acc, const TFHEpp::TLWE<P>& x)
{
    for (size_t i = 0; i < acc.size(); i++) acc[i] += x[i];
}

template<class P>
static void tlwe_sub_inplace(TFHEpp::TLWE<P>& acc, const TFHEpp::TLWE<P>& x)
{
    for (size_t i = 0; i < acc.size(); i++) acc[i] -= x[i];
}

template<class TorusT>
static inline uint64_t torus_to_u64(TorusT x)
{
    return (uint64_t)x; // works for uint32/uint64 torus
}

// Compare phase with expected torus coefficient by checking top bits.
template<class TorusT>
static inline bool equal_highbits(TorusT phase, TorusT expected, int drop_low_bits)
{
    // Compare (phase >> drop) == (expected >> drop)
    constexpr int digits = std::numeric_limits<TorusT>::digits;
    assert(drop_low_bits >= 0 && drop_low_bits < digits);
    return (torus_to_u64(phase) >> drop_low_bits) == (torus_to_u64(expected) >> drop_low_bits);
}

int main()
{
    using P = TFHEpp::lvl1param;  // accumulator ring level in PBS

    constexpr int N = P::n;
    constexpr int digits = std::numeric_limits<typename P::T>::digits;

    const int a = -3, b = 3;
    const int B = 4;
    assert(std::max(std::abs(a), std::abs(b)) * B < N);
    assert((b - a + 1) * B <= N);

    // How many low bits to ignore when comparing phase vs expected.
    // 16 is usually very safe; if you want stricter, try 12.
    const int DROP = 16;

    std::cout << "[metapbs_truncrepeat_cipher_coeffs] building keys...\n";
    TFHEpp::SecretKey sk;

    std::mt19937_64 rng(0);
    std::uniform_int_distribution<uint64_t> dist(0, (1ULL << 20) - 1); // small coefficients

    TruncRepeatPlan plan(N, a, b, B);

    for (int trial = 0; trial < 50; trial++) {
        // 1) Random plaintext polynomial in torus domain (coeffs are small ints embedded in torus)
        TFHEpp::Polynomial<P> M_poly{};
        for (int i = 0; i < N; i++) {
            M_poly[i] = (typename P::T)dist(rng);
        }

        // 2) Encrypt Cin = Enc(M_poly)
        TFHEpp::TRLWE<P> Cin{};
        TFHEpp::trlweSymEncrypt<P>(Cin, M_poly, sk.key.get<P>());

        // 3) For each out coefficient, compute TLWE(acc) = Σ ± SampleExtractIndex(Cin, in_idx)
        //    and compare its phase high-bits to expected plaintext coefficient.
        for (int out_idx = 0; out_idx < N; out_idx++) {
            TFHEpp::TLWE<P> acc{};
            tlwe_set_zero<P>(acc);

            // expected plaintext torus coefficient
            int64_t exp_acc = 0;

            for (auto [in_idx, sgn] : plan.out_terms[out_idx]) {
                TFHEpp::TLWE<P> tmp{};
                TFHEpp::SampleExtractIndex<P>(tmp, Cin, in_idx);

                if (sgn > 0) {
                    tlwe_add_inplace<P>(acc, tmp);
                    exp_acc += (int64_t)(uint64_t)M_poly[in_idx];
                } else {
                    tlwe_sub_inplace<P>(acc, tmp);
                    exp_acc -= (int64_t)(uint64_t)M_poly[in_idx];
                }
            }

            // Expected torus is exp_acc mod 2^digits
            typename P::T expected = (typename P::T)(uint64_t)exp_acc;

            // Phase = expected + noise
            auto phase = TFHEpp::tlweSymPhase<P>(acc, sk.key.get<P>());

            if (!equal_highbits<typename P::T>(phase, expected, DROP)) {
                std::cerr << "FAIL trial=" << trial << " out_idx=" << out_idx
                          << " DROP=" << DROP
                          << " phase=" << (uint64_t)phase
                          << " expected=" << (uint64_t)expected << "\n";
                // Optional: also show the highbits values
                std::cerr << " phase>>DROP=" << (torus_to_u64(phase) >> DROP)
                          << " expected>>DROP=" << (torus_to_u64(expected) >> DROP) << "\n";
                return 1;
            }
        }
    }

    std::cout << "PASS: truncRepeat coeff TLWE linear-combo matches expected plaintext (50 trials, DROP=" << DROP << ").\n";
    return 0;
}
