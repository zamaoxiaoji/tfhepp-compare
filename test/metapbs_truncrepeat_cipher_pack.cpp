// TFHEpp/test/metapbs_truncrepeat_cipher_pack.cpp
// Minimal Unit F3b (int-plaintext): build truncRepeat output coefficient TLWEs
// (via SampleExtract+plan), pack them to TRLWE, and verify extracted coeffs
// integer-decrypt to expected mod M exactly.

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <tfhe++.hpp>
#include <vector>
static inline uint32_t norm_mod_M(uint32_t x_u32, uint32_t M)
{
    int32_t x = (int32_t)x_u32;  // interpret as signed
    int32_t r = x % (int32_t)M;
    if (r < 0) r += (int32_t)M;
    return (uint32_t)r;
}

// ------------------ TruncRepeatPlan ------------------

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
    std::vector<std::vector<std::pair<int, int>>> out_terms;  // (in_idx, sign)

    static int wrap_exp_to_index(int N, int64_t exp, int& sign)
    {
        int64_t e = exp % (2LL * N);
        if (e < 0) e += 2LL * N;
        if (e >= N) {
            sign = -1;
            return (int)(e - N);
        }
        sign = +1;
        return (int)e;
    }

    TruncRepeatPlan(int N_, int a_, int b_, int B_)
        : N(N_), a(a_), b(b_), B(B_), out_terms(N_)
    {
        auto sym = B_sym(B);

        for (int j = a; j <= b; j++) {
            int in_idx = (j >= 0) ? j : (N + j);
            int64_t base_exp =
                (j >= 0) ? (int64_t)j * B : (int64_t)N + (int64_t)j * B;

            for (int i : sym) {
                int sgn = +1;
                int out_idx = wrap_exp_to_index(N, base_exp + i, sgn);
                out_terms[out_idx].push_back({in_idx, sgn});
            }
        }
    }
};

// ------------------ TLWE add/sub helpers ------------------

template <class P>
static void tlwe_set_zero(TFHEpp::TLWE<P>& c)
{
    for (auto& x : c) x = 0;
}

template <class P>
static void tlwe_add_inplace(TFHEpp::TLWE<P>& acc, const TFHEpp::TLWE<P>& x)
{
    for (size_t i = 0; i < acc.size(); i++) acc[i] += x[i];
}

template <class P>
static void tlwe_sub_inplace(TFHEpp::TLWE<P>& acc, const TFHEpp::TLWE<P>& x)
{
    for (size_t i = 0; i < acc.size(); i++) acc[i] -= x[i];
}

static inline uint32_t mod_int(int64_t x, uint32_t M)
{
    int64_t r = x % (int64_t)M;
    if (r < 0) r += M;
    return (uint32_t)r;
}

int main()
{
    using P = TFHEpp::lvl1param;
    constexpr int N = P::n;

    // integer plaintext modulus
    constexpr uint32_t M = 256;

    const int a = -3, b = 3;
    const int B = 4;
    assert(std::max(std::abs(a), std::abs(b)) * B < N);
    assert((b - a + 1) * B <= N);

    std::cout << "[metapbs_truncrepeat_cipher_pack:int] N=" << N << " M=" << M
              << "\n";
    std::cout << "[metapbs_truncrepeat_cipher_pack:int] building keys...\n";

    TFHEpp::SecretKey sk;

    // packing key
    TFHEpp::AnnihilateKey<P> ahk;
    TFHEpp::annihilatekeygen<P>(ahk, sk);

    TruncRepeatPlan plan(N, a, b, B);

    std::mt19937_64 rng(0);
    std::uniform_int_distribution<uint32_t> dist(0, M - 1);

    for (int trial = 0; trial < 5; trial++) {
        // 1) plaintext integer polynomial coefficients mod M
        std::vector<uint32_t> m_plain(N);
        std::vector<TFHEpp::TLWE<P>> in_coeffs(N);

        for (int i = 0; i < N; i++) {
            m_plain[i] = dist(rng);
            TFHEpp::tlweSymIntEncrypt<P, M>(in_coeffs[i],
                                            (typename P::T)m_plain[i], sk);
        }

        // 2) Pack input coeff-TLWEs to a TRLWE Cin
        TFHEpp::TRLWE<P> Cin{};
        TFHEpp::TLWE2TRLWEPacking<P>(Cin, in_coeffs, ahk);

        // 3) Build output coeff TLWEs by extraction+plan linear combination,
        // and expected integers
        std::vector<TFHEpp::TLWE<P>> out_coeffs(N);
        std::vector<uint32_t> expected(N);

        for (int out_idx = 0; out_idx < N; out_idx++) {
            TFHEpp::TLWE<P> acc{};
            tlwe_set_zero<P>(acc);

            int64_t exp_acc = 0;
            for (auto [in_idx, sgn] : plan.out_terms[out_idx]) {
                TFHEpp::TLWE<P> tmp{};
                TFHEpp::SampleExtractIndex<P>(tmp, Cin, in_idx);
                if (sgn > 0) {
                    tlwe_add_inplace<P>(acc, tmp);
                    exp_acc += (int64_t)m_plain[in_idx];
                }
                else {
                    tlwe_sub_inplace<P>(acc, tmp);
                    exp_acc -= (int64_t)m_plain[in_idx];
                }
            }
            out_coeffs[out_idx] = acc;
            expected[out_idx] = mod_int(exp_acc, M);
        }

        // 4) Pack output coeff TLWEs to TRLWE Cout
        TFHEpp::TRLWE<P> Cout{};
        TFHEpp::TLWE2TRLWEPacking<P>(Cout, out_coeffs, ahk);

        // 5) Verify by extracting and integer-decrypting each coeff
        for (int i = 0; i < N; i++) {
            TFHEpp::TLWE<P> ci{};
            TFHEpp::SampleExtractIndex<P>(ci, Cout, i);
            uint32_t got_raw =
                (uint32_t)TFHEpp::tlweSymIntDecrypt<P, M>(ci, sk);
            uint32_t got = norm_mod_M(got_raw, M);
            if (got != expected[i]) {
                std::cerr << "FAIL trial=" << trial << " coeff=" << i
                          << " got=" << got << " expect=" << expected[i]
                          << "\n";
                return 1;
            }
        }
    }

    std::cout << "PASS: truncRepeat coeff-TLWE -> packing -> TRLWE exact in "
                 "int domain (M="
              << M << ").\n";
    return 0;
}
