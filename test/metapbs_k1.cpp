// TFHEpp/test/metapbs_k1.cpp
// Minimal Unit G: K=1 Meta-PBS end-to-end in TFHEpp (adaptive to torus
// bit-width, no MulByXai dependency)
//
// Build & Run:
//   make -j
//   ./test/metapbs_k1

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <metapbs.hpp>
#include <tfhe++.hpp>
#include <type_traits>
#include <vector>

// ----------------- LUT: sign -----------------
static inline int f_sign(int m, int t) { return (m < (t / 2)) ? +1 : -1; }

// ----------------- TV construction (Algorithm 1 line 1 form) -----------------
template <class TargetP>
static inline void add_monomial_TV(TFHEpp::Polynomial<TargetP>& poly,
                                   typename TargetP::T coeff, int64_t exp)
{
    constexpr int64_t N = TargetP::n;
    int64_t e = exp % (2 * N);
    if (e < 0) e += 2 * N;
    if (e >= N)
        poly[e - N] -= coeff;  // X^{N+k}=-X^k
    else
        poly[e] += coeff;
}

template <class TargetP>
static TFHEpp::Polynomial<TargetP> build_tv_signlut(int t)
{
    constexpr int64_t N = TargetP::n;
    assert((2 * N) % t == 0);
    assert((N % t) == 0);
    assert((t & (t - 1)) == 0);  // power-of-two

    TFHEpp::Polynomial<TargetP> tv{};
    tv.fill(0);

    constexpr int digits = std::numeric_limits<typename TargetP::T>::digits;
    const int logt = __builtin_ctz((unsigned)t);
    // scale = 2^digits / t
    const typename TargetP::T scale = (typename TargetP::T)1 << (digits - logt);

    const int64_t r0 = (2 * N) / t;
    const int64_t shift0 = -(N / t);

    for (int i = 0; i < t / 2; i++) {
        const typename TargetP::T vm = scale;  // +1 half
        const int64_t base_exp = (2 * N * (int64_t)i) / t;
        for (int64_t j = 0; j < r0; j++) {
            add_monomial_TV<TargetP>(tv, vm, shift0 + base_exp + j);
        }
    }
    return tv;
}

// ----------------- centered mod -----------------
static inline int64_t centered_mod(int64_t x, int64_t mod)
{
    int64_t r = x % mod;
    if (r < 0) r += mod;
    if (r > mod / 2) r -= mod;
    return r;
}

// interpret torus word as signed int in Z
template <class TorusT>
static inline int64_t torus_to_signed_int(TorusT x)
{
    using ST = std::make_signed_t<TorusT>;
    return (int64_t)(ST)x;
}

// ----------------- HomDivRem on TLWE (integer domain) -----------------
template <class P>
static std::pair<std::vector<int64_t>, std::vector<int64_t>>
hom_div_rem_tlwe_int(const TFHEpp::TLWE<P>& c, int64_t qprime)
{
    constexpr int digits = std::numeric_limits<typename P::T>::digits;
    const int64_t q = (int64_t)1 << digits;  // 2^digits
    assert((q % qprime) == 0);
    const int64_t base = q / qprime;

    std::vector<int64_t> crem(c.size());
    std::vector<int64_t> cquo(c.size());

    for (size_t i = 0; i < c.size(); i++) {
        int64_t xi = torus_to_signed_int<typename P::T>(c[i]);
        int64_t ri = centered_mod(xi, base);
        crem[i] = ri;
        cquo[i] = (xi - ri) / base;
    }
    return {cquo, crem};
}

template <class P>
static TFHEpp::TLWE<P> tlwe_from_intvec_embedded(const std::vector<int64_t>& v,
                                                 int64_t embed)
{
    // embed should be q/qprime
    TFHEpp::TLWE<P> out{};
    assert(out.size() == v.size());

    using T = typename P::T;
    using ST = std::make_signed_t<T>;

    for (size_t i = 0; i < v.size(); i++) {
        // value in Z_{q'} -> embed into Z_q by multiplying embed factor
        int64_t x = v[i] * embed;
        out[i] = (T)(ST)x;  // two's complement cast
    }
    return out;
}

// ----------------- TruncRepeatPlan + Cipher implementation
// (extract+lincomb+packing) -----------------
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

// Δ_{B0,B} = min[BB0]_sym - B*min[B0]_sym - min[B]_sym
static inline int delta_B0_B(int B0, int B)
{
    auto min_sym = [](int X) { return -(X / 2); };
    return min_sym(B * B0) - B * min_sym(B0) - min_sym(B);
}

// Multiply TRLWE by X^k in negacyclic ring by coefficient rotation (no TFHEpp
// helper needed)
template <class P>
static TFHEpp::TRLWE<P> trlwe_mul_Xk(const TFHEpp::TRLWE<P>& in, int64_t k)
{
    constexpr int64_t N = P::n;
    // PolynomialMulByXai expects a in [0, 2N)
    int64_t a = k % (2 * N);
    if (a < 0) a += 2 * N;

    TFHEpp::TRLWE<P> out{};
    TFHEpp::PolynomialMulByXai<P>(out[0], in[0], (typename P::T)a);
    TFHEpp::PolynomialMulByXai<P>(out[1], in[1], (typename P::T)a);
    return out;
}

template <class P>
static TFHEpp::TRLWE<P> TruncRepeatCipher(const TFHEpp::TRLWE<P>& Cin, int T,
                                          int B,
                                          const TFHEpp::AnnihilateKey<P>& ahk)
{
    constexpr int N = P::n;
    const int a = -T, b = +T;
    TruncRepeatPlan plan(N, a, b, B);

    std::vector<TFHEpp::TLWE<P>> out_coeffs(N);

    for (int out_idx = 0; out_idx < N; out_idx++) {
        TFHEpp::TLWE<P> acc{};
        tlwe_set_zero<P>(acc);
        for (auto [in_idx, sgn] : plan.out_terms[out_idx]) {
            TFHEpp::TLWE<P> tmp{};
            TFHEpp::SampleExtractIndex<P>(tmp, Cin, in_idx);
            if (sgn > 0)
                tlwe_add_inplace<P>(acc, tmp);
            else
                tlwe_sub_inplace<P>(acc, tmp);
        }
        out_coeffs[out_idx] = acc;
    }

    TFHEpp::TRLWE<P> Cout{};
    TFHEpp::TLWE2TRLWEPacking<P>(Cout, out_coeffs, ahk);
    return Cout;
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
    using BRP = TFHEpp::lvl01param;

    constexpr int t = 8;

    // K=1 parameters
    const int beta0 = 4;
    const int T0 = 2;

    constexpr int N = TargetP::n;
    const int r0 = (2 * N) / t;

    std::cout << "[metapbs_k1] building keys...\n";
    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<BRP>(sk);

    TFHEpp::AnnihilateKey<TargetP> ahk;
    TFHEpp::annihilatekeygen<TargetP>(ahk, sk);

    const auto TV = build_tv_signlut<TargetP>(t);

    constexpr int digits0 = std::numeric_limits<typename DomainP::T>::digits;
    const int64_t q0 = (int64_t)1 << digits0;
    const int64_t embed0 = q0 / (2 * N);  // q/(2N)
    const int logt = __builtin_ctz((unsigned)t);
    const typename DomainP::T scale_in = (typename DomainP::T)1
                                         << (digits0 - logt);  // 2^digits / t

    bool ok = true;

    for (int m = 0; m < t; m++) {
        TFHEpp::TLWE<DomainP> cin{};
        const typename DomainP::T p = (typename DomainP::T)(m * scale_in);
        TFHEpp::tlweSymEncrypt<DomainP>(cin, p, sk);

        // HomDivRem(cin, 2N)
        auto [cquo0_int, crem0_int] = hom_div_rem_tlwe_int<DomainP>(cin, 2 * N);
        TFHEpp::TLWE<DomainP> cquo0 =
            tlwe_from_intvec_embedded<DomainP>(cquo0_int, embed0);

        // C0 = BlindRotate(cquo0, TV)
        TFHEpp::TRLWE<TargetP> C0{};
        TFHEpp::BlindRotate<BRP>(C0, cquo0, ek.getbkfft<BRP>(), TV);

        // C' = TruncRepeatCipher(C0, [-T0,T0], beta0) * X^{Δ}
        TFHEpp::TRLWE<TargetP> Cprime =
            TruncRepeatCipher<TargetP>(C0, T0, beta0, ahk);
        const int Delta = delta_B0_B(r0, beta0);
        Cprime = trlwe_mul_Xk<TargetP>(Cprime, Delta);

        // HomDivRem(crem0, beta0)
        TFHEpp::TLWE<DomainP> crem0 =
            tlwe_from_intvec_embedded<DomainP>(crem0_int, 1);
        auto [cquo1_int, crem1_int] =
            hom_div_rem_tlwe_int<DomainP>(crem0, beta0);
        const int64_t embed1 = q0 / beta0;
        TFHEpp::TLWE<DomainP> cquo1 =
            tlwe_from_intvec_embedded<DomainP>(cquo1_int, embed1);

        // C1 = BlindRotate(cquo1, C')
        TFHEpp::TRLWE<TargetP> C1{};
        TFHEpp::BlindRotate<BRP>(C1, cquo1, ek.getbkfft<BRP>(), Cprime);

        // SampleExtract
        TFHEpp::TLWE<TargetP> cout{};
        TFHEpp::SampleExtractIndex<TargetP>(cout, C1, 0);

        auto phase = TFHEpp::tlweSymPhase<TargetP>(cout, sk.key.get<TargetP>());
        int got = sign_of_torus<typename TargetP::T>(phase);
        int expect = f_sign(m, t);

        if (got != expect) {
            ok = false;
            std::cerr << "FAIL m=" << m << " expect=" << expect
                      << " got=" << got << "\n";
            break;
        }
    }

    if (ok) {
        std::cout << "PASS: K=1 Meta-PBS end-to-end (sign LUT) works.\n";
        return 0;
    }
    else {
        std::cout << "FAIL\n";
        return 1;
    }
}
