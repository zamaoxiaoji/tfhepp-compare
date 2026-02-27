// TFHEpp/test/metapbs_truncrepeat_plan.cpp
// Minimal Unit F1: build coefficient mapping plan for truncRepeat and test vs algebra reference.
// Build&Run: same as other tests: ./test/metapbs_truncrepeat_plan

#include <iostream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <random>

static inline int64_t mod_q(int64_t x, int64_t q)
{
    x %= q;
    if (x < 0) x += q;
    return x;
}

// Negacyclic multiply by X^k in R_q = Z_q[X]/(X^N+1)
static std::vector<int64_t> mul_Xk_negacyclic(const std::vector<int64_t>& p, int64_t k, int64_t q)
{
    const int64_t N = (int64_t)p.size();
    k %= (2 * N);
    if (k < 0) k += 2 * N;
    std::vector<int64_t> out(N, 0);
    for (int64_t i = 0; i < N; i++) {
        int64_t coeff = mod_q(p[i], q);
        if (coeff == 0) continue;
        int64_t deg = (i + k) % (2 * N);
        if (deg >= N) out[deg - N] = mod_q(out[deg - N] - coeff, q);
        else out[deg] = mod_q(out[deg] + coeff, q);
    }
    return out;
}

static std::vector<int> B_sym(int B)
{
    int lo = -(B / 2);
    int hi = (B + 1) / 2 - 1;
    std::vector<int> v;
    for (int i = lo; i <= hi; i++) v.push_back(i);
    return v;
}

// Lemma 5 exact truncPad
static std::vector<int64_t> truncPad(const std::vector<int64_t>& M, int a, int b, int B, int64_t q)
{
    const int N = (int)M.size();
    std::vector<int64_t> out(N, 0);
    std::vector<int64_t> base(N, 0);

    for (int j = a; j <= b; j++) {
        if (j >= 0) {
            int idx = j;
            int64_t coeff = mod_q(M[idx], q);
            if (!coeff) continue;
            base[0] = coeff;
            auto term = mul_Xk_negacyclic(base, (int64_t)j * B, q);
            for (int k = 0; k < N; k++) out[k] = mod_q(out[k] + term[k], q);
            base[0] = 0;
        } else {
            int idx = N + j;
            int64_t coeff = mod_q(M[idx], q);
            if (!coeff) continue;
            base[0] = coeff;
            auto term = mul_Xk_negacyclic(base, (int64_t)N + (int64_t)j * B, q);
            for (int k = 0; k < N; k++) out[k] = mod_q(out[k] + term[k], q);
            base[0] = 0;
        }
    }
    return out;
}

static std::vector<int64_t> truncRepeat_ref(const std::vector<int64_t>& M, int a, int b, int B, int64_t q)
{
    auto tp = truncPad(M, a, b, B, q);
    std::vector<int64_t> out(M.size(), 0);
    for (int i : B_sym(B)) {
        auto shifted = mul_Xk_negacyclic(tp, i, q);
        for (size_t k = 0; k < out.size(); k++) out[k] = mod_q(out[k] + shifted[k], q);
    }
    return out;
}

// ---------- TruncRepeatPlan ----------
// Plan produces mapping: each output coefficient index gets sum of some input coefficient indices.
struct TruncRepeatPlan {
    int N;
    int a, b;
    int B;
    std::vector<std::vector<std::pair<int,int>>> out_terms;
    // out_terms[out_idx] contains list of (in_idx, sign) where sign is +1 or -1.

    static int wrap_exp_to_index(int N, int64_t exp, int& sign)
    {
        // Map X^exp into coefficient index [0..N-1] with sign for exp>=N.
        int64_t e = exp % (2LL * N);
        if (e < 0) e += 2LL * N;
        if (e >= N) { sign = -1; return (int)(e - N); }
        sign = +1; return (int)e;
    }

    TruncRepeatPlan(int N_, int a_, int b_, int B_)
        : N(N_), a(a_), b(b_), B(B_), out_terms(N_)
    {
        // truncRepeat = truncPad * sum_{i in [B]_sym} X^i
        // truncPad places each input coeff:
        //  j>=0: M[j] at exp = jB
        //  j<0 : M[N+j] at exp = N + jB
        // then multiply by X^i: exp += i
        auto sym = B_sym(B);

        for (int j = a; j <= b; j++) {
            int in_idx = (j >= 0) ? j : (N + j);
            int64_t base_exp = (j >= 0) ? (int64_t)j * B : (int64_t)N + (int64_t)j * B;

            for (int i : sym) {
                int sgn1 = +1;
                int out_idx = wrap_exp_to_index(N, base_exp + i, sgn1);
                // term is (sgn1) * M[in_idx]
                out_terms[out_idx].push_back({in_idx, sgn1});
            }
        }
    }

    std::vector<int64_t> apply(const std::vector<int64_t>& M, int64_t q) const
    {
        std::vector<int64_t> out(N, 0);
        for (int out_idx = 0; out_idx < N; out_idx++) {
            int64_t acc = 0;
            for (auto [in_idx, sgn] : out_terms[out_idx]) {
                acc += (int64_t)sgn * mod_q(M[in_idx], q);
            }
            out[out_idx] = mod_q(acc, q);
        }
        return out;
    }
};

int main()
{
    const int64_t q = 12289;
    const int N = 64;
    const int a = -3, b = 3;
    const int B = 4;

    // constraints (for intended Lemma5/Def5 behavior in our setting)
    assert(std::max(std::abs(a), std::abs(b)) * B < N);
    assert((b - a + 1) * B <= N);

    std::mt19937_64 rng(0);
    std::uniform_int_distribution<int64_t> dist(0, q - 1);

    for (int trial = 0; trial < 200; trial++) {
        std::vector<int64_t> M(N);
        for (int i = 0; i < N; i++) M[i] = dist(rng);

        auto ref = truncRepeat_ref(M, a, b, B, q);

        TruncRepeatPlan plan(N, a, b, B);
        auto got = plan.apply(M, q);

        if (ref != got) {
            std::cerr << "FAIL at trial " << trial << "\n";
            // print first mismatch
            for (int i = 0; i < N; i++) {
                if (ref[i] != got[i]) {
                    std::cerr << "idx=" << i << " ref=" << ref[i] << " got=" << got[i] << "\n";
                    break;
                }
            }
            return 1;
        }
    }

    std::cout << "PASS: TruncRepeatPlan matches algebra reference for 200 random trials.\n";
    return 0;
}
