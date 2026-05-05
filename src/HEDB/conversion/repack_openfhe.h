#pragma once
// TFHE→CKKS Repack using OpenFHE
// Converts TLWE lvl1 ciphertexts to CKKS ciphertext without decryption.
//
// Algorithm (following OpenFHE's EvalFHEWtoCKKS / Bossuat et al.):
// 1. KeyGen: Encrypt TFHE sk in CKKS → Enc_CKKS(s_tfhe)
// 2. Scale A and b by prescale = 1/(q*K) where q=2^32, K=128
//    This maps the real-valued inner product into [-1, 1]
// 3. LinearTransform: BSGS diagonal method → Enc_CKKS(A·s * prescale)
// 4. Add b * prescale → Enc_CKKS(phase/(q*K))
// 5. HomMod: Chebyshev sin approximation to extract fractional part
// 6. Double-angle iterations to sharpen the approximation
// 7. Post-scale by 2π (for binary p≤4) to recover message

#include "openfhe.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <cloudkey.hpp>
#include "HEDB/comparison/tfhepp_utils.h"

namespace HEDB
{
    using namespace lbcrypto;

    template <typename T>
    static inline T CeilSqrt(T val) {
        return static_cast<T>(std::ceil(std::sqrt(1.0 * val)));
    }
    template <typename T>
    static inline T CeilDiv(T a, T b) {
        return (a + b - 1) / b;
    }

    // --- Repack pre-key: CKKS encryptions of TFHE secret key ---
    struct RepackKey {
        std::vector<Ciphertext<DCRTPoly>> rotated_sk;
    };

    // === KeyGen ===
    // Encrypts TFHE lvl1 secret key into CKKS with BSGS baby-step rotations
    inline void RepackKeyGen(
        RepackKey &rk,
        const CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        const TFHESecretKey &tfhe_sk,
        size_t tfhe_n)
    {
        std::vector<double> sk_vec(tfhe_n);
        for (size_t i = 0; i < tfhe_n; i++) {
            auto val = tfhe_sk.key.get<Lvl1>()[i];
            sk_vec[i] = (val > 1) ? -1.0 : static_cast<double>(val);
        }

        size_t g = CeilSqrt(tfhe_n);
        rk.rotated_sk.resize(g);

        for (size_t j = 0; j < g; j++) {
            auto pt = cc->MakeCKKSPackedPlaintext(sk_vec);
            rk.rotated_sk[j] = cc->Encrypt(keys.publicKey, pt);
            std::rotate(sk_vec.begin(), sk_vec.begin() + 1, sk_vec.end());
        }
    }

    // === BSGS Linear Transform ===
    // Computes Enc(A · s) where A is (num_lwes × tfhe_n) plaintext matrix
    inline Ciphertext<DCRTPoly> LinearTransformBSGS(
        const CryptoContext<DCRTPoly> &cc,
        const std::vector<std::vector<double>> &A,
        const RepackKey &rk,
        size_t tfhe_n)
    {
        size_t rows = A.size();
        size_t cols = tfhe_n;
        size_t max_len = std::max(rows, cols);
        size_t min_len = std::min(rows, cols);
        size_t g_tilde = CeilSqrt(min_len);
        size_t b_tilde = CeilDiv(min_len, g_tilde);

        Ciphertext<DCRTPoly> result;
        bool result_init = false;

        for (size_t b = 0; b < b_tilde && g_tilde * b < min_len; b++) {
            Ciphertext<DCRTPoly> sum;
            bool sum_init = false;

            for (size_t g = 0; g < g_tilde && b * g_tilde + g < min_len; g++) {
                size_t j = b * g_tilde + g;
                std::vector<double> diag(max_len);
                for (size_t r = 0; r < max_len; r++) {
                    diag[r] = A[r % rows][(r + j) % cols];
                }
                std::rotate(diag.rbegin(), diag.rbegin() + b * g_tilde, diag.rend());

                auto pt_diag = cc->MakeCKKSPackedPlaintext(diag);
                auto term = cc->EvalMult(rk.rotated_sk[g], pt_diag);

                if (!sum_init) { sum = term; sum_init = true; }
                else { sum = cc->EvalAdd(sum, term); }
            }

            if (b > 0) {
                sum = cc->EvalRotate(sum, (int32_t)(b * g_tilde));
            }

            if (!result_init) { result = sum; result_init = true; }
            else { result = cc->EvalAdd(result, sum); }
        }

        // Fold if rows < cols: rotate-and-add for power-of-2 folding
        if (rows < cols) {
            size_t gama = static_cast<size_t>(std::log2(1.0 * cols / rows));
            for (size_t j = 0; j < gama; j++) {
                auto temp = cc->EvalRotate(result, (int32_t)((1U << j) * rows));
                result = cc->EvalAdd(result, temp);
            }
        }

        return result;
    }

    // === HomRound ===
    // Maps values near 0 → 0 and near 1 → 1
    // Uses iterated sign polynomial: f(x) = 1.5x - 0.5x³
    // on the range [-1,1] (after mapping [0,1] → [-1,1])
    inline Ciphertext<DCRTPoly> HomRound(
        const CryptoContext<DCRTPoly> &cc,
        Ciphertext<DCRTPoly> ct)
    {
        // Map [0,1] → [-1,1]: y = 2x - 1
        ct = cc->EvalMult(ct, 2.0);
        ct = cc->EvalAdd(ct, -1.0);

        // Iterate sign approximation: g(y) = 1.5y - 0.5y³
        for (int iter = 0; iter < 2; iter++) {
            auto y2 = cc->EvalMult(ct, ct);    // y²
            auto y3 = cc->EvalMult(y2, ct);    // y³
            auto t1 = cc->EvalMult(ct, 1.5);   // 1.5y
            auto t2 = cc->EvalMult(y3, -0.5);  // -0.5y³
            ct = cc->EvalAdd(t1, t2);
        }

        // Map [-1,1] → [0,1]: x = (y+1)/2
        ct = cc->EvalAdd(ct, 1.0);
        ct = cc->EvalMult(ct, 0.5);

        return ct;
    }

    // === Full repack pipeline (no decryption) ===
    // TLWE lvl1 → CKKS, following OpenFHE's EvalFHEWtoCKKS approach:
    //
    // Key insight: prescale = 1/(q * K) maps the real-valued inner product
    // (b - a·s) into the range [-1, 1], because |(b - a·s)| ≤ q/2 and K ≈ q/2.
    // Then Chebyshev approximation of sin(2πx) on [-1, 1] extracts the
    // fractional part (message), followed by double-angle iterations.
    //
    // Parameters:
    //   K = 128: scaling constant, K > max(|noise|/q) ensures noise is absorbed
    //   BT_ITER = 3: number of double-angle iterations
    //   Chebyshev degree: determined by coefficients (119 for g_coefficientsFHEW128_8)
    inline Ciphertext<DCRTPoly> LWEsToOpenFHE(
        const CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        std::vector<TLWELvl1> &lwes,
        const RepackKey &rk,
        uint32_t scale_bits = 29,
        double K = 128.0)
    {
        size_t num_lwes = lwes.size();
        size_t tfhe_n = Lvl1::n;

        // prescale = 1/(q * K) where q = 2^32 (LWE modulus for Lvl1)
        // This maps (-a·s + b) which is O(2^31) into O(1/K) ∈ [-1, 1]
        double q_lwe = std::pow(2.0, 32);
        double prescale = 1.0 / (q_lwe * K);

        // Step 1: Build matrix A and vector b, scaled by prescale
        std::vector<std::vector<double>> A(num_lwes);
        std::vector<double> bvec(num_lwes);

        for (size_t i = 0; i < num_lwes; i++) {
            A[i].resize(tfhe_n);
            for (size_t j = 0; j < tfhe_n; j++) {
                // a values (negate for decryption formula: phase = b - a·s)
                A[i][j] = static_cast<double>(
                    static_cast<int32_t>(lwes[i][j])) * prescale;
            }
            bvec[i] = static_cast<double>(
                static_cast<int32_t>(lwes[i][tfhe_n])) * prescale;
        }

        // Step 2: Linear transform → Enc(A·s * prescale)
        // Note: A contains the original 'a' values (not negated).
        // The result is Enc(sum(a_j * s_j) * prescale).
        // We want Enc((b - a·s) * prescale), so in step 3 we compute b - A·s.
        auto AdotS = LinearTransformBSGS(cc, A, rk, tfhe_n);

        // Step 3: Compute b - A·s → Enc(phase * prescale)
        auto pt_b = cc->MakeCKKSPackedPlaintext(bvec);
        auto BminusAdotS = cc->EvalAdd(cc->EvalNegate(AdotS), pt_b);

        // Step 4: HomMod — Chebyshev approximation of sin(2πx)/(2π) on [-1, 1]
        // The input is phase/(q*K) which is in [-1/(2K), 1/(2K)] ⊂ [-1, 1]
        // plus integer multiples of 1/K (from the wrap-around).
        // sin(2πx) extracts the fractional part.
        //
        // Using OpenFHE's EvalChebyshevFunction to approximate
        // f(x) = (2π)^(-1/8) * cos(2π/8 * (x - 0.25))
        // which is equivalent to the double-angle formulation with BT_ITER=3.
        auto sinFunc = [](double x) -> double {
            // This is the base function before double-angle iterations.
            // With BT_ITER=3, the effective function applied is:
            // f(x) = sin(2πx)/(2π) after 3 double-angle steps.
            // The base Chebyshev approximates:
            //   (2π)^(-1/2^3) * cos(2π/2^3 * (x - 0.25))
            return std::pow(2.0 * M_PI, -1.0/8.0) *
                   std::cos(2.0 * M_PI / 8.0 * (x - 0.25));
        };

        uint32_t chebyDegree = 119;
        double a_cheby = -1.0;
        double b_cheby = 1.0;
        auto ct = cc->EvalChebyshevFunction(sinFunc, BminusAdotS, a_cheby, b_cheby, chebyDegree);

        // Step 5: Double-angle iterations (BT_ITER = 3)
        // Each iteration: y → 2y² - scalar
        const int32_t BT_ITER = 3;
        for (int32_t j = 1; j <= BT_ITER; j++) {
            ct = cc->EvalMult(ct, ct);     // y²
            ct = cc->EvalAdd(ct, ct);       // 2y²
            double scalar = 1.0 / std::pow(2.0 * M_PI, std::pow(2.0, j - BT_ITER));
            ct = cc->EvalSub(ct, scalar);   // 2y² - scalar
        }

        // Step 6: Post-scale
        // After HomMod, the output ≈ sin(2πx)/(2π) ≈ x for small x,
        // where x = m*Δ/(q*K). So output ≈ m*Δ/(q*K).
        // To recover m, multiply by q*K/Δ.
        double delta = std::pow(2.0, scale_bits);
        double postScale = q_lwe * K / delta;
        ct = cc->EvalMult(ct, postScale);

        return ct;
    }

    // === Simulated repack (for validation) ===
    inline Ciphertext<DCRTPoly> SimulatedRepack(
        const CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        std::vector<TLWELvl1> &lwes,
        const TFHESecretKey &sk,
        uint32_t scale_bits)
    {
        size_t n = lwes.size();
        std::vector<double> values(n);
        for (size_t i = 0; i < n; i++) {
            values[i] = static_cast<double>(
                TFHEpp::tlweSymInt32Decrypt<Lvl1>(
                    lwes[i], std::pow(2., scale_bits),
                    sk.key.get<Lvl1>()));
        }
        auto pt = cc->MakeCKKSPackedPlaintext(values);
        return cc->Encrypt(keys.publicKey, pt);
    }

    // === Utility: RotateAndSum ===
    inline Ciphertext<DCRTPoly> EvalSumSlots(
        const CryptoContext<DCRTPoly> &cc,
        Ciphertext<DCRTPoly> ct, size_t nSlots)
    {
        for (size_t step = 1; step < nSlots; step <<= 1) {
            auto rotated = cc->EvalRotate(ct, (int32_t)step);
            ct = cc->EvalAdd(ct, rotated);
        }
        return ct;
    }

} // namespace HEDB
