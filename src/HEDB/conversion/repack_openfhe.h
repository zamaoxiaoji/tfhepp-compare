#pragma once
// TFHE→CKKS Repack using OpenFHE
// Converts TLWE lvl1 ciphertexts to CKKS ciphertext without decryption.
//
// Algorithm (following HE3DB / Bossuat et al.):
// 1. KeyGen: Encrypt TFHE sk in CKKS → Enc_CKKS(s_tfhe)
// 2. LinearTransform: BSGS diagonal method → Enc_CKKS(-a·s/K)
// 3. Add b/K: Enc_CKKS(phase/K) where phase = (b - a·s) mod q
// 4. HomMod: Chebyshev approximation of sin(2πx)/(2π) → extract mod-q fractional part
// 5. Multiply by K → Enc_CKKS(m)
// 6. HomRound: polynomial sign approximation to snap to {0, 1}

#include "openfhe.h"
#include "math/chebyshev.h"
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

    // === HomMod ===
    // Performs approximate modular reduction in CKKS domain.
    // Input: Enc(phase/K) where phase = (b - a·s) is the real-valued inner product
    //   and K is a scaling constant (typically 25 for Lvl1).
    // Output: Enc(m/K) where m is the message extracted via mod reduction.
    //
    // Uses Chebyshev approximation of f(x) = sin(2πx) / (2π·K)
    // followed by r iterations of double-angle: y → 2y² - (1/(2π))^(2^i)
    //
    // Following Algorithm 7 in Bossuat et al. (2021).
    inline Ciphertext<DCRTPoly> HomMod(
        const CryptoContext<DCRTPoly> &cc,
        Ciphertext<DCRTPoly> ct,
        double K,
        uint32_t polyDegree = 59,
        uint32_t r = 2)
    {
        // The input is phase/K where phase ∈ [-q/2, q/2] and K = 25.
        // After HomMod we want to extract (phase mod 1) / K.
        // sin(2πx)/(2π) maps x ∈ R to its fractional part (near integers).
        //
        // Step 1: Scale the Chebyshev approximation coefficients
        // Following HE3DB: cnst_scale = (0.5/π)^(1/2^r)
        double cnst_scale = std::pow(0.5 / M_PI, 1.0 / (1 << r));

        // The function to approximate: cnst_scale * sin(2πx) / (2π)
        // Domain: [-1, 1] (after Chebyshev mapping from the input range)
        // The input range for phase/K:
        //   phase = sum(-a_i * s_i) + b ∈ [-q/2, q/2] where q = 2^32
        //   After /K: phase/K ∈ [-q/(2K), q/(2K)]
        //   But the actual values cluster around m/K + integer multiples.

        // For Lvl1 (32-bit), the range of phase/K is huge.
        // The sin function is periodic, so we evaluate on [-0.5, 0.5]
        // which captures one period. The offset -0.25/K shifts the encoding.

        // Step 1: Subtract 0.25/K (encoding offset for binary {0,1} → {0, Δ/2})
        ct = cc->EvalAdd(ct, -0.25 / K);

        // Step 2: Evaluate Chebyshev approximation of cnst_scale * sin(2πx)/(2π)
        // The input should already be mapped to the correct domain by the caller.
        // We evaluate on [-0.5, 0.5] which is one period of sin(2πx).
        auto sinFunc = [cnst_scale](double x) -> double {
            return cnst_scale * std::sin(2.0 * M_PI * x) / (2.0 * M_PI);
        };

        ct = cc->EvalChebyshevFunction(sinFunc, ct, -0.5, 0.5, polyDegree);

        // Step 3: Double-angle iterations: y → 2y² - θ
        // Each iteration squares the approximation, improving accuracy.
        double theta = std::pow(0.5 / M_PI, 1.0 / std::pow(2.0, r));
        for (uint32_t i = 0; i < r; i++) {
            theta *= theta;
            ct = cc->EvalMult(ct, ct);   // y²
            ct = cc->EvalAdd(ct, ct);     // 2y²
            ct = cc->EvalAdd(ct, -theta); // 2y² - θ
        }

        return ct;
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
    // TLWE lvl1 → CKKS, following HE3DB's approach:
    // 1. Scale by 1/K (K=25), compute A·s/K + b/K via BSGS linear transform
    // 2. HomMod to extract mod-q fractional part
    // 3. Multiply by K to recover message
    // 4. HomRound to snap to {0, 1}
    inline Ciphertext<DCRTPoly> LWEsToOpenFHE(
        const CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        std::vector<TLWELvl1> &lwes,
        const RepackKey &rk,
        uint32_t scale_bits = 29,
        double K = 25.0,
        uint32_t modPolyDegree = 59)
    {
        size_t num_lwes = lwes.size();
        size_t tfhe_n = Lvl1::n;

        // Step 1: Build matrix A/K and vector b/K
        // TLWE: (a, b) where phase = b - a·s ≡ Δ·m + e (mod q)
        // We compute: Enc_CKKS((-a·s + b) / K)
        std::vector<std::vector<double>> A(num_lwes);
        std::vector<double> bvec(num_lwes);

        for (size_t i = 0; i < num_lwes; i++) {
            A[i].resize(tfhe_n);
            for (size_t j = 0; j < tfhe_n; j++) {
                // Negate a, scale by 1/K (integer division like HE3DB)
                A[i][j] = static_cast<double>(
                    -static_cast<int32_t>(lwes[i][j])) / K;
            }
            bvec[i] = static_cast<double>(
                static_cast<int32_t>(lwes[i][tfhe_n])) / K;
        }

        // Step 2: Linear transform → Enc(-a·s / K)
        auto result = LinearTransformBSGS(cc, A, rk, tfhe_n);

        // Step 3: Add b/K → Enc(phase/K)
        auto pt_b = cc->MakeCKKSPackedPlaintext(bvec);
        result = cc->EvalAdd(result, pt_b);

        // Step 4: HomMod → Enc(m/K) (extract message via mod reduction)
        result = HomMod(cc, result, K, modPolyDegree);

        // Step 5: Multiply by K → Enc(m)
        result = cc->EvalMult(result, K);

        // Step 6: HomRound to snap to {0, 1}
        result = HomRound(cc, result);

        return result;
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
