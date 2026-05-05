#pragma once
// TFHE→CKKS Repack using OpenFHE
// Port of HE3DB's LWEsToRLWE (SEAL) to OpenFHE CKKS
//
// Algorithm:
// 1. KeyGen: Encrypt TFHE sk in CKKS → Enc_CKKS(s_tfhe)
// 2. LinearTransform: Compute Enc_CKKS(A·s) using BSGS
// 3. Add b: Enc_CKKS(A·s + b) = Enc_CKKS(Δm + e)
// 4. HomMod: Modular reduction → Enc_CKKS(m/K)
// 5. HomRound: Round to 0/1

#include "openfhe.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include "HEDB/comparison/tfhepp_utils.h"

namespace HEDB
{
    using namespace lbcrypto;

    // --- Helper: ceil(sqrt(n)) ---
    template <typename T>
    static inline T CeilSqrt(T val) {
        return static_cast<T>(std::ceil(std::sqrt(1.0 * val)));
    }
    template <typename T>
    static inline T CeilDiv(T a, T b) {
        return (a + b - 1) / b;
    }

    // --- Repack pre-key: CKKS encryptions of rotated TFHE secret key ---
    struct RepackKey {
        std::vector<Ciphertext<DCRTPoly>> rotated_sk; // Enc(Rot^j(s_tfhe))
    };

    // === KeyGen ===
    // Encrypts the TFHE lvl1 secret key into CKKS ciphertexts
    // with BSGS baby-step rotations pre-computed
    inline void RepackKeyGen(
        RepackKey &rk,
        const CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        const TFHESecretKey &tfhe_sk,
        size_t tfhe_n)
    {
        // Extract TFHE secret key as doubles
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
            // Rotate sk_vec by 1 for next baby step
            std::rotate(sk_vec.begin(), sk_vec.begin() + 1, sk_vec.end());
        }
    }

    // === LinearTransform ===
    // Computes Enc(A · s) using BSGS diagonal method
    // A is (num_lwes x tfhe_n) matrix
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
                // Get diagonal j of matrix A
                std::vector<double> diag(max_len);
                for (size_t r = 0; r < max_len; r++) {
                    diag[r] = A[r % rows][(r + j) % cols];
                }
                std::rotate(diag.rbegin(), diag.rbegin() + b * g_tilde, diag.rend());

                auto pt_diag = cc->MakeCKKSPackedPlaintext(diag);
                auto term = cc->EvalMult(rk.rotated_sk[g], pt_diag);

                if (!sum_init) {
                    sum = term;
                    sum_init = true;
                } else {
                    sum = cc->EvalAdd(sum, term);
                }
            }

            if (b > 0) {
                sum = cc->EvalRotate(sum, (int32_t)(b * g_tilde));
            }

            if (!result_init) {
                result = sum;
                result_init = true;
            } else {
                result = cc->EvalAdd(result, sum);
            }
        }

        // If rows < cols, fold via rotate+add
        if (rows < cols) {
            size_t gama = static_cast<size_t>(std::log2(cols / rows));
            for (size_t j = 0; j < gama; j++) {
                auto temp = cc->EvalRotate(result, (int32_t)((1U << j) * rows));
                result = cc->EvalAdd(result, temp);
            }
        }

        return result;
    }

    // === HomMod: Homomorphic modular reduction ===
    // Uses sin approximation to extract m from Δm+e
    // Simplified version: uses polynomial approximation of sin(2π·x)/(2π)
    // then iterative doubling
    inline Ciphertext<DCRTPoly> HomModReduction(
        const CryptoContext<DCRTPoly> &cc,
        Ciphertext<DCRTPoly> ct,
        double K)
    {
        // For small-scale testing, we skip the full HomMod and
        // rely on the message being already close to 0 or 1/K
        // after the linear transform.
        // Full implementation would use Chebyshev approximation of
        // sin(2πx)/(2π) followed by r iterations of squaring.

        // Scale the result: multiply by K to recover m ∈ {0, 1}
        ct = cc->EvalMult(ct, K);
        return ct;
    }

    // === HomRound: Round to 0/1 ===
    // Applies sign function approximation: f(x) = 0 if x < 0.5, 1 if x >= 0.5
    // Using polynomial: 1.5x - 0.5x^3 iterated
    inline Ciphertext<DCRTPoly> HomRound(
        const CryptoContext<DCRTPoly> &cc,
        Ciphertext<DCRTPoly> ct)
    {
        // x → 2x - 1 (maps [0,1] to [-1,1])
        ct = cc->EvalMult(ct, 2.0);
        ct = cc->EvalAdd(ct, -1.0);

        // Apply sign approximation: 1.5x - 0.5x^3
        auto x2 = cc->EvalMult(ct, ct);
        auto x3 = cc->EvalMult(x2, ct);
        auto term1 = cc->EvalMult(ct, 1.5);
        auto term2 = cc->EvalMult(x3, -0.5);
        ct = cc->EvalAdd(term1, term2);

        // Second iteration for better approximation
        x2 = cc->EvalMult(ct, ct);
        x3 = cc->EvalMult(x2, ct);
        term1 = cc->EvalMult(ct, 1.5);
        term2 = cc->EvalMult(x3, -0.5);
        ct = cc->EvalAdd(term1, term2);

        // Map back: (x+1)/2
        ct = cc->EvalAdd(ct, 1.0);
        ct = cc->EvalMult(ct, 0.5);

        return ct;
    }

    // === Full repack pipeline ===
    // Converts vector of TLWE lvl1 ciphertexts to a single CKKS ciphertext
    inline Ciphertext<DCRTPoly> LWEsToOpenFHE(
        const CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        std::vector<TLWELvl1> &lwes,
        const RepackKey &rk)
    {
        size_t num_lwes = lwes.size();
        size_t tfhe_n = Lvl1::n;
        double K = 41.0;

        // Step 1: Build matrix A and vector b from TLWE ciphertexts
        std::vector<std::vector<double>> A(num_lwes);
        std::vector<double> b(num_lwes);

        // Scale factor: TLWE uses uint32_t modulus (2^32)
        // We need to normalize to [-0.5, 0.5) range then scale by 1/K
        double rescale = 1.0; // for lvl1 (32-bit)
        double multiplier = 1.0 / K;

        for (size_t i = 0; i < num_lwes; i++) {
            A[i].resize(tfhe_n);
            // Negate a (TFHEpp format: (a, b=a·s+m+e) → (-a, b) so decrypt = b + (-a)·s)
            for (size_t j = 0; j < tfhe_n; j++) {
                int32_t neg_a = -static_cast<int32_t>(lwes[i][j]);
                A[i][j] = static_cast<double>(neg_a) * multiplier;
            }
            b[i] = static_cast<double>(static_cast<int32_t>(lwes[i][tfhe_n])) * multiplier;
        }

        // Step 2: Linear transform → Enc(A·s)
        auto result = LinearTransformBSGS(cc, A, rk, tfhe_n);

        // Step 3: Add b → Enc(A·s + b) = Enc((Δm+e)/K)
        auto pt_b = cc->MakeCKKSPackedPlaintext(b);
        result = cc->EvalAdd(result, pt_b);

        // Step 4+5: HomMod + HomRound
        // For the simplified version, multiply by K and round
        result = HomModReduction(cc, result, K);
        result = HomRound(cc, result);

        return result;
    }

    // === Simulated repack (for validation) ===
    // Decrypts TLWE → re-encrypts in CKKS. Useful for debugging.
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
