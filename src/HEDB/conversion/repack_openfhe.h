#pragma once
// TFHE→CKKS Repack using OpenFHE's native scheme switching
// Converts TFHEpp TLWE lvl1 ciphertexts to CKKS ciphertext.
//
// Approach: Use OpenFHE's EvalFHEWtoCKKS which handles the full pipeline
// (BSGS linear transform + Chebyshev HomMod + double-angle).
// The key step is converting TFHEpp format to OpenFHE format.
//
// For experimental evaluation, a SimulatedRepack (decrypt + re-encrypt)
// is provided as a fallback when the full pipeline is not needed.

#include "openfhe.h"
#include "binfhecontext.h"
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

    // === Convert TFHEpp secret key to OpenFHE LWE private key ===
    inline LWEPrivateKey TFHEppToOpenFHEKey(
        const TFHESecretKey &tfhe_sk,
        size_t n,
        NativeInteger modulus)
    {
        NativeVector s_vec(n, modulus);
        for (size_t i = 0; i < n; i++) {
            auto val = tfhe_sk.key.get<Lvl1>()[i];
            if (val > 1) {
                // -1 in TFHEpp binary → modulus - 1 in OpenFHE
                s_vec[i] = modulus - 1;
            } else {
                s_vec[i] = NativeInteger(val);
            }
        }
        return std::make_shared<LWEPrivateKeyImpl>(std::move(s_vec));
    }

    // === Convert TFHEpp TLWE ciphertexts to OpenFHE LWE ciphertexts ===
    // Performs modulus switch from q=2^32 to q_target
    inline std::vector<LWECiphertext> TFHEppToOpenFHECiphertexts(
        const std::vector<TLWELvl1> &lwes,
        size_t n,
        uint64_t q_target)
    {
        double q_ratio = static_cast<double>(q_target) / std::pow(2.0, 32);
        std::vector<LWECiphertext> result(lwes.size());

        for (size_t i = 0; i < lwes.size(); i++) {
            NativeVector a_vec(n, NativeInteger(q_target));
            for (size_t j = 0; j < n; j++) {
                // Modulus switch: round(a * q_target / q_orig)
                double a_d = static_cast<double>(lwes[i][j]) * q_ratio;
                int64_t a_round = static_cast<int64_t>(std::round(a_d));
                // Map to [0, q_target)
                a_round = ((a_round % static_cast<int64_t>(q_target)) +
                           static_cast<int64_t>(q_target)) % static_cast<int64_t>(q_target);
                a_vec[j] = NativeInteger(static_cast<uint64_t>(a_round));
            }
            // Same for b
            double b_d = static_cast<double>(lwes[i][n]) * q_ratio;
            int64_t b_round = static_cast<int64_t>(std::round(b_d));
            b_round = ((b_round % static_cast<int64_t>(q_target)) +
                       static_cast<int64_t>(q_target)) % static_cast<int64_t>(q_target);
            NativeInteger b_native(static_cast<uint64_t>(b_round));

            result[i] = std::make_shared<LWECiphertextImpl>(std::move(a_vec), b_native);
        }

        return result;
    }

    // === Full repack pipeline using OpenFHE scheme switching ===
    // Uses OpenFHE's EvalFHEWtoCKKS which handles:
    // 1. BSGS linear transform (homomorphic partial decryption)
    // 2. Chebyshev sin approximation (HomMod)
    // 3. Double-angle iterations
    // 4. Post-scaling
    //
    // Returns CKKS ciphertext with message values in the slots.
    inline Ciphertext<DCRTPoly> LWEsToOpenFHE(
        CryptoContext<DCRTPoly> &cc,
        const KeyPair<DCRTPoly> &keys,
        std::vector<TLWELvl1> &lwes,
        const TFHESecretKey &tfhe_sk,
        uint32_t numSlots,
        uint32_t scale_bits = 29,
        uint32_t logQ_LWE = 28)
    {
        size_t n = Lvl1::n;

        // Step 1: Create OpenFHE FHEW context with matching parameters
        auto ccLWE = std::make_shared<BinFHEContext>();
        ccLWE->BinFHEContext::GenerateBinFHEContext(
            TOY, false, logQ_LWE, 0, GINX, false);

        // Step 2: Convert TFHEpp secret key to OpenFHE format
        uint64_t q_target = 1ULL << logQ_LWE;
        auto lwesk = TFHEppToOpenFHEKey(tfhe_sk, n, NativeInteger(q_target));

        // Step 3: Setup scheme switching
        cc->EvalFHEWtoCKKSSetup(ccLWE, numSlots, logQ_LWE);
        cc->SetBinCCForSchemeSwitch(ccLWE);
        cc->EvalFHEWtoCKKSKeyGen(keys, lwesk);

        // Step 4: Convert TFHEpp ciphertexts to OpenFHE format
        auto lwes_openfhe = TFHEppToOpenFHECiphertexts(lwes, n, q_target);

        // Step 5: Perform the scheme switching
        // With default p=4, m=1 maps to ~0.707 = sin(π/4) = 1/sqrt(2)
        // This is because TFHEpp encodes at Δ=2^29 = q/8, while p=4 expects q/4.
        // Apply correction factor sqrt(2) to map to {0, 1}.
        auto result = cc->EvalFHEWtoCKKS(lwes_openfhe, numSlots, numSlots);
        result = cc->EvalMult(result, std::sqrt(2.0));

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
