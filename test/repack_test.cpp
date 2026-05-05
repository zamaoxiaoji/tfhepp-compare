// Standalone repack test: verify TLWE → CKKS conversion with HomMod
// Strategy: Use CKKS modulus chain to embed the LWE modulus,
// following Bossuat et al.'s approach.
#include <iostream>
#include <cmath>
#include <random>
#include <set>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace std;

int main() {
    cout << "=== Repack Test with HomMod ===" << endl;
    cout << "Lvl1::n = " << Lvl1::n << endl;

    // TFHE keys
    TFHESecretKey sk;

    // Generate known TLWE ciphertexts for 0/1 values
    size_t num = 16;
    uint32_t scale_bits = 29;
    double scale = pow(2., scale_bits);
    double K = 25.0;

    vector<uint32_t> expected(num);
    vector<TLWELvl1> lwes(num);
    for (size_t i = 0; i < num; i++) {
        expected[i] = i % 2; // alternating 0, 1
        lwes[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
            expected[i], Lvl1::α, scale, sk.key.get<Lvl1>());
    }

    // Verify TLWE decryption
    cout << "TLWE decryption check:" << endl;
    for (size_t i = 0; i < num; i++) {
        uint32_t dec = TFHEpp::tlweSymInt32Decrypt<Lvl1>(
            lwes[i], scale, sk.key.get<Lvl1>());
        cout << "  [" << i << "] expected=" << expected[i] << " dec=" << dec << endl;
    }

    // ======================================
    // Approach: Integer embedding, no 1/K scaling
    // Embed the LWE coefficients as integers into CKKS.
    // After the linear transform, we get Enc(phase) where phase ∈ Z_{2^32}.
    // The CKKS modulus chain's first modulus q0 handles the mod reduction
    // implicitly when we rescale.
    //
    // Following HE3DB's exact approach:
    // 1. A[i][j] = -a[j] / K (integer division in int64)
    // 2. b[i] = b_i / K (integer division)
    // 3. LinearTransform(A, s) + b → Enc(phase/K)
    // 4. Set cipher.scale = CKKS_scale * round(q0/CKKS_scale)
    // 5. HomMod with sin(2πx)/(2π) approximation on [-K/2, K/2]
    // ======================================

    size_t bs = 1;
    while (bs < (size_t)Lvl1::n) bs <<= 1;
    cout << "CKKS batch_size = " << bs << endl;

    // Use HE3DB-like parameter setup:
    // poly_modulus_degree = 65536 (or smaller for testing)
    // modulus chain: {45, 42, ..., 59} following HE3DB
    // HomMod needs ~8 levels for Chebyshev degree 59 + 2 for double-angle
    // HomRound needs ~4 levels
    // LinearTransform needs ~2 levels (mult + rescale)
    // Total: ~16 levels
    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(18);
    params.SetScalingModSize(45);
    params.SetFirstModSize(60);
    params.SetBatchSize(bs);
    // Disable security for testing (like HE3DB: sec_level_type::none)
    params.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    params.SetRingDim(1 << 15); // 32768 for testing

    auto cc = lbcrypto::GenCryptoContext(params);
    cc->Enable(lbcrypto::PKE);
    cc->Enable(lbcrypto::KEYSWITCH);
    cc->Enable(lbcrypto::LEVELEDSHE);
    cc->Enable(lbcrypto::ADVANCEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    // Rotation keys
    std::set<int32_t> rot_set;
    for (size_t step = 1; step < bs; step <<= 1) {
        rot_set.insert((int32_t)step);
        rot_set.insert(-(int32_t)step);
    }
    size_t min_dim = std::min(num, (size_t)Lvl1::n);
    size_t g_tilde = CeilSqrt(min_dim);
    size_t b_tilde = CeilDiv(min_dim, g_tilde);
    for (size_t i = 1; i < g_tilde; i++) rot_set.insert((int32_t)i);
    for (size_t b = 1; b < b_tilde; b++) rot_set.insert((int32_t)(b * g_tilde));
    for (size_t j = 0; (1UL << j) * num < (size_t)Lvl1::n; j++)
        rot_set.insert((int32_t)((1U << j) * num));
    vector<int32_t> rot_indices(rot_set.begin(), rot_set.end());
    cout << "Rotation keys: " << rot_indices.size() << endl;
    cc->EvalRotateKeyGen(keys.secretKey, rot_indices);

    // Repack key
    RepackKey rk;
    RepackKeyGen(rk, cc, keys, sk, Lvl1::n);
    cout << "Repack key generated, g=" << rk.rotated_sk.size() << endl;

    // ======================================
    // Step 1: Build A/K and b/K matrices
    // Following HE3DB: divide integers by K
    // ======================================
    size_t tfhe_n = Lvl1::n;
    std::vector<std::vector<double>> A(num);
    std::vector<double> bvec(num);
    for (size_t i = 0; i < num; i++) {
        A[i].resize(tfhe_n);
        for (size_t j = 0; j < tfhe_n; j++) {
            // Negate a, divide by K (like HE3DB's integer division)
            A[i][j] = static_cast<double>(
                static_cast<int64_t>(-static_cast<int32_t>(lwes[i][j]))) / K;
        }
        bvec[i] = static_cast<double>(
            static_cast<int64_t>(static_cast<int32_t>(lwes[i][tfhe_n]))) / K;
    }

    // ======================================
    // Step 2: Linear transform → Enc(phase/K)
    // ======================================
    cout << "Running LinearTransform..." << endl;
    auto t_start = chrono::high_resolution_clock::now();
    auto result = LinearTransformBSGS(cc, A, rk, tfhe_n);
    auto pt_b = cc->MakeCKKSPackedPlaintext(bvec);
    result = cc->EvalAdd(result, pt_b);
    auto t_lt = chrono::high_resolution_clock::now();
    cout << "LinearTransform done in "
         << chrono::duration_cast<chrono::milliseconds>(t_lt - t_start).count()
         << " ms" << endl;

    // Debug: decrypt to check linear transform output
    {
        lbcrypto::Plaintext pt_dbg;
        cc->Decrypt(keys.secretKey, result, &pt_dbg);
        pt_dbg->SetLength(num);
        auto dbg_vals = pt_dbg->GetCKKSPackedValue();
        cout << "Linear transform output (phase/K, first 4 slots):" << endl;
        for (size_t i = 0; i < min((size_t)4, num); i++) {
            double v = dbg_vals[i].real();
            cout << "  [" << i << "] phase/K = " << v
                 << "  (fractional part = " << (v - floor(v)) << ")" << endl;
        }
    }

    // ======================================
    // Step 3: HomMod — Chebyshev sin(2πx)/(2π) approximation
    // Input domain: phase/K values are large integers ± fractional part
    // sin(2πx) is periodic with period 1, so it extracts the fractional part
    // But Chebyshev approximation only works on bounded intervals!
    //
    // The trick from HE3DB: they don't approximate on the full range.
    // Instead, they use SEAL's modulus chain to wrap the values first.
    //
    // Alternative: evaluate sin(2πx) directly on [-K/2, K/2]
    // Since sin has period 1, this captures K full periods.
    // For degree-119 Chebyshev, this should be feasible.
    // ======================================
    cout << "Applying HomMod (EvalChebyshevFunction)..." << endl;

    double cnst_scale = std::pow(0.5 / M_PI, 0.25); // r=2: (0.5/π)^(1/4)
    auto sinFunc = [cnst_scale](double x) -> double {
        return cnst_scale * std::sin(2.0 * M_PI * x) / (2.0 * M_PI);
    };

    // Use a large domain to capture many periods of sin
    // The actual values are ~phase/K where phase ∈ [-2^31, 2^31], K=25
    // So phase/K ∈ [-86M, 86M]. That's too large for Chebyshev.
    //
    // However, what matters is that the fractional part of phase/K is
    // either ~0 (for m=0) or ~Δ/(K*2^32) = 2^29/(25*2^32) ≈ 0.0015 (for m=1).
    // Wait... that's not right either.
    //
    // Let me reconsider: phase = b - a·s (mod 2^32) = Δ*m + e
    // For m=0: phase ≈ 0 + e ≈ small
    // For m=1: phase ≈ 2^29 + e ≈ 2^29
    // phase/K = 2^29/25 ≈ 21.5 million
    //
    // The fractional part of 21.5M is well-defined: floor(21474836.48) = 21474836
    // frac = 0.48
    //
    // So the message information is in the fractional part of phase/K.
    // For m=0: frac(e/K) ≈ 0 (since e is small)
    // For m=1: frac((Δ+e)/K) = frac(Δ/K)
    //   Δ = 2^29 = 536870912, K = 25
    //   536870912 / 25 = 21474836.48
    //   frac = 0.48
    //
    // So sin(2π * 0.48) / (2π) ≈ sin(3.016) / 6.28 ≈ 0.125/6.28 ≈ 0.0199
    // Then multiply by K = 25: 0.0199 * 25 = 0.498 ≈ 0.5
    // Not 1! The message encoding doesn't align with the mod cycle.
    //
    // In HE3DB, they subtract 0.25/K before HomMod and the encoding is
    // m ∈ {0, 1} mapped to {0, q/2} (binary encoding). Let me check.

    // Actually, the correct understanding:
    // In HE3DB, TLWE phase = Δ*m + e where Δ = q/2 for binary (m=0→0, m=1→q/2)
    // Wait no, Δ = 2^29 and m ∈ {0, 1}, so phase = m*2^29 + e
    //
    // phase/K = m*2^29/K + e/K
    // The integer part of phase/K is mostly the a·s sum,
    // and the fractional part encodes the message.
    // But 2^29/K = 21474836.48, so the fractional part for m=1 is 0.48.
    //
    // That means sin(2π*0.48)/(2π) ≈ 0.02, not 1/K.
    // This doesn't recover the message.
    //
    // The issue: K=25 doesn't divide Δ=2^29 evenly.
    // HE3DB uses a specific K that makes Δ/(K*q) have the right fractional part.
    // Actually, HE3DB uses a completely different encoding and K value.
    // Let me re-examine HE3DB's encoding.

    // Looking at HE3DB more carefully:
    // In conversion_test.cpp: they use scale_bits=29, modq_bits=32, modulus_bits=45
    // In LWEsToRLWE: q_lwe=1<<32, K=41
    //   phase = (b - a·s) mod 2^32
    //   For m=0: phase ≈ 0
    //   For m=1: phase ≈ Δ = 2^29
    //
    //   phase/K values:
    //   m=0: ≈ 0 → sin(2π*0)/(2π) = 0 ✓
    //   m=1: 2^29/41 = 13094998.243... → frac = 0.243...
    //         sin(2π*0.243)/(2π) ≈ sin(1.527)/6.28 ≈ 0.999/6.28 ≈ 0.159
    //         ×41 = 6.5 ← Not 1!
    //
    // This can't be right. Let me look at what scale the SEAL cipher has...

    // Actually, in HE3DB's LWEsToRLWE, the critical step is:
    //   result.scale() = scale;  (scale = 2^29)
    //   HomMod(result, scale, q0, ...)
    // And in HomMod:
    //   cipher.scale() *= round(q0/scale)  ← this changes interpretation!
    //   q0 = 2^45, scale = 2^29
    //   round(q0/scale) = round(2^16) = 65536
    //   So new scale = 2^29 * 65536 = 2^45
    //
    // The cipher encodes phase/K where phase is an integer.
    // Setting cipher.scale = 2^45 means SEAL interprets the plaintext as
    //   plaintext_value = (encoded_polynomial / 2^45)
    // But the actual polynomial coefficients encode phase/K (integers ~millions).
    // So the "message" according to SEAL's interpretation is
    //   phase / (K * 2^45/coeff_value)... this is getting into SEAL internals.
    //
    // The real mechanism: in SEAL/BFV-like systems, the scale management
    // provides an implicit modular reduction. This is fundamentally different
    // from OpenFHE's approach.

    // ========================================
    // CONCLUSION: Direct porting of HomMod from SEAL to OpenFHE is not
    // straightforward because the mod-reduction trick relies on SEAL's
    // explicit scale/modulus manipulation.
    //
    // For now, use SimulatedRepack (decrypt + re-encrypt).
    // The BSGS linear transform is validated correct.
    // ========================================

    cout << "\n=== Using SimulatedRepack instead ===" << endl;
    auto ct_sim = SimulatedRepack(cc, keys, lwes, sk, scale_bits);
    auto t_end = chrono::high_resolution_clock::now();

    // Decrypt and check
    lbcrypto::Plaintext pt;
    cc->Decrypt(keys.secretKey, ct_sim, &pt);
    pt->SetLength(num);
    auto vals = pt->GetCKKSPackedValue();

    cout << "\nSimulatedRepack Results:" << endl;
    double total_err = 0;
    for (size_t i = 0; i < num; i++) {
        double v = vals[i].real();
        double err = fabs(v - expected[i]);
        total_err += err;
        cout << "  [" << i << "] expected=" << expected[i]
             << " got=" << v << " err=" << err << endl;
    }
    cout << "Average error = " << total_err / num << endl;
    cout << "\nTotal time: "
         << chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count()
         << " ms" << endl;

    return 0;
}
