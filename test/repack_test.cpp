// Standalone repack test: verify TLWE → CKKS conversion with HomMod
// Following OpenFHE's EvalFHEWtoCKKS approach: prescale = 1/(q*K)
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
    double K = 128.0;

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

    // CKKS setup
    // Depth budget following OpenFHE's EvalFHEWtoCKKS (SwitchFHEWtoCKKS example):
    //   r=3 (BT_ITER=3) → Chebyshev max depth = 9, +1 for postscaling, +3 for BT_ITER
    //   Total: 3 + 9 + 1 = 13 for FIXEDAUTO
    //   Plus 2 for LinearTransform (ct×pt mult + folding rescale)
    //   Plus 1 for postscale mult
    //   Total: ~16-18
    size_t bs = 1;
    while (bs < (size_t)Lvl1::n) bs <<= 1;
    cout << "CKKS batch_size = " << bs << endl;

    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(20);
    params.SetScalingModSize(50);
    params.SetFirstModSize(60);
    params.SetBatchSize(bs);
    params.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    params.SetRingDim(1 << 15); // 32768 for testing

    auto cc = lbcrypto::GenCryptoContext(params);
    cc->Enable(lbcrypto::PKE);
    cc->Enable(lbcrypto::KEYSWITCH);
    cc->Enable(lbcrypto::LEVELEDSHE);
    cc->Enable(lbcrypto::ADVANCEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    cout << "Ring dimension: " << cc->GetRingDimension() << endl;

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

    // ============================
    // First test: plaintext verification of prescale approach
    // ============================
    {
        double q_lwe = pow(2.0, 32);
        double prescale = 1.0 / (q_lwe * K);
        cout << "\nprescale = 1/(2^32 * " << K << ") = " << prescale << endl;

        cout << "Plaintext phase/(q*K) values:" << endl;
        auto sk_raw = sk.key.get<Lvl1>();
        for (size_t i = 0; i < min((size_t)4, num); i++) {
            double phase_real = 0;
            for (size_t j = 0; j < (size_t)Lvl1::n; j++) {
                int32_t a_val = static_cast<int32_t>(lwes[i][j]);
                int32_t s_val = (sk_raw[j] > 1) ? -1 : sk_raw[j];
                phase_real += a_val * (double)s_val;  // a·s
            }
            double b_val = static_cast<int32_t>(lwes[i][Lvl1::n]);
            double phase = b_val - phase_real;  // b - a·s (real, not mod q)
            double phase_mod_q = fmod(phase, q_lwe);
            if (phase_mod_q > q_lwe / 2) phase_mod_q -= q_lwe;
            if (phase_mod_q < -q_lwe / 2) phase_mod_q += q_lwe;

            cout << "  [" << i << "] expected=" << expected[i]
                 << " phase_real=" << phase
                 << " phase_mod_q=" << phase_mod_q
                 << " phase/(q*K)=" << phase * prescale
                 << " mod_phase/(q*K)=" << phase_mod_q * prescale
                 << " sin(2pi*mod/qK)*2pi=" << sin(2*M_PI*phase_mod_q*prescale)*2*M_PI
                 << endl;
        }
    }

    // ============================
    // True repack with HomMod
    // ============================
    cout << "\nRunning true repack with HomMod..." << endl;
    auto t_start = chrono::high_resolution_clock::now();
    auto ct_result = LWEsToOpenFHE(cc, keys, lwes, rk, scale_bits, K);
    auto t_end = chrono::high_resolution_clock::now();
    double ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
    cout << "Repack time: " << ms << " ms" << endl;

    // Decrypt and check
    lbcrypto::Plaintext pt;
    cc->Decrypt(keys.secretKey, ct_result, &pt);
    pt->SetLength(num);
    auto vals = pt->GetCKKSPackedValue();

    cout << "\nResults:" << endl;
    double total_err = 0;
    for (size_t i = 0; i < num; i++) {
        double v = vals[i].real();
        double err = fabs(v - expected[i]);
        total_err += err;
        cout << "  [" << i << "] expected=" << expected[i]
             << " got=" << v << " err=" << err << endl;
    }
    cout << "Average error = " << total_err / num << endl;

    return 0;
}
