// Standalone repack test: verify TLWE → CKKS conversion
#include <iostream>
#include <cmath>
#include <random>
#include <set>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace std;

int main() {
    cout << "=== Repack Test ===" << endl;
    cout << "Lvl1::n = " << Lvl1::n << endl;

    // TFHE keys
    TFHESecretKey sk;

    // Generate known TLWE ciphertexts for 0/1 values
    size_t num = 16;
    uint32_t scale_bits = 29;
    double scale = pow(2., scale_bits);

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
    size_t bs = 1;
    while (bs < (size_t)Lvl1::n) bs <<= 1;
    cout << "CKKS batch_size = " << bs << endl;

    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(10);
    params.SetScalingModSize(50);
    params.SetBatchSize(bs);
    params.SetSecurityLevel(lbcrypto::HEStd_128_classic);

    auto cc = lbcrypto::GenCryptoContext(params);
    cc->Enable(lbcrypto::PKE);
    cc->Enable(lbcrypto::KEYSWITCH);
    cc->Enable(lbcrypto::LEVELEDSHE);

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

    // True repack
    cout << "Running true repack..." << endl;
    auto ct_result = LWEsToOpenFHE(cc, keys, lwes, rk, scale_bits);

    // Decrypt and check
    lbcrypto::Plaintext pt;
    cc->Decrypt(keys.secretKey, ct_result, &pt);
    pt->SetLength(num);
    auto vals = pt->GetCKKSPackedValue();

    // Also compute expected via plaintext inner product
    double delta_inv = 1.0 / pow(2.0, scale_bits);
    auto sk_raw = sk.key.get<Lvl1>();

    cout << "\nResults:" << endl;
    double total_err = 0;
    for (size_t i = 0; i < num; i++) {
        // Compute plaintext phase: (b - a·s) / Δ
        double phase = 0;
        for (size_t j = 0; j < (size_t)Lvl1::n; j++) {
            int32_t a_val = static_cast<int32_t>(lwes[i][j]);
            int32_t s_val = (sk_raw[j] > 1) ? -1 : sk_raw[j];
            phase -= a_val * s_val;
        }
        phase += static_cast<int32_t>(lwes[i][Lvl1::n]);
        double expected_ckks = phase * delta_inv;

        double v = vals[i].real();
        double err = fabs(v - expected[i]);
        total_err += err;
        cout << "  [" << i << "] expected=" << expected[i]
             << " phase/Δ=" << expected_ckks
             << " got=" << v << " err=" << err << endl;
    }
    cout << "Average error = " << total_err / num << endl;

    return 0;
}
