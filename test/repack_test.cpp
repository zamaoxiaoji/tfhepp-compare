// Standalone repack test: verify TLWE → CKKS conversion
// Uses OpenFHE's native EvalFHEWtoCKKS scheme switching
#include <iostream>
#include <cmath>
#include <set>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace std;

int main() {
    cout << "=== Repack Test via OpenFHE Scheme Switching ===" << endl;
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

    // Verify modulus switch preserves message
    cout << "\nModulus switch verification (q=2^32 → q'=2^28):" << endl;
    {
        uint64_t q_target = 1ULL << 28;
        auto lwes_ofhe = TFHEppToOpenFHECiphertexts(lwes, Lvl1::n, q_target);
        auto lwesk = TFHEppToOpenFHEKey(sk, Lvl1::n, NativeInteger(q_target));

        for (size_t i = 0; i < min((size_t)4, num); i++) {
            // Manual decryption: phase = b - a·s mod q'
            auto& a_vec = lwes_ofhe[i]->GetA();
            auto b_val = lwes_ofhe[i]->GetB();
            auto& s_vec = lwesk->GetElement();

            // Compute a·s mod q'
            NativeInteger as_sum(0);
            for (size_t j = 0; j < (size_t)Lvl1::n; j++) {
                as_sum = as_sum.ModAdd(
                    a_vec[j].ModMul(s_vec[j], NativeInteger(q_target)),
                    NativeInteger(q_target));
            }
            auto phase = b_val.ModSub(as_sum, NativeInteger(q_target));
            // Map phase to signed: if phase > q'/2, phase -= q'
            int64_t phase_signed = phase.ConvertToInt();
            if (phase_signed > (int64_t)(q_target / 2))
                phase_signed -= q_target;

            // With scale_bits=29, Δ=2^29. After modulus switch:
            // Δ' = round(2^29 * 2^28 / 2^32) = round(2^25) = 2^25
            double delta_new = std::pow(2.0, scale_bits) * q_target / std::pow(2.0, 32);
            int32_t msg = (int32_t)std::round(phase_signed / delta_new);

            cout << "  [" << i << "] expected=" << expected[i]
                 << " phase=" << phase_signed
                 << " Δ'=" << delta_new
                 << " decoded=" << msg << endl;
        }
    }

    // CKKS setup for scheme switching
    // Following OpenFHE SwitchFHEWtoCKKS example:
    // For r=3 in FHEWtoCKKS, Chebyshev max depth = 9, +1 for postscaling
    // Total: 3 + 9 + 1 = 13
    uint32_t multDepth = 3 + 9 + 1;
    uint32_t scaleModSize = 50;
    uint32_t ringDim = 8192;
    uint32_t logQ_LWE = 28;
    uint32_t slots = num;
    uint32_t batchSize = slots;

    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(multDepth);
    params.SetScalingModSize(scaleModSize);
    params.SetScalingTechnique(lbcrypto::FIXEDAUTO);
    params.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    params.SetRingDim(ringDim);
    params.SetBatchSize(batchSize);

    auto cc = lbcrypto::GenCryptoContext(params);
    cc->Enable(lbcrypto::PKE);
    cc->Enable(lbcrypto::KEYSWITCH);
    cc->Enable(lbcrypto::LEVELEDSHE);
    cc->Enable(lbcrypto::ADVANCEDSHE);
    cc->Enable(lbcrypto::SCHEMESWITCH);

    auto keys = cc->KeyGen();

    cout << "\nCKKS ring dimension: " << cc->GetRingDimension() << endl;
    cout << "Multiplicative depth: " << multDepth << endl;
    cout << "Slots: " << slots << endl;

    // True repack via OpenFHE scheme switching
    cout << "\nRunning true repack via EvalFHEWtoCKKS..." << endl;
    auto t_start = chrono::high_resolution_clock::now();
    auto ct_result = LWEsToOpenFHE(cc, keys, lwes, sk, slots, logQ_LWE);
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
    int correct = 0;
    for (size_t i = 0; i < num; i++) {
        double v = vals[i].real();
        double err = fabs(v - expected[i]);
        total_err += err;
        bool ok = (err < 0.3);
        if (ok) correct++;
        cout << "  [" << i << "] expected=" << expected[i]
             << " got=" << v << " err=" << err
             << (ok ? " ✓" : " ✗") << endl;
    }
    cout << "Average error = " << total_err / num << endl;
    cout << "Correct: " << correct << "/" << num << endl;


    return 0;
}
