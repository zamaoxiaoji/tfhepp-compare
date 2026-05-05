// TPC-H Query 6 - End-to-end encrypted query
// WHERE filtering in TFHE (ETHMSB) → repack to CKKS → aggregation
//
// SELECT SUM(l_extendedprice * l_discount) AS revenue
// FROM lineitem
// WHERE l_shipdate >= 20101 AND l_shipdate < 21231
//   AND l_discount >= 8 AND l_discount <= 10
//   AND l_quantity < 32

#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace std;

void query_evaluation(size_t rows)
{
    cout << "=== TPC-H Q6 | rows=" << rows << " ===" << endl;

    // =================== Key Generation ===================
    TFHESecretKey sk;
    TFHEEvalKey ek;
    cout << "Generating TFHE keys..." << endl;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    // =================== Data Generation ===================
    cout << "Generating data..." << endl;
    random_device rd;
    default_random_engine eng(rd());

    uint32_t quantity_bits = 6, discount_bits = 4, ship_bits = 16;
    uint32_t quantity_scale = numeric_limits<Lvl1::T>::digits - quantity_bits - 1;
    uint32_t discount_scale = numeric_limits<Lvl1::T>::digits - discount_bits - 1;
    uint32_t ship_scale = numeric_limits<Lvl2::T>::digits - ship_bits - 1;

    uniform_int_distribution<Lvl1::T> qty_dist(0, (1 << quantity_bits) - 1);
    uniform_int_distribution<Lvl1::T> disc_dist(0, (1 << discount_bits) - 1);
    uniform_int_distribution<uint64_t> price_dist(1, 10);

    vector<uint32_t> quantity_data(rows), discount_data(rows);
    vector<uint64_t> ship_data(rows);
    vector<double> price_data(rows), discount_double(rows);

    auto generate_date = [&]() -> uint64_t {
        uniform_int_distribution<int> y(1, 2), m(1, 12), d(1, 28);
        return d(eng) + 100*m(eng) + 10000*(2020+y(eng));
    };

    for (size_t i = 0; i < rows; i++) {
        quantity_data[i] = qty_dist(eng);
        discount_data[i] = disc_dist(eng);
        ship_data[i] = generate_date();
        price_data[i] = (double)price_dist(eng);
        discount_double[i] = (double)discount_data[i];
    }
    // Force a known-good record
    quantity_data[0] = 1; discount_data[0] = 9; ship_data[0] = 21215;

    // =================== TFHE Encryption ===================
    cout << "Encrypting with TFHE..." << endl;
    vector<TLWELvl1> qty_ct(rows), disc_ct(rows);
    vector<TLWELvl2> ship_ct(rows);
    for (size_t i = 0; i < rows; i++) {
        qty_ct[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity_data[i], Lvl1::α, pow(2., quantity_scale), sk.key.get<Lvl1>());
        disc_ct[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(discount_data[i], Lvl1::α, pow(2., discount_scale), sk.key.get<Lvl1>());
        ship_ct[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(ship_data[i], Lvl2::α, pow(2., ship_scale), sk.key.get<Lvl2>());
    }

    // Predicate thresholds
    Lvl2::T pred1 = 20101, pred2 = 21231;
    Lvl1::T pred3 = 8, pred4 = 10, pred5 = 32;
    auto ct_pred1 = TFHEpp::tlweSymInt32Encrypt<Lvl2>(pred1, Lvl2::α, pow(2., ship_scale), sk.key.get<Lvl2>());
    auto ct_pred2 = TFHEpp::tlweSymInt32Encrypt<Lvl2>(pred2, Lvl2::α, pow(2., ship_scale), sk.key.get<Lvl2>());
    auto ct_pred3 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(pred3, Lvl1::α, pow(2., discount_scale), sk.key.get<Lvl1>());
    auto ct_pred4 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(pred4, Lvl1::α, pow(2., discount_scale), sk.key.get<Lvl1>());
    auto ct_pred5 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(pred5, Lvl1::α, pow(2., quantity_scale), sk.key.get<Lvl1>());

    // =================== Predicate Evaluation (ETHMSB) ===================
    cout << "Evaluating predicates (ETHMSB)..." << endl;
    vector<TLWELvl1> pred_cres(rows);
    vector<uint32_t> pred_plain(rows);

    // Compute plaintext predicate results
    for (size_t i = 0; i < rows; i++) {
        uint32_t r1 = (ship_data[i] >= pred1) ? 1 : 0;
        uint32_t r2 = (ship_data[i] < pred2) ? 1 : 0;
        uint32_t r3 = (discount_data[i] >= pred3) ? 1 : 0;
        uint32_t r4 = (discount_data[i] <= pred4) ? 1 : 0;
        uint32_t r5 = (quantity_data[i] < pred5) ? 1 : 0;
        pred_plain[i] = r1 * r2 * r3 * r4 * r5;
    }

    auto t_filter_start = chrono::high_resolution_clock::now();

    vector<TLWELvl1> cres1(rows), cres2(rows), cres3(rows), cres4(rows), cres5(rows);
    for (size_t i = 0; i < rows; i++) {
        // ship >= pred1
        ethmsb_greater_than_equal<Lvl2>(ship_ct[i], ct_pred1, cres1[i], ship_bits, ek, LOGIC);
        // ship < pred2
        ethmsb_less_than<Lvl2>(ship_ct[i], ct_pred2, cres2[i], ship_bits, ek, LOGIC);
        HomAND(pred_cres[i], cres1[i], cres2[i], ek, LOGIC);

        // disc >= pred3
        ethmsb_greater_than_equal<Lvl1>(disc_ct[i], ct_pred3, cres3[i], discount_bits, ek, LOGIC);
        HomAND(pred_cres[i], pred_cres[i], cres3[i], ek, LOGIC);

        // disc <= pred4
        ethmsb_less_than_equal<Lvl1>(disc_ct[i], ct_pred4, cres4[i], discount_bits, ek, LOGIC);
        HomAND(pred_cres[i], pred_cres[i], cres4[i], ek, LOGIC);

        // qty < pred5
        ethmsb_less_than<Lvl1>(qty_ct[i], ct_pred5, cres5[i], quantity_bits, ek, LOGIC);
        HomAND(pred_cres[i], pred_cres[i], cres5[i], ek, ARITHMETIC);
    }

    auto t_filter_end = chrono::high_resolution_clock::now();
    double filter_ms = chrono::duration_cast<chrono::milliseconds>(t_filter_end - t_filter_start).count();

    // Rescale to 29-bit for repack compatibility
    uint32_t rlwe_scale_bits = 29;
    for (size_t i = 0; i < rows; i++) {
        TFHEpp::ari_rescale(pred_cres[i], pred_cres[i], rlwe_scale_bits, ek);
    }

    // Verify predicate results
    size_t pred_errors = 0;
    for (size_t i = 0; i < rows; i++) {
        uint32_t dec = TFHEpp::tlweSymInt32Decrypt<Lvl1>(
            pred_cres[i], pow(2., rlwe_scale_bits), sk.key.get<Lvl1>());
        if (dec != pred_plain[i]) pred_errors++;
    }
    cout << "  Filter time: " << filter_ms << " ms" << endl;
    cout << "  Predicate errors: " << pred_errors << "/" << rows << endl;

    // =================== CKKS Setup ===================
    cout << "Setting up OpenFHE CKKS..." << endl;
    CCParams<CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(17); // 13 repack + 1 sqrt2 + 1 mask*data + 2 spare
    params.SetScalingModSize(50);
    params.SetScalingTechnique(FIXEDAUTO);
    params.SetSecurityLevel(HEStd_128_classic);
    size_t bs = 1; while (bs < rows) bs <<= 1;
    params.SetBatchSize(bs);

    auto cc = GenCryptoContext(params);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(SCHEMESWITCH);

    auto ckks_keys = cc->KeyGen();
    cc->EvalMultKeyGen(ckks_keys.secretKey);

    // Rotation keys for RotateAndSum
    vector<int32_t> rot_indices;
    for (size_t step = 1; step < bs; step <<= 1)
        rot_indices.push_back((int32_t)step);
    cc->EvalRotateKeyGen(ckks_keys.secretKey, rot_indices);

    // =================== Repack: TFHE → CKKS (true homomorphic) ===================
    cout << "Repacking TFHE→CKKS (EvalFHEWtoCKKS)..." << endl;
    auto t_repack_start = chrono::high_resolution_clock::now();

    auto ct_mask = LWEsToOpenFHE(cc, ckks_keys, pred_cres, sk, rows, rlwe_scale_bits);

    auto t_repack_end = chrono::high_resolution_clock::now();
    double repack_ms = chrono::duration_cast<chrono::milliseconds>(t_repack_end - t_repack_start).count();
    cout << "  Repack time: " << repack_ms << " ms" << endl;

    // =================== CKKS Aggregation ===================
    cout << "CKKS aggregation..." << endl;

    // Compute price * discount and encrypt
    vector<double> price_discount(rows);
    for (size_t i = 0; i < rows; i++) {
        price_discount[i] = price_data[i] * discount_double[i];
    }
    auto pt_pd = cc->MakeCKKSPackedPlaintext(price_discount);
    auto ct_pd = cc->Encrypt(ckks_keys.publicKey, pt_pd);

    auto t_agg_start = chrono::high_resolution_clock::now();

    // mask × price_discount
    auto ct_filtered = cc->EvalMult(ct_mask, ct_pd);

    // RotateAndSum
    auto ct_sum = EvalSumSlots(cc, ct_filtered, rows);

    auto t_agg_end = chrono::high_resolution_clock::now();
    double agg_ms = chrono::duration_cast<chrono::milliseconds>(t_agg_end - t_agg_start).count();

    // =================== Decrypt and verify ===================
    Plaintext pt_result;
    cc->Decrypt(ckks_keys.secretKey, ct_sum, &pt_result);
    pt_result->SetLength(1);
    double encrypted_result = pt_result->GetCKKSPackedValue()[0].real();

    // Compute plaintext result
    double plain_result = 0;
    for (size_t i = 0; i < rows; i++) {
        plain_result += price_data[i] * discount_double[i] * pred_plain[i];
    }

    cout << "\n=== Results ===" << endl;
    cout << "  Plain result:     " << plain_result << endl;
    cout << "  Encrypted result: " << round(encrypted_result) << endl;
    cout << "  Error:            " << abs(encrypted_result - plain_result) << endl;
    cout << "\n=== Timing ===" << endl;
    cout << "  Filter (TFHE):    " << filter_ms / 1000 << " s" << endl;
    cout << "  Repack:           " << repack_ms / 1000 << " s" << endl;
    cout << "  Aggregation:      " << agg_ms / 1000 << " s" << endl;
    cout << "  Total:            " << (filter_ms + repack_ms + agg_ms) / 1000 << " s" << endl;
}

int main()
{
    query_evaluation(16);
    return 0;
}
