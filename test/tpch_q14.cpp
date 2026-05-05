// TPC-H Query 14 - End-to-end encrypted query
// WHERE filtering in TFHE (ETHMSB) → repack → CKKS GROUP BY → SUM
//
// SELECT 100.00 * SUM(CASE WHEN p_type LIKE 'PROMO%'
//        THEN l_extendedprice * (1 - l_discount) ELSE 0 END)
//        / SUM(l_extendedprice * (1 - l_discount)) AS promo_revenue
// FROM lineitem, part
// WHERE l_partkey = p_partkey
//   AND l_shipdate >= date '1995-09-01'
//   AND l_shipdate < date '1995-10-01'
//
// Simplified: p_type is encoded as integer categories
// Category 0 = PROMO, Category 1..M-1 = non-PROMO

#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace lbcrypto;
using namespace std;

// Lagrange indicator polynomial coefficients for category matching
// For M categories, l_j(x) = product_{k!=j} (x - k) / (j - k)
vector<vector<double>> BuildLagrangeCoeffs(size_t M) {
    vector<vector<double>> coeffs(M, vector<double>(M, 0.0));
    for (size_t j = 0; j < M; j++) {
        vector<double> poly(1, 1.0);
        double denom = 1.0;
        for (size_t k = 0; k < M; k++) {
            if (k == j) continue;
            vector<double> newPoly(poly.size() + 1, 0.0);
            for (size_t d = 0; d < poly.size(); d++) {
                newPoly[d] += -double(k) * poly[d];
                newPoly[d + 1] += poly[d];
            }
            poly.swap(newPoly);
            denom *= double(j) - double(k);
        }
        for (size_t d = 0; d < poly.size(); d++) {
            coeffs[j][d] = poly[d] / denom;
        }
    }
    return coeffs;
}

void query_evaluation(size_t rows)
{
    cout << "=== TPC-H Q14 | rows=" << rows << " ===" << endl;

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

    uint32_t ship_bits = 16;
    uint32_t ship_scale = numeric_limits<Lvl2::T>::digits - ship_bits - 1;

    // p_type: 4 categories (0=PROMO, 1,2,3=non-PROMO)
    size_t num_types = 4;

    vector<uint64_t> ship_data(rows);
    vector<int> ptype_data(rows);
    vector<double> price_data(rows), discount_data(rows);

    uniform_int_distribution<int> type_dist(0, num_types - 1);
    uniform_real_distribution<double> price_ddist(1.0, 100.0);
    uniform_real_distribution<double> disc_ddist(0.0, 0.1);

    auto generate_date = [&]() -> uint64_t {
        uniform_int_distribution<int> m(1, 12), d(1, 28);
        return d(eng) + 100*m(eng) + 10000*1995;
    };

    for (size_t i = 0; i < rows; i++) {
        ship_data[i] = generate_date();
        ptype_data[i] = type_dist(eng);
        price_data[i] = price_ddist(eng);
        discount_data[i] = disc_ddist(eng);
    }

    // =================== TFHE Encryption ===================
    cout << "Encrypting with TFHE..." << endl;
    vector<TLWELvl2> ship_ct(rows);
    for (size_t i = 0; i < rows; i++) {
        ship_ct[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(
            ship_data[i], Lvl2::α, pow(2., ship_scale), sk.key.get<Lvl2>());
    }

    // Predicates: ship >= 19950901 AND ship < 19951001
    Lvl2::T pred1 = 19501, pred2 = 11001; // Simplified: use MMDD encoding
    // Actually for 16-bit encoding, use smaller date range
    pred1 = 901;   // September 1
    pred2 = 1001;  // October 1

    // Re-encode dates to fit 16-bit: use MMDD format
    for (size_t i = 0; i < rows; i++) {
        ship_data[i] = ship_data[i] % 10000; // extract MMDD
    }
    // Re-encrypt with MMDD encoding
    uint32_t date_bits = 16;
    uint32_t date_scale = numeric_limits<Lvl2::T>::digits - date_bits - 1;
    for (size_t i = 0; i < rows; i++) {
        ship_ct[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(
            ship_data[i], Lvl2::α, pow(2., date_scale), sk.key.get<Lvl2>());
    }
    auto ct_pred1 = TFHEpp::tlweSymInt32Encrypt<Lvl2>(pred1, Lvl2::α, pow(2., date_scale), sk.key.get<Lvl2>());
    auto ct_pred2 = TFHEpp::tlweSymInt32Encrypt<Lvl2>(pred2, Lvl2::α, pow(2., date_scale), sk.key.get<Lvl2>());

    // Plaintext predicate evaluation
    vector<uint32_t> pred_plain(rows);
    for (size_t i = 0; i < rows; i++) {
        pred_plain[i] = (ship_data[i] >= pred1 && ship_data[i] < pred2) ? 1 : 0;
    }

    // =================== Predicate Evaluation (ETHMSB) ===================
    cout << "Evaluating predicates (ETHMSB)..." << endl;
    vector<TLWELvl1> pred_cres(rows);
    auto t_filter_start = chrono::high_resolution_clock::now();

    vector<TLWELvl1> cres1(rows), cres2(rows);
    for (size_t i = 0; i < rows; i++) {
        ethmsb_greater_than_equal<Lvl2>(ship_ct[i], ct_pred1, cres1[i], date_bits, ek, LOGIC);
        ethmsb_less_than<Lvl2>(ship_ct[i], ct_pred2, cres2[i], date_bits, ek, LOGIC);
        HomAND(pred_cres[i], cres1[i], cres2[i], ek, ARITHMETIC);
    }

    auto t_filter_end = chrono::high_resolution_clock::now();
    double filter_ms = chrono::duration_cast<chrono::milliseconds>(t_filter_end - t_filter_start).count();

    // Rescale for repack
    uint32_t rlwe_scale_bits = 29;
    for (size_t i = 0; i < rows; i++) {
        TFHEpp::ari_rescale(pred_cres[i], pred_cres[i], rlwe_scale_bits, ek);
    }

    // Verify predicates
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
    params.SetMultiplicativeDepth(15);
    params.SetScalingModSize(50);
    params.SetBatchSize(rows);
    params.SetSecurityLevel(HEStd_128_classic);

    auto cc = GenCryptoContext(params);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    vector<int32_t> rot_indices;
    for (size_t step = 1; step < rows; step <<= 1) {
        rot_indices.push_back((int32_t)step);
        rot_indices.push_back(-(int32_t)step);
    }
    cc->EvalRotateKeyGen(keys.secretKey, rot_indices);

    // =================== Repack: TFHE → CKKS ===================
    cout << "Repacking TFHE→CKKS..." << endl;
    auto t_repack_start = chrono::high_resolution_clock::now();
    auto ct_mask = SimulatedRepack(cc, keys, pred_cres, sk, rlwe_scale_bits);
    auto t_repack_end = chrono::high_resolution_clock::now();
    double repack_ms = chrono::duration_cast<chrono::milliseconds>(t_repack_end - t_repack_start).count();

    // =================== CKKS GROUP BY (Lagrange) ===================
    cout << "CKKS GROUP BY (p_type)..." << endl;
    auto t_agg_start = chrono::high_resolution_clock::now();

    // Encrypt p_type column and price*(1-discount) column
    vector<double> ptype_double(rows), revenue(rows);
    for (size_t i = 0; i < rows; i++) {
        ptype_double[i] = (double)ptype_data[i];
        revenue[i] = price_data[i] * (1.0 - discount_data[i]);
    }
    auto ct_ptype = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(ptype_double));
    auto ct_revenue = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(revenue));

    // Build powers of ptype: ctB^0 = 1, ctB^1 = ptype, ctB^2, ...
    auto ones = vector<double>(rows, 1.0);
    auto ct_one = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(ones));

    size_t deg = num_types - 1;
    vector<Ciphertext<DCRTPoly>> powers(deg + 1);
    powers[0] = ct_one;
    if (deg >= 1) {
        powers[1] = ct_ptype;
        for (size_t d = 2; d <= deg; d++) {
            powers[d] = cc->EvalMult(powers[d - 1], ct_ptype);
        }
    }

    // Build Lagrange coefficients
    auto lagrange_coeffs = BuildLagrangeCoeffs(num_types);

    // For each group j: ctL_j = sum_d alpha[j][d] * powers[d]
    // Then ctMasked_j = ctL_j * ct_mask * ct_revenue
    // Then sum_j = RotateAndSum(ctMasked_j)
    vector<double> group_sums(num_types, 0.0);

    for (size_t j = 0; j < num_types; j++) {
        // Build mask for category j via Lagrange
        Ciphertext<DCRTPoly> ctL;
        bool first = true;
        for (size_t d = 0; d <= deg; d++) {
            double c = lagrange_coeffs[j][d];
            if (fabs(c) < 1e-12) continue;
            auto term = cc->EvalMult(powers[d], c);
            if (first) { ctL = term; first = false; }
            else { ctL = cc->EvalAdd(ctL, term); }
        }

        // Apply filter mask and revenue
        auto ctTmp = cc->EvalMult(ctL, ct_mask);
        ctTmp = cc->EvalMult(ctTmp, ct_revenue);

        // RotateAndSum
        auto ctSum = EvalSumSlots(cc, ctTmp, rows);

        Plaintext ptRes;
        cc->Decrypt(keys.secretKey, ctSum, &ptRes);
        ptRes->SetLength(1);
        group_sums[j] = ptRes->GetCKKSPackedValue()[0].real();
    }

    auto t_agg_end = chrono::high_resolution_clock::now();
    double agg_ms = chrono::duration_cast<chrono::milliseconds>(t_agg_end - t_agg_start).count();

    // =================== Compute plaintext result ===================
    double plain_promo_sum = 0, plain_total_sum = 0;
    for (size_t i = 0; i < rows; i++) {
        if (pred_plain[i]) {
            double r = price_data[i] * (1.0 - discount_data[i]);
            plain_total_sum += r;
            if (ptype_data[i] == 0) plain_promo_sum += r; // PROMO
        }
    }

    double promo_revenue = (plain_total_sum > 0) ? 100.0 * plain_promo_sum / plain_total_sum : 0;
    double enc_total = 0;
    for (size_t j = 0; j < num_types; j++) enc_total += group_sums[j];
    double enc_promo_revenue = (enc_total > 0) ? 100.0 * group_sums[0] / enc_total : 0;

    cout << "\n=== Results ===" << endl;
    for (size_t j = 0; j < num_types; j++) {
        double psum = 0;
        for (size_t i = 0; i < rows; i++)
            if (pred_plain[i] && ptype_data[i] == (int)j)
                psum += price_data[i] * (1.0 - discount_data[i]);
        cout << "  Group " << j << ": plain=" << psum << ", encrypted=" << group_sums[j]
             << ", err=" << abs(psum - group_sums[j]) << endl;
    }
    cout << "  Promo revenue (plain):     " << promo_revenue << "%" << endl;
    cout << "  Promo revenue (encrypted): " << enc_promo_revenue << "%" << endl;
    cout << "\n=== Timing ===" << endl;
    cout << "  Filter (TFHE):    " << filter_ms / 1000 << " s" << endl;
    cout << "  Repack:           " << repack_ms / 1000 << " s" << endl;
    cout << "  GROUP BY + AGG:   " << agg_ms / 1000 << " s" << endl;
    cout << "  Total:            " << (filter_ms + repack_ms + agg_ms) / 1000 << " s" << endl;
}

int main()
{
    query_evaluation(16);
    return 0;
}
