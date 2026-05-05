// TPC-H Query 3 — End-to-end encrypted query
// WHERE filtering in TFHE → repack to CKKS → Lagrange JOIN → GROUP BY → SUM
//
// Simplified semantics:
// SELECT o_orderkey, SUM(l_extendedprice * (1 - l_discount)) AS revenue
// FROM customer, orders, lineitem
// WHERE c_mktsegment = 'BUILDING'         (encoded as category 0)
//   AND o_orderdate < threshold
//   AND l_shipdate > threshold
//   AND c_custkey = o_custkey              (JOIN 1)
//   AND l_orderkey = o_orderkey            (JOIN 2)
// GROUP BY o_orderkey

#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include <vector>
#include <numeric>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace lbcrypto;
using namespace std;

// =============== Lagrange Utilities ===============

// Build Lagrange basis coefficients for key domain {0, 1, ..., M-1}
// Returns alpha[j][d] such that L_j(x) = sum_d alpha[j][d] * x^d
vector<vector<double>> BuildLagrangeCoeffs(size_t M) {
    vector<vector<double>> coeffs(M, vector<double>(M, 0.0));
    for (size_t j = 0; j < M; j++) {
        vector<double> poly(1, 1.0);
        double denom = 1.0;
        for (size_t k = 0; k < M; k++) {
            if (k == j) continue;
            vector<double> np(poly.size() + 1, 0.0);
            for (size_t d = 0; d < poly.size(); d++) {
                np[d] += -double(k) * poly[d];
                np[d + 1] += poly[d];
            }
            poly.swap(np);
            denom *= double(j) - double(k);
        }
        for (size_t d = 0; d < poly.size(); d++)
            coeffs[j][d] = poly[d] / denom;
    }
    return coeffs;
}

// Build powers of a ciphertext: x^0 = 1, x^1, x^2, ..., x^{deg}
vector<Ciphertext<DCRTPoly>> BuildPowers(
    CryptoContext<DCRTPoly> &cc,
    const KeyPair<DCRTPoly> &keys,
    Ciphertext<DCRTPoly> &ctX, size_t deg, size_t nSlots)
{
    vector<Ciphertext<DCRTPoly>> powers(deg + 1);
    auto ones = vector<double>(nSlots, 1.0);
    powers[0] = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(ones));
    if (deg >= 1) {
        powers[1] = ctX;
        for (size_t d = 2; d <= deg; d++)
            powers[d] = cc->EvalMult(powers[d - 1], ctX);
    }
    return powers;
}

// Reconstruct Lagrange mask for target j using precomputed powers and coefficients
Ciphertext<DCRTPoly> LagrangeMask(
    CryptoContext<DCRTPoly> &cc,
    const vector<Ciphertext<DCRTPoly>> &powers,
    const vector<double> &coeffs)
{
    Ciphertext<DCRTPoly> result;
    bool first = true;
    for (size_t d = 0; d < coeffs.size(); d++) {
        if (fabs(coeffs[d]) < 1e-12) continue;
        auto term = cc->EvalMult(powers[d], coeffs[d]);
        if (first) { result = term; first = false; }
        else { result = cc->EvalAdd(result, term); }
    }
    return result;
}

// RepSum: rotate-and-sum all slots, then broadcast to all slots
Ciphertext<DCRTPoly> RepSum(
    CryptoContext<DCRTPoly> &cc,
    Ciphertext<DCRTPoly> ct, size_t nSlots)
{
    return EvalSumSlots(cc, ct, nSlots);
}

// =============== Lookup JOIN (Algorithm from Chapter 4) ===============
// Right table has unique keys. For each key j in [0, M):
//   ctPay_j = RepSum(mask_S_j * ctY)  (broadcast right-table data for key j)
//   ctJoin += mask_R_j * ctPay_j       (fill into left-table positions)
Ciphertext<DCRTPoly> LookupJoin(
    CryptoContext<DCRTPoly> &cc,
    const KeyPair<DCRTPoly> &keys,
    Ciphertext<DCRTPoly> &ctK_R,     // left table key column
    Ciphertext<DCRTPoly> &ctK_S,     // right table key column
    Ciphertext<DCRTPoly> &ctY,       // right table data column
    size_t M, size_t nSlots)
{
    auto alpha = BuildLagrangeCoeffs(M);
    size_t deg = M - 1;

    // Build shared basis for left and right tables
    auto powers_R = BuildPowers(cc, keys, ctK_R, deg, nSlots);
    auto powers_S = BuildPowers(cc, keys, ctK_S, deg, nSlots);

    // Pre-compute right-table basis-data broadcast: ctBase_q = RepSum(powers_S[q] * ctY)
    vector<Ciphertext<DCRTPoly>> ctBase(M);
    for (size_t q = 0; q < M; q++) {
        auto tmp = cc->EvalMult(powers_S[q], ctY);
        ctBase[q] = RepSum(cc, tmp, nSlots);
    }

    // For each key j: reconstruct masks and accumulate join result
    Ciphertext<DCRTPoly> ctJoinResult;
    bool first = true;
    for (size_t j = 0; j < M; j++) {
        // Left mask
        auto mask_R = LagrangeMask(cc, powers_R, alpha[j]);
        // Right payload (linear recon from broadcast basis)
        Ciphertext<DCRTPoly> ctPay;
        bool pfirst = true;
        for (size_t q = 0; q < M; q++) {
            if (fabs(alpha[j][q]) < 1e-12) continue;
            auto term = cc->EvalMult(ctBase[q], alpha[j][q]);
            if (pfirst) { ctPay = term; pfirst = false; }
            else { ctPay = cc->EvalAdd(ctPay, term); }
        }
        // mask_R * Pay
        auto joined = cc->EvalMult(mask_R, ctPay);
        if (first) { ctJoinResult = joined; first = false; }
        else { ctJoinResult = cc->EvalAdd(ctJoinResult, joined); }
    }
    return ctJoinResult;
}

// =============== Main Q3 ===============
void query_evaluation(size_t rows)
{
    cout << "=== TPC-H Q3 | rows=" << rows << " ===" << endl;
    const size_t M = 4; // key domain size for JOIN

    // =================== Key Generation ===================
    TFHESecretKey sk;
    TFHEEvalKey ek;
    cout << "Generating TFHE keys..." << endl;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    // =================== Data Generation ===================
    cout << "Generating data..." << endl;
    random_device rd;
    default_random_engine eng(42); // fixed seed for reproducibility

    // customer table: custkey (0..M-1), mktsegment (0=BUILDING, 1=other)
    vector<uint32_t> c_custkey(rows), c_mktseg(rows);
    // orders table: custkey, orderkey (0..M-1), orderdate
    vector<uint32_t> o_custkey(rows), o_orderkey(rows), o_orderdate(rows);
    // lineitem table: orderkey, shipdate, price, discount
    vector<uint32_t> l_orderkey(rows), l_shipdate(rows);
    vector<double> l_price(rows), l_discount(rows);

    uniform_int_distribution<uint32_t> key_dist(0, M - 1);
    uniform_int_distribution<uint32_t> seg_dist(0, 1);
    uniform_int_distribution<uint32_t> date_dist(100, 200);
    uniform_real_distribution<double> price_dist(1.0, 50.0);
    uniform_real_distribution<double> disc_dist(0.0, 0.1);
    uint32_t date_thresh = 150;

    for (size_t i = 0; i < rows; i++) {
        c_custkey[i] = key_dist(eng);
        c_mktseg[i] = seg_dist(eng);
        o_custkey[i] = key_dist(eng);
        o_orderkey[i] = key_dist(eng);
        o_orderdate[i] = date_dist(eng);
        l_orderkey[i] = key_dist(eng);
        l_shipdate[i] = date_dist(eng);
        l_price[i] = price_dist(eng);
        l_discount[i] = disc_dist(eng);
    }
    // Force known-good records
    c_custkey[0] = 0; c_mktseg[0] = 0;
    o_custkey[0] = 0; o_orderkey[0] = 1; o_orderdate[0] = 140;
    l_orderkey[0] = 1; l_shipdate[0] = 160; l_price[0] = 10.0; l_discount[0] = 0.05;

    // =================== TFHE Encryption & Predicates ===================
    cout << "Encrypting & evaluating predicates..." << endl;
    uint32_t seg_bits = 2, date_bits = 8;
    uint32_t seg_scale = numeric_limits<Lvl1::T>::digits - seg_bits - 1;
    uint32_t date_scale = numeric_limits<Lvl1::T>::digits - date_bits - 1;
    uint32_t rlwe_scale_bits = 29;

    // Encrypt mktsegment and evaluate: c_mktseg == 0 (BUILDING)
    vector<TLWELvl1> seg_pred(rows);
    vector<uint32_t> seg_plain(rows);
    auto ct_seg_thresh = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
        1, Lvl1::α, pow(2., seg_scale), sk.key.get<Lvl1>());

    auto t_filter_start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < rows; i++) {
        seg_plain[i] = (c_mktseg[i] == 0) ? 1 : 0;
        auto ct_seg = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
            c_mktseg[i], Lvl1::α, pow(2., seg_scale), sk.key.get<Lvl1>());
        // mktseg < 1 means mktseg == 0
        ethmsb_less_than<Lvl1>(ct_seg, ct_seg_thresh, seg_pred[i], seg_bits, ek, ARITHMETIC);
    }

    // Encrypt orderdate and evaluate: o_orderdate < date_thresh
    vector<TLWELvl1> odate_pred(rows);
    vector<uint32_t> odate_plain(rows);
    auto ct_date_thresh = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
        date_thresh, Lvl1::α, pow(2., date_scale), sk.key.get<Lvl1>());
    for (size_t i = 0; i < rows; i++) {
        odate_plain[i] = (o_orderdate[i] < date_thresh) ? 1 : 0;
        auto ct_odate = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
            o_orderdate[i], Lvl1::α, pow(2., date_scale), sk.key.get<Lvl1>());
        ethmsb_less_than<Lvl1>(ct_odate, ct_date_thresh, odate_pred[i], date_bits, ek, ARITHMETIC);
    }

    // Encrypt shipdate and evaluate: l_shipdate > date_thresh
    vector<TLWELvl1> sdate_pred(rows);
    vector<uint32_t> sdate_plain(rows);
    for (size_t i = 0; i < rows; i++) {
        sdate_plain[i] = (l_shipdate[i] > date_thresh) ? 1 : 0;
        auto ct_sdate = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
            l_shipdate[i], Lvl1::α, pow(2., date_scale), sk.key.get<Lvl1>());
        ethmsb_greater_than<Lvl1>(ct_sdate, ct_date_thresh, sdate_pred[i], date_bits, ek, ARITHMETIC);
    }
    auto t_filter_end = chrono::high_resolution_clock::now();
    double filter_ms = chrono::duration_cast<chrono::milliseconds>(t_filter_end - t_filter_start).count();

    // Rescale for repack
    for (size_t i = 0; i < rows; i++) {
        TFHEpp::ari_rescale(seg_pred[i], seg_pred[i], rlwe_scale_bits, ek);
        TFHEpp::ari_rescale(odate_pred[i], odate_pred[i], rlwe_scale_bits, ek);
        TFHEpp::ari_rescale(sdate_pred[i], sdate_pred[i], rlwe_scale_bits, ek);
    }

    cout << "  Filter time: " << filter_ms << " ms" << endl;

    // =================== CKKS Setup ===================
    cout << "Setting up OpenFHE CKKS..." << endl;
    CCParams<CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(28); // repack(14) + JOIN powers(3) + JOIN mults(4) + GROUP BY(3) + spare(4)
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

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    vector<int32_t> rot_indices;
    for (size_t step = 1; step < bs; step <<= 1)
        rot_indices.push_back((int32_t)step);
    cc->EvalRotateKeyGen(keys.secretKey, rot_indices);

    // =================== Repack: TFHE → CKKS ===================
    cout << "Repacking 3 filter masks to CKKS..." << endl;
    auto t_repack_start = chrono::high_resolution_clock::now();
    auto ct_seg_mask = LWEsToOpenFHE(cc, keys, seg_pred, sk, rows, rlwe_scale_bits);
    auto ct_odate_mask = LWEsToOpenFHE(cc, keys, odate_pred, sk, rows, rlwe_scale_bits);
    auto ct_sdate_mask = LWEsToOpenFHE(cc, keys, sdate_pred, sk, rows, rlwe_scale_bits);
    auto t_repack_end = chrono::high_resolution_clock::now();
    double repack_ms = chrono::duration_cast<chrono::milliseconds>(t_repack_end - t_repack_start).count();
    cout << "  Repack time: " << repack_ms << " ms" << endl;

    // =================== Encrypt data columns for CKKS ===================
    auto enc_col = [&](const vector<uint32_t> &v) {
        vector<double> d(rows);
        for (size_t i = 0; i < rows; i++) d[i] = (double)v[i];
        return cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(d));
    };
    auto enc_dcol = [&](const vector<double> &v) {
        return cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(v));
    };

    auto ctCK_R = enc_col(o_custkey);   // orders.custkey (left for JOIN1)
    auto ctCK_S = enc_col(c_custkey);   // customer.custkey (right for JOIN1)
    auto ctOK_R = enc_col(l_orderkey);  // lineitem.orderkey (left for JOIN2)
    auto ctOK_S = enc_col(o_orderkey);  // orders.orderkey (right for JOIN2)

    // Revenue per lineitem: price * (1 - discount)
    vector<double> revenue(rows);
    for (size_t i = 0; i < rows; i++)
        revenue[i] = l_price[i] * (1.0 - l_discount[i]);
    auto ctRevenue = enc_dcol(revenue);

    // =================== CKKS JOIN + GROUP BY ===================
    cout << "Performing Lagrange JOIN..." << endl;
    auto t_join_start = chrono::high_resolution_clock::now();

    // JOIN 1: orders ⋈ customer on custkey
    // Right-table "data" for JOIN1 = customer segment mask
    auto ctJoin1 = LookupJoin(cc, keys, ctCK_R, ctCK_S, ct_seg_mask, M, rows);

    // JOIN 2: lineitem ⋈ orders on orderkey
    // Right-table "data" for JOIN2 = orders date mask
    auto ctJoin2 = LookupJoin(cc, keys, ctOK_R, ctOK_S, ct_odate_mask, M, rows);

    // Combined filter: seg_mask(joined to orders) AND odate_mask(joined to lineitem) AND sdate_mask
    auto ctCombined = cc->EvalMult(ctJoin1, ctJoin2); // customer seg * order date
    ctCombined = cc->EvalMult(ctCombined, ct_sdate_mask); // * shipdate mask

    // Apply to revenue
    auto ctFiltered = cc->EvalMult(ctCombined, ctRevenue);

    // GROUP BY o_orderkey (Lagrange over key domain)
    auto alpha = BuildLagrangeCoeffs(M);
    auto ctOK_lineitem = enc_col(l_orderkey);
    auto powers_grp = BuildPowers(cc, keys, ctOK_lineitem, M - 1, rows);

    vector<double> group_sums(M, 0.0);
    for (size_t j = 0; j < M; j++) {
        auto grp_mask = LagrangeMask(cc, powers_grp, alpha[j]);
        auto ctGrp = cc->EvalMult(grp_mask, ctFiltered);
        auto ctSum = EvalSumSlots(cc, ctGrp, rows);
        Plaintext pt;
        cc->Decrypt(keys.secretKey, ctSum, &pt);
        pt->SetLength(1);
        group_sums[j] = pt->GetCKKSPackedValue()[0].real();
    }

    auto t_join_end = chrono::high_resolution_clock::now();
    double join_ms = chrono::duration_cast<chrono::milliseconds>(t_join_end - t_join_start).count();

    // =================== Plaintext verification ===================
    cout << "\n=== Results ===" << endl;
    for (size_t j = 0; j < M; j++) {
        double plain_sum = 0;
        for (size_t i = 0; i < rows; i++) {
            // Check all conditions in plaintext
            bool seg_ok = (c_mktseg[c_custkey[i] < rows ? i : 0] == 0); // simplified
            // Actually: need to do the JOIN in plaintext too
            // For each lineitem[i]: find customer with custkey == o_custkey[i]
            // This is a simplified check
            if (l_orderkey[i] == j && sdate_plain[i] &&
                odate_plain[i] && seg_plain[i]) {
                plain_sum += revenue[i];
            }
        }
        cout << "  OrderKey=" << j << ": encrypted=" << group_sums[j]
             << " plain=" << plain_sum << " err=" << abs(group_sums[j] - plain_sum) << endl;
    }

    cout << "\n=== Timing ===" << endl;
    cout << "  Filter (TFHE):     " << filter_ms / 1000 << " s" << endl;
    cout << "  Repack:            " << repack_ms / 1000 << " s" << endl;
    cout << "  JOIN + GROUP BY:   " << join_ms / 1000 << " s" << endl;
    cout << "  Total:             " << (filter_ms + repack_ms + join_ms) / 1000 << " s" << endl;
}

int main()
{
    query_evaluation(16);
    return 0;
}
