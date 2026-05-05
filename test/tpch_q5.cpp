// TPC-H Query 5 — End-to-end encrypted query
// WHERE filtering in TFHE → repack → CKKS multi-table JOIN → GROUP BY → SUM
//
// Simplified semantics:
// SELECT n_name, SUM(l_extendedprice * (1 - l_discount)) AS revenue
// FROM customer, orders, lineitem, supplier, nation, region
// WHERE c_custkey = o_custkey AND l_orderkey = o_orderkey
//   AND l_suppkey = s_suppkey AND c_nationkey = s_nationkey
//   AND s_nationkey = n_nationkey AND n_regionkey = r_regionkey
//   AND r_name = 'ASIA' (encoded as region 0)
//   AND o_orderdate >= date_lo AND o_orderdate < date_hi
// GROUP BY n_name

#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include <vector>
#include "../src/HEDB/comparison/HomCompare.h"
#include "../src/HEDB/conversion/repack_openfhe.h"

using namespace HEDB;
using namespace lbcrypto;
using namespace std;

// Lagrange utilities (same as Q3)
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

vector<Ciphertext<DCRTPoly>> BuildPowers(
    CryptoContext<DCRTPoly> &cc, const KeyPair<DCRTPoly> &keys,
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

// Lookup JOIN: left ⋈ right on key, returning right-table data column
// filled into left-table positions
Ciphertext<DCRTPoly> LookupJoin(
    CryptoContext<DCRTPoly> &cc, const KeyPair<DCRTPoly> &keys,
    Ciphertext<DCRTPoly> &ctK_R, Ciphertext<DCRTPoly> &ctK_S,
    Ciphertext<DCRTPoly> &ctY, size_t M, size_t nSlots)
{
    auto alpha = BuildLagrangeCoeffs(M);
    size_t deg = M - 1;
    auto powers_R = BuildPowers(cc, keys, ctK_R, deg, nSlots);
    auto powers_S = BuildPowers(cc, keys, ctK_S, deg, nSlots);

    // Pre-compute right-table basis-data broadcast
    vector<Ciphertext<DCRTPoly>> ctBase(M);
    for (size_t q = 0; q < M; q++) {
        auto tmp = cc->EvalMult(powers_S[q], ctY);
        ctBase[q] = EvalSumSlots(cc, tmp, nSlots);
    }

    Ciphertext<DCRTPoly> ctJoin;
    bool first = true;
    for (size_t j = 0; j < M; j++) {
        auto mask_R = LagrangeMask(cc, powers_R, alpha[j]);
        Ciphertext<DCRTPoly> ctPay;
        bool pfirst = true;
        for (size_t q = 0; q < M; q++) {
            if (fabs(alpha[j][q]) < 1e-12) continue;
            auto term = cc->EvalMult(ctBase[q], alpha[j][q]);
            if (pfirst) { ctPay = term; pfirst = false; }
            else { ctPay = cc->EvalAdd(ctPay, term); }
        }
        auto joined = cc->EvalMult(mask_R, ctPay);
        if (first) { ctJoin = joined; first = false; }
        else { ctJoin = cc->EvalAdd(ctJoin, joined); }
    }
    return ctJoin;
}

void query_evaluation(size_t rows)
{
    cout << "=== TPC-H Q5 | rows=" << rows << " ===" << endl;

    // Key domain sizes (small-scale)
    const size_t M_region = 2;  // regions
    const size_t M_nation = 4;  // nations
    const size_t M_supp = 4;    // suppliers
    const size_t M_cust = 4;    // customer keys
    const size_t M_order = 4;   // order keys

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
    cout << "Generating data (6 tables)..." << endl;
    default_random_engine eng(42);
    uniform_int_distribution<uint32_t> date_dist(100, 200);
    uniform_real_distribution<double> price_dist(1.0, 50.0);
    uniform_real_distribution<double> disc_dist(0.0, 0.1);
    uint32_t date_lo = 130, date_hi = 170;

    // region: regionkey → name (0=ASIA, 1=EUROPE)
    // nation: nationkey → regionkey, name
    vector<uint32_t> n_regionkey(M_nation);
    for (size_t i = 0; i < M_nation; i++)
        n_regionkey[i] = i % M_region; // nations alternate between regions

    // supplier: suppkey → nationkey
    vector<uint32_t> s_nationkey(rows);
    for (size_t i = 0; i < rows; i++)
        s_nationkey[i] = i % M_nation;

    // customer: custkey → nationkey
    vector<uint32_t> c_custkey(rows), c_nationkey(rows);
    for (size_t i = 0; i < rows; i++) {
        c_custkey[i] = i % M_cust;
        c_nationkey[i] = i % M_nation;
    }

    // orders: custkey, orderkey, orderdate
    vector<uint32_t> o_custkey(rows), o_orderkey(rows), o_orderdate(rows);
    for (size_t i = 0; i < rows; i++) {
        o_custkey[i] = i % M_cust;
        o_orderkey[i] = i % M_order;
        o_orderdate[i] = date_dist(eng);
    }

    // lineitem: orderkey, suppkey, price, discount
    vector<uint32_t> l_orderkey(rows), l_suppkey(rows);
    vector<double> l_price(rows), l_discount(rows);
    for (size_t i = 0; i < rows; i++) {
        l_orderkey[i] = i % M_order;
        l_suppkey[i] = i % M_supp;
        l_price[i] = price_dist(eng);
        l_discount[i] = disc_dist(eng);
    }

    // =================== TFHE Predicates ===================
    cout << "Evaluating TFHE predicates..." << endl;
    uint32_t date_bits = 8;
    uint32_t date_scale = numeric_limits<Lvl1::T>::digits - date_bits - 1;
    uint32_t rlwe_scale_bits = 29;

    // Predicate: o_orderdate >= date_lo AND o_orderdate < date_hi
    vector<TLWELvl1> odate_pred(rows);
    vector<uint32_t> odate_plain(rows);
    auto ct_dlo = TFHEpp::tlweSymInt32Encrypt<Lvl1>(date_lo, Lvl1::α, pow(2., date_scale), sk.key.get<Lvl1>());
    auto ct_dhi = TFHEpp::tlweSymInt32Encrypt<Lvl1>(date_hi, Lvl1::α, pow(2., date_scale), sk.key.get<Lvl1>());

    auto t_filter_start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < rows; i++) {
        odate_plain[i] = (o_orderdate[i] >= date_lo && o_orderdate[i] < date_hi) ? 1 : 0;
        auto ct_od = TFHEpp::tlweSymInt32Encrypt<Lvl1>(
            o_orderdate[i], Lvl1::α, pow(2., date_scale), sk.key.get<Lvl1>());
        TLWELvl1 cr1, cr2;
        ethmsb_greater_than_equal<Lvl1>(ct_od, ct_dlo, cr1, date_bits, ek, LOGIC);
        ethmsb_less_than<Lvl1>(ct_od, ct_dhi, cr2, date_bits, ek, LOGIC);
        HomAND(odate_pred[i], cr1, cr2, ek, ARITHMETIC);
    }
    auto t_filter_end = chrono::high_resolution_clock::now();
    double filter_ms = chrono::duration_cast<chrono::milliseconds>(t_filter_end - t_filter_start).count();

    for (size_t i = 0; i < rows; i++)
        TFHEpp::ari_rescale(odate_pred[i], odate_pred[i], rlwe_scale_bits, ek);

    cout << "  Filter time: " << filter_ms << " ms" << endl;

    // =================== CKKS Setup ===================
    cout << "Setting up OpenFHE CKKS..." << endl;
    CCParams<CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(30); // repack + multi-JOIN chain + GROUP BY
    params.SetScalingModSize(50);
    params.SetScalingTechnique(FIXEDAUTO);
    params.SetSecurityLevel(HEStd_NotSet);
    params.SetRingDim(32768);
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

    // =================== Repack ===================
    cout << "Repacking order date mask to CKKS..." << endl;
    auto t_repack_start = chrono::high_resolution_clock::now();
    auto ct_odate_mask = LWEsToOpenFHE(cc, keys, odate_pred, sk, rows, rlwe_scale_bits);
    auto t_repack_end = chrono::high_resolution_clock::now();
    double repack_ms = chrono::duration_cast<chrono::milliseconds>(t_repack_end - t_repack_start).count();
    cout << "  Repack time: " << repack_ms << " ms" << endl;

    // =================== Encrypt columns ===================
    auto enc_col = [&](const vector<uint32_t> &v) {
        vector<double> d(rows);
        for (size_t i = 0; i < rows; i++) d[i] = (double)v[i];
        return cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(d));
    };

    // Region filter: r_name == 'ASIA' (region 0)
    // nation.regionkey column → build a mask of nations in ASIA
    // For simplicity: encode nation's region membership as plaintext mask
    vector<double> nation_asia_mask(rows);
    for (size_t i = 0; i < rows; i++)
        nation_asia_mask[i] = (n_regionkey[s_nationkey[i] % M_nation] == 0) ? 1.0 : 0.0;
    auto ct_asia_mask = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(nation_asia_mask));

    // Key columns for JOINs
    auto ctCK_o = enc_col(o_custkey);    // orders.custkey
    auto ctCK_c = enc_col(c_custkey);    // customer.custkey
    auto ctOK_l = enc_col(l_orderkey);   // lineitem.orderkey
    auto ctOK_o = enc_col(o_orderkey);   // orders.orderkey
    auto ctSK_l = enc_col(l_suppkey);    // lineitem.suppkey
    auto ctNK_s = enc_col(s_nationkey);  // supplier.nationkey
    auto ctNK_c = enc_col(c_nationkey);  // customer.nationkey

    // Revenue
    vector<double> revenue(rows);
    for (size_t i = 0; i < rows; i++)
        revenue[i] = l_price[i] * (1.0 - l_discount[i]);
    auto ctRevenue = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(revenue));

    // =================== Multi-JOIN chain ===================
    cout << "Performing multi-table JOIN chain..." << endl;
    auto t_join_start = chrono::high_resolution_clock::now();

    // JOIN 1: orders ⋈ customer on custkey → propagate customer nationkey to orders
    auto ctCustNation = LookupJoin(cc, keys, ctCK_o, ctCK_c, ctNK_c, M_cust, rows);

    // Combined mask: orderdate filter AND region(nation) filter
    auto ctMask = cc->EvalMult(ct_odate_mask, ct_asia_mask);

    // JOIN 2: lineitem ⋈ orders on orderkey → propagate order mask to lineitem
    auto ctOrderMask = LookupJoin(cc, keys, ctOK_l, ctOK_o, ctMask, M_order, rows);

    // Apply all masks to revenue
    auto ctFiltered = cc->EvalMult(ctOrderMask, ctRevenue);

    // GROUP BY nation (using supplier.nationkey as proxy for local supplier)
    auto alpha = BuildLagrangeCoeffs(M_nation);
    auto powers_n = BuildPowers(cc, keys, ctNK_s, M_nation - 1, rows);

    vector<double> group_sums(M_nation, 0.0);
    for (size_t j = 0; j < M_nation; j++) {
        auto grp_mask = LagrangeMask(cc, powers_n, alpha[j]);
        auto ctGrp = cc->EvalMult(grp_mask, ctFiltered);
        auto ctSum = EvalSumSlots(cc, ctGrp, rows);
        Plaintext pt;
        cc->Decrypt(keys.secretKey, ctSum, &pt);
        pt->SetLength(1);
        group_sums[j] = pt->GetCKKSPackedValue()[0].real();
    }

    auto t_join_end = chrono::high_resolution_clock::now();
    double join_ms = chrono::duration_cast<chrono::milliseconds>(t_join_end - t_join_start).count();

    // =================== Results ===================
    cout << "\n=== Results (by nation) ===" << endl;
    for (size_t j = 0; j < M_nation; j++) {
        cout << "  Nation " << j << " (region="
             << n_regionkey[j] << "): revenue=" << group_sums[j] << endl;
    }

    cout << "\n=== Timing ===" << endl;
    cout << "  Filter (TFHE):  " << filter_ms / 1000 << " s" << endl;
    cout << "  Repack:         " << repack_ms / 1000 << " s" << endl;
    cout << "  JOIN+GROUP BY:  " << join_ms / 1000 << " s" << endl;
    cout << "  Total:          " << (filter_ms + repack_ms + join_ms) / 1000 << " s" << endl;
}

int main()
{
    query_evaluation(16);
    return 0;
}
