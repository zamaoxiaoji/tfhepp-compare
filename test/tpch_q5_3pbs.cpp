#include "../comparison/comparison.h"
#include "ckks_relational.h"
#include "ckks_repack.h"
#include "gate.hpp"
#include <iomanip>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>

using namespace tfhepp_compare;
using namespace tfhepp_compare::three_pbs;
using namespace seal;

/***
 * TPC-H Query 5
 * select
 *      n_name,
 *      sum(l_extendedprice * (1 - l_discount)) as revenue
 *  from
 *      customer, orders, lineitem, supplier, nation, region
 *  where
 *      c_custkey = o_custkey
 *      and l_orderkey = o_orderkey
 *      and l_suppkey = s_suppkey
 *      and c_nationkey = s_nationkey
 *      and s_nationkey = n_nationkey
 *      and n_regionkey = r_regionkey
 *      and r_name = ':1'
 *      and o_orderdate >= date ':2'
 *      and o_orderdate < date ':2' + interval '1' year
 *  group by
 *      n_name
 *  order by
 *      revenue desc;
 *
 *  Consider the pre-joined customer x orders x lineitem x supplier x nation
 *  x region table (HE3DB-style: the region join is folded into a per-row
 *  is_asia_region flag derived from n_regionkey == r_regionkey('ASIA')).
 *  Filtering uses the 3-PBS pruned comparator (Algorithm 2, Chapter 3 of
 *  Tang Li's thesis); CKKS aggregation reuses this repository's repack
 *  pipeline (PackLWEsToCKKS + HomomorphicRound).
 */

void lift_and_and(TLWELvl1 &cipher1, TLWELvl1 &cipher2, TLWELvl1 &res,
                  uint32_t scale_bits, TFHEEvalKey &ek)
{
    using namespace TFHEpp;
    TLWELvl1 temp;
    for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        temp[i] = cipher1[i] + cipher2[i];
    temp[Lvl1::k * Lvl1::n] -= Lvl1::μ;
    Lvl1::T c = (Lvl1::T(1) << (scale_bits - 1));
    TLWELvl0 tlwelvl0;
    IdentityKeySwitch<Lvl10>(tlwelvl0, temp, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01,
                                          μ_polygen<Lvl1>(-c));
    res[Lvl1::k * Lvl1::n] += c;
}

template <class P>
void greater_than(const TFHEpp::TLWE<P> &cipher1,
                  const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                  uint32_t plain_bits, const TFHEEvalKey &ek,
                  const FastB2AEvalKeyPack &micro_pack, bool result_type)
{
    TFHEpp::TLWE<P> sub;
    for (int i = 0; i <= P::k * P::n; i++) sub[i] = cipher2[i] - cipher1[i];
    HomMSB(res, sub, plain_bits + 1, ek, micro_pack, result_type);
}

template <class P>
void less_than(const TFHEpp::TLWE<P> &cipher1,
               const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
               uint32_t plain_bits, const TFHEEvalKey &ek,
               const FastB2AEvalKeyPack &micro_pack, bool result_type)
{
    TFHEpp::TLWE<P> sub;
    for (int i = 0; i <= P::k * P::n; i++) sub[i] = cipher1[i] - cipher2[i];
    HomMSB(res, sub, plain_bits + 1, ek, micro_pack, result_type);
}

template <class P>
void greater_than_equal(const TFHEpp::TLWE<P> &cipher1,
                        const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        const FastB2AEvalKeyPack &micro_pack,
                        bool result_type)
{
    TFHEpp::TLWE<P> sub;
    for (int i = 0; i <= P::k * P::n; i++) sub[i] = cipher1[i] - cipher2[i];
    HomMSB(res, sub, plain_bits + 1, ek, micro_pack, LOGIC);
    HomNOT<Lvl1>(res, res);
    if (IS_ARITHMETIC(result_type)) LOG_to_ARI(res, res, ek);
}

template <class P>
void less_than_equal(const TFHEpp::TLWE<P> &cipher1,
                     const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                     uint32_t plain_bits, const TFHEEvalKey &ek,
                     const FastB2AEvalKeyPack &micro_pack, bool result_type)
{
    TFHEpp::TLWE<P> sub;
    for (int i = 0; i <= P::k * P::n; i++) sub[i] = cipher2[i] - cipher1[i];
    HomMSB(res, sub, plain_bits + 1, ek, micro_pack, LOGIC);
    HomNOT<Lvl1>(res, res);
    if (IS_ARITHMETIC(result_type)) LOG_to_ARI(res, res, ek);
}

template <class P>
void equal(const TFHEpp::TLWE<P> &cipher1, const TFHEpp::TLWE<P> &cipher2,
           TLWELvl1 &res, uint32_t plain_bits, const TFHEEvalKey &ek,
           const FastB2AEvalKeyPack &micro_pack, bool result_type)
{
    TLWELvl1 greater_tlwe, less_tlwe;
    greater_than_equal<P>(cipher1, cipher2, greater_tlwe, plain_bits, ek,
                          micro_pack, LOGIC);
    less_than_equal<P>(cipher1, cipher2, less_tlwe, plain_bits, ek, micro_pack,
                       LOGIC);
    HomAND(res, greater_tlwe, less_tlwe, ek, result_type);
}

double relational_query5(size_t num)
{
    std::cout << "Relational SQL Query5 Test (3-PBS pruned + CKKS repack): "
              << std::endl;
    std::cout << "--------------------------------------------------------"
              << std::endl;
    std::cout << "Records: " << num << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    using bkP = Lvl01;
    using iksP = Lvl10;
    std::uniform_int_distribution<uint32_t> orderdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> nationkey_message(0, 3);
    std::uniform_int_distribution<uint32_t> regionkey_message(0, 1);
    std::uniform_int_distribution<uint64_t> revenue_message(0, 100);
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
    const auto micro_pack = GenerateFastB2AEvalKeyPack(sk, true);

    // Filtering (operating on the pre-joined customer x orders x lineitem
    // x supplier x nation x region table; HE3DB style.)
    std::vector<uint64_t> order_date(num);
    std::vector<uint64_t> c_nationkey(num), s_nationkey(num);
    std::vector<uint64_t> is_asia(num);  // n_regionkey == ASIA_REGION_KEY
    std::vector<TLWELvl2> orderdate_ciphers(num);
    std::vector<TLWELvl2> c_nation_ciphers(num), s_nation_ciphers(num);
    std::vector<TLWELvl2> is_asia_ciphers(num);
    std::vector<TLWELvl2> nation_group_ciphers(num);

    uint32_t num_bits = 16;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - num_bits - 1;

    TLWELvl2 predicate1_cipher, predicate2_cipher, asia_predicate;
    TLWELvl2 nation_predicate[4];
    uint64_t predicate1_value = 12000;  // date '1994-01-01' lower window
    uint64_t predicate2_value = 13000;  // date '1995-01-01' upper window
    uint64_t asia_predicate_value = 1;  // r_name = 'ASIA' → region key 1
    predicate1_cipher = tlweSymInt32Encrypt<Lvl2>(predicate1_value, Lvl2::α,
                                                  pow(2., scale_bits),
                                                  sk.key.get<Lvl2>());
    predicate2_cipher = tlweSymInt32Encrypt<Lvl2>(predicate2_value, Lvl2::α,
                                                  pow(2., scale_bits),
                                                  sk.key.get<Lvl2>());
    asia_predicate = tlweSymInt32Encrypt<Lvl2>(asia_predicate_value, Lvl2::α,
                                               pow(2., scale_bits),
                                               sk.key.get<Lvl2>());
    for (uint32_t g = 0; g < 4; g++) {
        nation_predicate[g] = tlweSymInt32Encrypt<Lvl2>(
            g, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    }

    // Start sql evaluation: one filter cipher per nation group.
    std::vector<TLWELvl1> filter_res_n0(num), filter_res_n1(num),
        filter_res_n2(num), filter_res_n3(num);

    std::vector<double> revenue(num);

    for (size_t i = 0; i < num; i++) {
        revenue[i] = revenue_message(engine);
    }

    for (size_t i = 0; i < num; i++) {
        // Generate data
        order_date[i] = orderdate_message(engine);
        c_nationkey[i] = nationkey_message(engine);
        s_nationkey[i] = nationkey_message(engine);
        is_asia[i] = regionkey_message(engine);
        orderdate_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(order_date[i],
                                                         Lvl2::α,
                                                         pow(2., scale_bits),
                                                         sk.key.get<Lvl2>());
        c_nation_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(c_nationkey[i],
                                                        Lvl2::α,
                                                        pow(2., scale_bits),
                                                        sk.key.get<Lvl2>());
        s_nation_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(s_nationkey[i],
                                                        Lvl2::α,
                                                        pow(2., scale_bits),
                                                        sk.key.get<Lvl2>());
        is_asia_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(is_asia[i], Lvl2::α,
                                                       pow(2., scale_bits),
                                                       sk.key.get<Lvl2>());
        // For grouping we use the customer-side nation key (which equals the
        // supplier-side under the join predicate).
        nation_group_ciphers[i] = c_nation_ciphers[i];
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

    for (size_t i = 0; i < num; i++) {
        TLWELvl1 base_filter, pre_res, group_eq;
        // base_filter = (o_orderdate > lo) AND (o_orderdate < hi)
        //               AND (c_nationkey == s_nationkey)
        //               AND (n_regionkey == ASIA)
        greater_than<Lvl2>(orderdate_ciphers[i], predicate1_cipher, base_filter,
                           num_bits, ek, micro_pack, LOGIC);
        less_than<Lvl2>(orderdate_ciphers[i], predicate2_cipher, pre_res,
                        num_bits, ek, micro_pack, LOGIC);
        TFHEpp::HomAND(base_filter, pre_res, base_filter, ek);
        equal<Lvl2>(c_nation_ciphers[i], s_nation_ciphers[i], pre_res, num_bits,
                    ek, micro_pack, LOGIC);
        TFHEpp::HomAND(base_filter, pre_res, base_filter, ek);
        equal<Lvl2>(is_asia_ciphers[i], asia_predicate, pre_res, num_bits, ek,
                    micro_pack, LOGIC);
        TFHEpp::HomAND(base_filter, pre_res, base_filter, ek);

        // group filter = base_filter AND (nationkey == g)
        TLWELvl1 *group_outputs[4] = {&filter_res_n0[i], &filter_res_n1[i],
                                       &filter_res_n2[i], &filter_res_n3[i]};
        for (uint32_t g = 0; g < 4; g++) {
            equal<Lvl2>(nation_group_ciphers[i], nation_predicate[g], group_eq,
                        num_bits, ek, micro_pack, LOGIC);
            lift_and_and(base_filter, group_eq, *group_outputs[g], 29, ek);
        }
    }
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start)
                         .count();

    std::vector<std::vector<uint64_t>> plain_filter_res(
        4, std::vector<uint64_t>(num, 0));
    std::vector<uint64_t> plain_agg_res(4, 0);
    for (size_t i = 0; i < num; i++) {
        if (order_date[i] > predicate1_value &&
            order_date[i] < predicate2_value &&
            c_nationkey[i] == s_nationkey[i] &&
            is_asia[i] == asia_predicate_value) {
            const uint64_t g = c_nationkey[i];
            if (g < 4) {
                plain_filter_res[g][i] = 1;
                plain_agg_res[g] += revenue[i];
            }
        }
    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = num;
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(
        poly_modulus_degree,
        {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45,
         45, 45, 45, 59}));
    double scale = std::pow(2.0, scale_bits);

    // context instance
    seal::SEALContext context(parms, true, seal::sec_level_type::none);

    // key generation
    seal::KeyGenerator keygen(context);
    seal::SecretKey seal_secret_key = keygen.secret_key();
    seal::PublicKey seal_public_key;
    keygen.create_public_key(seal_public_key);
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    std::vector<int> rotation_steps;
    for (size_t step = 1; step < poly_modulus_degree / 2; step <<= 1) {
        rotation_steps.push_back(static_cast<int>(step));
        rotation_steps.push_back(-static_cast<int>(step));
    }
    seal::GaloisKeys galois_keys;
    keygen.create_galois_keys(rotation_steps, galois_keys);

    // utils
    seal::Encryptor encryptor(context, seal_public_key);
    seal::Encryptor symmetric_encryptor(context, seal_secret_key);
    seal::Evaluator evaluator(context);
    seal::Decryptor decryptor(context, seal_secret_key);

    // encoder
    seal::CKKSEncoder ckks_encoder(context);

    // generate evaluation key
    std::cout << "Generating Conversion Key..." << std::endl;
    auto repack_config =
        tfhepp_ckks::DefaultRepackConfig<Lvl1>(scale_bits, modulus_bits);
    tfhepp_ckks::RepackEvaluationKey repack_key;
    tfhepp_ckks::GenerateRepackKey<Lvl1>(repack_key, sk, repack_config.key_scale,
                                          ckks_encoder, symmetric_encryptor,
                                          context);

    // conversion
    std::cout << "Starting Conversion..." << std::endl;
    seal::Ciphertext result_n0, result_n1, result_n2, result_n3;
    start = std::chrono::system_clock::now();
    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(result_n0, filter_res_n0, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(result_n0, result_n0.scale(), ckks_encoder,
                                   relin_keys, evaluator, context);

    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(result_n1, filter_res_n1, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(result_n1, result_n1.scale(), ckks_encoder,
                                   relin_keys, evaluator, context);

    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(result_n2, filter_res_n2, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(result_n2, result_n2.scale(), ckks_encoder,
                                   relin_keys, evaluator, context);

    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(result_n3, filter_res_n3, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(result_n3, result_n3.scale(), ckks_encoder,
                                   relin_keys, evaluator, context);
    end = std::chrono::system_clock::now();
    aggregation_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                           end - start)
                           .count();
    printf("Repack time = %f\n", aggregation_time);
    seal::Plaintext plain;
    std::vector<double> computed_n0(slots_count), computed_n1(slots_count),
        computed_n2(slots_count), computed_n3(slots_count);
    decryptor.decrypt(result_n0, plain);
    ckks_encoder.decode(plain, computed_n0);
    decryptor.decrypt(result_n1, plain);
    ckks_encoder.decode(plain, computed_n1);
    decryptor.decrypt(result_n2, plain);
    ckks_encoder.decode(plain, computed_n2);
    decryptor.decrypt(result_n3, plain);
    ckks_encoder.decode(plain, computed_n3);

    double err0 = 0., err1 = 0., err2 = 0., err3 = 0.;
    for (size_t i = 0; i < slots_count; ++i) {
        err0 += std::abs(computed_n0[i] - plain_filter_res[0][i]);
        err1 += std::abs(computed_n1[i] - plain_filter_res[1][i]);
        err2 += std::abs(computed_n2[i] - plain_filter_res[2][i]);
        err3 += std::abs(computed_n3[i] - plain_filter_res[3][i]);
    }
    printf("Repack nation0 average error = %f ~ 2^%.1f\n", err0 / slots_count,
           std::log2(err0 / slots_count));
    printf("Repack nation1 average error = %f ~ 2^%.1f\n", err1 / slots_count,
           std::log2(err1 / slots_count));
    printf("Repack nation2 average error = %f ~ 2^%.1f\n", err2 / slots_count,
           std::log2(err2 / slots_count));
    printf("Repack nation3 average error = %f ~ 2^%.1f\n", err3 / slots_count,
           std::log2(err3 / slots_count));

    // Filter result * data
    seal::Ciphertext revenue_cipher;
    double qd = parms.coeff_modulus()[result_n0.coeff_modulus_size() - 1]
                    .value();
    std::vector<double> revenue_slots(ckks_encoder.slot_count(), 0.);
    for (size_t i = 0; i < num; i++) revenue_slots[i] = revenue[i];
    ckks_encoder.encode(revenue_slots, result_n0.parms_id(), qd, plain);
    symmetric_encryptor.encrypt_symmetric(plain, revenue_cipher);

    std::cout << "Aggregating price and discount .." << std::endl;
    start = std::chrono::system_clock::now();
    seal::Ciphertext sum_revenue_n0, sum_revenue_n1, sum_revenue_n2,
        sum_revenue_n3;
    evaluator.multiply(result_n0, revenue_cipher, sum_revenue_n0);
    evaluator.relinearize_inplace(sum_revenue_n0, relin_keys);
    evaluator.rescale_to_next_inplace(sum_revenue_n0);
    evaluator.multiply(result_n1, revenue_cipher, sum_revenue_n1);
    evaluator.relinearize_inplace(sum_revenue_n1, relin_keys);
    evaluator.rescale_to_next_inplace(sum_revenue_n1);
    evaluator.multiply(result_n2, revenue_cipher, sum_revenue_n2);
    evaluator.relinearize_inplace(sum_revenue_n2, relin_keys);
    evaluator.rescale_to_next_inplace(sum_revenue_n2);
    evaluator.multiply(result_n3, revenue_cipher, sum_revenue_n3);
    evaluator.relinearize_inplace(sum_revenue_n3, relin_keys);
    evaluator.rescale_to_next_inplace(sum_revenue_n3);

    int logrow = log2(num);
    seal::Ciphertext temp;
    for (int i = 0; i < logrow; i++) {
        size_t step = 1 << (logrow - i - 1);
        temp = sum_revenue_n0;
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_revenue_n0, temp);
        temp = sum_revenue_n1;
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_revenue_n1, temp);
        temp = sum_revenue_n2;
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_revenue_n2, temp);
        temp = sum_revenue_n3;
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(sum_revenue_n3, temp);
    }
    end = std::chrono::system_clock::now();
    double ta =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count() /
        1000.0;
    printf("Agg time = %f ms\n", ta);
    aggregation_time +=
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::vector<double> agg_n0(slots_count), agg_n1(slots_count),
        agg_n2(slots_count), agg_n3(slots_count);
    decryptor.decrypt(sum_revenue_n0, plain);
    ckks_encoder.decode(plain, agg_n0);
    decryptor.decrypt(sum_revenue_n1, plain);
    ckks_encoder.decode(plain, agg_n1);
    decryptor.decrypt(sum_revenue_n2, plain);
    ckks_encoder.decode(plain, agg_n2);
    decryptor.decrypt(sum_revenue_n3, plain);
    ckks_encoder.decode(plain, agg_n3);

    std::cout << "--------------------------------------------------------"
              << std::endl;
    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time
              << " ms" << std::endl;
    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(12) << "n_name" << "|" << std::setw(12) << "revenue"
              << std::endl;
    std::cout << std::setw(12) << 0 << "|" << std::setw(12)
              << std::round(agg_n0[0]) << std::endl;
    std::cout << std::setw(12) << 1 << "|" << std::setw(12)
              << std::round(agg_n1[0]) << std::endl;
    std::cout << std::setw(12) << 2 << "|" << std::setw(12)
              << std::round(agg_n2[0]) << std::endl;
    std::cout << std::setw(12) << 3 << "|" << std::setw(12)
              << std::round(agg_n3[0]) << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) << "n_name" << "|" << std::setw(12) << "revenue"
              << std::endl;
    std::cout << std::setw(12) << 0 << "|" << std::setw(12) << plain_agg_res[0]
              << std::endl;
    std::cout << std::setw(12) << 1 << "|" << std::setw(12) << plain_agg_res[1]
              << std::endl;
    std::cout << std::setw(12) << 2 << "|" << std::setw(12) << plain_agg_res[2]
              << std::endl;
    std::cout << std::setw(12) << 3 << "|" << std::setw(12) << plain_agg_res[3]
              << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    return (filtering_time + aggregation_time) / 1000;
}

int main(int argc, char **argv)
{
    size_t num = 16;
    if (argc > 1) num = static_cast<size_t>(std::stoull(argv[1]));
    std::cout << "-------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "TPC-H Q5 evaluated via 3-PBS pruned compare + repack-based "
              << "CKKS aggregation" << std::endl;
    std::cout << std::endl;
    relational_query5(num);
    return 0;
}
