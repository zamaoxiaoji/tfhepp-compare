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
 * TPC-H Query 14
 * select
 *      100.00 * sum(case
 *          when p_type like 'PROMO%'
 *              then l_extendedprice * (1 - l_discount)
 *          else 0
 *      end) / sum(l_extendedprice * (1 - l_discount)) as promo_revenue
 *  from
 *      lineitem, part
 *  where
 *      l_partkey = p_partkey
 *      and l_shipdate >= date ':1'
 *      and l_shipdate < date ':1' + interval '1' month;
 *
 *  Filtering is rebuilt with the 3-PBS pruned comparator
 *  (Algorithm 2, Chapter 3 of Tang Li's thesis); CKKS aggregation reuses
 *  this repository's repack pipeline (PackLWEsToCKKS + HomomorphicRound).
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

double relational_query14(size_t num)
{
    std::cout << "Relational SQL Query14 Test (3-PBS pruned + CKKS repack): "
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
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> revenue_message(0, 100);
    std::uniform_int_distribution<uint32_t> ptype_message(0, 100);
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
    const auto micro_pack = GenerateFastB2AEvalKeyPack(sk, true);

    // Filtering
    std::vector<uint64_t> ship_date(num);
    std::vector<uint64_t> ptype(num);
    std::vector<TLWELvl2> shipdate_ciphers(num);
    std::vector<TLWELvl2> ptype_ciphers(num);

    uint32_t num_bits = 16;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - num_bits - 1;

    TLWELvl2 predicate1_cipher, predicate2_cipher;
    TLWELvl2 predicate3_cipher, predicate4_cipher;
    uint64_t predicate1_value = 10592, predicate2_value = 10957;
    uint64_t predicate3_value = 30, predicate4_value = 70;
    predicate1_cipher = tlweSymInt32Encrypt<Lvl2>(predicate1_value, Lvl2::α,
                                                  pow(2., scale_bits),
                                                  sk.key.get<Lvl2>());
    predicate2_cipher = tlweSymInt32Encrypt<Lvl2>(predicate2_value, Lvl2::α,
                                                  pow(2., scale_bits),
                                                  sk.key.get<Lvl2>());
    predicate3_cipher = tlweSymInt32Encrypt<Lvl2>(predicate3_value, Lvl2::α,
                                                  pow(2., scale_bits),
                                                  sk.key.get<Lvl2>());
    predicate4_cipher = tlweSymInt32Encrypt<Lvl2>(predicate4_value, Lvl2::α,
                                                  pow(2., scale_bits),
                                                  sk.key.get<Lvl2>());

    // Start sql evaluation
    std::vector<TLWELvl1> filter_res(num), filter_case_res(num);
    std::vector<TLWELvl2> aggregation_res(num);
    TLWELvl2 count_res;

    std::vector<double> revenue(num);

    for (size_t i = 0; i < num; i++) {
        revenue[i] = revenue_message(engine);
    }

    for (size_t i = 0; i < num; i++) {
        // Generate data
        ship_date[i] = shipdate_message(engine);
        ptype[i] = ptype_message(engine);
        shipdate_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(ship_date[i], Lvl2::α,
                                                        pow(2., scale_bits),
                                                        sk.key.get<Lvl2>());
        ptype_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(ptype[i], Lvl2::α,
                                                     pow(2., scale_bits),
                                                     sk.key.get<Lvl2>());
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

    for (size_t i = 0; i < num; i++) {
        TLWELvl1 pre_res;
        greater_than<Lvl2>(shipdate_ciphers[i], predicate1_cipher,
                           filter_res[i], num_bits, ek, micro_pack, LOGIC);
        less_than<Lvl2>(shipdate_ciphers[i], predicate2_cipher, pre_res,
                        num_bits, ek, micro_pack, LOGIC);
        TFHEpp::HomAND(filter_res[i], pre_res, filter_res[i], ek);
        greater_than<Lvl2>(ptype_ciphers[i], predicate3_cipher, pre_res,
                           num_bits, ek, micro_pack, LOGIC);
        TFHEpp::HomAND(filter_case_res[i], pre_res, filter_res[i], ek);
        less_than<Lvl2>(ptype_ciphers[i], predicate4_cipher, pre_res, num_bits,
                        ek, micro_pack, LOGIC);
        lift_and_and(filter_case_res[i], pre_res, filter_case_res[i], 29, ek);
        lift_and_and(filter_res[i], filter_res[i], filter_res[i], 29, ek);
    }
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start)
                         .count();

    std::vector<uint64_t> plain_filter_res(num), plain_filter_case_res(num);
    uint64_t plain_agg_res = 0, plain_agg_case_res = 0;
    for (size_t i = 0; i < num; i++) {
        if (ship_date[i] > predicate1_value &&
            ship_date[i] < predicate2_value) {
            plain_filter_res[i] = 1;
            plain_agg_res += revenue[i];
            if (ptype[i] > predicate3_value && ptype[i] < predicate4_value) {
                plain_filter_case_res[i] = 1;
                plain_agg_case_res += revenue[i];
            }
            else {
                plain_filter_case_res[i] = 0;
            }
        }
        else {
            plain_filter_res[i] = 0;
            plain_filter_case_res[i] = 0;
        }
    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = filter_res.size();
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
    seal::Ciphertext result, result_case;
    start = std::chrono::system_clock::now();
    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(result, filter_res, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(result, result.scale(), ckks_encoder,
                                   relin_keys, evaluator, context);

    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(result_case, filter_case_res, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(result_case, result_case.scale(),
                                   ckks_encoder, relin_keys, evaluator,
                                   context);
    end = std::chrono::system_clock::now();
    aggregation_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                           end - start)
                           .count();
    seal::Plaintext plain;
    std::vector<double> computed(slots_count), computed_case(slots_count);
    decryptor.decrypt(result, plain);
    ckks_encoder.decode(plain, computed);

    decryptor.decrypt(result_case, plain);
    ckks_encoder.decode(plain, computed_case);

    double err1 = 0., err2 = 0.;

    for (size_t i = 0; i < slots_count; ++i) {
        err1 += std::abs(computed[i] - plain_filter_res[i]);
        err2 += std::abs(computed_case[i] - plain_filter_case_res[i]);
    }

    printf("Repack average error = %f ~ 2^%.1f\n", err1 / slots_count,
           std::log2(err1 / slots_count));
    printf("Repack case average error = %f ~ 2^%.1f\n", err2 / slots_count,
           std::log2(err2 / slots_count));

    // Filter result * data
    seal::Ciphertext revenue_cipher;
    double qd = parms.coeff_modulus()[result.coeff_modulus_size() - 1].value();
    std::vector<double> revenue_slots(ckks_encoder.slot_count(), 0.);
    for (size_t i = 0; i < num; i++) revenue_slots[i] = revenue[i];
    ckks_encoder.encode(revenue_slots, result.parms_id(), qd, plain);
    symmetric_encryptor.encrypt_symmetric(plain, revenue_cipher);

    std::cout << "Aggregating price and discount .." << std::endl;
    start = std::chrono::system_clock::now();
    seal::Ciphertext aggregated, aggregated_case;
    evaluator.multiply(result, revenue_cipher, aggregated);
    evaluator.relinearize_inplace(aggregated, relin_keys);
    evaluator.rescale_to_next_inplace(aggregated);

    evaluator.multiply(result_case, revenue_cipher, aggregated_case);
    evaluator.relinearize_inplace(aggregated_case, relin_keys);
    evaluator.rescale_to_next_inplace(aggregated_case);
    std::cout << "Remian modulus: " << aggregated.coeff_modulus_size()
              << std::endl;
    int logrow = log2(num);

    seal::Ciphertext temp;
    for (int i = 0; i < logrow; i++) {
        temp = aggregated;
        size_t step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(aggregated, temp);

        temp = aggregated_case;
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(aggregated_case, temp);
    }
    end = std::chrono::system_clock::now();
    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(
                            end - start)
                            .count();
    std::vector<double> agg_result(slots_count), agg_case_result(slots_count);
    decryptor.decrypt(aggregated, plain);
    ckks_encoder.decode(plain, agg_result);

    decryptor.decrypt(aggregated_case, plain);
    ckks_encoder.decode(plain, agg_case_result);

    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time
              << " ms" << std::endl;

    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(16) << "promo_revenue" << std::endl;
    std::cout << std::setw(16)
              << (std::abs(agg_result[0]) < 1e-9
                      ? 0.
                      : 100.0 * agg_case_result[0] / agg_result[0])
              << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(16) << "promo_revenue" << std::endl;
    std::cout << std::setw(16)
              << (plain_agg_res == 0
                      ? 0.
                      : 100.0 * (plain_agg_case_res + 0.) / plain_agg_res)
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
    std::cout << "TPC-H Q14 evaluated via 3-PBS pruned compare + repack-based "
              << "CKKS aggregation" << std::endl;
    std::cout << std::endl;
    relational_query14(num);
    return 0;
}
