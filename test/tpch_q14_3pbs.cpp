#include "../comparison/comparison.h"
#include "ckks_relational.h"
#include "ckks_repack.h"
#include "gate.hpp"
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>

using namespace tfhepp_compare;
using namespace tfhepp_compare::three_pbs;
using namespace seal;

/***
 * TPC-H Query 14 -- compliant pipeline
 * select
 *     100.00 * sum(case when p_type like 'PROMO%'
 *                       then l_extendedprice * (1 - l_discount)
 *                       else 0
 *                  end) / sum(l_extendedprice * (1 - l_discount)) as promo_revenue
 *  from lineitem, part
 *  where l_partkey = p_partkey
 *    and l_shipdate >= date ':1'
 *    and l_shipdate <  date ':1' + interval '1' month;
 *
 *  Strict Chapter 4 contract:
 *    - WHERE predicate (l_shipdate range) is the only TFHE 3-PBS comparison.
 *    - The client encrypts raw lineitem/part keys and the part promo payload.
 *    - The server derives CKKS Lagrange masks, performs the encrypted lookup
 *      join, then aggregates the filtered revenue.
 *    - Plain data below is used only by the standalone driver to build the
 *      reference answer and to simulate client-side encryption.
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

// -------------------- CKKS level/scale helpers --------------------
namespace {

double LastCoeffModulus(const seal::Ciphertext &cipher,
                        const seal::SEALContext &context);

void ApplyActiveSlotMaskInPlace(seal::Ciphertext &cipher,
                                std::size_t active_slots,
                                seal::CKKSEncoder &encoder,
                                const seal::SEALContext &context,
                                seal::Evaluator &evaluator)
{
    std::vector<double> slots(encoder.slot_count(), 0.0);
    std::fill(slots.begin(), slots.begin() + active_slots, 1.0);
    seal::Plaintext plain;
    encoder.encode(slots, cipher.parms_id(),
                   std::min(std::ldexp(1.0, 45),
                            LastCoeffModulus(cipher, context) / 2.0),
                   plain);
    evaluator.multiply_plain_inplace(cipher, plain);
    evaluator.rescale_to_next_inplace(cipher);
}

void ModSwitchToCommonLevel(seal::Ciphertext &lhs, seal::Ciphertext &rhs,
                            seal::Evaluator &evaluator)
{
    if (lhs.parms_id() == rhs.parms_id()) return;
    if (lhs.coeff_modulus_size() > rhs.coeff_modulus_size())
        evaluator.mod_switch_to_inplace(lhs, rhs.parms_id());
    else if (rhs.coeff_modulus_size() > lhs.coeff_modulus_size())
        evaluator.mod_switch_to_inplace(rhs, lhs.parms_id());
}

seal::Ciphertext MultiplyAndRescale(const seal::Ciphertext &lhs_in,
                                    const seal::Ciphertext &rhs_in,
                                    const seal::RelinKeys &relin_keys,
                                    seal::Evaluator &evaluator)
{
    seal::Ciphertext lhs = lhs_in, rhs = rhs_in;
    ModSwitchToCommonLevel(lhs, rhs, evaluator);
    seal::Ciphertext result;
    evaluator.multiply(lhs, rhs, result);
    evaluator.relinearize_inplace(result, relin_keys);
    evaluator.rescale_to_next_inplace(result);
    return result;
}

double LastCoeffModulus(const seal::Ciphertext &cipher,
                        const seal::SEALContext &context)
{
    const auto context_data = context.get_context_data(cipher.parms_id());
    const auto &moduli = context_data->parms().coeff_modulus();
    return static_cast<double>(moduli[cipher.coeff_modulus_size() - 1].value());
}

seal::Ciphertext EncryptAtLevel(const std::vector<double> &slots,
                                seal::parms_id_type parms_id, double scale,
                                seal::CKKSEncoder &encoder,
                                seal::Encryptor &encryptor)
{
    seal::Plaintext plain;
    encoder.encode(slots, parms_id, scale, plain);
    seal::Ciphertext cipher;
    encryptor.encrypt(plain, cipher);
    return cipher;
}

std::vector<double> SlotsFrom(const std::vector<double> &values,
                              std::size_t slot_count)
{
    std::vector<double> s(slot_count, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) s[i] = values[i];
    return s;
}

template <typename T>
std::vector<double> SlotsFromIntegral(const std::vector<T> &values,
                                      std::size_t slot_count)
{
    std::vector<double> s(slot_count, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
        s[i] = static_cast<double>(values[i]);
    return s;
}

} // namespace

double relational_query14(size_t num)
{
    std::cout << "Relational SQL Query14 Test (WHERE->TFHE 3-PBS pruned, "
              << "JOIN+aggregation->CKKS Lagrange):\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "Records: " << num << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    TFHESecretKey sk;
    TFHEEvalKey ek;

    // -------- data domains --------
    // Lineitem: shipdate (WHERE), partkey (JOIN key)
    // Part:     partkey  (JOIN key), promo (payload, 0/1)
    constexpr std::size_t kPartDomain = 4; // part rows, partkey ∈ {0,1,2,3}
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> partkey_message(0, kPartDomain - 1);
    std::uniform_int_distribution<uint32_t> revenue_message(0, 100);
    std::uniform_int_distribution<uint32_t> promo_message(0, 1);

    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
    const auto micro_pack = GenerateFastB2AEvalKeyPack(sk, true);

    // -------- generate plaintext source tables (the "client side") --------
    std::vector<uint64_t> ship_date(num), line_partkey(num);
    std::vector<double> revenue(num);
    for (size_t i = 0; i < num; i++) {
        ship_date[i] = shipdate_message(engine);
        line_partkey[i] = partkey_message(engine);
        revenue[i] = revenue_message(engine);
    }
    std::vector<uint32_t> part_keys(kPartDomain), part_promo(kPartDomain);
    for (size_t j = 0; j < kPartDomain; j++) {
        part_keys[j] = static_cast<uint32_t>(j);
        part_promo[j] = promo_message(engine);
    }
    uint64_t pred_lo = 10592, pred_hi = 10957;
    if (num > 0) {
        ship_date[0] = pred_lo;
        line_partkey[0] = 0;
        revenue[0] = std::max<double>(revenue[0], 1.0);
        part_promo[0] = 1;
    }

    // -------- TFHE encryption of WHERE-relevant column only --------
    uint32_t num_bits = 16;
    uint32_t compprecision = 32;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - num_bits - 1;
    std::vector<TLWELvl2> shipdate_ciphers(num);
    for (size_t i = 0; i < num; i++) {
        shipdate_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(
            ship_date[i], Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    }
    TLWELvl2 pred_lo_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_lo - 1, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    TLWELvl2 pred_hi_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_hi, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());

    // -------- WHERE in TFHE 3-PBS --------
    std::vector<TLWELvl1> filter_res(num);
    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, ckks_time = 0;
    start = std::chrono::system_clock::now();
    for (size_t i = 0; i < num; i++) {
        TLWELvl1 pre_res;
        greater_than<Lvl2>(shipdate_ciphers[i], pred_lo_ct, filter_res[i],
                           compprecision, ek, micro_pack, LOGIC);
        less_than<Lvl2>(shipdate_ciphers[i], pred_hi_ct, pre_res,
                        compprecision, ek, micro_pack, LOGIC);
        lift_and_and(filter_res[i], pre_res, filter_res[i], 29, ek);
    }
    end = std::chrono::system_clock::now();
    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start).count();

    // -------- plain reference --------
    std::vector<uint64_t> plain_filter(num, 0);
    double plain_total = 0.0, plain_promo = 0.0;
    for (size_t i = 0; i < num; i++) {
        const bool date_ok =
            ship_date[i] >= pred_lo && ship_date[i] < pred_hi;
        plain_filter[i] = date_ok ? 1 : 0;
        if (date_ok) {
            plain_total += revenue[i];
            plain_promo += revenue[i] * static_cast<double>(
                                            part_promo[line_partkey[i]]);
        }
    }
    std::cout << "Filtering finish" << std::endl;

    // -------- CKKS setup --------
    std::cout << "Aggregation :" << std::endl;
    uint32_t mask_scale_bits = 29;
    uint64_t modulus_bits = 45;
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    // Coeff chain is intentionally generous: repack consumes several levels
    // and the encrypted Lagrange lookup consumes more before aggregation.
    parms.set_coeff_modulus(seal::CoeffModulus::Create(
        poly_modulus_degree,
        {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45,
         45, 45, 45, 45, 45, 45, 45, 45, 59}));
    seal::SEALContext context(parms, true, seal::sec_level_type::none);
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
    seal::Encryptor encryptor(context, seal_public_key);
    seal::Encryptor symmetric_encryptor(context, seal_secret_key);
    seal::Evaluator evaluator(context);
    seal::Decryptor decryptor(context, seal_secret_key);
    seal::CKKSEncoder ckks_encoder(context);
    auto repack_config =
        tfhepp_ckks::DefaultRepackConfig<Lvl1>(mask_scale_bits, modulus_bits);
    tfhepp_ckks::RepackEvaluationKey repack_key;
    std::cout << "Generating Conversion Key..." << std::endl;
    tfhepp_ckks::GenerateRepackKey<Lvl1>(repack_key, sk, repack_config.key_scale,
                                          ckks_encoder, symmetric_encryptor,
                                          context);

    // -------- 1. Repack TFHE WHERE mask to CKKS --------
    std::cout << "Starting Conversion..." << std::endl;
    start = std::chrono::system_clock::now();
    seal::Ciphertext mask_ckks;
    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(mask_ckks, filter_res, repack_key,
                                       repack_config, ckks_encoder, galois_keys,
                                       relin_keys, evaluator, context);
    tfhepp_ckks::HomomorphicRound(mask_ckks, mask_ckks.scale(), ckks_encoder,
                                   relin_keys, evaluator, context);
    {
        seal::Plaintext pt;
        std::vector<double> decoded;
        decryptor.decrypt(mask_ckks, pt);
        ckks_encoder.decode(pt, decoded);
        double e = 0.0;
        for (size_t i = 0; i < num; ++i)
            e += std::abs(decoded[i] - static_cast<double>(plain_filter[i]));
        printf("Repack mask average error = %g ~ 2^%.1f  "
               "(active slots = %zu)\n",
               e / num, std::log2(e / std::max<double>(num, 1)), num);
    }

    // -------- 2. Encrypt JOIN keys and payload as CKKS (client-side) --------
    // (We encrypt here in the test driver since the test plays both roles.)
    std::vector<double> line_partkey_slots =
        SlotsFromIntegral(line_partkey, ckks_encoder.slot_count());
    std::vector<double> part_partkey_slots =
        SlotsFromIntegral(part_keys, ckks_encoder.slot_count());
    std::vector<double> part_promo_slots =
        SlotsFromIntegral(part_promo, ckks_encoder.slot_count());

    const double scale = std::pow(2.0, 40);
    auto line_partkey_ct =
        EncryptAtLevel(line_partkey_slots, context.first_parms_id(), scale,
                       ckks_encoder, encryptor);
    auto part_partkey_ct =
        EncryptAtLevel(part_partkey_slots, context.first_parms_id(), scale,
                       ckks_encoder, encryptor);
    auto part_promo_ct =
        EncryptAtLevel(part_promo_slots, context.first_parms_id(), scale,
                       ckks_encoder, encryptor);

    // -------- 3. CKKS payload JOIN: broadcast part_promo to each lineitem --------
    std::vector<double> partkey_domain(kPartDomain);
    for (size_t j = 0; j < kPartDomain; ++j)
        partkey_domain[j] = static_cast<double>(j);
    auto line_partkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        line_partkey_ct, partkey_domain, relin_keys, ckks_encoder, evaluator);
    auto part_partkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        part_partkey_ct, partkey_domain, relin_keys, ckks_encoder, evaluator);
    auto promo_on_line = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        line_partkey_masks, part_partkey_masks, part_promo_ct, kPartDomain,
        num, relin_keys, galois_keys, ckks_encoder, context, evaluator);
    ApplyActiveSlotMaskInPlace(promo_on_line, num, ckks_encoder, context,
                                evaluator);

    // -------- 4. Encrypt revenue at matching level for multiplication --------
    double qd_mask = LastCoeffModulus(mask_ckks, context);
    auto revenue_ct = EncryptAtLevel(
        SlotsFrom(revenue, ckks_encoder.slot_count()), mask_ckks.parms_id(),
        qd_mask, ckks_encoder, encryptor);

    auto filtered_revenue =
        MultiplyAndRescale(mask_ckks, revenue_ct, relin_keys, evaluator);
    // Restrict to the active region (public num).
    ApplyActiveSlotMaskInPlace(filtered_revenue, num, ckks_encoder, context,
                                evaluator);
    auto promo_revenue =
        MultiplyAndRescale(filtered_revenue, promo_on_line, relin_keys,
                           evaluator);
    ApplyActiveSlotMaskInPlace(promo_revenue, num, ckks_encoder, context,
                                evaluator);

    int logrow = static_cast<int>(std::ceil(std::log2(num)));
    seal::Ciphertext tot_sum = filtered_revenue;
    seal::Ciphertext prm_sum = promo_revenue;
    for (int i = 0; i < logrow; ++i) {
        seal::Ciphertext temp = tot_sum;
        size_t step = 1ULL << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, static_cast<int>(step),
                                         galois_keys);
        evaluator.add_inplace(tot_sum, temp);
        temp = prm_sum;
        evaluator.rotate_vector_inplace(temp, static_cast<int>(step),
                                         galois_keys);
        evaluator.add_inplace(prm_sum, temp);
    }
    end = std::chrono::system_clock::now();
    ckks_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                    end - start).count();

    // -------- Verification --------
    std::vector<double> total_decoded, promo_decoded;
    seal::Plaintext plain_out;
    decryptor.decrypt(tot_sum, plain_out);
    ckks_encoder.decode(plain_out, total_decoded);
    decryptor.decrypt(prm_sum, plain_out);
    ckks_encoder.decode(plain_out, promo_decoded);
    const double encrypted_total = total_decoded[0];
    const double encrypted_promo = promo_decoded[0];
    const double encrypted_ratio =
        std::abs(encrypted_total) < 0.5 ? 0.0
                                          : 100.0 * encrypted_promo /
                                                encrypted_total;
    const double plain_ratio =
        plain_total < 1e-9 ? 0.0 : 100.0 * plain_promo / plain_total;

    std::cout << "Filtering Time: " << filtering_time << " ms" << std::endl;
    std::cout << "CKKS (Repack + Join + Aggregation) Time: " << ckks_time
              << " ms" << std::endl;
    std::cout << "Query Evaluation Time: " << filtering_time + ckks_time
              << " ms" << std::endl;
    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(16) << "promo_revenue" << std::endl;
    std::cout << std::setw(16) << encrypted_ratio << std::endl;
    std::cout << std::setw(16) << "(total=" << encrypted_total
              << ", promo=" << encrypted_promo << ")" << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(16) << "promo_revenue" << std::endl;
    std::cout << std::setw(16) << plain_ratio << std::endl;
    std::cout << std::setw(16) << "(total=" << plain_total
              << ", promo=" << plain_promo << ")" << std::endl;
    std::cout << std::endl;
    return (filtering_time + ckks_time) / 1000.0;
}

int main(int argc, char **argv)
{
    size_t num = 16;
    if (argc > 1) num = static_cast<size_t>(std::stoull(argv[1]));
    std::cout << "----------------------------------------------------"
              << std::endl;
    std::cout << "TPC-H Q14 (strict Chapter 4): WHERE = TFHE 3-PBS, "
              << "JOIN + aggregation = CKKS Lagrange" << std::endl;
    std::cout << std::endl;
    relational_query14(num);
    return 0;
}
