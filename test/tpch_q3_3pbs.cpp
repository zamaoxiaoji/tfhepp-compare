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
 * TPC-H Query 3 -- compliant pipeline
 * select
 *     l_orderkey,
 *     sum(l_extendedprice * (1 - l_discount)) as revenue,
 *     o_orderdate,
 *     o_shippriority
 *  from customer, orders, lineitem
 *  where c_mktsegment = ':1'
 *    and c_custkey   = o_custkey
 *    and l_orderkey  = o_orderkey
 *    and o_orderdate < date ':2'
 *    and l_shipdate  > date ':2'
 *  group by l_orderkey, o_orderdate, o_shippriority;
 *
 *  Compliance contract (per Chapter 4):
 *    - The three WHERE predicates (c_mktsegment == 'BUILDING',
 *      o_orderdate < D, l_shipdate > D) are the only TFHE 3-PBS comparisons.
 *    - All foreign keys (l_orderkey, o_orderkey, o_custkey, c_custkey) and
 *      every group-by payload (o_shippriority) are sent as CKKS ciphertexts.
 *    - customer ⋈ orders is done with CKKS Lagrange masks built from the
 *      encrypted custkey columns and LookupJoinFromEncryptedMasks
 *      (Algorithm 4.3).
 *    - orders ⋈ lineitem is done the same way on orderkey.
 *    - GROUP BY orderkey is realised with CKKS Lagrange indicator masks on
 *      the same encrypted lineitem.orderkey column.
 *    - The server sees no plaintext join key, no plaintext selector, no
 *      plaintext group id; it only sees public domain sizes and the test
 *      driver's plaintext tables (which simulate the client).
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
    TLWELvl1 ge, le;
    greater_than_equal<P>(cipher1, cipher2, ge, plain_bits, ek, micro_pack,
                          LOGIC);
    less_than_equal<P>(cipher1, cipher2, le, plain_bits, ek, micro_pack,
                       LOGIC);
    HomAND(res, ge, le, ek, result_type);
}

// -------------------- CKKS level/scale helpers --------------------
namespace {

void ApplyActiveSlotMaskInPlace(seal::Ciphertext &cipher,
                                std::size_t active_slots,
                                seal::CKKSEncoder &encoder,
                                const seal::SEALContext &context,
                                seal::Evaluator &evaluator)
{
    std::vector<double> slots(encoder.slot_count(), 0.0);
    std::fill(slots.begin(), slots.begin() + active_slots, 1.0);
    seal::Plaintext plain;
    const auto context_data = context.get_context_data(cipher.parms_id());
    const auto &moduli = context_data->parms().coeff_modulus();
    const double qd =
        static_cast<double>(moduli[cipher.coeff_modulus_size() - 1].value());
    const double target_scale = cipher.scale();
    encoder.encode(slots, cipher.parms_id(), qd, plain);
    evaluator.multiply_plain_inplace(cipher, plain);
    evaluator.rescale_to_next_inplace(cipher);
    cipher.scale() = target_scale;
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

template <typename T>
std::vector<double> SlotsFromIntegral(const std::vector<T> &values,
                                      std::size_t slot_count)
{
    std::vector<double> s(slot_count, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
        s[i] = static_cast<double>(values[i]);
    return s;
}

std::vector<double> SlotsFromDouble(const std::vector<double> &values,
                                    std::size_t slot_count)
{
    std::vector<double> s(slot_count, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) s[i] = values[i];
    return s;
}

seal::Ciphertext PackTfheMask(std::vector<TLWELvl1> masks, std::size_t active,
                              tfhepp_ckks::RepackEvaluationKey &repack_key,
                              tfhepp_ckks::RepackConfig &repack_config,
                              seal::CKKSEncoder &encoder,
                              seal::GaloisKeys &galois_keys,
                              seal::RelinKeys &relin_keys,
                              seal::Evaluator &evaluator,
                              seal::SEALContext &context)
{
    // Pad to a power-of-two of at least 2*active to satisfy PackLWEsToCKKS.
    std::size_t target = 1;
    while (target < std::max<std::size_t>(2 * active, 4)) target <<= 1;
    masks.resize(target);
    seal::Ciphertext packed;
    tfhepp_ckks::PackLWEsToCKKS<Lvl1>(packed, masks, repack_key, repack_config,
                                       encoder, galois_keys, relin_keys,
                                       evaluator, context);
    tfhepp_ckks::HomomorphicRound(packed, packed.scale(), encoder, relin_keys,
                                   evaluator, context);
    return packed;
}

} // namespace

double relational_query3(size_t num)
{
    std::cout << "Relational SQL Query3 Test (WHERE->TFHE 3-PBS pruned, "
              << "JOIN+GROUP BY->CKKS Lagrange):\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "Records: " << num << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    TFHESecretKey sk;
    TFHEEvalKey ek;
    constexpr std::size_t kKeyDomain = 4;
    constexpr std::size_t kSegmentDomain = 8; // mktsegment 0..7
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> orderdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> seg_message(0, 4);
    std::uniform_int_distribution<uint32_t> key_message(0, kKeyDomain - 1);
    std::uniform_int_distribution<uint32_t> shippri_message(0, 1);
    std::uniform_int_distribution<uint32_t> revenue_message(0, 100);

    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
    const auto micro_pack = GenerateFastB2AEvalKeyPack(sk, true);

    // -------- generate plaintext source tables --------
    // Lineitem (num rows): orderkey, shipdate, revenue
    std::vector<uint64_t> ship_date(num), line_orderkey(num);
    std::vector<double> revenue(num);
    for (size_t i = 0; i < num; i++) {
        ship_date[i] = shipdate_message(engine);
        line_orderkey[i] = key_message(engine);
        revenue[i] = revenue_message(engine);
    }
    // Orders (kKeyDomain rows): orderkey, custkey, orderdate, shippriority
    std::vector<uint64_t> order_orderkey(kKeyDomain), order_custkey(kKeyDomain),
        order_date(kKeyDomain), order_shippriority(kKeyDomain);
    for (size_t j = 0; j < kKeyDomain; j++) {
        order_orderkey[j] = j;
        order_custkey[j] = key_message(engine);
        order_date[j] = orderdate_message(engine);
        order_shippriority[j] = shippri_message(engine);
    }
    // Customer (kKeyDomain rows): custkey, mktsegment
    std::vector<uint64_t> cust_custkey(kKeyDomain), cust_mktsegment(kKeyDomain);
    for (size_t k = 0; k < kKeyDomain; k++) {
        cust_custkey[k] = k;
        cust_mktsegment[k] = seg_message(engine);
    }

    // -------- TFHE encryption of WHERE-relevant columns only --------
    uint32_t num_bits = 16;
    uint32_t compprecision = 32;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - num_bits - 1;
    std::vector<TLWELvl2> shipdate_ciphers(num), orderdate_ciphers(kKeyDomain),
        mktsegment_ciphers(kKeyDomain);
    for (size_t i = 0; i < num; i++)
        shipdate_ciphers[i] = tlweSymInt32Encrypt<Lvl2>(
            ship_date[i], Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    for (size_t j = 0; j < kKeyDomain; j++) {
        orderdate_ciphers[j] = tlweSymInt32Encrypt<Lvl2>(
            order_date[j], Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
        mktsegment_ciphers[j] = tlweSymInt32Encrypt<Lvl2>(
            cust_mktsegment[j], Lvl2::α, pow(2., scale_bits),
            sk.key.get<Lvl2>());
    }
    uint64_t pred_date = 15000;
    uint64_t pred_segment = 1; // BUILDING
    TLWELvl2 pred_date_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_date, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    TLWELvl2 pred_segment_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_segment, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());

    // -------- WHERE in TFHE 3-PBS (3 columns, one mask per row in each table) --------
    std::vector<TLWELvl1> line_mask(num), order_mask(kKeyDomain),
        cust_mask(kKeyDomain);
    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, ckks_time = 0;
    start = std::chrono::system_clock::now();
    for (size_t i = 0; i < num; i++) {
        TLWELvl1 raw;
        greater_than<Lvl2>(shipdate_ciphers[i], pred_date_ct, raw,
                           compprecision, ek, micro_pack, LOGIC);
        lift_and_and(raw, raw, line_mask[i], 29, ek);
    }
    for (size_t j = 0; j < kKeyDomain; j++) {
        TLWELvl1 raw;
        less_than<Lvl2>(orderdate_ciphers[j], pred_date_ct, raw, compprecision,
                        ek, micro_pack, LOGIC);
        lift_and_and(raw, raw, order_mask[j], 29, ek);
    }
    for (size_t k = 0; k < kKeyDomain; k++) {
        TLWELvl1 raw;
        equal<Lvl2>(mktsegment_ciphers[k], pred_segment_ct, raw, compprecision,
                    ek, micro_pack, LOGIC);
        lift_and_and(raw, raw, cust_mask[k], 29, ek);
    }
    end = std::chrono::system_clock::now();
    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start).count();

    // -------- plain reference (per-table masks + final aggregate) --------
    std::vector<uint64_t> plain_line_mask(num, 0);
    std::vector<uint64_t> plain_order_mask(kKeyDomain, 0);
    std::vector<uint64_t> plain_cust_mask(kKeyDomain, 0);
    for (size_t i = 0; i < num; i++)
        plain_line_mask[i] = (ship_date[i] > pred_date) ? 1 : 0;
    for (size_t j = 0; j < kKeyDomain; j++)
        plain_order_mask[j] = (order_date[j] < pred_date) ? 1 : 0;
    for (size_t k = 0; k < kKeyDomain; k++)
        plain_cust_mask[k] = (cust_mktsegment[k] == pred_segment) ? 1 : 0;

    std::vector<double> plain_agg(kKeyDomain, 0.0);
    for (size_t i = 0; i < num; i++) {
        size_t oj = 0;
        for (; oj < kKeyDomain; oj++)
            if (order_orderkey[oj] == line_orderkey[i]) break;
        if (oj >= kKeyDomain) continue;
        size_t ck = 0;
        for (; ck < kKeyDomain; ck++)
            if (cust_custkey[ck] == order_custkey[oj]) break;
        if (ck >= kKeyDomain) continue;
        if (plain_line_mask[i] && plain_order_mask[oj] && plain_cust_mask[ck])
            plain_agg[line_orderkey[i]] += revenue[i];
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
    parms.set_coeff_modulus(seal::CoeffModulus::Create(
        poly_modulus_degree,
        {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45,
         45, 45, 45, 45, 45, 45, 45, 59}));
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

    // -------- 1. Repack the three TFHE WHERE masks to CKKS --------
    std::cout << "Starting Conversion..." << std::endl;
    start = std::chrono::system_clock::now();
    auto line_mask_ckks = PackTfheMask(line_mask, num, repack_key,
                                        repack_config, ckks_encoder,
                                        galois_keys, relin_keys, evaluator,
                                        context);
    auto order_mask_ckks = PackTfheMask(order_mask, kKeyDomain, repack_key,
                                         repack_config, ckks_encoder,
                                         galois_keys, relin_keys, evaluator,
                                         context);
    auto cust_mask_ckks = PackTfheMask(cust_mask, kKeyDomain, repack_key,
                                        repack_config, ckks_encoder,
                                        galois_keys, relin_keys, evaluator,
                                        context);
    auto report_repack = [&](const seal::Ciphertext &ct,
                              const std::vector<uint64_t> &plain,
                              const char *name) {
        seal::Plaintext pt;
        std::vector<double> dec;
        decryptor.decrypt(ct, pt);
        ckks_encoder.decode(pt, dec);
        double e = 0.0;
        for (size_t i = 0; i < plain.size(); ++i)
            e += std::abs(dec[i] - static_cast<double>(plain[i]));
        printf("Repack %s average error = %g ~ 2^%.1f  (active=%zu)\n", name,
               e / plain.size(),
               std::log2(e / std::max<double>(plain.size(), 1)),
               plain.size());
    };
    report_repack(line_mask_ckks, plain_line_mask, "line_mask");
    report_repack(order_mask_ckks, plain_order_mask, "order_mask");
    report_repack(cust_mask_ckks, plain_cust_mask, "cust_mask");

    // -------- 2. Encrypt JOIN keys and group-by payload as CKKS --------
    double scale = std::pow(2.0, 40);
    seal::Plaintext plain_tmp;
    ckks_encoder.encode(SlotsFromIntegral(order_custkey,
                                          ckks_encoder.slot_count()),
                        scale, plain_tmp);
    seal::Ciphertext order_custkey_ct;
    encryptor.encrypt(plain_tmp, order_custkey_ct);
    ckks_encoder.encode(SlotsFromIntegral(cust_custkey,
                                          ckks_encoder.slot_count()),
                        scale, plain_tmp);
    seal::Ciphertext cust_custkey_ct;
    encryptor.encrypt(plain_tmp, cust_custkey_ct);
    ckks_encoder.encode(SlotsFromIntegral(line_orderkey,
                                          ckks_encoder.slot_count()),
                        scale, plain_tmp);
    seal::Ciphertext line_orderkey_ct;
    encryptor.encrypt(plain_tmp, line_orderkey_ct);
    ckks_encoder.encode(SlotsFromIntegral(order_orderkey,
                                          ckks_encoder.slot_count()),
                        scale, plain_tmp);
    seal::Ciphertext order_orderkey_ct;
    encryptor.encrypt(plain_tmp, order_orderkey_ct);

    // -------- 3. JOIN customer→orders on custkey via CKKS Lagrange --------
    std::vector<double> key_domain(kKeyDomain);
    for (size_t j = 0; j < kKeyDomain; j++) key_domain[j] = double(j);

    auto order_custkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        order_custkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto cust_custkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        cust_custkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto cust_mask_on_orders = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        order_custkey_masks, cust_custkey_masks, cust_mask_ckks,
        kKeyDomain, relin_keys, galois_keys, evaluator);

    // AND with order_mask → order_filter (still indexed by orders rows)
    auto order_filter =
        MultiplyAndRescale(cust_mask_on_orders, order_mask_ckks, relin_keys,
                           evaluator);

    // -------- 4. JOIN orders→lineitem on orderkey via CKKS Lagrange --------
    auto line_orderkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        line_orderkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto order_orderkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        order_orderkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto order_filter_on_line = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        line_orderkey_masks, order_orderkey_masks, order_filter,
        kKeyDomain, relin_keys, galois_keys, evaluator);

    // AND with l_mask (lineitem-level shipdate filter)
    auto row_filter =
        MultiplyAndRescale(order_filter_on_line, line_mask_ckks, relin_keys,
                           evaluator);

    // -------- 5. Multiply by revenue and group by lineitem.orderkey --------
    double qd = LastCoeffModulus(row_filter, context);
    auto revenue_ct = EncryptAtLevel(
        SlotsFromDouble(revenue, ckks_encoder.slot_count()),
        row_filter.parms_id(), qd, ckks_encoder, encryptor);
    auto filtered_revenue =
        MultiplyAndRescale(row_filter, revenue_ct, relin_keys, evaluator);
    // Zero out padding slots so noise outside the active region does not leak
    // into the rotate-and-sum aggregate.
    ApplyActiveSlotMaskInPlace(filtered_revenue, num, ckks_encoder, context,
                                evaluator);

    // GROUP BY: reuse line_orderkey_masks but they may be at a higher level;
    // MultiplyAndRescale will mod-switch automatically.
    std::vector<seal::Ciphertext> grouped(kKeyDomain);
    int logrow = static_cast<int>(std::ceil(std::log2(num)));
    for (size_t g = 0; g < kKeyDomain; g++) {
        auto weighted = MultiplyAndRescale(filtered_revenue,
                                           line_orderkey_masks[g],
                                           relin_keys, evaluator);
        for (int b = 0; b < logrow; ++b) {
            seal::Ciphertext temp = weighted;
            size_t step = 1ULL << (logrow - b - 1);
            evaluator.rotate_vector_inplace(temp, static_cast<int>(step),
                                             galois_keys);
            evaluator.add_inplace(weighted, temp);
        }
        grouped[g] = std::move(weighted);
    }
    end = std::chrono::system_clock::now();
    ckks_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                    end - start).count();

    // -------- Verification --------
    std::cout << "Filtering Time: " << filtering_time << " ms" << std::endl;
    std::cout << "CKKS (Repack + Join + Group + Agg) Time: " << ckks_time
              << " ms" << std::endl;
    std::cout << "Query Evaluation Time: " << filtering_time + ckks_time
              << " ms" << std::endl;
    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(12) << "orderkey" << "|" << std::setw(12)
              << "revenue" << std::endl;
    for (size_t g = 0; g < kKeyDomain; g++) {
        std::vector<double> dec;
        seal::Plaintext plain_g;
        decryptor.decrypt(grouped[g], plain_g);
        ckks_encoder.decode(plain_g, dec);
        std::cout << std::setw(12) << g << "|" << std::setw(12)
                  << std::round(dec[0]) << std::endl;
    }
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) << "orderkey" << "|" << std::setw(12)
              << "revenue" << std::endl;
    for (size_t g = 0; g < kKeyDomain; g++) {
        std::cout << std::setw(12) << g << "|" << std::setw(12)
                  << plain_agg[g] << std::endl;
    }
    std::cout << std::endl;
    return (filtering_time + ckks_time) / 1000.0;
}

int main(int argc, char **argv)
{
    size_t num = 16;
    if (argc > 1) num = static_cast<size_t>(std::stoull(argv[1]));
    std::cout << "----------------------------------------------------"
              << std::endl;
    std::cout << "TPC-H Q3 (compliant): WHERE = TFHE 3-PBS, "
              << "JOIN+GROUP BY = CKKS Lagrange (Algorithm 4.3)" << std::endl;
    std::cout << std::endl;
    relational_query3(num);
    return 0;
}
