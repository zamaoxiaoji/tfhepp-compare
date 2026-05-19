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
 * TPC-H Query 5 -- compliant pipeline
 * select
 *     n_name,
 *     sum(l_extendedprice * (1 - l_discount)) as revenue
 *  from customer, orders, lineitem, supplier, nation, region
 *  where c_custkey   = o_custkey
 *    and l_orderkey  = o_orderkey
 *    and l_suppkey   = s_suppkey
 *    and c_nationkey = s_nationkey
 *    and s_nationkey = n_nationkey
 *    and n_regionkey = r_regionkey
 *    and r_name      = ':1'
 *    and o_orderdate >= date ':2'
 *    and o_orderdate <  date ':2' + interval '1' year
 *  group by n_name;
 *
 *  Strict Chapter 4 contract:
 *    - WHERE predicates (o_orderdate range, r_name == 'ASIA') are the only
 *      TFHE 3-PBS comparisons.
 *    - The client encrypts raw key/payload columns; the server derives
 *      Lagrange masks, lookup joins, equality masks, and group-by masks over
 *      CKKS.
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

// -------------------- CKKS helpers --------------------
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

void KeepFirstSlotInPlace(seal::Ciphertext &cipher,
                          seal::CKKSEncoder &encoder,
                          const seal::SEALContext &context,
                          seal::Evaluator &evaluator)
{
    std::vector<double> slots(encoder.slot_count(), 0.0);
    slots[0] = 1.0;
    seal::Plaintext plain;
    encoder.encode(slots, cipher.parms_id(),
                   std::min(std::ldexp(1.0, 45),
                            LastCoeffModulus(cipher, context) / 2.0),
                   plain);
    evaluator.multiply_plain_inplace(cipher, plain);
    evaluator.rescale_to_next_inplace(cipher);
}

void ReduceToFirstSlotInPlace(seal::Ciphertext &cipher,
                              std::size_t active_slots,
                              const seal::GaloisKeys &galois_keys,
                              seal::CKKSEncoder &encoder,
                              const seal::SEALContext &context,
                              seal::Evaluator &evaluator)
{
    tfhepp_ckks::RotateAndSumInPlace(cipher, active_slots, galois_keys,
                                     evaluator);
    KeepFirstSlotInPlace(cipher, encoder, context, evaluator);
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

std::vector<double> SlotsFromDouble(const std::vector<double> &values,
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

seal::Ciphertext PackTfheMask(std::vector<TLWELvl1> masks, std::size_t active,
                              tfhepp_ckks::RepackEvaluationKey &repack_key,
                              tfhepp_ckks::RepackConfig &repack_config,
                              seal::CKKSEncoder &encoder,
                              seal::GaloisKeys &galois_keys,
                              seal::RelinKeys &relin_keys,
                              seal::Evaluator &evaluator,
                              seal::SEALContext &context)
{
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

double relational_query5(size_t num)
{
    std::cout << "Relational SQL Query5 Test (WHERE->TFHE 3-PBS, "
              << "JOIN+GROUP BY->CKKS Lagrange):\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "Records: " << num << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    TFHESecretKey sk;
    TFHEEvalKey ek;
    constexpr std::size_t kKeyDomain    = 4; // # of orders/customers/suppliers
    constexpr std::size_t kNationDomain = 4;
    constexpr std::size_t kRegionDomain = 2; // ASIA = 1, OTHER = 0
    std::uniform_int_distribution<uint32_t> orderdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> key_message(0, kKeyDomain - 1);
    std::uniform_int_distribution<uint32_t> nat_message(0, kNationDomain - 1);
    std::uniform_int_distribution<uint32_t> reg_message(0, kRegionDomain - 1);
    std::uniform_int_distribution<uint32_t> revenue_message(0, 100);

    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
    const auto micro_pack = GenerateFastB2AEvalKeyPack(sk, true);

    // -------- generate plaintext source tables (the "client") --------
    // Lineitem: orderkey, suppkey, revenue
    std::vector<uint64_t> line_orderkey(num), line_suppkey(num);
    std::vector<double> revenue(num);
    for (size_t i = 0; i < num; i++) {
        line_orderkey[i] = key_message(engine);
        line_suppkey[i] = key_message(engine);
        revenue[i] = revenue_message(engine);
    }
    // Orders: orderkey, custkey, orderdate
    std::vector<uint64_t> order_orderkey(kKeyDomain), order_custkey(kKeyDomain),
        order_date(kKeyDomain);
    for (size_t j = 0; j < kKeyDomain; j++) {
        order_orderkey[j] = j;
        order_custkey[j] = key_message(engine);
        order_date[j] = orderdate_message(engine);
    }
    // Customer: custkey, nationkey
    std::vector<uint64_t> cust_custkey(kKeyDomain), cust_nationkey(kKeyDomain);
    for (size_t k = 0; k < kKeyDomain; k++) {
        cust_custkey[k] = k;
        cust_nationkey[k] = nat_message(engine);
    }
    // Supplier: suppkey, nationkey
    std::vector<uint64_t> supp_suppkey(kKeyDomain), supp_nationkey(kKeyDomain);
    for (size_t s = 0; s < kKeyDomain; s++) {
        supp_suppkey[s] = s;
        supp_nationkey[s] = nat_message(engine);
    }
    // Nation: nationkey, regionkey
    std::vector<uint64_t> nation_nationkey(kNationDomain),
        nation_regionkey(kNationDomain);
    for (size_t n = 0; n < kNationDomain; n++) {
        nation_nationkey[n] = n;
        nation_regionkey[n] = reg_message(engine);
    }
    // Region: regionkey, name (binary: 1 = ASIA, 0 = other)
    std::vector<uint64_t> region_regionkey(kRegionDomain),
        region_name(kRegionDomain);
    for (size_t r = 0; r < kRegionDomain; r++) {
        region_regionkey[r] = r;
        region_name[r] = (r == 1) ? 1 : 0; // pretend regionkey 1 is "ASIA"
    }
    uint64_t pred_lo = 12000, pred_hi = 13000;  // orderdate range
    uint64_t pred_asia = 1;                     // r_name == 'ASIA' flag
    if (num > 0) {
        line_orderkey[0] = 0;
        line_suppkey[0] = 0;
        revenue[0] = std::max<double>(revenue[0], 1.0);
        order_custkey[0] = 0;
        order_date[0] = pred_lo;
        cust_nationkey[0] = 0;
        supp_nationkey[0] = 0;
        nation_regionkey[0] = pred_asia;
        for (size_t j = 1; j < kKeyDomain; ++j)
            order_date[j] = pred_hi + 1000 + j;
    }

    // -------- TFHE encryption of WHERE-relevant columns only --------
    uint32_t num_bits = 16;
    uint32_t compprecision = 32;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - num_bits - 1;
    uint32_t region_bits = 2;
    uint32_t region_scale_bits =
        std::numeric_limits<Lvl1::T>::digits - region_bits - 1;
    std::vector<TLWELvl2> orderdate_ciphers(kKeyDomain);
    std::vector<TLWELvl1> region_name_ciphers(kRegionDomain);
    for (size_t j = 0; j < kKeyDomain; j++)
        orderdate_ciphers[j] = tlweSymInt32Encrypt<Lvl2>(
            order_date[j], Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    for (size_t r = 0; r < kRegionDomain; r++)
        region_name_ciphers[r] = tlweSymInt32Encrypt<Lvl1>(
            region_name[r], Lvl1::α, pow(2., region_scale_bits),
            sk.key.get<Lvl1>());

    TLWELvl2 pred_lo_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_lo - 1, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    TLWELvl2 pred_hi_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_hi, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    TLWELvl1 pred_asia_ct = tlweSymInt32Encrypt<Lvl1>(
        pred_asia, Lvl1::α, pow(2., region_scale_bits),
        sk.key.get<Lvl1>());

    // -------- WHERE in TFHE 3-PBS --------
    std::vector<TLWELvl1> order_mask(kKeyDomain), region_mask(kRegionDomain);
    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, ckks_time = 0;
    start = std::chrono::system_clock::now();
    for (size_t j = 0; j < kKeyDomain; j++) {
        TLWELvl1 gt, lt;
        greater_than<Lvl2>(orderdate_ciphers[j], pred_lo_ct, gt, compprecision,
                           ek, micro_pack, LOGIC);
        less_than<Lvl2>(orderdate_ciphers[j], pred_hi_ct, lt, compprecision,
                        ek, micro_pack, LOGIC);
        lift_and_and(gt, lt, order_mask[j], 29, ek);
    }
    for (size_t r = 0; r < kRegionDomain; r++) {
        TLWELvl1 raw;
        equal<Lvl1>(region_name_ciphers[r], pred_asia_ct, raw, region_bits,
                    ek, micro_pack, LOGIC);
        lift_and_and(raw, raw, region_mask[r], 29, ek);
    }
    end = std::chrono::system_clock::now();
    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                         end - start).count();

    // -------- plain reference (per-table masks + final aggregate) --------
    std::vector<uint64_t> plain_order_mask(kKeyDomain, 0),
        plain_region_mask(kRegionDomain, 0);
    for (size_t j = 0; j < kKeyDomain; j++)
        plain_order_mask[j] =
            (order_date[j] >= pred_lo && order_date[j] < pred_hi) ? 1 : 0;
    for (size_t r = 0; r < kRegionDomain; r++)
        plain_region_mask[r] = (region_name[r] == pred_asia) ? 1 : 0;

    std::vector<double> plain_agg(kNationDomain, 0.0);
    auto find_idx = [](const std::vector<uint64_t> &v, uint64_t key) {
        for (size_t i = 0; i < v.size(); i++) if (v[i] == key) return i;
        return v.size();
    };
    for (size_t i = 0; i < num; i++) {
        size_t oj = find_idx(order_orderkey, line_orderkey[i]);
        if (oj == kKeyDomain) continue;
        size_t cj = find_idx(cust_custkey, order_custkey[oj]);
        if (cj == kKeyDomain) continue;
        size_t sj = find_idx(supp_suppkey, line_suppkey[i]);
        if (sj == kKeyDomain) continue;
        if (cust_nationkey[cj] != supp_nationkey[sj]) continue;
        size_t nj = find_idx(nation_nationkey, cust_nationkey[cj]);
        if (nj == kNationDomain) continue;
        size_t rj = find_idx(region_regionkey, nation_regionkey[nj]);
        if (rj == kRegionDomain) continue;
        if (region_name[rj] != pred_asia) continue;
        if (!(order_date[oj] >= pred_lo && order_date[oj] < pred_hi)) continue;
        plain_agg[cust_nationkey[cj]] += revenue[i];
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
    // Q5 needs a deep chain: two repacked WHERE masks, several encrypted
    // lookup joins, nation equality, then encrypted group-by aggregation.
    std::vector<int> coeff_modulus_bits{
        59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59};
    parms.set_coeff_modulus(
        seal::CoeffModulus::Create(poly_modulus_degree, coeff_modulus_bits));
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

    // -------- 1. Repack TFHE WHERE masks --------
    std::cout << "Starting Conversion..." << std::endl;
    start = std::chrono::system_clock::now();
    auto order_mask_ckks = PackTfheMask(order_mask, kKeyDomain, repack_key,
                                         repack_config, ckks_encoder,
                                         galois_keys, relin_keys, evaluator,
                                         context);
    auto region_mask_ckks = PackTfheMask(region_mask, kRegionDomain, repack_key,
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
    report_repack(order_mask_ckks, plain_order_mask, "order_mask");
    report_repack(region_mask_ckks, plain_region_mask, "region_mask");

    // -------- 2. Encrypted-key joins and equality masks --------
    const double key_scale = std::pow(2.0, 40);
    auto line_orderkey_ct = EncryptAtLevel(
        SlotsFromIntegral(line_orderkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto line_suppkey_ct = EncryptAtLevel(
        SlotsFromIntegral(line_suppkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto order_orderkey_ct = EncryptAtLevel(
        SlotsFromIntegral(order_orderkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto order_custkey_ct = EncryptAtLevel(
        SlotsFromIntegral(order_custkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto cust_custkey_ct = EncryptAtLevel(
        SlotsFromIntegral(cust_custkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto cust_nationkey_ct = EncryptAtLevel(
        SlotsFromIntegral(cust_nationkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto supp_suppkey_ct = EncryptAtLevel(
        SlotsFromIntegral(supp_suppkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto supp_nationkey_ct = EncryptAtLevel(
        SlotsFromIntegral(supp_nationkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto nation_nationkey_ct = EncryptAtLevel(
        SlotsFromIntegral(nation_nationkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto nation_regionkey_ct = EncryptAtLevel(
        SlotsFromIntegral(nation_regionkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto region_regionkey_ct = EncryptAtLevel(
        SlotsFromIntegral(region_regionkey, ckks_encoder.slot_count()),
        context.first_parms_id(), key_scale, ckks_encoder, encryptor);
    auto revenue_ct = EncryptAtLevel(
        SlotsFromDouble(revenue, ckks_encoder.slot_count()),
        context.first_parms_id(), std::pow(2.0, 45), ckks_encoder, encryptor);

    std::vector<double> key_domain(kKeyDomain), nation_domain(kNationDomain),
        region_domain(kRegionDomain);
    for (size_t j = 0; j < kKeyDomain; ++j)
        key_domain[j] = static_cast<double>(j);
    for (size_t n = 0; n < kNationDomain; ++n)
        nation_domain[n] = static_cast<double>(n);
    for (size_t r = 0; r < kRegionDomain; ++r)
        region_domain[r] = static_cast<double>(r);

    auto line_orderkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        line_orderkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto line_suppkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        line_suppkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto order_orderkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        order_orderkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto order_custkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        order_custkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto cust_custkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        cust_custkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto cust_nation_masks = tfhepp_ckks::BuildLagrangeMasks(
        cust_nationkey_ct, nation_domain, relin_keys, ckks_encoder, evaluator);
    auto supp_suppkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        supp_suppkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto supp_nation_masks = tfhepp_ckks::BuildLagrangeMasks(
        supp_nationkey_ct, nation_domain, relin_keys, ckks_encoder, evaluator);
    auto nation_nationkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        nation_nationkey_ct, nation_domain, relin_keys, ckks_encoder,
        evaluator);
    auto nation_regionkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        nation_regionkey_ct, region_domain, relin_keys, ckks_encoder,
        evaluator);
    auto region_regionkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        region_regionkey_ct, region_domain, relin_keys, ckks_encoder,
        evaluator);

    auto asia_on_nation = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        nation_regionkey_masks, region_regionkey_masks, region_mask_ckks,
        kRegionDomain, kNationDomain, relin_keys, galois_keys, ckks_encoder,
        context, evaluator);

    // -------- 3. Low-depth per-nation encrypted filters and aggregation --------
    std::vector<seal::Ciphertext> grouped;
    grouped.reserve(kNationDomain);
    for (size_t n = 0; n < kNationDomain; ++n) {
        auto cust_nat_on_orders = tfhepp_ckks::LookupJoinFromEncryptedMasks(
            order_custkey_masks, cust_custkey_masks, cust_nation_masks[n],
            kKeyDomain, kKeyDomain, relin_keys, galois_keys, ckks_encoder,
            context, evaluator);
        auto order_nat_filter = MultiplyAndRescale(
            order_mask_ckks, cust_nat_on_orders, relin_keys, evaluator);
        auto order_nat_filter_on_line =
            tfhepp_ckks::LookupJoinFromEncryptedMasks(
                line_orderkey_masks, order_orderkey_masks, order_nat_filter,
                kKeyDomain, num, relin_keys, galois_keys, ckks_encoder,
                context, evaluator);

        auto supp_nat_on_line = tfhepp_ckks::LookupJoinFromEncryptedMasks(
            line_suppkey_masks, supp_suppkey_masks, supp_nation_masks[n],
            kKeyDomain, num, relin_keys, galois_keys, ckks_encoder, context,
            evaluator);

        auto group_filter = MultiplyAndRescale(
            order_nat_filter_on_line, supp_nat_on_line, relin_keys,
            evaluator);
        ApplyActiveSlotMaskInPlace(group_filter, num, ckks_encoder, context,
                                    evaluator);

        auto selected_revenue =
            MultiplyAndRescale(group_filter, revenue_ct, relin_keys,
                               evaluator);
        ApplyActiveSlotMaskInPlace(selected_revenue, num, ckks_encoder,
                                    context, evaluator);
        tfhepp_ckks::RotateAndSumInPlace(selected_revenue, num, galois_keys,
                                         evaluator);
        auto asia_scalar = MultiplyAndRescale(nation_nationkey_masks[n],
                                              asia_on_nation, relin_keys,
                                              evaluator);
        ReduceToFirstSlotInPlace(asia_scalar, kNationDomain, galois_keys,
                                  ckks_encoder, context, evaluator);
        selected_revenue =
            MultiplyAndRescale(selected_revenue, asia_scalar, relin_keys,
                               evaluator);
        grouped.push_back(std::move(selected_revenue));
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
    std::cout << std::setw(12) << "n_nationkey" << "|" << std::setw(12)
              << "revenue" << std::endl;
    for (size_t g = 0; g < kNationDomain; g++) {
        std::vector<double> dec;
        seal::Plaintext plain_g;
        decryptor.decrypt(grouped[g], plain_g);
        ckks_encoder.decode(plain_g, dec);
        std::cout << std::setw(12) << g << "|" << std::setw(12)
                  << std::round(dec[0]) << std::endl;
    }
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) << "n_nationkey" << "|" << std::setw(12)
              << "revenue" << std::endl;
    for (size_t g = 0; g < kNationDomain; g++) {
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
    std::cout << "TPC-H Q5 (strict Chapter 4): WHERE = TFHE 3-PBS, "
              << "JOIN+GROUP BY = CKKS Lagrange"
              << std::endl;
    std::cout << std::endl;
    relational_query5(num);
    return 0;
}
