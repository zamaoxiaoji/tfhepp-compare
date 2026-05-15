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
 *  Compliance contract (per Chapter 4):
 *    - WHERE predicates (o_orderdate range, r_name == 'ASIA') are the only
 *      TFHE 3-PBS comparisons.
 *    - All foreign keys (l_orderkey, l_suppkey, o_orderkey, o_custkey,
 *      c_custkey, c_nationkey, s_suppkey, s_nationkey, n_nationkey,
 *      n_regionkey, r_regionkey) are CKKS ciphertexts.
 *    - JOINs are realised by CKKS Lagrange masks +
 *      LookupJoinFromEncryptedMasks (Algorithm 4.3):
 *        nation⋈region on regionkey  → broadcasts r_mask onto nations
 *        customer⋈nation on nationkey → broadcasts asia flag onto customers
 *        customer⋈orders on custkey   → broadcasts onto orders
 *        orders⋈lineitem on orderkey  → broadcasts onto lineitem
 *        supplier⋈lineitem on suppkey → broadcasts s_nationkey onto lineitem
 *    - The c_nationkey == s_nationkey check is the inner product
 *        Σ_n mask_c_nat[n] ⊙ mask_s_nat[n]
 *      computed in CKKS, not in TFHE.
 *    - GROUP BY n_nationkey uses CKKS Lagrange indicator masks on the
 *      encrypted lineitem-side c_nationkey column.
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

seal::Ciphertext EncryptIntegralCt(const std::vector<uint64_t> &v, double scale,
                                   seal::CKKSEncoder &encoder,
                                   seal::Encryptor &encryptor)
{
    seal::Plaintext plain;
    encoder.encode(SlotsFromIntegral(v, encoder.slot_count()), scale, plain);
    seal::Ciphertext c;
    encryptor.encrypt(plain, c);
    return c;
}

} // namespace

double relational_query5(size_t num)
{
    std::cout << "Relational SQL Query5 Test (WHERE->TFHE 3-PBS pruned, "
              << "JOINs+GROUP BY->CKKS Lagrange):\n";
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

    // -------- TFHE encryption of WHERE-relevant columns only --------
    uint32_t num_bits = 16;
    uint32_t compprecision = 32;
    uint32_t scale_bits = std::numeric_limits<Lvl2::T>::digits - num_bits - 1;
    std::vector<TLWELvl2> orderdate_ciphers(kKeyDomain),
        region_name_ciphers(kRegionDomain);
    for (size_t j = 0; j < kKeyDomain; j++)
        orderdate_ciphers[j] = tlweSymInt32Encrypt<Lvl2>(
            order_date[j], Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    for (size_t r = 0; r < kRegionDomain; r++)
        region_name_ciphers[r] = tlweSymInt32Encrypt<Lvl2>(
            region_name[r], Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());

    uint64_t pred_lo = 12000, pred_hi = 13000;  // orderdate range
    uint64_t pred_asia = 1;                     // r_name == 'ASIA' flag
    TLWELvl2 pred_lo_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_lo, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    TLWELvl2 pred_hi_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_hi, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
    TLWELvl2 pred_asia_ct = tlweSymInt32Encrypt<Lvl2>(
        pred_asia, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());

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
        equal<Lvl2>(region_name_ciphers[r], pred_asia_ct, raw, compprecision,
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
            (order_date[j] > pred_lo && order_date[j] < pred_hi) ? 1 : 0;
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
        if (!(order_date[oj] > pred_lo && order_date[oj] < pred_hi)) continue;
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
    // Q5 needs a much deeper chain than Q3:
    //   repack (≈7) + 4 joins (≈8) + 4 ANDs (≈4) + group-by mul (1) ≈ 20+
    parms.set_coeff_modulus(seal::CoeffModulus::Create(
        poly_modulus_degree,
        {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45,
         45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
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

    // -------- 2. Encrypt all the join keys and the nation-grouping key --------
    double scale = std::pow(2.0, 40);
    auto line_orderkey_ct =
        EncryptIntegralCt(line_orderkey, scale, ckks_encoder, encryptor);
    auto line_suppkey_ct =
        EncryptIntegralCt(line_suppkey, scale, ckks_encoder, encryptor);
    auto order_orderkey_ct =
        EncryptIntegralCt(order_orderkey, scale, ckks_encoder, encryptor);
    auto order_custkey_ct =
        EncryptIntegralCt(order_custkey, scale, ckks_encoder, encryptor);
    auto cust_custkey_ct =
        EncryptIntegralCt(cust_custkey, scale, ckks_encoder, encryptor);
    auto cust_nationkey_ct =
        EncryptIntegralCt(cust_nationkey, scale, ckks_encoder, encryptor);
    auto supp_suppkey_ct =
        EncryptIntegralCt(supp_suppkey, scale, ckks_encoder, encryptor);
    auto supp_nationkey_ct =
        EncryptIntegralCt(supp_nationkey, scale, ckks_encoder, encryptor);
    auto nation_nationkey_ct =
        EncryptIntegralCt(nation_nationkey, scale, ckks_encoder, encryptor);
    auto nation_regionkey_ct =
        EncryptIntegralCt(nation_regionkey, scale, ckks_encoder, encryptor);
    auto region_regionkey_ct =
        EncryptIntegralCt(region_regionkey, scale, ckks_encoder, encryptor);

    // -------- 3. JOIN nation⋈region on regionkey → asia flag onto nation --------
    std::vector<double> region_domain(kRegionDomain);
    for (size_t r = 0; r < kRegionDomain; r++) region_domain[r] = double(r);
    auto nation_regionkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        nation_regionkey_ct, region_domain, relin_keys, ckks_encoder,
        evaluator);
    auto region_regionkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        region_regionkey_ct, region_domain, relin_keys, ckks_encoder,
        evaluator);
    auto asia_on_nation = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        nation_regionkey_masks, region_regionkey_masks, region_mask_ckks,
        kRegionDomain, relin_keys, galois_keys, evaluator);

    // -------- 4. JOIN customer⋈nation on nationkey → asia flag onto customer --------
    std::vector<double> nation_domain(kNationDomain);
    for (size_t n = 0; n < kNationDomain; n++) nation_domain[n] = double(n);
    auto cust_nationkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        cust_nationkey_ct, nation_domain, relin_keys, ckks_encoder, evaluator);
    auto nation_nationkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        nation_nationkey_ct, nation_domain, relin_keys, ckks_encoder,
        evaluator);
    auto asia_on_cust = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        cust_nationkey_masks, nation_nationkey_masks, asia_on_nation,
        kNationDomain, relin_keys, galois_keys, evaluator);

    // Also broadcast c_nationkey onto orders later via custkey join; first put
    // it on customer rows -- that's already in cust_nationkey_ct.

    // -------- 5. JOIN customer⋈orders on custkey: broadcast (asia, c_nat) onto orders --------
    std::vector<double> key_domain(kKeyDomain);
    for (size_t j = 0; j < kKeyDomain; j++) key_domain[j] = double(j);
    auto order_custkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        order_custkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto cust_custkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        cust_custkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto asia_on_orders = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        order_custkey_masks, cust_custkey_masks, asia_on_cust,
        kKeyDomain, relin_keys, galois_keys, evaluator);
    auto cust_nat_on_orders = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        order_custkey_masks, cust_custkey_masks, cust_nationkey_ct,
        kKeyDomain, relin_keys, galois_keys, evaluator);

    // AND order_mask × asia_on_orders → orders-level filter
    auto orders_filter = MultiplyAndRescale(order_mask_ckks, asia_on_orders,
                                             relin_keys, evaluator);

    // -------- 6. JOIN orders⋈lineitem on orderkey: broadcast onto lineitem --------
    auto line_orderkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        line_orderkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto order_orderkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        order_orderkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto orders_filter_on_line = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        line_orderkey_masks, order_orderkey_masks, orders_filter,
        kKeyDomain, relin_keys, galois_keys, evaluator);
    auto cust_nat_on_line = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        line_orderkey_masks, order_orderkey_masks, cust_nat_on_orders,
        kKeyDomain, relin_keys, galois_keys, evaluator);

    // -------- 7. JOIN supplier⋈lineitem on suppkey: broadcast s_nat onto lineitem --------
    auto line_suppkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        line_suppkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto supp_suppkey_masks = tfhepp_ckks::BuildLagrangeMasks(
        supp_suppkey_ct, key_domain, relin_keys, ckks_encoder, evaluator);
    auto supp_nat_on_line = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        line_suppkey_masks, supp_suppkey_masks, supp_nationkey_ct,
        kKeyDomain, relin_keys, galois_keys, evaluator);

    // -------- 8. c_nationkey == s_nationkey check, in CKKS --------
    //   Build Lagrange masks for cust_nat_on_line and supp_nat_on_line, then
    //   sum_n  mask_c[n] ⊙ mask_s[n].
    auto cust_nat_on_line_masks = tfhepp_ckks::BuildLagrangeMasks(
        cust_nat_on_line, nation_domain, relin_keys, ckks_encoder, evaluator);
    auto supp_nat_on_line_masks = tfhepp_ckks::BuildLagrangeMasks(
        supp_nat_on_line, nation_domain, relin_keys, ckks_encoder, evaluator);
    seal::Ciphertext nat_eq_mask;
    {
        auto first = MultiplyAndRescale(cust_nat_on_line_masks[0],
                                        supp_nat_on_line_masks[0], relin_keys,
                                        evaluator);
        nat_eq_mask = std::move(first);
        for (size_t n = 1; n < kNationDomain; n++) {
            auto term = MultiplyAndRescale(cust_nat_on_line_masks[n],
                                           supp_nat_on_line_masks[n],
                                           relin_keys, evaluator);
            ModSwitchToCommonLevel(nat_eq_mask, term, evaluator);
            term.scale() = nat_eq_mask.scale();
            evaluator.add_inplace(nat_eq_mask, term);
        }
    }

    // -------- 9. Combine into row_filter = orders_filter_on_line × nat_eq_mask --------
    auto row_filter = MultiplyAndRescale(orders_filter_on_line, nat_eq_mask,
                                          relin_keys, evaluator);

    // -------- 10. Multiply by revenue and group-by c_nat_on_line --------
    double qd = LastCoeffModulus(row_filter, context);
    auto revenue_ct = EncryptAtLevel(
        SlotsFromDouble(revenue, ckks_encoder.slot_count()),
        row_filter.parms_id(), qd, ckks_encoder, encryptor);
    auto filtered_revenue =
        MultiplyAndRescale(row_filter, revenue_ct, relin_keys, evaluator);
    // Zero out padding slots so noise outside the active region does not leak.
    ApplyActiveSlotMaskInPlace(filtered_revenue, num, ckks_encoder, context,
                                evaluator);

    int logrow = static_cast<int>(std::ceil(std::log2(num)));
    std::vector<seal::Ciphertext> grouped(kNationDomain);
    for (size_t g = 0; g < kNationDomain; g++) {
        auto weighted = MultiplyAndRescale(filtered_revenue,
                                            cust_nat_on_line_masks[g],
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
    std::cout << "CKKS (Repack + 5 Joins + Group + Agg) Time: " << ckks_time
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
    std::cout << "TPC-H Q5 (compliant): WHERE = TFHE 3-PBS, "
              << "JOINs + GROUP BY = CKKS Lagrange (Algorithm 4.3)"
              << std::endl;
    std::cout << std::endl;
    relational_query5(num);
    return 0;
}
