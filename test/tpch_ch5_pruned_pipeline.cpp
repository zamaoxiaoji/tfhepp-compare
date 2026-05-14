#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../comparison/comparison.h"
#include "ckks_relational.h"
#include "ckks_repack.h"

namespace
{
    using namespace tfhepp_compare;

    constexpr uint32_t kQuantityBits = 6;
    constexpr uint32_t kShipDateBits = 16;
    constexpr uint32_t kDiscountBits = 4;
    constexpr uint32_t kSmallKeyBits = 4;
    constexpr uint32_t kRepackScaleBits = 29;
    constexpr uint64_t kQ6ShipDateLo = 20101;
    constexpr uint64_t kQ6ShipDateHi = 21101;
    constexpr uint32_t kQ6DiscountLo = 8;
    constexpr uint32_t kQ6DiscountHi = 10;
    constexpr uint32_t kQ6QuantityHi = 10;
    constexpr uint64_t kQ14ShipDateLo = 20101;
    constexpr uint64_t kQ14ShipDateHi = 20201;
    constexpr uint64_t kQ3Date = 20101;
    constexpr uint64_t kQ5OrderDateLo = 20101;
    constexpr uint64_t kQ5OrderDateHi = 21231;

    struct Options {
        std::size_t rows = 16;
        std::string query = "all";
        uint32_t seed = 0x5eed5000u;
        bool encrypted = false;
    };

    struct LineitemRow {
        uint32_t orderkey = 0;
        uint32_t partkey = 0;
        uint32_t suppkey = 0;
        uint32_t quantity = 0;
        uint32_t discount = 0;
        uint64_t shipdate = 0;
        double extendedprice = 0.0;
    };

    struct OrdersRow {
        uint32_t orderkey = 0;
        uint32_t custkey = 0;
        uint32_t shippriority = 0;
        uint64_t orderdate = 0;
    };

    struct CustomerRow {
        uint32_t custkey = 0;
        uint32_t nationkey = 0;
        uint32_t mktsegment = 0;
    };

    struct PartRow {
        uint32_t partkey = 0;
        uint32_t promo = 0;
    };

    struct SupplierRow {
        uint32_t suppkey = 0;
        uint32_t nationkey = 0;
    };

    struct NationRow {
        uint32_t nationkey = 0;
        uint32_t regionkey = 0;
    };

    struct RegionRow {
        uint32_t regionkey = 0;
        uint32_t name = 0;
    };

    struct TpchData {
        std::vector<LineitemRow> lineitem;
        std::vector<OrdersRow> orders;
        std::vector<CustomerRow> customers;
        std::vector<PartRow> parts;
        std::vector<SupplierRow> suppliers;
        std::vector<NationRow> nations;
        std::vector<RegionRow> regions;
        std::size_t key_domain = 4;
        std::size_t nation_domain = 4;
        std::size_t priority_domain = 2;
    };

    struct MaskResult {
        std::vector<TLWELvl1> masks;
        std::vector<uint32_t> expected;
        double compare_ms = 0.0;
        std::size_t compare_errors = 0;
        std::size_t repack_errors = 0;
    };

    struct QueryPlainResult {
        std::string query;
        std::vector<double> values;
    };

    enum class CkksProfile {
        He3dbQ6,
        Relational,
        JoinSmoke,
        DeepJoin,
    };

    const char *ProfileName(CkksProfile profile)
    {
        switch (profile) {
        case CkksProfile::He3dbQ6:
            return "he3db_q6";
        case CkksProfile::Relational:
            return "relational";
        case CkksProfile::JoinSmoke:
            return "join_smoke";
        case CkksProfile::DeepJoin:
            return "deep_join";
        }
        return "unknown";
    }

    std::vector<int> RotationSteps(std::size_t slot_count)
    {
        std::vector<int> steps;
        for (std::size_t step = 1; step < slot_count; step <<= 1)
            steps.push_back(static_cast<int>(step));
        return steps;
    }

    std::size_t NextPowerOfTwo(std::size_t value)
    {
        std::size_t result = 1;
        while (result < value) result <<= 1;
        return result;
    }

    std::size_t PaddedPackSize(std::size_t active_slots)
    {
        return NextPowerOfTwo(std::max<std::size_t>(1, active_slots * 2));
    }

    std::vector<double> Domain(std::size_t size)
    {
        std::vector<double> domain(size);
        for (std::size_t i = 0; i < size; ++i)
            domain[i] = static_cast<double>(i);
        return domain;
    }

    uint64_t GenerateDate(std::default_random_engine &engine, uint64_t down,
                          uint64_t up)
    {
        const uint64_t dyear = down / 10000;
        const uint64_t dmonth = (down / 100) % 100;
        const uint64_t dday = down % 100;
        const uint64_t uyear = up / 10000;
        const uint64_t umonth = (up / 100) % 100;
        const uint64_t uday = up % 100;
        std::uniform_int_distribution<uint64_t> day_message(dday, uday);
        std::uniform_int_distribution<uint64_t> month_message(dmonth, umonth);
        std::uniform_int_distribution<uint64_t> year_message(dyear, uyear);
        return day_message(engine) + 100 * month_message(engine) +
               10000 * year_message(engine);
    }

    TpchData GenerateTpchLikeData(const Options &opts)
    {
        TpchData data;
        data.key_domain = std::max<std::size_t>(2, std::min<std::size_t>(8, opts.rows));
        data.nation_domain = std::max<std::size_t>(2, std::min<std::size_t>(4, data.key_domain));
        data.priority_domain = 2;

        std::default_random_engine engine(opts.seed);
        std::uniform_int_distribution<uint32_t> quantity_message(
            0, (uint32_t{1} << kQuantityBits) - 1);
        std::uniform_int_distribution<uint32_t> discount_message(
            0, (uint32_t{1} << kDiscountBits) - 1);
        std::uniform_int_distribution<uint32_t> key_message(
            0, static_cast<uint32_t>(data.key_domain - 1));
        std::uniform_int_distribution<uint32_t> nation_message(
            0, static_cast<uint32_t>(data.nation_domain - 1));
        std::uniform_int_distribution<uint32_t> bool_message(0, 1);
        std::uniform_int_distribution<uint32_t> segment_message(0, 3);
        std::uniform_int_distribution<uint64_t> extendedprice_message(1, 10);

        data.lineitem.resize(opts.rows);
        for (std::size_t i = 0; i < opts.rows; ++i) {
            auto &row = data.lineitem[i];
            row.orderkey = key_message(engine);
            row.partkey = key_message(engine);
            row.suppkey = key_message(engine);
            row.quantity = quantity_message(engine);
            row.discount = discount_message(engine);
            row.shipdate = GenerateDate(engine, 10101, 21230);
            row.extendedprice = static_cast<double>(extendedprice_message(engine));
        }

        // Match HE3DB's Q6 fixed row.
        if (!data.lineitem.empty()) {
            data.lineitem[0].quantity = 1;
            data.lineitem[0].discount = 9;
            data.lineitem[0].shipdate = 21215;
        }

        data.orders.resize(data.key_domain);
        for (std::size_t i = 0; i < data.orders.size(); ++i) {
            data.orders[i].orderkey = static_cast<uint32_t>(i);
            data.orders[i].custkey = key_message(engine);
            data.orders[i].shippriority = bool_message(engine);
            data.orders[i].orderdate = GenerateDate(engine, 10101, 21230);
        }

        data.customers.resize(data.key_domain);
        for (std::size_t i = 0; i < data.customers.size(); ++i) {
            data.customers[i].custkey = static_cast<uint32_t>(i);
            data.customers[i].nationkey = nation_message(engine);
            data.customers[i].mktsegment = segment_message(engine);
        }

        data.parts.resize(data.key_domain);
        for (std::size_t i = 0; i < data.parts.size(); ++i) {
            data.parts[i].partkey = static_cast<uint32_t>(i);
            data.parts[i].promo = bool_message(engine);
        }

        data.suppliers.resize(data.key_domain);
        for (std::size_t i = 0; i < data.suppliers.size(); ++i) {
            data.suppliers[i].suppkey = static_cast<uint32_t>(i);
            data.suppliers[i].nationkey = nation_message(engine);
        }

        data.nations.resize(data.nation_domain);
        for (std::size_t i = 0; i < data.nations.size(); ++i) {
            data.nations[i].nationkey = static_cast<uint32_t>(i);
            data.nations[i].regionkey = static_cast<uint32_t>(i % 2);
        }

        data.regions.resize(2);
        for (std::size_t i = 0; i < data.regions.size(); ++i) {
            data.regions[i].regionkey = static_cast<uint32_t>(i);
            data.regions[i].name = static_cast<uint32_t>(i);
        }

        // Keep HE3DB's Q6 fixed row, then add a second witness row for the
        // one-year Q6 window and the join/group-by queries when possible.
        if (!data.orders.empty() && !data.customers.empty() &&
            !data.parts.empty() && !data.suppliers.empty()) {
            data.orders[0].orderdate = 11215;
            data.orders[0].shippriority = 1;
            data.orders[0].custkey = 0;
            data.customers[0].mktsegment = 1;
            data.customers[0].nationkey = 0;
            data.parts[0].promo = 1;
            data.suppliers[0].nationkey = 0;
            data.nations[0].regionkey = 1;
            data.lineitem[0].orderkey = 0;
            data.lineitem[0].partkey = 0;
            data.lineitem[0].suppkey = 0;
        }
        if (data.lineitem.size() > 1 && data.orders.size() > 1 &&
            data.customers.size() > 1 && data.suppliers.size() > 1 &&
            data.nations.size() > 1 && !data.parts.empty()) {
            data.lineitem[1].quantity = 1;
            data.lineitem[1].discount = 9;
            data.lineitem[1].shipdate = 20115;
            data.lineitem[1].orderkey = 1;
            data.lineitem[1].partkey = 0;
            data.lineitem[1].suppkey = 1;
            data.orders[1].orderkey = 1;
            data.orders[1].custkey = 1;
            data.orders[1].orderdate = 20115;
            data.customers[1].nationkey = 1;
            data.suppliers[1].nationkey = 1;
            data.nations[1].regionkey = 1;
        }
        return data;
    }

    template <class P>
    uint32_t ScaleBits(uint32_t plain_bits)
    {
        return std::numeric_limits<typename P::T>::digits - plain_bits - 1;
    }

    template <class P>
    TFHEpp::TLWE<P> EncryptInteger(uint64_t value, uint32_t plain_bits,
                                   const TFHESecretKey &sk)
    {
        return tlweSymInt32Encrypt<P>(
            static_cast<typename P::T>(value), P::α,
            std::pow(2.0, ScaleBits<P>(plain_bits)), sk.key.get<P>());
    }

    void AriRescale(TLWELvl1 &res, const TLWELvl1 &tlwe,
                    uint32_t scale_bits, const TFHEEvalKey &ek)
    {
        const Lvl1::T mu = Lvl1::T{1} << (scale_bits - 1);
        constexpr uint64_t offset =
            uint64_t{1} << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlwe_offset = tlwe;
        tlwe_offset[Lvl1::k * Lvl1::n] += offset;
        TLWELvl0 tlwe_lvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwe_lvl0, tlwe_offset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwe_lvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(mu));
        res[Lvl1::k * Lvl1::n] += mu;
    }

    TLWELvl1 AndPredicatesToArithmetic(const std::vector<TLWELvl1> &predicates,
                                       const TFHEEvalKey &ek)
    {
        if (predicates.empty())
            throw std::invalid_argument("empty predicate list");
        if (predicates.size() == 1) {
            TLWELvl1 result;
            LOG_to_ARI(result, predicates.front(), ek);
            return result;
        }

        TLWELvl1 acc = predicates[0];
        for (std::size_t i = 1; i + 1 < predicates.size(); ++i) {
            TLWELvl1 next;
            HomAND(next, acc, predicates[i], ek, LOGIC);
            acc = next;
        }
        TLWELvl1 result;
        HomAND(result, acc, predicates.back(), ek, ARITHMETIC);
        return result;
    }

    std::size_t DecodeMask(const TLWELvl1 &mask, const TFHESecretKey &sk,
                           double scale)
    {
        return static_cast<std::size_t>(
            tlweSymInt32Decrypt<Lvl1>(mask, scale, sk.key.get<Lvl1>()));
    }

    template <typename T>
    std::vector<double> SlotsFromIntegral(const std::vector<T> &values,
                                          std::size_t slot_count)
    {
        std::vector<double> slots(slot_count, 0.0);
        for (std::size_t i = 0; i < values.size(); ++i)
            slots[i] = static_cast<double>(values[i]);
        return slots;
    }

    template <typename T>
    std::vector<double> SlotsFromIntegralWithInactive(
        const std::vector<T> &values, std::size_t slot_count,
        double inactive_value)
    {
        std::vector<double> slots(slot_count, inactive_value);
        for (std::size_t i = 0; i < values.size(); ++i)
            slots[i] = static_cast<double>(values[i]);
        return slots;
    }

    std::vector<double> SlotsFromDouble(const std::vector<double> &values,
                                        std::size_t slot_count)
    {
        std::vector<double> slots(slot_count, 0.0);
        for (std::size_t i = 0; i < values.size(); ++i) slots[i] = values[i];
        return slots;
    }

    seal::Ciphertext EncryptSlots(const std::vector<double> &slots, double scale,
                                  seal::CKKSEncoder &encoder,
                                  seal::Encryptor &encryptor)
    {
        seal::Plaintext plain;
        encoder.encode(slots, scale, plain);
        seal::Ciphertext cipher;
        encryptor.encrypt(plain, cipher);
        return cipher;
    }

    seal::Ciphertext EncryptSlotsAtLevel(const std::vector<double> &slots,
                                         seal::parms_id_type parms_id,
                                         double scale,
                                         seal::CKKSEncoder &encoder,
                                         seal::Encryptor &encryptor)
    {
        seal::Plaintext plain;
        encoder.encode(slots, parms_id, scale, plain);
        seal::Ciphertext cipher;
        encryptor.encrypt(plain, cipher);
        return cipher;
    }

    std::vector<double> DecryptSlots(const seal::Ciphertext &cipher,
                                     seal::Decryptor &decryptor,
                                     seal::CKKSEncoder &encoder)
    {
        seal::Plaintext plain;
        decryptor.decrypt(cipher, plain);
        std::vector<double> decoded;
        encoder.decode(plain, decoded);
        return decoded;
    }

    double LastCoeffModulus(const seal::Ciphertext &cipher,
                            const seal::SEALContext &context)
    {
        const auto context_data = context.get_context_data(cipher.parms_id());
        if (!context_data)
            throw std::invalid_argument("invalid ciphertext parms_id");
        const auto &moduli = context_data->parms().coeff_modulus();
        return static_cast<double>(moduli[cipher.coeff_modulus_size() - 1].value());
    }

    void ModSwitchToCommonLevel(seal::Ciphertext &lhs, seal::Ciphertext &rhs,
                                seal::Evaluator &evaluator)
    {
        if (lhs.parms_id() == rhs.parms_id()) return;
        if (lhs.coeff_modulus_size() > rhs.coeff_modulus_size())
            evaluator.mod_switch_to_inplace(lhs, rhs.parms_id());
        else if (rhs.coeff_modulus_size() > lhs.coeff_modulus_size())
            evaluator.mod_switch_to_inplace(rhs, lhs.parms_id());
        else
            throw std::invalid_argument("same-level parms_id mismatch");
    }

    seal::Ciphertext MultiplyAndRescale(const seal::Ciphertext &lhs_in,
                                        const seal::Ciphertext &rhs_in,
                                        const seal::RelinKeys &relin_keys,
                                        seal::Evaluator &evaluator)
    {
        seal::Ciphertext lhs = lhs_in;
        seal::Ciphertext rhs = rhs_in;
        ModSwitchToCommonLevel(lhs, rhs, evaluator);
        seal::Ciphertext result;
        evaluator.multiply(lhs, rhs, result);
        evaluator.relinearize_inplace(result, relin_keys);
        evaluator.rescale_to_next_inplace(result);
        return result;
    }

    void AddAlignedInPlace(seal::Ciphertext &acc,
                           const seal::Ciphertext &term_in,
                           seal::Evaluator &evaluator)
    {
        seal::Ciphertext term = term_in;
        ModSwitchToCommonLevel(acc, term, evaluator);
        term.scale() = acc.scale();
        evaluator.add_inplace(acc, term);
    }

    std::vector<double> Q6RevenueColumn(const TpchData &data)
    {
        std::vector<double> revenue(data.lineitem.size(), 0.0);
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            revenue[i] = data.lineitem[i].extendedprice *
                         static_cast<double>(data.lineitem[i].discount);
        }
        return revenue;
    }

    std::vector<double> DiscountedRevenueColumn(const TpchData &data)
    {
        std::vector<double> revenue(data.lineitem.size(), 0.0);
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            revenue[i] = data.lineitem[i].extendedprice *
                         (100.0 - static_cast<double>(data.lineitem[i].discount));
        }
        return revenue;
    }

    std::vector<uint32_t> LineitemOrderKeys(const TpchData &data)
    {
        std::vector<uint32_t> out(data.lineitem.size());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = data.lineitem[i].orderkey;
        return out;
    }

    std::vector<uint32_t> LineitemPartKeys(const TpchData &data)
    {
        std::vector<uint32_t> out(data.lineitem.size());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = data.lineitem[i].partkey;
        return out;
    }

    std::vector<uint32_t> LineitemSuppKeys(const TpchData &data)
    {
        std::vector<uint32_t> out(data.lineitem.size());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = data.lineitem[i].suppkey;
        return out;
    }

    template <typename Row, typename Getter>
    std::vector<uint32_t> ExtractU32(const std::vector<Row> &rows, Getter getter)
    {
        std::vector<uint32_t> out(rows.size());
        for (std::size_t i = 0; i < rows.size(); ++i) out[i] = getter(rows[i]);
        return out;
    }

    template <typename P>
    void PrunedGreaterThan(const TFHEpp::TLWE<P> &cipher1,
                           const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                           uint32_t plain_bits, const TFHEEvalKey &ek,
                           const three_pbs::FastB2AEvalKeyPack &micro_pack,
                           bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek, micro_pack,
                          result_type);
    }

    template <typename P>
    void PrunedGreaterThanEqual(
        const TFHEpp::TLWE<P> &cipher1, const TFHEpp::TLWE<P> &cipher2,
        TLWELvl1 &res, uint32_t plain_bits, const TFHEEvalKey &ek,
        const three_pbs::FastB2AEvalKeyPack &micro_pack, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek, micro_pack,
                          LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void PrunedLessThan(const TFHEpp::TLWE<P> &cipher1,
                        const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                        uint32_t plain_bits, const TFHEEvalKey &ek,
                        const three_pbs::FastB2AEvalKeyPack &micro_pack,
                        bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek, micro_pack,
                          result_type);
    }

    template <typename P>
    void PrunedLessThanEqual(
        const TFHEpp::TLWE<P> &cipher1, const TFHEpp::TLWE<P> &cipher2,
        TLWELvl1 &res, uint32_t plain_bits, const TFHEEvalKey &ek,
        const three_pbs::FastB2AEvalKeyPack &micro_pack, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe;
        for (size_t i = 0; i <= P::k * P::n; i++)
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek, micro_pack,
                          LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void PrunedEqual(const TFHEpp::TLWE<P> &cipher1,
                     const TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                     uint32_t plain_bits, const TFHEEvalKey &ek,
                     const three_pbs::FastB2AEvalKeyPack &micro_pack,
                     bool result_type)
    {
        TLWELvl1 ge_tlwe, le_tlwe;
        PrunedGreaterThanEqual<P>(cipher1, cipher2, ge_tlwe, plain_bits, ek,
                                  micro_pack, LOGIC);
        PrunedLessThanEqual<P>(cipher1, cipher2, le_tlwe, plain_bits, ek,
                               micro_pack, LOGIC);
        HomAND(res, ge_tlwe, le_tlwe, ek, result_type);
    }

    MaskResult EvalQ6Mask(const TpchData &data, const TFHESecretKey &sk,
                          TFHEEvalKey &ek,
                          const three_pbs::FastB2AEvalKeyPack &micro_pack)
    {
        auto pred_cipher1 =
            EncryptInteger<Lvl2>(kQ6ShipDateLo, kShipDateBits, sk);
        auto pred_cipher2 =
            EncryptInteger<Lvl2>(kQ6ShipDateHi - 1, kShipDateBits, sk);
        auto pred_cipher3 =
            EncryptInteger<Lvl1>(kQ6DiscountLo, kDiscountBits, sk);
        auto pred_cipher4 =
            EncryptInteger<Lvl1>(kQ6DiscountHi, kDiscountBits, sk);
        auto pred_cipher5 =
            EncryptInteger<Lvl1>(kQ6QuantityHi, kQuantityBits, sk);

        MaskResult result;
        result.masks.resize(data.lineitem.size());
        result.expected.resize(data.lineitem.size(), 0);

        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            const auto &row = data.lineitem[i];
            auto quantity =
                EncryptInteger<Lvl1>(row.quantity, kQuantityBits, sk);
            auto discount =
                EncryptInteger<Lvl1>(row.discount, kDiscountBits, sk);
            auto shipdate =
                EncryptInteger<Lvl2>(row.shipdate, kShipDateBits, sk);

            TLWELvl1 p1, p2, p3, p4, p5;
            three_pbs::greater_than_equal<Lvl2>(shipdate, pred_cipher1,
                                                p1, kShipDateBits, ek, LOGIC);
            three_pbs::greater_than_equal<Lvl2>(pred_cipher2, shipdate,
                                                p2, kShipDateBits, ek, LOGIC);
            three_pbs::greater_than_equal<Lvl1>(discount, pred_cipher3,
                                                p3, kDiscountBits, ek, LOGIC);
            three_pbs::less_than_equal<Lvl1>(discount, pred_cipher4,
                                             p4, kDiscountBits, ek, LOGIC);
            three_pbs::less_than<Lvl1>(quantity, pred_cipher5,
                                       p5, kQuantityBits, ek, LOGIC);
            result.masks[i] = AndPredicatesToArithmetic({p1, p2, p3, p4, p5}, ek);
            result.expected[i] =
                (row.shipdate >= kQ6ShipDateLo &&
                 row.shipdate < kQ6ShipDateHi &&
                 row.discount >= kQ6DiscountLo &&
                 row.discount <= kQ6DiscountHi &&
                 row.quantity < kQ6QuantityHi)
                    ? 1
                    : 0;
            const auto decoded_mask =
                DecodeMask(result.masks[i], sk, std::pow(2.0, 31));
            if (decoded_mask != result.expected[i]) {
                ++result.compare_errors;
                std::cerr << "debug=q6_compare_mismatch,row=" << i
                          << ",shipdate=" << row.shipdate
                          << ",discount=" << row.discount
                          << ",quantity=" << row.quantity
                          << ",expected=" << result.expected[i]
                          << ",got=" << decoded_mask
                          << ",p_ship_ge="
                          << TFHEpp::tlweSymDecrypt<Lvl1>(p1, sk.key.lvl1)
                          << ",p_ship_lt="
                          << TFHEpp::tlweSymDecrypt<Lvl1>(p2, sk.key.lvl1)
                          << ",p_discount_ge="
                          << TFHEpp::tlweSymDecrypt<Lvl1>(p3, sk.key.lvl1)
                          << ",p_discount_le="
                          << TFHEpp::tlweSymDecrypt<Lvl1>(p4, sk.key.lvl1)
                          << ",p_quantity_lt="
                          << TFHEpp::tlweSymDecrypt<Lvl1>(p5, sk.key.lvl1)
                          << "\n";
            }
        }
        const auto end = std::chrono::steady_clock::now();
        result.compare_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        return result;
    }

    MaskResult EvalQ14DateMask(const TpchData &data, const TFHESecretKey &sk,
                               TFHEEvalKey &ek)
    {
        auto clo = EncryptInteger<Lvl2>(kQ14ShipDateLo, kShipDateBits, sk);
        auto chi = EncryptInteger<Lvl2>(kQ14ShipDateHi - 1, kShipDateBits, sk);

        MaskResult result;
        result.masks.resize(data.lineitem.size());
        result.expected.resize(data.lineitem.size(), 0);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            auto shipdate =
                EncryptInteger<Lvl2>(data.lineitem[i].shipdate, kShipDateBits, sk);
            TLWELvl1 ge, lt;
            three_pbs::greater_than_equal<Lvl2>(shipdate, clo,
                                                ge, kShipDateBits, ek, LOGIC);
            three_pbs::greater_than_equal<Lvl2>(chi, shipdate,
                                                lt, kShipDateBits, ek, LOGIC);
            result.masks[i] = AndPredicatesToArithmetic({ge, lt}, ek);
            result.expected[i] =
                (data.lineitem[i].shipdate >= kQ14ShipDateLo &&
                 data.lineitem[i].shipdate < kQ14ShipDateHi)
                    ? 1
                    : 0;
            if (DecodeMask(result.masks[i], sk, std::pow(2.0, 31)) !=
                result.expected[i])
                ++result.compare_errors;
        }
        const auto end = std::chrono::steady_clock::now();
        result.compare_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        return result;
    }

    MaskResult EvalQ3LineitemMask(const TpchData &data, const TFHESecretKey &sk,
                                  TFHEEvalKey &ek)
    {
        auto cdate = EncryptInteger<Lvl2>(kQ3Date - 1, kShipDateBits, sk);
        MaskResult result;
        result.masks.resize(data.lineitem.size());
        result.expected.resize(data.lineitem.size(), 0);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            auto shipdate =
                EncryptInteger<Lvl2>(data.lineitem[i].shipdate, kShipDateBits, sk);
            TLWELvl1 gt;
            three_pbs::greater_than<Lvl2>(shipdate, cdate,
                                          gt, kShipDateBits, ek, LOGIC);
            result.masks[i] = AndPredicatesToArithmetic({gt}, ek);
            result.expected[i] = data.lineitem[i].shipdate > kQ3Date ? 1 : 0;
            if (DecodeMask(result.masks[i], sk, std::pow(2.0, 31)) !=
                result.expected[i])
                ++result.compare_errors;
        }
        const auto end = std::chrono::steady_clock::now();
        result.compare_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        return result;
    }

    MaskResult EvalQ3OrderDateMask(const TpchData &data,
                                   const TFHESecretKey &sk, TFHEEvalKey &ek)
    {
        auto cdate = EncryptInteger<Lvl2>(kQ3Date, kShipDateBits, sk);
        MaskResult result;
        result.masks.resize(data.orders.size());
        result.expected.resize(data.orders.size(), 0);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < data.orders.size(); ++i) {
            const auto &order = data.orders[i];
            auto orderdate = EncryptInteger<Lvl2>(order.orderdate, kShipDateBits, sk);
            TLWELvl1 date_ok;
            three_pbs::greater_than_equal<Lvl2>(cdate, orderdate, date_ok,
                                                kShipDateBits, ek, LOGIC);
            result.masks[i] = AndPredicatesToArithmetic({date_ok}, ek);
            result.expected[i] = (order.orderdate < kQ3Date) ? 1 : 0;
            if (DecodeMask(result.masks[i], sk, std::pow(2.0, 31)) !=
                result.expected[i])
                ++result.compare_errors;
        }
        const auto end = std::chrono::steady_clock::now();
        result.compare_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        return result;
    }

    MaskResult EvalQ3CustomerSegmentMask(const TpchData &data,
                                         const TFHESecretKey &sk,
                                         TFHEEvalKey &ek)
    {
        const uint32_t segment = 1;
        auto csegment = EncryptInteger<Lvl1>(segment, kSmallKeyBits, sk);
        MaskResult result;
        result.masks.resize(data.customers.size());
        result.expected.resize(data.customers.size(), 0);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < data.customers.size(); ++i) {
            auto mktsegment =
                EncryptInteger<Lvl1>(data.customers[i].mktsegment, kSmallKeyBits, sk);
            TLWELvl1 segment_ok;
            three_pbs::equal<Lvl1>(mktsegment, csegment, segment_ok,
                                   kSmallKeyBits, ek, LOGIC);
            result.masks[i] = AndPredicatesToArithmetic({segment_ok}, ek);
            result.expected[i] =
                (data.customers[i].mktsegment == segment) ? 1 : 0;
            if (DecodeMask(result.masks[i], sk, std::pow(2.0, 31)) !=
                result.expected[i])
                ++result.compare_errors;
        }
        const auto end = std::chrono::steady_clock::now();
        result.compare_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        return result;
    }

    MaskResult EvalQ5OrderMask(const TpchData &data, const TFHESecretKey &sk,
                               TFHEEvalKey &ek)
    {
        auto clo = EncryptInteger<Lvl2>(kQ5OrderDateLo, kShipDateBits, sk);
        auto chi = EncryptInteger<Lvl2>(kQ5OrderDateHi - 1, kShipDateBits, sk);
        MaskResult result;
        result.masks.resize(data.orders.size());
        result.expected.resize(data.orders.size(), 0);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < data.orders.size(); ++i) {
            auto orderdate =
                EncryptInteger<Lvl2>(data.orders[i].orderdate, kShipDateBits, sk);
            TLWELvl1 ge, lt;
            three_pbs::greater_than_equal<Lvl2>(orderdate, clo,
                                                ge, kShipDateBits, ek, LOGIC);
            three_pbs::greater_than_equal<Lvl2>(chi, orderdate,
                                                lt, kShipDateBits, ek, LOGIC);
            result.masks[i] = AndPredicatesToArithmetic({ge, lt}, ek);
            result.expected[i] =
                (data.orders[i].orderdate >= kQ5OrderDateLo &&
                 data.orders[i].orderdate < kQ5OrderDateHi)
                    ? 1
                    : 0;
            if (DecodeMask(result.masks[i], sk, std::pow(2.0, 31)) !=
                result.expected[i])
                ++result.compare_errors;
        }
        const auto end = std::chrono::steady_clock::now();
        result.compare_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        return result;
    }

    QueryPlainResult PlainQ6(const TpchData &data)
    {
        QueryPlainResult out{"q6", {0.0}};
        const auto revenue = Q6RevenueColumn(data);
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            const auto &row = data.lineitem[i];
            if (row.shipdate >= kQ6ShipDateLo &&
                row.shipdate < kQ6ShipDateHi &&
                row.discount >= kQ6DiscountLo &&
                row.discount <= kQ6DiscountHi &&
                row.quantity < kQ6QuantityHi)
                out.values[0] += revenue[i];
        }
        return out;
    }

    QueryPlainResult PlainQ14Groups(const TpchData &data)
    {
        QueryPlainResult out{"q14_groups", std::vector<double>(2, 0.0)};
        const auto revenue = DiscountedRevenueColumn(data);
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            const auto &line = data.lineitem[i];
            if (!(line.shipdate >= kQ14ShipDateLo &&
                  line.shipdate < kQ14ShipDateHi))
                continue;
            const auto promo = data.parts[line.partkey % data.parts.size()].promo;
            out.values[promo] += revenue[i];
        }
        return out;
    }

    double PromoRevenueRatio(const std::vector<double> &promo_groups)
    {
        if (promo_groups.size() < 2) return 0.0;
        const double total = promo_groups[0] + promo_groups[1];
        return std::abs(total) < 1e-9 ? 0.0 : 100.0 * promo_groups[1] / total;
    }

    QueryPlainResult PlainQ14(const TpchData &data)
    {
        const auto groups = PlainQ14Groups(data);
        return QueryPlainResult{"q14", {PromoRevenueRatio(groups.values)}};
    }

    QueryPlainResult PlainQ3(const TpchData &data)
    {
        QueryPlainResult out{
            "q3", std::vector<double>(data.key_domain * data.priority_domain, 0.0)};
        const auto revenue = DiscountedRevenueColumn(data);
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            const auto &line = data.lineitem[i];
            const auto &order = data.orders[line.orderkey % data.orders.size()];
            const auto &customer =
                data.customers[order.custkey % data.customers.size()];
            if (customer.mktsegment == 1 && order.orderdate < kQ3Date &&
                line.shipdate > kQ3Date) {
                const std::size_t group =
                    order.orderkey * data.priority_domain + order.shippriority;
                out.values[group] += revenue[i];
            }
        }
        return out;
    }

    QueryPlainResult PlainQ5(const TpchData &data)
    {
        QueryPlainResult out{"q5", std::vector<double>(data.nation_domain, 0.0)};
        const auto revenue = DiscountedRevenueColumn(data);
        for (std::size_t i = 0; i < data.lineitem.size(); ++i) {
            const auto &line = data.lineitem[i];
            const auto &order = data.orders[line.orderkey % data.orders.size()];
            const auto &customer =
                data.customers[order.custkey % data.customers.size()];
            const auto &supplier =
                data.suppliers[line.suppkey % data.suppliers.size()];
            const auto &nation =
                data.nations[customer.nationkey % data.nations.size()];
            if (order.orderdate >= kQ5OrderDateLo &&
                order.orderdate < kQ5OrderDateHi &&
                customer.nationkey == supplier.nationkey && nation.regionkey == 1) {
                out.values[customer.nationkey] += revenue[i];
            }
        }
        return out;
    }

    void PrintResult(const std::string &mode, const QueryPlainResult &expected,
                     const std::vector<double> &got, double compare_ms,
                     std::size_t compare_errors, std::size_t repack_errors)
    {
        double max_abs_error = 0.0;
        for (std::size_t i = 0; i < expected.values.size(); ++i) {
            max_abs_error =
                std::max(max_abs_error, std::abs(expected.values[i] - got[i]));
        }
        std::cout << "query=" << expected.query << ",mode=" << mode
                  << ",groups=" << expected.values.size()
                  << ",compare_ms=" << compare_ms
                  << ",compare_errors=" << compare_errors
                  << ",repack_errors=" << repack_errors
                  << ",max_abs_error=" << std::fixed << std::setprecision(6)
                  << max_abs_error << "\n";
    }

    void PrintPlainOnly(const QueryPlainResult &plain)
    {
        PrintResult("plain", plain, plain.values, 0.0, 0, 0);
    }

    struct CkksEnv {
        seal::EncryptionParameters parms;
        seal::SEALContext context;
        seal::KeyGenerator keygen;
        seal::SecretKey secret_key;
        seal::PublicKey public_key;
        seal::RelinKeys relin_keys;
        seal::GaloisKeys galois_keys;
        seal::Encryptor encryptor;
        seal::Encryptor symmetric_encryptor;
        seal::Decryptor decryptor;
        seal::Evaluator evaluator;
        seal::CKKSEncoder encoder;
        double scale = std::pow(2.0, 45);

        static seal::PublicKey MakePublicKey(seal::KeyGenerator &keygen)
        {
            seal::PublicKey public_key;
            keygen.create_public_key(public_key);
            return public_key;
        }

        static seal::EncryptionParameters MakeParameters()
        {
            return MakeParameters(CkksProfile::DeepJoin);
        }

        static seal::EncryptionParameters MakeParameters(CkksProfile profile)
        {
            if (profile == CkksProfile::He3dbQ6) {
                return tfhepp_ckks::MakeDefaultCKKSParameters(
                    65536, {59, 42, 42, 42, 42, 42, 42, 42, 45, 45,
                            45, 45, 45, 45, 45, 45, 45, 45, 45, 59});
            }
            if (profile == CkksProfile::Relational) {
                return tfhepp_ckks::MakeDefaultCKKSParameters(
                    65536, {59, 42, 42, 42, 42, 42, 42, 42, 45, 45,
                            45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
                            45, 45, 45, 45, 45, 59});
            }
            if (profile == CkksProfile::JoinSmoke) {
                return tfhepp_ckks::MakeDefaultCKKSParameters(
                    65536, {59, 42, 42, 42, 42, 42, 42, 42, 45, 45,
                            45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
                            45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
                            45, 59});
            }
            std::vector<int> coeff_modulus_bits{60};
            coeff_modulus_bits.insert(coeff_modulus_bits.end(), 35, 45);
            return tfhepp_ckks::MakeDefaultCKKSParameters(
                65536, coeff_modulus_bits);
        }

        static double MakeScale(CkksProfile profile)
        {
            return profile == CkksProfile::He3dbQ6 ? std::pow(2.0, 45)
                                                   : std::pow(2.0, 40);
        }

        explicit CkksEnv(CkksProfile profile = CkksProfile::DeepJoin)
            : parms(MakeParameters(profile)),
              context(parms, true, seal::sec_level_type::none),
              keygen(context),
              secret_key(keygen.secret_key()),
              public_key(MakePublicKey(keygen)),
              encryptor(context, public_key),
              symmetric_encryptor(context, secret_key),
              decryptor(context, secret_key),
              evaluator(context),
              encoder(context),
              scale(MakeScale(profile))
        {
            keygen.create_relin_keys(relin_keys);
            keygen.create_galois_keys(RotationSteps(encoder.slot_count()),
                                      galois_keys);
        }

        seal::Ciphertext encrypt(const std::vector<double> &active)
        {
            return EncryptSlots(SlotsFromDouble(active, encoder.slot_count()),
                                scale, encoder, encryptor);
        }

        seal::Ciphertext encrypt_for_multiply(
            const std::vector<double> &active,
            const seal::Ciphertext &reference)
        {
            const double modulus_bound = LastCoeffModulus(reference, context);
            const double encoded_scale =
                std::min(reference.scale(), modulus_bound / 4.0);
            return EncryptSlotsAtLevel(
                SlotsFromDouble(active, encoder.slot_count()),
                reference.parms_id(), encoded_scale, encoder, encryptor);
        }

        template <typename T>
        seal::Ciphertext encrypt_integral(const std::vector<T> &active)
        {
            return EncryptSlots(SlotsFromIntegral(active, encoder.slot_count()),
                                scale, encoder, encryptor);
        }

        template <typename T>
        seal::Ciphertext encrypt_integral_with_inactive(
            const std::vector<T> &active, double inactive_value)
        {
            return EncryptSlots(
                SlotsFromIntegralWithInactive(active, encoder.slot_count(),
                                              inactive_value),
                scale, encoder, encryptor);
        }
    };

    void DebugSlots(const char *label, const seal::Ciphertext &cipher,
                    std::size_t count, CkksEnv &ckks)
    {
        const auto decoded = DecryptSlots(cipher, ckks.decryptor, ckks.encoder);
        std::cerr << "debug=" << label << ",slots=";
        for (std::size_t i = 0; i < count; ++i) {
            if (i != 0) std::cerr << "/";
            std::cerr << decoded[i];
        }
        std::cerr << "\n";
    }

    void ApplyActiveSlotMaskInPlace(seal::Ciphertext &cipher,
                                    std::size_t active_slots,
                                    CkksEnv &ckks)
    {
        if (active_slots > ckks.encoder.slot_count())
            throw std::invalid_argument("active slot count exceeds CKKS slots");
        std::vector<double> slots(ckks.encoder.slot_count(), 0.0);
        std::fill(slots.begin(), slots.begin() + active_slots, 1.0);
        const double target_scale = cipher.scale();
        seal::Plaintext plain;
        ckks.encoder.encode(slots, cipher.parms_id(),
                            LastCoeffModulus(cipher, ckks.context), plain);
        ckks.evaluator.multiply_plain_inplace(cipher, plain);
        ckks.evaluator.rescale_to_next_inplace(cipher);
        cipher.scale() = target_scale;
    }

    seal::Ciphertext PackMaskToCkks(
        MaskResult &mask, const TFHESecretKey &tfhe_sk, TFHEEvalKey &tfhe_ek,
        tfhepp_ckks::RepackEvaluationKey &repack_key,
        tfhepp_ckks::RepackConfig &repack_config, CkksEnv &ckks)
    {
        std::vector<TLWELvl1> rescaled(mask.masks.size());
        for (std::size_t i = 0; i < mask.masks.size(); ++i) {
            AriRescale(rescaled[i], mask.masks[i], kRepackScaleBits, tfhe_ek);
            const auto decoded =
                DecodeMask(rescaled[i], tfhe_sk, std::pow(2.0, kRepackScaleBits));
            if (decoded != mask.expected[i]) ++mask.repack_errors;
        }

        seal::Ciphertext packed;
        std::vector<TLWELvl1> padded_rescaled =
            std::move(rescaled);
        padded_rescaled.resize(PaddedPackSize(mask.masks.size()));
        try {
            tfhepp_ckks::PackLWEsToCKKS<Lvl1>(
                packed, padded_rescaled, repack_key, repack_config, ckks.encoder,
                ckks.galois_keys, ckks.relin_keys, ckks.evaluator,
                ckks.context);
        }
        catch (const std::exception &e) {
            throw std::runtime_error(std::string("PackLWEsToCKKS failed: ") +
                                     e.what());
        }
        try {
            tfhepp_ckks::HomomorphicRound(packed, packed.scale(), ckks.encoder,
                                          ckks.relin_keys, ckks.evaluator,
                                          ckks.context);
        }
        catch (const std::exception &e) {
            throw std::runtime_error(std::string("HomomorphicRound failed: ") +
                                     e.what());
        }
        const auto decoded = DecryptSlots(packed, ckks.decryptor, ckks.encoder);
        double max_mask_error = 0.0;
        for (std::size_t i = 0; i < mask.expected.size(); ++i) {
            max_mask_error = std::max(
                max_mask_error,
                std::abs(decoded[i] - static_cast<double>(mask.expected[i])));
        }
        std::cerr << "debug=packed_mask,max_abs_error=" << max_mask_error
                  << ",scale_log2=" << std::log2(packed.scale())
                  << ",level=" << packed.coeff_modulus_size() << "\n";
        return packed;
    }

    std::vector<seal::Ciphertext> BuildMasksForColumn(
        const seal::Ciphertext &column, std::size_t domain_size, CkksEnv &ckks)
    {
        return tfhepp_ckks::BuildLagrangeMasks(
            column, Domain(domain_size), ckks.relin_keys, ckks.encoder,
            ckks.evaluator);
    }

    std::vector<seal::Ciphertext> GroupByFilteredRevenue(
        const seal::Ciphertext &row_filter,
        const std::vector<seal::Ciphertext> &group_masks,
        const std::vector<double> &revenue, std::size_t active_slots,
        CkksEnv &ckks)
    {
        std::vector<seal::Ciphertext> sums;
        sums.reserve(group_masks.size());
        for (const auto &group_mask : group_masks) {
            auto revenue_ct =
                ckks.encrypt_for_multiply(revenue, group_mask);
            auto weighted_group =
                MultiplyAndRescale(group_mask, revenue_ct, ckks.relin_keys,
                                   ckks.evaluator);
            auto selected =
                MultiplyAndRescale(row_filter, weighted_group,
                                   ckks.relin_keys, ckks.evaluator);
            tfhepp_ckks::RotateAndSumInPlace(
                selected, active_slots, ckks.galois_keys, ckks.evaluator);
            sums.push_back(std::move(selected));
        }
        return sums;
    }

    seal::Ciphertext LookupJoinPayload(
        const std::vector<uint32_t> &left_keys,
        const std::vector<uint32_t> &right_keys,
        const std::vector<double> &right_payload,
        std::size_t domain_size, CkksEnv &ckks)
    {
        auto left_key_ct = ckks.encrypt_integral(left_keys);
        auto right_key_ct = ckks.encrypt_integral(right_keys);
        auto payload_ct = ckks.encrypt(right_payload);
        auto left_masks = BuildMasksForColumn(left_key_ct, domain_size, ckks);
        auto right_masks = BuildMasksForColumn(right_key_ct, domain_size, ckks);
        auto joined = tfhepp_ckks::LookupJoinFromEncryptedMasks(
            left_masks, right_masks, payload_ct, right_keys.size(),
            ckks.relin_keys, ckks.galois_keys, ckks.evaluator);
        ApplyActiveSlotMaskInPlace(joined, left_keys.size(), ckks);
        return joined;
    }

    seal::Ciphertext LookupJoinPayloadCipher(
        const seal::Ciphertext &left_key_ct,
        const seal::Ciphertext &right_key_ct,
        const seal::Ciphertext &right_payload_ct, std::size_t left_active_slots,
        std::size_t right_active_slots, std::size_t domain_size, CkksEnv &ckks)
    {
        auto left_masks = BuildMasksForColumn(left_key_ct, domain_size, ckks);
        auto right_masks = BuildMasksForColumn(right_key_ct, domain_size, ckks);
        auto joined = tfhepp_ckks::LookupJoinFromEncryptedMasks(
            left_masks, right_masks, right_payload_ct, right_active_slots,
            ckks.relin_keys, ckks.galois_keys, ckks.evaluator);
        ApplyActiveSlotMaskInPlace(joined, left_active_slots, ckks);
        return joined;
    }

    seal::Ciphertext LookupJoinPayloadCipher(
        const std::vector<uint32_t> &left_keys,
        const std::vector<uint32_t> &right_keys,
        const seal::Ciphertext &right_payload_ct,
        std::size_t domain_size, CkksEnv &ckks)
    {
        auto left_key_ct = ckks.encrypt_integral(left_keys);
        auto right_key_ct = ckks.encrypt_integral(right_keys);
        return LookupJoinPayloadCipher(left_key_ct, right_key_ct,
                                       right_payload_ct, left_keys.size(),
                                       right_keys.size(), domain_size, ckks);
    }

    seal::Ciphertext LookupJoinPayloadFromCipherKey(
        const seal::Ciphertext &left_key_ct,
        const std::vector<uint32_t> &right_keys,
        const std::vector<double> &right_payload,
        std::size_t left_active_slots, std::size_t domain_size, CkksEnv &ckks)
    {
        auto right_key_ct = ckks.encrypt_integral(right_keys);
        auto payload_ct = ckks.encrypt(right_payload);
        return LookupJoinPayloadCipher(left_key_ct, right_key_ct, payload_ct,
                                       left_active_slots, right_keys.size(),
                                       domain_size, ckks);
    }

    QueryPlainResult EncryptedQ6(const TpchData &data, const TFHESecretKey &tfhe_sk,
                                 TFHEEvalKey &tfhe_ek,
                                 tfhepp_ckks::RepackEvaluationKey &repack_key,
                                 tfhepp_ckks::RepackConfig &repack_config,
                                 CkksEnv &ckks)
    {
        std::cerr << "stage=q6,start\n";
        auto mask = EvalQ6Mask(data, tfhe_sk, tfhe_ek);
        auto mask_ct =
            PackMaskToCkks(mask, tfhe_sk, tfhe_ek, repack_key, repack_config, ckks);
        std::cerr << "debug=q6_mask,scale_log2=" << std::log2(mask_ct.scale())
                  << ",level=" << mask_ct.coeff_modulus_size() << "\n";
        auto revenue_ct =
            ckks.encrypt_for_multiply(Q6RevenueColumn(data), mask_ct);
        std::cerr << "debug=q6_revenue,scale_log2=" << std::log2(revenue_ct.scale())
                  << ",level=" << revenue_ct.coeff_modulus_size() << "\n";
        auto filtered =
            MultiplyAndRescale(mask_ct, revenue_ct, ckks.relin_keys, ckks.evaluator);
        auto sum_ct = tfhepp_ckks::RotateAndSum(
            filtered, data.lineitem.size(), ckks.galois_keys, ckks.evaluator);
        const auto decoded = DecryptSlots(sum_ct, ckks.decryptor, ckks.encoder);
        QueryPlainResult got{"q6", {decoded[0]}};
        const auto expected = PlainQ6(data);
        PrintResult("encrypted", expected, got.values, mask.compare_ms,
                    mask.compare_errors, mask.repack_errors);
        return got;
    }

    QueryPlainResult EncryptedQ14(const TpchData &data,
                                  const TFHESecretKey &tfhe_sk,
                                  TFHEEvalKey &tfhe_ek,
                                  tfhepp_ckks::RepackEvaluationKey &repack_key,
                                  tfhepp_ckks::RepackConfig &repack_config,
                                  CkksEnv &ckks)
    {
        std::cerr << "stage=q14,start\n";
        auto date_mask = EvalQ14DateMask(data, tfhe_sk, tfhe_ek);
        auto date_mask_ct = PackMaskToCkks(date_mask, tfhe_sk, tfhe_ek,
                                           repack_key, repack_config, ckks);
        std::cerr << "stage=q14,date_mask_packed,scale_log2="
                  << std::log2(date_mask_ct.scale())
                  << ",level=" << date_mask_ct.coeff_modulus_size() << "\n";

        const auto line_part_keys = LineitemPartKeys(data);
        const auto part_keys = ExtractU32(data.parts, [](const PartRow &r) {
            return r.partkey;
        });
        std::vector<double> promo_payload(data.parts.size());
        for (std::size_t i = 0; i < data.parts.size(); ++i)
            promo_payload[i] = static_cast<double>(data.parts[i].promo);
        std::cerr << "stage=q14,promo_join,start\n";
        auto promo_on_line = LookupJoinPayload(line_part_keys, part_keys,
                                               promo_payload, data.key_domain, ckks);
        std::cerr << "stage=q14,promo_join,done,scale_log2="
                  << std::log2(promo_on_line.scale())
                  << ",level=" << promo_on_line.coeff_modulus_size() << "\n";

        std::cerr << "stage=q14,revenue_filter,start\n";
        auto revenue_ct =
            ckks.encrypt_for_multiply(DiscountedRevenueColumn(data), date_mask_ct);
        auto filtered_revenue =
            MultiplyAndRescale(date_mask_ct, revenue_ct, ckks.relin_keys,
                               ckks.evaluator);
        std::cerr << "stage=q14,revenue_filter,done,scale_log2="
                  << std::log2(filtered_revenue.scale())
                  << ",level=" << filtered_revenue.coeff_modulus_size()
                  << "\n";
        std::cerr << "stage=q14,promo_masks,start\n";
        auto promo_masks = BuildMasksForColumn(promo_on_line, 2, ckks);
        std::cerr << "stage=q14,promo_masks,done\n";
        std::cerr << "stage=q14,groupby,start\n";
        auto sums = tfhepp_ckks::GroupBySumFromEncryptedMasks(
            filtered_revenue, promo_masks, data.lineitem.size(), ckks.relin_keys,
            ckks.galois_keys, ckks.evaluator);
        std::cerr << "stage=q14,groupby,done\n";
        std::vector<double> group_values(2, 0.0);
        for (std::size_t i = 0; i < sums.size(); ++i)
            group_values[i] =
                DecryptSlots(sums[i], ckks.decryptor, ckks.encoder)[0];
        QueryPlainResult got{"q14", {PromoRevenueRatio(group_values)}};
        const auto expected = PlainQ14(data);
        PrintResult("encrypted", expected, got.values, date_mask.compare_ms,
                    date_mask.compare_errors, date_mask.repack_errors);
        return got;
    }

    QueryPlainResult EncryptedQ3(const TpchData &data, const TFHESecretKey &tfhe_sk,
                                 TFHEEvalKey &tfhe_ek,
                                 tfhepp_ckks::RepackEvaluationKey &repack_key,
                                 tfhepp_ckks::RepackConfig &repack_config,
                                 CkksEnv &ckks)
    {
        std::cerr << "stage=q3,start\n";
        auto line_mask = EvalQ3LineitemMask(data, tfhe_sk, tfhe_ek);
        auto order_date_mask = EvalQ3OrderDateMask(data, tfhe_sk, tfhe_ek);
        auto customer_segment_mask =
            EvalQ3CustomerSegmentMask(data, tfhe_sk, tfhe_ek);
        auto line_mask_ct =
            PackMaskToCkks(line_mask, tfhe_sk, tfhe_ek, repack_key, repack_config, ckks);
        auto order_date_mask_ct = PackMaskToCkks(order_date_mask, tfhe_sk,
                                                 tfhe_ek, repack_key,
                                                 repack_config, ckks);
        auto customer_segment_mask_ct =
            PackMaskToCkks(customer_segment_mask, tfhe_sk, tfhe_ek,
                           repack_key, repack_config, ckks);
        DebugSlots("q3_line_mask", line_mask_ct, data.lineitem.size(), ckks);
        DebugSlots("q3_order_date_mask", order_date_mask_ct,
                   data.orders.size(), ckks);
        DebugSlots("q3_customer_segment_mask", customer_segment_mask_ct,
                   data.customers.size(), ckks);

        const auto line_order_keys = LineitemOrderKeys(data);
        const auto order_keys = ExtractU32(data.orders, [](const OrdersRow &r) {
            return r.orderkey;
        });
        const auto order_customer_keys =
            ExtractU32(data.orders, [](const OrdersRow &r) {
                return r.custkey;
            });
        const auto customer_keys =
            ExtractU32(data.customers, [](const CustomerRow &r) {
                return r.custkey;
            });
        const auto order_priorities =
            ExtractU32(data.orders, [](const OrdersRow &r) {
                return r.shippriority;
            });
        std::vector<double> priority_payload(data.orders.size());
        for (std::size_t i = 0; i < data.orders.size(); ++i) {
            priority_payload[i] = static_cast<double>(order_priorities[i]);
        }

        auto segment_on_order = LookupJoinPayloadCipher(
            order_customer_keys, customer_keys, customer_segment_mask_ct,
            data.key_domain, ckks);
        std::cerr << "debug=q3_segment_on_order,scale_log2="
                  << std::log2(segment_on_order.scale())
                  << ",level=" << segment_on_order.coeff_modulus_size()
                  << "\n";
        DebugSlots("q3_segment_on_order", segment_on_order,
                   data.orders.size(), ckks);
        auto order_mask_ct =
            MultiplyAndRescale(order_date_mask_ct, segment_on_order,
                               ckks.relin_keys, ckks.evaluator);
        std::cerr << "debug=q3_order_mask,scale_log2="
                  << std::log2(order_mask_ct.scale())
                  << ",level=" << order_mask_ct.coeff_modulus_size() << "\n";
        auto order_mask_on_line = LookupJoinPayloadCipher(
            line_order_keys, order_keys, order_mask_ct, data.key_domain, ckks);
        std::cerr << "debug=q3_order_mask_on_line,scale_log2="
                  << std::log2(order_mask_on_line.scale())
                  << ",level=" << order_mask_on_line.coeff_modulus_size()
                  << "\n";
        DebugSlots("q3_order_mask_on_line", order_mask_on_line,
                   data.lineitem.size(), ckks);
        auto priority_on_line = LookupJoinPayload(
            line_order_keys, order_keys, priority_payload, data.key_domain, ckks);
        std::cerr << "debug=q3_priority_on_line,scale_log2="
                  << std::log2(priority_on_line.scale())
                  << ",level=" << priority_on_line.coeff_modulus_size()
                  << "\n";
        DebugSlots("q3_priority_on_line", priority_on_line,
                   data.lineitem.size(), ckks);

        auto filtered =
            MultiplyAndRescale(line_mask_ct, order_mask_on_line, ckks.relin_keys,
                               ckks.evaluator);
        std::cerr << "debug=q3_filtered_mask,scale_log2="
                  << std::log2(filtered.scale())
                  << ",level=" << filtered.coeff_modulus_size() << "\n";
        DebugSlots("q3_filtered_mask", filtered, data.lineitem.size(), ckks);
        auto order_key_ct = ckks.encrypt_integral(line_order_keys);
        std::vector<seal::Ciphertext> group_columns{order_key_ct,
                                                    priority_on_line};
        std::vector<std::vector<double>> group_domains{
            Domain(data.key_domain), Domain(data.priority_domain)};
        auto order_priority_masks = tfhepp_ckks::BuildTensorLagrangeMasks(
            group_columns, group_domains,
            ckks.relin_keys, ckks.encoder, ckks.evaluator);
        auto sums = GroupByFilteredRevenue(
            filtered, order_priority_masks, DiscountedRevenueColumn(data),
            data.lineitem.size(), ckks);
        QueryPlainResult got{
            "q3", std::vector<double>(data.key_domain * data.priority_domain)};
        for (std::size_t i = 0; i < sums.size(); ++i)
            got.values[i] = DecryptSlots(sums[i], ckks.decryptor, ckks.encoder)[0];
        const auto expected = PlainQ3(data);
        PrintResult("encrypted", expected, got.values,
                    line_mask.compare_ms + order_date_mask.compare_ms +
                        customer_segment_mask.compare_ms,
                    line_mask.compare_errors + order_date_mask.compare_errors +
                        customer_segment_mask.compare_errors,
                    line_mask.repack_errors + order_date_mask.repack_errors +
                        customer_segment_mask.repack_errors);
        return got;
    }

    QueryPlainResult EncryptedQ5(const TpchData &data, const TFHESecretKey &tfhe_sk,
                                 TFHEEvalKey &tfhe_ek,
                                 tfhepp_ckks::RepackEvaluationKey &repack_key,
                                 tfhepp_ckks::RepackConfig &repack_config,
                                 CkksEnv &ckks)
    {
        std::cerr << "stage=q5,start\n";
        auto order_mask = EvalQ5OrderMask(data, tfhe_sk, tfhe_ek);
        auto order_mask_ct = PackMaskToCkks(order_mask, tfhe_sk, tfhe_ek,
                                            repack_key, repack_config, ckks);

        const auto line_order_keys = LineitemOrderKeys(data);
        const auto line_supp_keys = LineitemSuppKeys(data);
        const auto order_keys = ExtractU32(data.orders, [](const OrdersRow &r) {
            return r.orderkey;
        });
        const auto order_customer_keys =
            ExtractU32(data.orders, [](const OrdersRow &r) {
                return r.custkey;
            });
        const auto customer_keys =
            ExtractU32(data.customers, [](const CustomerRow &r) {
                return r.custkey;
            });
        const auto supplier_keys =
            ExtractU32(data.suppliers, [](const SupplierRow &r) {
                return r.suppkey;
            });
        std::vector<double> cust_nation_payload(data.customers.size());
        std::vector<double> order_cust_payload(data.orders.size());
        for (std::size_t i = 0; i < data.orders.size(); ++i) {
            order_cust_payload[i] = static_cast<double>(data.orders[i].custkey);
        }
        for (std::size_t i = 0; i < data.customers.size(); ++i) {
            cust_nation_payload[i] =
                static_cast<double>(data.customers[i].nationkey);
        }
        std::vector<double> supp_nation_payload(data.suppliers.size());
        for (std::size_t i = 0; i < data.suppliers.size(); ++i)
            supp_nation_payload[i] =
                static_cast<double>(data.suppliers[i].nationkey);

        auto order_mask_on_line = LookupJoinPayloadCipher(
            line_order_keys, order_keys, order_mask_ct, data.key_domain, ckks);
        auto order_customer_on_line = LookupJoinPayload(
            line_order_keys, order_keys, order_cust_payload, data.key_domain, ckks);
        auto customer_nation_on_line = LookupJoinPayloadFromCipherKey(
            order_customer_on_line, customer_keys, cust_nation_payload,
            data.lineitem.size(), data.key_domain, ckks);
        auto supplier_nation_on_line = LookupJoinPayload(
            line_supp_keys, supplier_keys, supp_nation_payload, data.key_domain, ckks);

        auto nation_masks = BuildMasksForColumn(customer_nation_on_line,
                                                data.nation_domain, ckks);
        auto supplier_nation_masks = BuildMasksForColumn(supplier_nation_on_line,
                                                         data.nation_domain, ckks);
        seal::Ciphertext same_nation;
        bool initialized = false;
        for (std::size_t i = 0; i < data.nation_domain; ++i) {
            auto term = MultiplyAndRescale(nation_masks[i], supplier_nation_masks[i],
                                           ckks.relin_keys, ckks.evaluator);
            if (!initialized) {
                same_nation = term;
                initialized = true;
            }
            else {
                AddAlignedInPlace(same_nation, term, ckks.evaluator);
            }
        }

        std::vector<double> region_ok_by_nation(data.nation_domain, 0.0);
        for (const auto &nation : data.nations)
            region_ok_by_nation[nation.nationkey] =
                (nation.regionkey == 1) ? 1.0 : 0.0;
        auto nation_keys = ExtractU32(data.nations, [](const NationRow &r) {
            return r.nationkey;
        });
        auto region_on_line = LookupJoinPayloadFromCipherKey(
            customer_nation_on_line, nation_keys, region_ok_by_nation,
            data.lineitem.size(), data.nation_domain, ckks);

        auto filtered =
            MultiplyAndRescale(order_mask_on_line, same_nation, ckks.relin_keys,
                               ckks.evaluator);
        filtered =
            MultiplyAndRescale(filtered, region_on_line, ckks.relin_keys,
                               ckks.evaluator);
        auto sums = GroupByFilteredRevenue(
            filtered, nation_masks, DiscountedRevenueColumn(data),
            data.lineitem.size(), ckks);
        QueryPlainResult got{"q5", std::vector<double>(data.nation_domain)};
        for (std::size_t i = 0; i < sums.size(); ++i)
            got.values[i] = DecryptSlots(sums[i], ckks.decryptor, ckks.encoder)[0];
        const auto expected = PlainQ5(data);
        PrintResult("encrypted", expected, got.values, order_mask.compare_ms,
                    order_mask.compare_errors, order_mask.repack_errors);
        return got;
    }

    bool WantsQuery(const Options &opts, const std::string &query)
    {
        return opts.query == "all" || opts.query == query;
    }

    CkksProfile SelectCkksProfile(const Options &opts)
    {
        if (opts.query == "q6")
            return CkksProfile::He3dbQ6;
        if (opts.query == "q14")
            return CkksProfile::Relational;
        if (opts.query == "q3" || opts.query == "q5")
            return CkksProfile::JoinSmoke;
        return CkksProfile::DeepJoin;
    }

    Options ParseOptions(int argc, char **argv)
    {
        Options opts;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            const auto read_value = [&](const char *name) -> std::string {
                if (i + 1 >= argc)
                    throw std::invalid_argument(std::string("missing value for ") +
                                                name);
                return argv[++i];
            };
            if (arg == "--rows") {
                opts.rows = static_cast<std::size_t>(
                    std::stoull(read_value("--rows")));
            }
            else if (arg == "--query") {
                opts.query = read_value("--query");
            }
            else if (arg == "--seed") {
                opts.seed = static_cast<uint32_t>(
                    std::stoul(read_value("--seed")));
            }
            else if (arg == "--encrypted") {
                opts.encrypted = true;
            }
            else if (arg == "--plain-only") {
                opts.encrypted = false;
            }
            else {
                throw std::invalid_argument("unknown option: " + arg);
            }
        }
        if (opts.rows == 0) throw std::invalid_argument("--rows must be positive");
        if (opts.rows > 16 && opts.encrypted)
            throw std::invalid_argument(
                "encrypted smoke is intentionally capped at 16 rows");
        return opts;
    }

} // namespace

int main(int argc, char **argv)
{
    const Options opts = ParseOptions(argc, argv);
    const TpchData data = GenerateTpchLikeData(opts);

    std::cout << "tpch_3pbs_q6_q14_q3_q5"
              << ",rows=" << opts.rows
              << ",seed=" << opts.seed
              << ",encrypted=" << (opts.encrypted ? "true" : "false")
              << ",quantity_bits=" << kQuantityBits
              << ",ship_bits=" << kShipDateBits
              << ",discount_bits=" << kDiscountBits
              << ",compare=pruned_3pbs\n";

    if (!opts.encrypted) {
        if (WantsQuery(opts, "q6")) PrintPlainOnly(PlainQ6(data));
        if (WantsQuery(opts, "q14")) PrintPlainOnly(PlainQ14(data));
        if (WantsQuery(opts, "q3")) PrintPlainOnly(PlainQ3(data));
        if (WantsQuery(opts, "q5")) PrintPlainOnly(PlainQ5(data));
        return 0;
    }

    TFHESecretKey tfhe_sk;
    TFHEEvalKey tfhe_ek;
    std::cerr << "stage=tfhe_eval_keygen,start\n";
    tfhe_ek.emplacebkfft<Lvl01>(tfhe_sk);
    std::cerr << "stage=tfhe_eval_keygen,bkfft_lvl01_done\n";
    tfhe_ek.emplacebkfft<Lvl02>(tfhe_sk);
    std::cerr << "stage=tfhe_eval_keygen,bkfft_lvl02_done\n";
    tfhe_ek.emplaceiksk<Lvl10>(tfhe_sk);
    std::cerr << "stage=tfhe_eval_keygen,iksk_lvl10_done\n";
    tfhe_ek.emplaceiksk<Lvl20>(tfhe_sk);
    std::cerr << "stage=tfhe_eval_keygen,iksk_lvl20_done\n";
    tfhe_ek.emplaceiksk<Lvl21>(tfhe_sk);
    std::cerr << "stage=tfhe_eval_keygen,done\n";

    try {
        const CkksProfile ckks_profile = SelectCkksProfile(opts);
        std::cerr << "stage=ckks_context,start,profile="
                  << ProfileName(ckks_profile) << "\n";
        CkksEnv ckks(ckks_profile);
        std::cerr << "stage=ckks_context,done\n";
        auto repack_config =
            tfhepp_ckks::DefaultRepackConfig<Lvl1>(kRepackScaleBits, 45);
        tfhepp_ckks::RepackEvaluationKey repack_key;
        try {
            std::cerr << "stage=repack_keygen,start\n";
            tfhepp_ckks::GenerateRepackKey<Lvl1>(
                repack_key, tfhe_sk, repack_config.key_scale, ckks.encoder,
                ckks.symmetric_encryptor, ckks.context);
            std::cerr << "stage=repack_keygen,done\n";
        }
        catch (const std::exception &e) {
            throw std::runtime_error(
                std::string("GenerateRepackKey failed: ") + e.what());
        }

        if (WantsQuery(opts, "q6"))
            try {
                (void)EncryptedQ6(data, tfhe_sk, tfhe_ek, repack_key,
                                  repack_config, ckks);
            }
            catch (const std::exception &e) {
                throw std::runtime_error(std::string("EncryptedQ6 failed: ") +
                                         e.what());
            }
        if (WantsQuery(opts, "q14"))
            try {
                (void)EncryptedQ14(data, tfhe_sk, tfhe_ek, repack_key,
                                   repack_config, ckks);
            }
            catch (const std::exception &e) {
                throw std::runtime_error(std::string("EncryptedQ14 failed: ") +
                                         e.what());
            }
        if (WantsQuery(opts, "q3"))
            try {
                (void)EncryptedQ3(data, tfhe_sk, tfhe_ek, repack_key,
                                  repack_config, ckks);
            }
            catch (const std::exception &e) {
                throw std::runtime_error(std::string("EncryptedQ3 failed: ") +
                                         e.what());
            }
        if (WantsQuery(opts, "q5"))
            try {
                (void)EncryptedQ5(data, tfhe_sk, tfhe_ek, repack_key,
                                  repack_config, ckks);
            }
            catch (const std::exception &e) {
                throw std::runtime_error(std::string("EncryptedQ5 failed: ") +
                                         e.what());
            }
    }
    catch (const std::exception &e) {
        std::cerr << "encrypted pipeline failed: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
