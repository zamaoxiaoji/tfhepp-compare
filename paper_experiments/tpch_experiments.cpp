#include "tpch_experiments.hpp"

#include "algorithms.hpp"

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>

namespace PaperReview {

namespace {

using Clock = std::chrono::steady_clock;

struct br_lvl22param {
    using domainP = TFHEpp::lvl2param;
    using targetP = TFHEpp::lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

struct br_lvl11param : public TFHEpp::lvl11param {
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

using brP_meta = br_lvl22param;
using brP_logari = TFHEpp::lvl02param;
using brP_base = TFHEpp::lvl01param;
using brP_mask = br_lvl11param;
using iksP_t = TFHEpp::lvl20param;
using P_in = TFHEpp::lvl2param;
using P_out = TFHEpp::lvl1param;

template <class P>
typename P::T EncodeInteger(std::uint64_t message, int comparison_bits) {
    const int precision_bits = comparison_bits + 1;
    return static_cast<typename P::T>(message)
           << (std::numeric_limits<typename P::T>::digits - precision_bits);
}

std::uint64_t MaxComparableValue(int comparison_bits) {
    if (comparison_bits <= 1 ||
        comparison_bits >= std::numeric_limits<std::uint64_t>::digits)
        throw std::invalid_argument("comparison bit width is outside the supported range");
    return (std::uint64_t{1} << comparison_bits) - 1;
}

void RequireComparable(std::uint64_t value, int comparison_bits,
                       const std::string& name) {
    if (value > MaxComparableValue(comparison_bits))
        throw std::invalid_argument(
            name + " exceeds the configured [0, 2^bits-1] comparison range");
}

int DecodeComparisonBit(const TFHEpp::TLWE<P_out>& ct,
                        const TFHEpp::SecretKey& sk) {
    const auto phase = TFHEpp::tlweSymPhase<P_out>(ct, sk.key.get<P_out>());
    return static_cast<std::make_signed_t<typename P_out::T>>(phase) < 0 ? 1 : 0;
}

std::size_t NextPowerOfTwo(std::size_t x) {
    std::size_t out = 1;
    while (out < x) out <<= 1;
    return out;
}

template <typename Fn>
double TimeMs(Fn&& fn) {
    const auto start = Clock::now();
    fn();
    const auto stop = Clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
}

void AddTiming(
    std::vector<StageTiming>& timings,
    const std::string& category,
    const std::string& name,
    double milliseconds) {
    timings.push_back({category, name, milliseconds});
}

double SumTimings(
    const std::vector<StageTiming>& timings,
    const std::string& category) {
    double total = 0.0;
    for (const auto& t : timings) {
        if (t.category == category) total += t.milliseconds;
    }
    return total;
}

std::vector<double> AsDouble(const std::vector<std::uint32_t>& values, std::size_t slots) {
    std::vector<double> out(slots, 0.0);
    for (std::size_t i = 0; i < values.size(); i++)
        out[i] = static_cast<double>(values[i]);
    return out;
}

std::vector<double> PadDouble(const std::vector<double>& values, std::size_t slots) {
    std::vector<double> out(slots, 0.0);
    std::copy(values.begin(), values.end(), out.begin());
    return out;
}

double DecryptSlot0(const CkksRuntime& runtime, const CkksCiphertext& ct) {
    lbcrypto::Plaintext plain;
    runtime.cc->Decrypt(runtime.keys.secretKey, ct, &plain);
    plain->SetLength(runtime.slots);
    return plain->GetRealPackedValue().front();
}

std::vector<double> DecryptGroupSums(
    const CkksRuntime& runtime,
    const std::vector<CkksCiphertext>& sums) {
    std::vector<double> out;
    out.reserve(sums.size());
    for (const auto& ct : sums) out.push_back(DecryptSlot0(runtime, ct));
    return out;
}

double MaxAbsError(
    const std::vector<double>& expected,
    const std::vector<double>& actual) {
    if (expected.size() != actual.size())
        throw std::invalid_argument("group size mismatch");
    double err = 0.0;
    for (std::size_t i = 0; i < expected.size(); i++)
        err = std::max(err, std::abs(expected[i] - actual[i]));
    return err;
}

struct DomainArtifacts {
    std::vector<AttributeSpec> attrs;
    std::vector<DomainPoint> domain;
    std::vector<BasisFunction> basis;
    BasisMatrix matrix;
    AlphaTable alpha;
};

DomainArtifacts MakeDomain(const std::vector<std::size_t>& domain_sizes) {
    DomainArtifacts d;
    if (domain_sizes.empty())
        throw std::invalid_argument("domain must contain at least one attribute");
    d.attrs.reserve(domain_sizes.size());
    for (const auto size : domain_sizes) {
        if (size == 0)
            throw std::invalid_argument("domain attribute size must be nonzero");
        d.attrs.push_back(AttributeSpec{size});
    }
    d.domain = EnumerateDomainTuples(d.attrs);
    d.basis = EnumerateMonomialBasis(d.attrs);
    d.matrix = BuildBasisMatrix(d.domain, d.basis);
    d.alpha = SolveCoeffTableForAllTargets(d.matrix);
    if (!VerifyAlphaTable(d.matrix, d.alpha))
        throw std::runtime_error("AlphaTable verification failed");
    return d;
}

DomainArtifacts MakeDomain(std::size_t domain_size) {
    return MakeDomain(std::vector<std::size_t>{domain_size});
}

std::vector<CkksCiphertext> OneHotMasks(
    const CkksRuntime& runtime,
    const CkksCiphertext& key_column,
    const DomainArtifacts& domain) {
    auto basis_ct = BuildEncryptedBasis(
        runtime.cc, runtime.keys.publicKey, domain.attrs,
        {key_column}, domain.basis, runtime.slots);
    std::vector<CkksCiphertext> masks;
    masks.reserve(domain.alpha.coefficients.size());
    for (const auto& alpha : domain.alpha.coefficients)
        masks.push_back(LinearCombination(runtime.cc, basis_ct, alpha));
    return masks;
}

CkksCiphertext EqualityMask(
    const CkksRuntime& runtime,
    const CkksCiphertext& lhs,
    const CkksCiphertext& rhs,
    const DomainArtifacts& domain) {
    const auto lhs_masks = OneHotMasks(runtime, lhs, domain);
    const auto rhs_masks = OneHotMasks(runtime, rhs, domain);
    auto result = runtime.cc->EvalMult(lhs_masks.front(), 0.0);
    for (std::size_t j = 0; j < lhs_masks.size(); j++)
        result = runtime.cc->EvalAdd(
            result, runtime.cc->EvalMult(lhs_masks[j], rhs_masks[j]));
    return result;
}

// HE3DB-style WHERE stage: encrypt predicate operands as TFHE TLWEs, evaluate
// each comparison with Chapter-3 HomCompare/GapMSB, convert the sign bit to a
// 29-bit arithmetic TLWE mask, then repack that mask into CKKS slots.
class WhereEvaluator {
public:
    WhereEvaluator(CkksRuntime& runtime, const ExperimentOptions& options,
                   std::vector<StageTiming>& timings)
        : runtime_(runtime), options_(options), timings_(timings) {
        hom_options_.kappa = options.compare_kappa;
    }

    CkksCiphertext CompareConstant(
        const std::vector<std::uint32_t>& values,
        std::uint32_t constant,
        int comparison_bits,
        ComparePredicate predicate,
        const std::string& name,
        std::size_t* error_count = nullptr) {
        if (values.size() > runtime_.slots)
            throw std::invalid_argument("comparison vector exceeds slot count");
        RequireComparable(constant, comparison_bits, name + " constant");
        for (const auto value : values)
            RequireComparable(value, comparison_bits, name + " value");

        EnsureTfheKeys();
        EnsureRepackContext();

        std::vector<TFHEpp::TLWE<P_out>> signs(runtime_.slots);
        std::vector<TFHEpp::TLWE<P_out>> arithmetic_masks(runtime_.slots);
        std::size_t local_errors = 0;
        const double compare_ms = TimeMs([&] {
            TFHEpp::TLWE<P_in> rhs_ct{};
            TFHEpp::tlweSymEncrypt<P_in>(
                rhs_ct, EncodeInteger<P_in>(constant, comparison_bits),
                P_in::α, sk_.key.get<P_in>());
            for (std::size_t i = 0; i < values.size(); i++) {
                TFHEpp::TLWE<P_in> lhs_ct{};
                TFHEpp::tlweSymEncrypt<P_in>(
                    lhs_ct, EncodeInteger<P_in>(values[i], comparison_bits),
                    P_in::α, sk_.key.get<P_in>());
                signs[i] = Chapter3HomCompare<
                    brP_meta, brP_logari, brP_base, iksP_t>(
                        lhs_ct, rhs_ct, comparison_bits, predicate,
                        *bk_meta_, trkeys_, cfg_, *bk_logari_, *bk_base_,
                        *iksk_, hom_options_, nullptr);
                const int got = DecodeComparisonBit(signs[i], sk_);
                const int expected = ExpectedPredicate(values[i], constant, predicate);
                if (got != expected) local_errors++;
            }
            for (std::size_t i = values.size(); i < runtime_.slots; i++) {
                TFHEpp::tlweSymEncrypt<P_out>(
                    signs[i], static_cast<typename P_out::T>(P_out::μ),
                    P_out::α, sk_.key.get<P_out>());
            }
        });
        AddTiming(timings_, "filter", "where_gapmsb_compare_" + name, compare_ms);
        if (error_count) *error_count += local_errors;

        const double mask_ms = TimeMs([&] {
            for (std::size_t i = 0; i < values.size(); i++)
                arithmetic_masks[i] = SignToArithmeticMask(signs[i]);
            for (std::size_t i = values.size(); i < runtime_.slots; i++) {
                TFHEpp::tlweSymEncrypt<P_out>(
                    arithmetic_masks[i], typename P_out::T(0),
                    P_out::α, sk_.key.get<P_out>());
            }
        });
        AddTiming(timings_, "aggregation", "where_sign_to_tfhe_mask_" + name, mask_ms);

        CkksCiphertext binary_mask;
        const double repack_ms = TimeMs([&] {
            binary_mask = RepackExecute<P_out>(*repack_context_, arithmetic_masks);
        });
        AddTiming(timings_, "aggregation", "where_repack_mask_" + name, repack_ms);
        return binary_mask;
    }

private:
    static int ExpectedPredicate(
        std::uint32_t lhs,
        std::uint32_t rhs,
        ComparePredicate predicate) {
        switch (predicate) {
        case ComparePredicate::LessThan:
            return lhs < rhs ? 1 : 0;
        case ComparePredicate::GreaterThan:
            return lhs > rhs ? 1 : 0;
        case ComparePredicate::LessEqual:
            return lhs <= rhs ? 1 : 0;
        case ComparePredicate::GreaterEqual:
            return lhs >= rhs ? 1 : 0;
        }
        throw std::invalid_argument("unknown predicate");
    }

    void EnsureTfheKeys() {
        if (tfhe_ready_) return;
        const double keygen_ms = TimeMs([&] {
            bk_meta_ = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_meta>>();
            bk_logari_ =
                std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_logari>>();
            bk_base_ = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_base>>();
            bk_mask_ = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_mask>>();
            iksk_ = std::make_unique<TFHEpp::KeySwitchingKey<iksP_t>>();

            TFHEpp::bkfftgen<brP_meta>(*bk_meta_, sk_);
            TFHEpp::bkfftgen<brP_logari>(*bk_logari_, sk_);
            TFHEpp::bkfftgen<brP_base>(*bk_base_, sk_);
            TFHEpp::bkfftgen<brP_mask>(*bk_mask_, sk_);
            TFHEpp::ikskgen<iksP_t>(*iksk_, sk_);
            trkeys_.push_back(MetaPBS2::GenerateTruncRepeatKey<P_in>(
                sk_.key.get<P_in>(), cfg_.rounds[0].beta));
            trkeys_.push_back(MetaPBS2::GenerateTruncRepeatKey<P_in>(
                sk_.key.get<P_in>(), cfg_.rounds[1].beta));
        });
        AddTiming(timings_, "setup", "where_gapmsb_keygen", keygen_ms);
        tfhe_ready_ = true;
    }

    TFHEpp::TLWE<P_out> SignToArithmeticMask(
        const TFHEpp::TLWE<P_out>& sign) const {
        const auto half = static_cast<typename P_out::T>(P_out::μ >> 1);
        TFHEpp::Polynomial<P_out> tv{};
        tv.fill(static_cast<typename P_out::T>(-half));

        TFHEpp::TLWE<P_out> out{};
        TFHEpp::GateBootstrappingTLWE2TLWE<brP_mask>(
            out, sign, *bk_mask_, tv);
        out[P_out::k * P_out::n] += half;
        return out;
    }

    void EnsureRepackContext() {
        if (repack_context_) return;
        double setup_ms = TimeMs([&] {
            repack_context_.emplace(RepackSetup<P_out>(
                runtime_.cc, runtime_.keys,
                sk_.key.get<P_out>(),
                static_cast<std::uint32_t>(runtime_.slots),
                kLogQLwe));
        });
        AddTiming(timings_, "setup", "repack_setup", setup_ms);
    }

    static constexpr std::uint32_t kLogQLwe = 28;

    CkksRuntime& runtime_;
    const ExperimentOptions& options_;
    std::vector<StageTiming>& timings_;
    TFHEpp::SecretKey sk_;
    MetaPBS2::Algorithm1Config cfg_ = Chapter3MetaPBSConfig();
    MetaPBS2::HomMSBOptions hom_options_{};
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_meta>> bk_meta_;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_logari>> bk_logari_;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_base>> bk_base_;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<brP_mask>> bk_mask_;
    std::unique_ptr<TFHEpp::KeySwitchingKey<iksP_t>> iksk_;
    std::vector<MetaPBS2::TruncRepeatKey<P_in>> trkeys_;
    bool tfhe_ready_ = false;
    std::optional<RepackContext<P_out>> repack_context_;
};

CkksCiphertext EncryptOneMask(
    const CkksRuntime& runtime,
    std::vector<StageTiming>& timings,
    const std::string& name) {
    CkksCiphertext one;
    AddTiming(timings, "aggregation", "where_constant_one_" + name, TimeMs([&] {
        one = EncryptConstant(runtime.cc, runtime.keys.publicKey, 1.0, runtime.slots);
    }));
    return one;
}

CkksCiphertext GreaterEqualViaStrict(
    WhereEvaluator& where,
    const CkksRuntime& runtime,
    std::vector<StageTiming>& timings,
    const std::vector<std::uint32_t>& values,
    std::uint32_t lower,
    int bits,
    const std::string& name,
    std::size_t& predicate_errors) {
    if (lower == 0) return EncryptOneMask(runtime, timings, name + "_ge_zero");
    return where.CompareConstant(
        values, lower - 1, bits, ComparePredicate::GreaterThan,
        name + "_gt_lower_minus_one", &predicate_errors);
}

CkksCiphertext LessEqualViaStrict(
    WhereEvaluator& where,
    const std::vector<std::uint32_t>& values,
    std::uint32_t upper,
    int bits,
    const std::string& name,
    std::size_t& predicate_errors) {
    RequireComparable(static_cast<std::uint64_t>(upper) + 1, bits, name + " upper+1");
    return where.CompareConstant(
        values, upper + 1, bits, ComparePredicate::LessThan,
        name + "_lt_upper_plus_one", &predicate_errors);
}

CkksCiphertext EqualViaStrictRange(
    WhereEvaluator& where,
    const CkksRuntime& runtime,
    std::vector<StageTiming>& timings,
    const std::vector<std::uint32_t>& values,
    std::uint32_t target,
    int bits,
    const std::string& name,
    std::size_t& predicate_errors) {
    std::vector<CkksCiphertext> masks;
    if (target > 0) {
        masks.push_back(where.CompareConstant(
            values, target - 1, bits, ComparePredicate::GreaterThan,
            name + "_gt_target_minus_one", &predicate_errors));
    }
    RequireComparable(static_cast<std::uint64_t>(target) + 1, bits, name + " target+1");
    masks.push_back(where.CompareConstant(
        values, target + 1, bits, ComparePredicate::LessThan,
        name + "_lt_target_plus_one", &predicate_errors));
    if (masks.empty()) return EncryptOneMask(runtime, timings, name + "_all");
    if (masks.size() == 1) return masks.front();
    return MultiplyBalanced(runtime.cc, masks);
}

CkksCiphertext EncryptColumn(
    const CkksRuntime& runtime,
    const std::vector<double>& values) {
    return EncryptVector(
        runtime.cc, runtime.keys.publicKey,
        PadDouble(values, runtime.slots), runtime.slots);
}

CkksCiphertext EncryptColumn(
    const CkksRuntime& runtime,
    const std::vector<std::uint32_t>& values) {
    return EncryptVector(
        runtime.cc, runtime.keys.publicKey,
        AsDouble(values, runtime.slots), runtime.slots);
}

std::vector<std::uint32_t> PadU32(
    const std::vector<std::uint32_t>& values,
    std::size_t slots) {
    std::vector<std::uint32_t> out(slots, 0);
    std::copy(values.begin(), values.end(), out.begin());
    return out;
}

std::vector<double> PadRevenue(
    const std::vector<double>& values,
    std::size_t slots) {
    return PadDouble(values, slots);
}

void AddParam(ExperimentResult& result, const std::string& name,
              std::uint64_t value) {
    result.parameters.push_back({name, std::to_string(value)});
}

void FinalizeHe3dbTiming(ExperimentResult& result) {
    result.filter_time_ms = SumTimings(result.timings, "filter");
    result.aggregation_time_ms = SumTimings(result.timings, "aggregation");
    result.total_query_time_ms = result.filter_time_ms + result.aggregation_time_ms;
}

std::uint32_t CkksDepth(
    const ExperimentOptions& options,
    std::uint32_t default_depth) {
    return options.ckks_depth_override == 0
               ? default_depth
               : options.ckks_depth_override;
}

}  // namespace

ExperimentOptions ParseExperimentOptions(int argc, char** argv) {
    ExperimentOptions options;
    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        if (arg == "--rows" && i + 1 < argc) {
            options.rows = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--row-exp" && i + 1 < argc) {
            const auto exp = std::stoul(argv[++i]);
            if (exp >= std::numeric_limits<std::size_t>::digits)
                throw std::invalid_argument("row exponent is too large");
            options.rows = std::size_t{1} << exp;
        } else if (arg == "--seed" && i + 1 < argc) {
            options.seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        } else if ((arg == "--output" || arg == "--output-file") && i + 1 < argc) {
            options.output_path = argv[++i];
        } else if (arg == "--system" && i + 1 < argc) {
            options.system = argv[++i];
        } else if (arg == "--ckks-depth" && i + 1 < argc) {
            options.ckks_depth_override = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--kappa" && i + 1 < argc) {
            options.compare_kappa = std::stoi(argv[++i]);
        } else if (arg == "--date-bits" && i + 1 < argc) {
            options.date_bits = std::stoi(argv[++i]);
        } else if (arg == "--quantity-bits" && i + 1 < argc) {
            options.quantity_bits = std::stoi(argv[++i]);
        } else if (arg == "--discount-bits" && i + 1 < argc) {
            options.discount_bits = std::stoi(argv[++i]);
        } else if (arg == "--key-bits" && i + 1 < argc) {
            options.key_bits = std::stoi(argv[++i]);
        } else if (arg == "--q6-shipdate-data-min" && i + 1 < argc) {
            options.q6_shipdate_data_min = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q6-shipdate-data-max" && i + 1 < argc) {
            options.q6_shipdate_data_max = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q6-shipdate-lower" && i + 1 < argc) {
            options.q6_shipdate_lower = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q6-shipdate-upper" && i + 1 < argc) {
            options.q6_shipdate_upper = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q6-discount-lower" && i + 1 < argc) {
            options.q6_discount_lower = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q6-discount-upper" && i + 1 < argc) {
            options.q6_discount_upper = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q6-quantity-upper" && i + 1 < argc) {
            options.q6_quantity_upper = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q14-shipdate-lower" && i + 1 < argc) {
            options.q14_shipdate_lower = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q14-shipdate-upper" && i + 1 < argc) {
            options.q14_shipdate_upper = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q14-shipdate-data-min" && i + 1 < argc) {
            options.q14_shipdate_data_min = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q14-shipdate-data-max" && i + 1 < argc) {
            options.q14_shipdate_data_max = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q14-type-domain" && i + 1 < argc) {
            options.q14_type_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q14-promo-type" && i + 1 < argc) {
            options.q14_promo_type = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-key-domain" && i + 1 < argc) {
            options.q3_key_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q3-priority-domain" && i + 1 < argc) {
            options.q3_priority_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q3-segment-domain" && i + 1 < argc) {
            options.q3_segment_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q3-segment" && i + 1 < argc) {
            options.q3_customer_segment = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-orderdate-base" && i + 1 < argc) {
            options.q3_orderdate_base = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-orderdate-step" && i + 1 < argc) {
            options.q3_orderdate_step = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-orderdate-upper" && i + 1 < argc) {
            options.q3_orderdate_upper = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-shipdate-data-min" && i + 1 < argc) {
            options.q3_shipdate_data_min = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-shipdate-data-max" && i + 1 < argc) {
            options.q3_shipdate_data_max = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q3-shipdate-lower" && i + 1 < argc) {
            options.q3_shipdate_lower = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q5-key-domain" && i + 1 < argc) {
            options.q5_key_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q5-nation-domain" && i + 1 < argc) {
            options.q5_nation_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q5-region-domain" && i + 1 < argc) {
            options.q5_region_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
        } else if (arg == "--q5-region" && i + 1 < argc) {
            options.q5_region = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q5-orderdate-base" && i + 1 < argc) {
            options.q5_orderdate_base = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q5-orderdate-step" && i + 1 < argc) {
            options.q5_orderdate_step = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q5-orderdate-lower" && i + 1 < argc) {
            options.q5_orderdate_lower = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--q5-orderdate-upper" && i + 1 < argc) {
            options.q5_orderdate_upper = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else {
            throw std::invalid_argument(
                "usage: --rows N|--row-exp E --seed S "
                "[--system ours|he3db|arcedb|all] [--output PATH] [query parameters]");
        }
    }
    if (options.rows == 0) throw std::invalid_argument("rows must be nonzero");
    if (options.compare_kappa <= 0)
        throw std::invalid_argument("kappa must be positive");
    return options;
}

void PrintExperimentResult(const ExperimentResult& result, std::ostream& os) {
    os << "=== " << result.query_name << " ===\n";
    os << "rows=" << result.rows << " slots=" << result.slots
       << " mask_mode=" << result.mask_mode << "\n";
    os << "system=" << result.system << "\n";
    if (!result.parameters.empty()) {
        os << "parameter,value\n";
        for (const auto& [name, value] : result.parameters)
            os << name << "," << value << "\n";
    }
    os << "predicate_errors=" << result.predicate_errors << "\n";
    os << "he3db_timing_policy=filter excludes keygen/setup/repack; aggregation includes TFHE-to-CKKS repack and CKKS online aggregation\n";
    os << "filter_time_ms=" << result.filter_time_ms << "\n";
    os << "aggregation_time_ms=" << result.aggregation_time_ms << "\n";
    os << "total_query_time_ms=" << result.total_query_time_ms << "\n";
    os << "total_query_time_s=" << result.total_query_time_ms / 1000.0 << "\n";
    os << std::fixed << std::setprecision(6);
    if (!result.plain_groups.empty()) {
        os << "group,plain,encrypted,abs_error\n";
        for (std::size_t i = 0; i < result.plain_groups.size(); i++) {
            const double err =
                std::abs(result.plain_groups[i] - result.encrypted_groups[i]);
            os << i << "," << result.plain_groups[i] << ","
               << result.encrypted_groups[i] << "," << err << "\n";
        }
        os << "max_abs_error=" << result.abs_error << "\n";
        os << "scalar_plain=" << result.plain_scalar
           << " scalar_encrypted=" << result.encrypted_scalar
           << " scalar_abs_error="
           << std::abs(result.plain_scalar - result.encrypted_scalar) << "\n";
    } else {
        os << "plain=" << result.plain_scalar
           << " encrypted=" << result.encrypted_scalar
           << " abs_error=" << result.abs_error << "\n";
    }
    os << "stage_category,stage,ms\n";
    double total = 0.0;
    for (const auto& t : result.timings) {
        os << t.category << "," << t.name << "," << t.milliseconds << "\n";
        total += t.milliseconds;
    }
    os << "all,total_measured_ms," << total << "\n";
}

void WriteExperimentResultIfRequested(
    const ExperimentResult& result,
    const ExperimentOptions& options) {
    if (options.output_path.empty()) return;
    std::ofstream out(options.output_path);
    if (!out)
        throw std::runtime_error("failed to open output file: " + options.output_path);
    PrintExperimentResult(result, out);
}

ExperimentResult RunTpchQ6Experiment(const ExperimentOptions& options) {
    const std::size_t slots = NextPowerOfTwo(options.rows);
    ExperimentResult result;
    result.query_name = "TPC-H Q6 synthetic: TFHE GapMSB WHERE + CKKS SUM";
    result.rows = options.rows;
    result.slots = slots;
    result.mask_mode = "tfhe-gapmsb-compare-repacked";
    result.system = options.system;
    AddParam(result, "date_bits", options.date_bits);
    AddParam(result, "discount_bits", options.discount_bits);
    AddParam(result, "quantity_bits", options.quantity_bits);
    AddParam(result, "compare_kappa", options.compare_kappa);
    AddParam(result, "shipdate_lower", options.q6_shipdate_lower);
    AddParam(result, "shipdate_upper", options.q6_shipdate_upper);
    AddParam(result, "discount_lower", options.q6_discount_lower);
    AddParam(result, "discount_upper", options.q6_discount_upper);
    AddParam(result, "quantity_upper", options.q6_quantity_upper);
    if (options.q6_shipdate_lower >= options.q6_shipdate_upper)
        throw std::invalid_argument("q6 shipdate lower must be smaller than upper");
    if (options.q6_discount_lower > options.q6_discount_upper)
        throw std::invalid_argument("q6 discount lower must be no larger than upper");
    if (options.q6_quantity_upper == 0)
        throw std::invalid_argument("q6 quantity upper must be positive");
    if (options.q6_shipdate_lower < options.q6_shipdate_data_min ||
        options.q6_shipdate_upper < options.q6_shipdate_data_min ||
        options.q6_shipdate_data_max < options.q6_shipdate_data_min)
        throw std::invalid_argument("q6 shipdate offsets require min <= thresholds and max");
    const auto q6_shipdate_lower_code =
        options.q6_shipdate_lower - options.q6_shipdate_data_min;
    const auto q6_shipdate_upper_code =
        options.q6_shipdate_upper - options.q6_shipdate_data_min;
    const auto q6_shipdate_max_code =
        options.q6_shipdate_data_max - options.q6_shipdate_data_min;
    AddParam(result, "shipdate_offset_min", options.q6_shipdate_data_min);
    AddParam(result, "shipdate_lower_offset", q6_shipdate_lower_code);
    AddParam(result, "shipdate_upper_offset", q6_shipdate_upper_code);
    RequireComparable(q6_shipdate_max_code, options.date_bits, "q6 shipdate offset max");
    RequireComparable(q6_shipdate_lower_code, options.date_bits, "q6 shipdate lower offset");
    RequireComparable(q6_shipdate_upper_code, options.date_bits, "q6 shipdate upper offset");
    RequireComparable(options.q6_discount_upper, options.discount_bits, "q6 discount upper");
    RequireComparable(options.q6_quantity_upper, options.quantity_bits, "q6 quantity upper");

    std::mt19937_64 rng(options.seed);
    const auto quantity_max =
        static_cast<int>(std::min<std::uint64_t>(50, MaxComparableValue(options.quantity_bits)));
    const auto discount_max =
        static_cast<int>(std::min<std::uint64_t>(15, MaxComparableValue(options.discount_bits)));
    std::uniform_int_distribution<int> quantity_dist(1, quantity_max);
    std::uniform_int_distribution<int> discount_dist(0, discount_max);
    std::uniform_int_distribution<std::uint32_t> date_dist(
        options.q6_shipdate_data_min, options.q6_shipdate_data_max);
    std::uniform_real_distribution<double> price_dist(1.0, 100.0);

    std::vector<int> mask(options.rows);
    std::vector<std::uint32_t> quantity(options.rows), discount(options.rows);
    std::vector<std::uint32_t> shipdate(options.rows), shipdate_code(options.rows);
    std::vector<double> revenue(options.rows);
    double plain = 0.0;

    const double data_ms = TimeMs([&] {
        for (std::size_t i = 0; i < options.rows; i++) {
            quantity[i] = static_cast<std::uint32_t>(quantity_dist(rng));
            discount[i] = static_cast<std::uint32_t>(discount_dist(rng));
            shipdate[i] = date_dist(rng);
            shipdate_code[i] = shipdate[i] - options.q6_shipdate_data_min;
            const double price = price_dist(rng);
            const bool pass =
                shipdate[i] >= options.q6_shipdate_lower &&
                shipdate[i] < options.q6_shipdate_upper &&
                discount[i] >= options.q6_discount_lower &&
                discount[i] <= options.q6_discount_upper &&
                quantity[i] < options.q6_quantity_upper;
            mask[i] = pass ? 1 : 0;
            revenue[i] = price * static_cast<double>(discount[i]);
            plain += mask[i] ? revenue[i] : 0.0;
        }
        if (!mask.empty()) {
            shipdate[0] = options.q6_shipdate_lower +
                          (options.q6_shipdate_upper -
                           options.q6_shipdate_lower) /
                              2;
            shipdate_code[0] = shipdate[0] - options.q6_shipdate_data_min;
            discount[0] = options.q6_discount_lower +
                          (options.q6_discount_upper -
                           options.q6_discount_lower) /
                              2;
            quantity[0] = std::max<std::uint32_t>(
                1, options.q6_quantity_upper / 2);
            mask[0] = 1;
            revenue[0] = 90.0;
            plain = 0.0;
            for (std::size_t i = 0; i < options.rows; i++)
                plain += mask[i] ? revenue[i] : 0.0;
        }
    });
    AddTiming(result.timings, "setup", "generate_q6_data_and_plain_baseline", data_ms);

    CkksRuntime runtime = MakeCkksRuntime(
        slots,
        CkksDepth(options, 20),
        true);
    WhereEvaluator where(runtime, options, result.timings);

    CkksCiphertext ct_filter;
    std::vector<CkksCiphertext> q6_masks;
    q6_masks.push_back(GreaterEqualViaStrict(
        where, runtime, result.timings, shipdate_code,
        q6_shipdate_lower_code, options.date_bits, "q6_shipdate",
        result.predicate_errors));
    q6_masks.push_back(where.CompareConstant(
        shipdate_code, q6_shipdate_upper_code, options.date_bits,
        ComparePredicate::LessThan, "q6_shipdate_lt",
        &result.predicate_errors));
    q6_masks.push_back(GreaterEqualViaStrict(
        where, runtime, result.timings, discount,
        options.q6_discount_lower, options.discount_bits, "q6_discount",
        result.predicate_errors));
    q6_masks.push_back(LessEqualViaStrict(
        where, discount, options.q6_discount_upper, options.discount_bits,
        "q6_discount", result.predicate_errors));
    q6_masks.push_back(where.CompareConstant(
        quantity, options.q6_quantity_upper, options.quantity_bits,
        ComparePredicate::LessThan, "q6_quantity_lt",
        &result.predicate_errors));
    AddTiming(result.timings, "aggregation", "where_ckks_mask_product_q6", TimeMs([&] {
        ct_filter = MultiplyBalanced(runtime.cc, q6_masks);
    }));
    CkksCiphertext ct_revenue;
    AddTiming(result.timings, "setup", "encrypt_revenue_ckks", TimeMs([&] {
        ct_revenue = EncryptColumn(runtime, revenue);
    }));

    CkksCiphertext ct_sum;
    AddTiming(result.timings, "aggregation", "ckks_mask_mul_rotate_sum", TimeMs([&] {
        auto filtered = runtime.cc->EvalMult(ct_filter, ct_revenue);
        ct_sum = EvalRotateAndSum(runtime.cc, filtered, runtime.slots);
    }));

    result.plain_scalar = plain;
    result.encrypted_scalar = DecryptSlot0(runtime, ct_sum);
    result.abs_error = std::abs(result.plain_scalar - result.encrypted_scalar);
    FinalizeHe3dbTiming(result);
    return result;
}

ExperimentResult RunTpchQ14Experiment(const ExperimentOptions& options) {
    const std::size_t type_domain = options.q14_type_domain;
    const std::size_t slots = NextPowerOfTwo(std::max(options.rows, type_domain));
    ExperimentResult result;
    result.query_name = "TPC-H Q14 synthetic: TFHE GapMSB WHERE + CKKS promo GROUP BY";
    result.rows = options.rows;
    result.slots = slots;
    result.mask_mode = "tfhe-gapmsb-compare-repacked";
    result.system = options.system;
    AddParam(result, "date_bits", options.date_bits);
    AddParam(result, "compare_kappa", options.compare_kappa);
    AddParam(result, "shipdate_lower", options.q14_shipdate_lower);
    AddParam(result, "shipdate_upper", options.q14_shipdate_upper);
    AddParam(result, "shipdate_data_min", options.q14_shipdate_data_min);
    AddParam(result, "shipdate_data_max", options.q14_shipdate_data_max);
    AddParam(result, "type_domain", type_domain);
    AddParam(result, "promo_type", options.q14_promo_type);
    if (type_domain == 0) throw std::invalid_argument("q14 type domain must be nonzero");
    if (options.q14_promo_type >= type_domain)
        throw std::invalid_argument("q14 promo type must be inside the type domain");
    if (options.q14_shipdate_data_min > options.q14_shipdate_data_max)
        throw std::invalid_argument("q14 shipdate data range is invalid");
    if (options.q14_shipdate_lower >= options.q14_shipdate_upper)
        throw std::invalid_argument("q14 shipdate lower must be smaller than upper");
    RequireComparable(options.q14_shipdate_lower, options.date_bits, "q14 shipdate lower");
    RequireComparable(options.q14_shipdate_upper, options.date_bits, "q14 shipdate upper");
    RequireComparable(options.q14_shipdate_data_max, options.date_bits, "q14 shipdate data max");

    std::mt19937_64 rng(options.seed);
    std::uniform_int_distribution<int> type_dist(0, static_cast<int>(type_domain - 1));
    std::uniform_int_distribution<std::uint32_t> date_dist(
        options.q14_shipdate_data_min, options.q14_shipdate_data_max);
    std::uniform_real_distribution<double> price_dist(1.0, 100.0);
    std::uniform_real_distribution<double> discount_dist(0.0, 0.10);

    std::vector<int> mask(options.rows);
    std::vector<std::uint32_t> shipdate(options.rows);
    std::vector<std::uint32_t> part_type(options.rows);
    std::vector<double> revenue(options.rows);
    std::vector<double> plain_groups(type_domain, 0.0);

    AddTiming(result.timings, "setup", "generate_q14_data_and_plain_baseline", TimeMs([&] {
        for (std::size_t i = 0; i < options.rows; i++) {
            shipdate[i] = date_dist(rng);
            part_type[i] = static_cast<std::uint32_t>(type_dist(rng));
            revenue[i] = price_dist(rng) * (1.0 - discount_dist(rng));
            mask[i] = shipdate[i] >= options.q14_shipdate_lower &&
                              shipdate[i] < options.q14_shipdate_upper
                          ? 1
                          : 0;
            if (mask[i]) plain_groups[part_type[i]] += revenue[i];
        }
        if (!mask.empty()) {
            shipdate[0] = options.q14_shipdate_lower +
                          (options.q14_shipdate_upper -
                           options.q14_shipdate_lower) /
                              2;
            mask[0] = 1;
            part_type[0] = options.q14_promo_type;
            revenue[0] = 50.0;
            std::fill(plain_groups.begin(), plain_groups.end(), 0.0);
            for (std::size_t i = 0; i < options.rows; i++)
                if (mask[i]) plain_groups[part_type[i]] += revenue[i];
        }
    }));

    CkksRuntime runtime = MakeCkksRuntime(
        slots,
        CkksDepth(options, 24),
        true);
    WhereEvaluator where(runtime, options, result.timings);
    const auto domain = MakeDomain(type_domain);

    CkksCiphertext ct_mask;
    std::vector<CkksCiphertext> q14_masks;
    q14_masks.push_back(GreaterEqualViaStrict(
        where, runtime, result.timings, shipdate,
        options.q14_shipdate_lower, options.date_bits, "q14_shipdate",
        result.predicate_errors));
    q14_masks.push_back(where.CompareConstant(
        shipdate, options.q14_shipdate_upper, options.date_bits,
        ComparePredicate::LessThan, "q14_shipdate_lt",
        &result.predicate_errors));
    AddTiming(result.timings, "aggregation", "where_ckks_mask_product_q14", TimeMs([&] {
        ct_mask = MultiplyBalanced(runtime.cc, q14_masks);
    }));
    CkksCiphertext ct_type, ct_revenue, ct_masked_revenue;
    AddTiming(result.timings, "setup", "encrypt_type_and_revenue_ckks", TimeMs([&] {
        ct_type = EncryptColumn(runtime, part_type);
        ct_revenue = EncryptColumn(runtime, revenue);
    }));
    AddTiming(result.timings, "aggregation", "apply_filter_mask", TimeMs([&] {
        ct_masked_revenue = runtime.cc->EvalMult(ct_mask, ct_revenue);
    }));

    std::vector<CkksCiphertext> group_sums;
    AddTiming(result.timings, "aggregation", "matrix_group_by_sum", TimeMs([&] {
        group_sums = MatrixGroupBySum(
            runtime.cc, runtime.keys.publicKey, domain.attrs, {ct_type},
            ct_masked_revenue, domain.basis, domain.alpha, runtime.slots);
    }));

    result.plain_groups = plain_groups;
    result.encrypted_groups = DecryptGroupSums(runtime, group_sums);
    result.abs_error = MaxAbsError(result.plain_groups, result.encrypted_groups);
    const double plain_total = std::accumulate(plain_groups.begin(), plain_groups.end(), 0.0);
    const double enc_total = std::accumulate(
        result.encrypted_groups.begin(), result.encrypted_groups.end(), 0.0);
    const auto promo = static_cast<std::size_t>(options.q14_promo_type);
    result.plain_scalar =
        plain_total == 0.0 ? 0.0 : 100.0 * plain_groups[promo] / plain_total;
    result.encrypted_scalar =
        enc_total == 0.0 ? 0.0 : 100.0 * result.encrypted_groups[promo] / enc_total;
    FinalizeHe3dbTiming(result);
    return result;
}

ExperimentResult RunTpchQ3Experiment(const ExperimentOptions& options) {
    const std::size_t key_domain = options.q3_key_domain;
    const std::size_t priority_domain = options.q3_priority_domain;
    const std::size_t slots = NextPowerOfTwo(std::max(options.rows, key_domain));
    ExperimentResult result;
    result.query_name = "TPC-H Q3 synthetic: TFHE GapMSB WHERE + CKKS joins + order/priority GROUP BY";
    result.rows = options.rows;
    result.slots = slots;
    result.mask_mode = "tfhe-gapmsb-compare-repacked";
    result.system = options.system;
    AddParam(result, "date_bits", options.date_bits);
    AddParam(result, "key_bits", options.key_bits);
    AddParam(result, "compare_kappa", options.compare_kappa);
    AddParam(result, "key_domain", key_domain);
    AddParam(result, "priority_domain", priority_domain);
    AddParam(result, "segment_domain", options.q3_segment_domain);
    AddParam(result, "customer_segment", options.q3_customer_segment);
    AddParam(result, "orderdate_base", options.q3_orderdate_base);
    AddParam(result, "orderdate_step", options.q3_orderdate_step);
    AddParam(result, "orderdate_upper", options.q3_orderdate_upper);
    AddParam(result, "shipdate_data_min", options.q3_shipdate_data_min);
    AddParam(result, "shipdate_data_max", options.q3_shipdate_data_max);
    AddParam(result, "shipdate_lower", options.q3_shipdate_lower);
    if (key_domain == 0) throw std::invalid_argument("q3 key domain must be nonzero");
    if (priority_domain == 0) throw std::invalid_argument("q3 priority domain must be nonzero");
    if (options.q3_segment_domain == 0)
        throw std::invalid_argument("q3 segment domain must be nonzero");
    RequireComparable(key_domain - 1, options.key_bits, "q3 key domain");
    RequireComparable(priority_domain - 1, options.key_bits, "q3 priority domain");
    RequireComparable(options.q3_segment_domain - 1, options.key_bits, "q3 segment domain");
    RequireComparable(options.q3_customer_segment, options.key_bits, "q3 segment");
    if (options.q3_customer_segment >= options.q3_segment_domain)
        throw std::invalid_argument("q3 segment must be inside the segment domain");
    RequireComparable(options.q3_orderdate_upper, options.date_bits, "q3 orderdate upper");
    RequireComparable(options.q3_orderdate_base +
                          options.q3_orderdate_step *
                              static_cast<std::uint32_t>(key_domain - 1),
                      options.date_bits, "q3 generated orderdate max");
    if (options.q3_shipdate_data_min > options.q3_shipdate_data_max)
        throw std::invalid_argument("q3 shipdate data range is invalid");
    if (options.q3_shipdate_lower >= options.q3_shipdate_data_max)
        throw std::invalid_argument("q3 shipdate lower must leave at least one passing generated date");
    RequireComparable(options.q3_shipdate_data_max, options.date_bits, "q3 shipdate data max");
    RequireComparable(options.q3_shipdate_lower, options.date_bits, "q3 shipdate lower");

    std::mt19937_64 rng(options.seed);
    std::uniform_int_distribution<int> key_dist(0, static_cast<int>(key_domain - 1));
    std::uniform_int_distribution<std::uint32_t> date_dist(
        options.q3_shipdate_data_min, options.q3_shipdate_data_max);
    std::uniform_real_distribution<double> price_dist(1.0, 50.0);
    std::uniform_real_distribution<double> discount_dist(0.0, 0.10);

    std::vector<std::uint32_t> customer_key(slots), customer_seg(slots);
    std::vector<std::uint32_t> order_key(slots), order_customer_key(slots), order_date(slots);
    std::vector<std::uint32_t> order_priority(slots);
    std::vector<std::uint32_t> line_order_key(options.rows), line_shipdate(options.rows);
    std::vector<double> revenue(options.rows);
    std::vector<int> customer_segment_mask(slots), order_date_mask(slots), line_ship_mask(options.rows);
    std::vector<double> plain_groups(key_domain * priority_domain, 0.0);

    AddTiming(result.timings, "setup", "generate_q3_data_and_plain_baseline", TimeMs([&] {
        for (std::size_t j = 0; j < key_domain; j++) {
            customer_key[j] = static_cast<std::uint32_t>(j);
            customer_seg[j] =
                static_cast<std::uint32_t>(j % options.q3_segment_domain);
            customer_segment_mask[j] =
                customer_seg[j] == options.q3_customer_segment ? 1 : 0;
            order_key[j] = static_cast<std::uint32_t>(j);
            order_customer_key[j] = static_cast<std::uint32_t>(j % key_domain);
            order_priority[j] = static_cast<std::uint32_t>(j % priority_domain);
            order_date[j] = options.q3_orderdate_base +
                            options.q3_orderdate_step *
                                static_cast<std::uint32_t>(j);
            order_date_mask[j] =
                order_date[j] < options.q3_orderdate_upper ? 1 : 0;
        }
        for (std::size_t j = key_domain; j < slots; j++) {
            customer_key[j] = 0;
            customer_seg[j] = options.q3_customer_segment == 0 ? 1 : 0;
            customer_segment_mask[j] = 0;
            order_key[j] = 0;
            order_customer_key[j] = 0;
            order_priority[j] = 0;
            order_date[j] = options.q3_orderdate_upper;
            order_date_mask[j] = 0;
        }
        for (std::size_t i = 0; i < options.rows; i++) {
            line_order_key[i] = static_cast<std::uint32_t>(key_dist(rng));
            line_shipdate[i] = date_dist(rng);
            revenue[i] = price_dist(rng) * (1.0 - discount_dist(rng));
            line_ship_mask[i] =
                line_shipdate[i] > options.q3_shipdate_lower ? 1 : 0;
            const auto ok = line_order_key[i];
            const auto ck = order_customer_key[ok];
            const auto pr = order_priority[ok];
            if (customer_segment_mask[ck] && order_date_mask[ok] && line_ship_mask[i])
                plain_groups[ok + key_domain * pr] += revenue[i];
        }
        if (!line_order_key.empty()) {
            line_order_key[0] = 0;
            line_shipdate[0] = options.q3_shipdate_lower +
                               std::max<std::uint32_t>(
                                   1,
                                   (options.q3_shipdate_data_max -
                                    options.q3_shipdate_lower) /
                                       2);
            line_ship_mask[0] = 1;
            revenue[0] = 30.0;
            std::fill(plain_groups.begin(), plain_groups.end(), 0.0);
            for (std::size_t i = 0; i < options.rows; i++) {
                const auto ok = line_order_key[i];
                const auto ck = order_customer_key[ok];
                const auto pr = order_priority[ok];
                if (customer_segment_mask[ck] && order_date_mask[ok] && line_ship_mask[i])
                    plain_groups[ok + key_domain * pr] += revenue[i];
            }
        }
    }));

    CkksRuntime runtime = MakeCkksRuntime(
        slots,
        CkksDepth(options, 32),
        true);
    WhereEvaluator where(runtime, options, result.timings);
    const auto domain = MakeDomain(key_domain);
    const auto group_domain = MakeDomain({key_domain, priority_domain});

    CkksCiphertext ct_customer_key, ct_order_customer_key, ct_order_key, ct_order_priority;
    CkksCiphertext ct_line_order_key, ct_revenue;
    AddTiming(result.timings, "setup", "encrypt_join_and_value_columns_ckks", TimeMs([&] {
        ct_customer_key = EncryptColumn(runtime, customer_key);
        ct_order_customer_key = EncryptColumn(runtime, order_customer_key);
        ct_order_key = EncryptColumn(runtime, order_key);
        ct_order_priority = EncryptColumn(runtime, order_priority);
        ct_line_order_key = EncryptColumn(runtime, PadU32(line_order_key, slots));
        ct_revenue = EncryptColumn(runtime, PadRevenue(revenue, slots));
    }));

    CkksCiphertext ct_customer_segment, ct_order_date, ct_line_ship;
    ct_customer_segment = EqualViaStrictRange(
        where, runtime, result.timings, customer_seg,
        options.q3_customer_segment, options.key_bits,
        "q3_customer_segment", result.predicate_errors);
    ct_order_date = where.CompareConstant(
        order_date, options.q3_orderdate_upper, options.date_bits,
        ComparePredicate::LessThan, "q3_orderdate_lt",
        &result.predicate_errors);
    ct_line_ship = where.CompareConstant(
        line_shipdate, options.q3_shipdate_lower, options.date_bits,
        ComparePredicate::GreaterThan, "q3_line_shipdate_gt",
        &result.predicate_errors);

    std::vector<CkksCiphertext> group_sums;
    AddTiming(result.timings, "aggregation", "ckks_lookup_join_join_groupby", TimeMs([&] {
        auto order_segment = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_order_customer_key, {},
            ct_customer_key, {ct_customer_segment}, domain.basis, domain.alpha,
            runtime.slots).front();
        auto order_pass = runtime.cc->EvalMult(order_segment, ct_order_date);
        auto line_order_pass = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_line_order_key, {},
            ct_order_key, {order_pass}, domain.basis, domain.alpha,
            runtime.slots).front();
        auto line_order_priority = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_line_order_key, {},
            ct_order_key, {ct_order_priority}, domain.basis, domain.alpha,
            runtime.slots).front();
        auto combined = runtime.cc->EvalMult(line_order_pass, ct_line_ship);
        auto filtered_revenue = runtime.cc->EvalMult(combined, ct_revenue);
        group_sums = MatrixGroupBySum(
            runtime.cc, runtime.keys.publicKey, group_domain.attrs,
            {ct_line_order_key, line_order_priority}, filtered_revenue,
            group_domain.basis, group_domain.alpha, runtime.slots);
    }));

    result.plain_groups = plain_groups;
    result.encrypted_groups = DecryptGroupSums(runtime, group_sums);
    result.abs_error = MaxAbsError(result.plain_groups, result.encrypted_groups);
    result.plain_scalar = std::accumulate(plain_groups.begin(), plain_groups.end(), 0.0);
    result.encrypted_scalar = std::accumulate(
        result.encrypted_groups.begin(), result.encrypted_groups.end(), 0.0);
    FinalizeHe3dbTiming(result);
    return result;
}

ExperimentResult RunTpchQ5Experiment(const ExperimentOptions& options) {
    const std::size_t key_domain = options.q5_key_domain;
    const std::size_t nation_domain = options.q5_nation_domain;
    const std::size_t slots =
        NextPowerOfTwo(std::max({options.rows, key_domain, nation_domain}));
    ExperimentResult result;
    result.query_name = "TPC-H Q5 synthetic: TFHE GapMSB WHERE + CKKS multijoin + nation GROUP BY";
    result.rows = options.rows;
    result.slots = slots;
    result.mask_mode = "tfhe-gapmsb-compare-repacked";
    result.system = options.system;
    AddParam(result, "date_bits", options.date_bits);
    AddParam(result, "key_bits", options.key_bits);
    AddParam(result, "compare_kappa", options.compare_kappa);
    AddParam(result, "key_domain", key_domain);
    AddParam(result, "nation_domain", nation_domain);
    AddParam(result, "region_domain", options.q5_region_domain);
    AddParam(result, "region", options.q5_region);
    AddParam(result, "orderdate_base", options.q5_orderdate_base);
    AddParam(result, "orderdate_step", options.q5_orderdate_step);
    AddParam(result, "orderdate_lower", options.q5_orderdate_lower);
    AddParam(result, "orderdate_upper", options.q5_orderdate_upper);
    if (key_domain == 0) throw std::invalid_argument("q5 key domain must be nonzero");
    if (nation_domain == 0) throw std::invalid_argument("q5 nation domain must be nonzero");
    if (options.q5_region_domain == 0) throw std::invalid_argument("q5 region domain must be nonzero");
    RequireComparable(key_domain - 1, options.key_bits, "q5 key domain");
    RequireComparable(nation_domain - 1, options.key_bits, "q5 nation domain");
    RequireComparable(options.q5_region_domain - 1, options.key_bits, "q5 region domain");
    RequireComparable(options.q5_region, options.key_bits, "q5 region");
    if (options.q5_region >= options.q5_region_domain)
        throw std::invalid_argument("q5 region must be inside the region domain");
    RequireComparable(options.q5_orderdate_base +
                          options.q5_orderdate_step *
                              static_cast<std::uint32_t>(key_domain - 1),
                      options.date_bits, "q5 generated orderdate max");
    RequireComparable(options.q5_orderdate_lower, options.date_bits, "q5 orderdate lower");
    RequireComparable(options.q5_orderdate_upper, options.date_bits, "q5 orderdate upper");
    if (options.q5_orderdate_lower >= options.q5_orderdate_upper)
        throw std::invalid_argument("q5 orderdate lower must be smaller than upper");

    std::mt19937_64 rng(options.seed);
    std::uniform_int_distribution<int> key_dist(0, static_cast<int>(key_domain - 1));
    std::uniform_real_distribution<double> price_dist(1.0, 50.0);
    std::uniform_real_distribution<double> discount_dist(0.0, 0.10);

    std::vector<std::uint32_t> nation_key(slots), nation_region(slots);
    std::vector<std::uint32_t> customer_key(slots), customer_nation(slots);
    std::vector<std::uint32_t> supplier_key(slots), supplier_nation(slots);
    std::vector<std::uint32_t> order_key(slots), order_customer_key(slots), order_date(slots);
    std::vector<std::uint32_t> line_order_key(options.rows), line_supp_key(options.rows);
    std::vector<double> revenue(options.rows);
    std::vector<int> nation_asia_mask(slots), order_date_mask(slots);
    std::vector<double> plain_groups(nation_domain, 0.0);

    AddTiming(result.timings, "setup", "generate_q5_data_and_plain_baseline", TimeMs([&] {
        for (std::size_t j = 0; j < nation_domain; j++) {
            nation_key[j] = static_cast<std::uint32_t>(j);
            nation_region[j] = static_cast<std::uint32_t>(j % options.q5_region_domain);
            nation_asia_mask[j] = nation_region[j] == options.q5_region ? 1 : 0;
        }
        for (std::size_t j = nation_domain; j < slots; j++) {
            nation_key[j] = 0;
            nation_region[j] = options.q5_region == 0 ? 1 : 0;
            nation_asia_mask[j] = 0;
        }
        for (std::size_t j = 0; j < key_domain; j++) {
            customer_key[j] = static_cast<std::uint32_t>(j);
            customer_nation[j] = static_cast<std::uint32_t>(j % nation_domain);
            supplier_key[j] = static_cast<std::uint32_t>(j);
            supplier_nation[j] = static_cast<std::uint32_t>(j % nation_domain);
            order_key[j] = static_cast<std::uint32_t>(j);
            order_customer_key[j] = static_cast<std::uint32_t>(j % key_domain);
            order_date[j] = options.q5_orderdate_base +
                            options.q5_orderdate_step *
                                static_cast<std::uint32_t>(j);
            order_date_mask[j] =
                order_date[j] >= options.q5_orderdate_lower &&
                        order_date[j] < options.q5_orderdate_upper
                    ? 1
                    : 0;
        }
        for (std::size_t j = key_domain; j < slots; j++) {
            customer_key[j] = 0;
            customer_nation[j] = 0;
            supplier_key[j] = 0;
            supplier_nation[j] = 0;
            order_key[j] = 0;
            order_customer_key[j] = 0;
            order_date[j] = options.q5_orderdate_upper;
            order_date_mask[j] = 0;
        }
        for (std::size_t i = 0; i < options.rows; i++) {
            line_order_key[i] = static_cast<std::uint32_t>(key_dist(rng));
            line_supp_key[i] = static_cast<std::uint32_t>(key_dist(rng));
            revenue[i] = price_dist(rng) * (1.0 - discount_dist(rng));
            const auto ok = line_order_key[i];
            const auto sk = line_supp_key[i];
            const auto cn = customer_nation[order_customer_key[ok]];
            const auto sn = supplier_nation[sk];
            if (order_date_mask[ok] && cn == sn && nation_asia_mask[sn])
                plain_groups[sn] += revenue[i];
        }
        if (!line_order_key.empty()) {
            const auto forced_key =
                static_cast<std::uint32_t>(std::min<std::size_t>(2, key_domain - 1));
            line_order_key[0] = forced_key;
            line_supp_key[0] = forced_key;
            order_date[forced_key] = options.q5_orderdate_lower +
                                     (options.q5_orderdate_upper -
                                      options.q5_orderdate_lower) /
                                         2;
            order_date_mask[forced_key] = 1;
            customer_nation[order_customer_key[forced_key]] =
                static_cast<std::uint32_t>(forced_key % nation_domain);
            supplier_nation[forced_key] =
                static_cast<std::uint32_t>(forced_key % nation_domain);
            nation_region[forced_key % nation_domain] = options.q5_region;
            nation_asia_mask[forced_key % nation_domain] = 1;
            revenue[0] = 40.0;
            std::fill(plain_groups.begin(), plain_groups.end(), 0.0);
            for (std::size_t i = 0; i < options.rows; i++) {
                const auto ok = line_order_key[i];
                const auto sk = line_supp_key[i];
                const auto cn = customer_nation[order_customer_key[ok]];
                const auto sn = supplier_nation[sk];
                if (order_date_mask[ok] && cn == sn && nation_asia_mask[sn])
                    plain_groups[sn] += revenue[i];
            }
        }
    }));

    CkksRuntime runtime = MakeCkksRuntime(
        slots,
        CkksDepth(options, 36),
        true);
    WhereEvaluator where(runtime, options, result.timings);
    const auto key_dom = MakeDomain(key_domain);
    const auto nation_dom = MakeDomain(nation_domain);

    CkksCiphertext ct_nation_key, ct_customer_key, ct_customer_nation;
    CkksCiphertext ct_supplier_key, ct_supplier_nation;
    CkksCiphertext ct_order_key, ct_order_customer_key;
    CkksCiphertext ct_line_order_key, ct_line_supp_key, ct_revenue;
    AddTiming(result.timings, "setup", "encrypt_six_table_columns_ckks", TimeMs([&] {
        ct_nation_key = EncryptColumn(runtime, nation_key);
        ct_customer_key = EncryptColumn(runtime, customer_key);
        ct_customer_nation = EncryptColumn(runtime, customer_nation);
        ct_supplier_key = EncryptColumn(runtime, supplier_key);
        ct_supplier_nation = EncryptColumn(runtime, supplier_nation);
        ct_order_key = EncryptColumn(runtime, order_key);
        ct_order_customer_key = EncryptColumn(runtime, order_customer_key);
        ct_line_order_key = EncryptColumn(runtime, PadU32(line_order_key, slots));
        ct_line_supp_key = EncryptColumn(runtime, PadU32(line_supp_key, slots));
        ct_revenue = EncryptColumn(runtime, PadRevenue(revenue, slots));
    }));

    CkksCiphertext ct_nation_asia, ct_order_date;
    ct_nation_asia = EqualViaStrictRange(
        where, runtime, result.timings, nation_region, options.q5_region,
        options.key_bits, "q5_region", result.predicate_errors);
    auto order_date_ge = GreaterEqualViaStrict(
        where, runtime, result.timings, order_date,
        options.q5_orderdate_lower, options.date_bits, "q5_orderdate",
        result.predicate_errors);
    auto order_date_lt = where.CompareConstant(
        order_date, options.q5_orderdate_upper, options.date_bits,
        ComparePredicate::LessThan, "q5_orderdate_lt",
        &result.predicate_errors);
    AddTiming(result.timings, "aggregation", "where_ckks_mask_product_q5", TimeMs([&] {
        ct_order_date = runtime.cc->EvalMult(order_date_ge, order_date_lt);
    }));

    std::vector<CkksCiphertext> group_sums;
    AddTiming(result.timings, "aggregation", "ckks_multijoin_groupby", TimeMs([&] {
        auto order_customer_nation = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_order_customer_key, {},
            ct_customer_key, {ct_customer_nation}, key_dom.basis, key_dom.alpha,
            runtime.slots).front();
        auto line_customer_nation = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_line_order_key, {},
            ct_order_key, {order_customer_nation}, key_dom.basis, key_dom.alpha,
            runtime.slots).front();
        auto line_order_date = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_line_order_key, {},
            ct_order_key, {ct_order_date}, key_dom.basis, key_dom.alpha,
            runtime.slots).front();
        auto line_supplier_nation = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_line_supp_key, {},
            ct_supplier_key, {ct_supplier_nation}, key_dom.basis, key_dom.alpha,
            runtime.slots).front();
        auto supplier_asia = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_supplier_nation, {},
            ct_nation_key, {ct_nation_asia}, nation_dom.basis, nation_dom.alpha,
            runtime.slots).front();
        auto line_asia = LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_line_supp_key, {},
            ct_supplier_key, {supplier_asia}, key_dom.basis, key_dom.alpha,
            runtime.slots).front();
        auto same_nation = EqualityMask(
            runtime, line_customer_nation, line_supplier_nation, nation_dom);
        auto mask = runtime.cc->EvalMult(line_order_date, line_asia);
        mask = runtime.cc->EvalMult(mask, same_nation);
        auto filtered_revenue = runtime.cc->EvalMult(mask, ct_revenue);
        group_sums = MatrixGroupBySum(
            runtime.cc, runtime.keys.publicKey, nation_dom.attrs,
            {line_supplier_nation}, filtered_revenue, nation_dom.basis,
            nation_dom.alpha, runtime.slots);
    }));

    result.plain_groups = plain_groups;
    result.encrypted_groups = DecryptGroupSums(runtime, group_sums);
    result.abs_error = MaxAbsError(result.plain_groups, result.encrypted_groups);
    result.plain_scalar = std::accumulate(plain_groups.begin(), plain_groups.end(), 0.0);
    result.encrypted_scalar = std::accumulate(
        result.encrypted_groups.begin(), result.encrypted_groups.end(), 0.0);
    FinalizeHe3dbTiming(result);
    return result;
}

}  // namespace PaperReview
