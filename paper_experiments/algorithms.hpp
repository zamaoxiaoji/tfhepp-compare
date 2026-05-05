#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "binfhecontext.h"
#include "openfhe.h"

#include "key.hpp"
#include "metapbs2/bit_extraction.hpp"
#include "metapbs2/gap_msb.hpp"
#include "metapbs2/hom_compare.hpp"
#include "metapbs2/paper_params.hpp"
#include "tlwe.hpp"

namespace PaperReview {

using CkksCiphertext = lbcrypto::Ciphertext<lbcrypto::DCRTPoly>;
using CkksContext = lbcrypto::CryptoContext<lbcrypto::DCRTPoly>;
using CkksKeyPair = lbcrypto::KeyPair<lbcrypto::DCRTPoly>;
using CkksPublicKey = lbcrypto::PublicKey<lbcrypto::DCRTPoly>;

constexpr int kRequiredSecurityBits = 128;

struct Chapter3Params {
    int p = 12;
    int k = 7;
    int kappa = 5;
    int rounds = 2;
    bool enable_periodic_pruning = true;
};

MetaPBS2::Algorithm1Config Chapter3MetaPBSConfig();
MetaPBS2::BitExtractOptions ToBitExtractOptions(const Chapter3Params& params);
MetaPBS2::GapMSBOptions ToGapMSBOptions(const Chapter3Params& params);
std::string Chapter3ParameterSummary(const Chapter3Params& params);
std::vector<std::string> Chapter3AlgorithmTrace();
std::string SecuritySummary();

template <class brP>
TFHEpp::TLWE<typename brP::targetP> Chapter3BitExtract(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<MetaPBS2::TruncRepeatKey<typename brP::targetP>>& trkeys,
    const MetaPBS2::Algorithm1Config& cfg,
    const Chapter3Params& params,
    MetaPBS2::BlindRotatePruneStats* prune_stats = nullptr) {
    return MetaPBS2::BitExtract<brP>(
        ct, bkfft, trkeys, cfg, ToBitExtractOptions(params), prune_stats);
}

template <class brP>
TFHEpp::TLWE<typename brP::targetP> Chapter3GapMSB(
    const TFHEpp::TLWE<typename brP::domainP>& ct,
    const TFHEpp::BootstrappingKeyFFT<brP>& bkfft,
    const std::vector<MetaPBS2::TruncRepeatKey<typename brP::targetP>>& trkeys,
    const MetaPBS2::Algorithm1Config& cfg,
    const Chapter3Params& params,
    MetaPBS2::BlindRotatePruneStats* prune_stats = nullptr) {
    return MetaPBS2::GapMSB<brP>(
        ct, bkfft, trkeys, cfg, ToGapMSBOptions(params), prune_stats);
}

enum class ComparePredicate { LessThan, GreaterThan, LessEqual, GreaterEqual };

template <class brP_metapbs, class brP_logari, class brP_base, class iksP>
TFHEpp::TLWE<typename brP_base::targetP> Chapter3HomCompare(
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& lhs,
    const TFHEpp::TLWE<typename brP_metapbs::domainP>& rhs,
    int plain_bits,
    ComparePredicate predicate,
    const TFHEpp::BootstrappingKeyFFT<brP_metapbs>& bkfft,
    const std::vector<MetaPBS2::TruncRepeatKey<typename brP_metapbs::targetP>>& trkeys,
    const MetaPBS2::Algorithm1Config& cfg,
    const TFHEpp::BootstrappingKeyFFT<brP_logari>& bkfft_logari,
    const TFHEpp::BootstrappingKeyFFT<brP_base>& bkfft_base,
    const TFHEpp::KeySwitchingKey<iksP>& iksk,
    const MetaPBS2::HomMSBOptions& options = MetaPBS2::HomMSBOptions{},
    MetaPBS2::BlindRotatePruneStats* prune_stats = nullptr) {
    switch (predicate) {
    case ComparePredicate::LessThan:
        return MetaPBS2::HomLessThan<brP_metapbs, brP_logari, brP_base, iksP>(
            lhs, rhs, plain_bits, bkfft, trkeys, cfg,
            bkfft_logari, bkfft_base, iksk, options, prune_stats);
    case ComparePredicate::GreaterThan:
        return MetaPBS2::HomGreaterThan<brP_metapbs, brP_logari, brP_base, iksP>(
            lhs, rhs, plain_bits, bkfft, trkeys, cfg,
            bkfft_logari, bkfft_base, iksk, options, prune_stats);
    case ComparePredicate::LessEqual:
        return MetaPBS2::HomLessThanEqual<brP_metapbs, brP_logari, brP_base, iksP>(
            lhs, rhs, plain_bits, bkfft, trkeys, cfg,
            bkfft_logari, bkfft_base, iksk, options, prune_stats);
    case ComparePredicate::GreaterEqual:
        return MetaPBS2::HomGreaterThanEqual<brP_metapbs, brP_logari, brP_base, iksP>(
            lhs, rhs, plain_bits, bkfft, trkeys, cfg,
            bkfft_logari, bkfft_base, iksk, options, prune_stats);
    }
    throw std::invalid_argument("unknown comparison predicate");
}

struct AttributeSpec {
    std::size_t num_categories = 0;
};

struct DomainPoint {
    std::vector<int> coordinates;
};

struct BasisFunction {
    std::vector<int> exponents;
};

struct BasisMatrix {
    std::vector<std::vector<double>> values;
};

struct AlphaTable {
    std::vector<std::vector<double>> coefficients;
    std::vector<std::size_t> target_indices;
};

std::size_t TotalDomainSize(const std::vector<AttributeSpec>& attrs);
std::vector<DomainPoint> EnumerateDomainTuples(const std::vector<AttributeSpec>& attrs);
std::vector<BasisFunction> EnumerateMonomialBasis(const std::vector<AttributeSpec>& attrs);
double EvaluateBasisAtPoint(const BasisFunction& basis, const DomainPoint& point);
BasisMatrix BuildBasisMatrix(
    const std::vector<DomainPoint>& domain,
    const std::vector<BasisFunction>& basis);
void BuildBasis(
    const std::vector<std::vector<int>>& domains,
    const std::vector<BasisFunction>& basis,
    BasisMatrix& matrix,
    std::vector<DomainPoint>& points);
AlphaTable SolveCoeffTable(
    const BasisMatrix& matrix,
    const std::vector<std::size_t>& target_indices);
AlphaTable SolveCoeffTableForAllTargets(const BasisMatrix& matrix);
bool VerifyAlphaTable(
    const BasisMatrix& matrix,
    const AlphaTable& alpha_table,
    double tolerance = 1e-8);

struct CkksRuntime {
    CkksContext cc;
    CkksKeyPair keys;
    std::size_t slots = 0;
};

CkksRuntime MakeCkksRuntime(
    std::size_t slots,
    std::uint32_t multiplicative_depth,
    bool enable_scheme_switch = false,
    std::uint32_t scaling_mod_size = 50);
std::vector<int32_t> RotationIndicesForPowerOfTwoSum(std::size_t slots);
CkksCiphertext EncryptVector(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<double>& values,
    std::size_t slots);
CkksCiphertext EncryptConstant(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    double value,
    std::size_t slots);
CkksCiphertext EvalRotateAndSum(
    const CkksContext& cc,
    CkksCiphertext ct,
    std::size_t slots);
std::vector<std::vector<CkksCiphertext>> PrecomputeAttributePowers(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<AttributeSpec>& attrs,
    const std::vector<CkksCiphertext>& encrypted_attributes,
    std::size_t slots);
CkksCiphertext MultiplyBalanced(
    const CkksContext& cc,
    const std::vector<CkksCiphertext>& factors);
std::vector<CkksCiphertext> BuildEncryptedBasis(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<AttributeSpec>& attrs,
    const std::vector<CkksCiphertext>& encrypted_attributes,
    const std::vector<BasisFunction>& basis,
    std::size_t slots);
CkksCiphertext LinearCombination(
    const CkksContext& cc,
    const std::vector<CkksCiphertext>& terms,
    const std::vector<double>& coefficients);
std::vector<CkksCiphertext> MatrixGroupBySum(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<AttributeSpec>& group_attributes,
    const std::vector<CkksCiphertext>& encrypted_group_columns,
    const CkksCiphertext& encrypted_value_column,
    const std::vector<BasisFunction>& basis,
    const AlphaTable& alpha_table,
    std::size_t slots);
std::vector<CkksCiphertext> LookupJoin(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const CkksCiphertext& left_key,
    const std::vector<CkksCiphertext>& left_payload_columns,
    const CkksCiphertext& right_key,
    const std::vector<CkksCiphertext>& right_payload_columns,
    const std::vector<BasisFunction>& basis,
    const AlphaTable& alpha_table,
    std::size_t slots);

std::vector<int> BaseDigits(int value, int base, int digit_count);
std::vector<std::vector<double>> PlainDigitColumns(
    const std::vector<int>& values,
    int base,
    int digit_count,
    std::size_t slots);
std::vector<CkksCiphertext> EncryptDigitColumns(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<int>& values,
    int base,
    int digit_count,
    std::size_t slots);
CkksCiphertext DigitMaskFromDigits(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<CkksCiphertext>& digit_ciphertexts,
    int target_value,
    int base,
    int digit_count,
    std::size_t slots);

template <class LweP>
lbcrypto::LWEPrivateKey TFHEppKeyToOpenFHE(
    const TFHEpp::Key<LweP>& secret_key,
    std::uint64_t modulus) {
    const std::size_t dimension = LweP::k * LweP::n;
    lbcrypto::NativeVector lwe_key(dimension, lbcrypto::NativeInteger(modulus));
    for (std::size_t i = 0; i < dimension; i++) {
        const auto raw = secret_key[i];
        if constexpr (LweP::key_value_min < 0) {
            if (raw > static_cast<typename LweP::T>(LweP::key_value_max))
                lwe_key[i] = lbcrypto::NativeInteger(modulus - 1);
            else
                lwe_key[i] = lbcrypto::NativeInteger(static_cast<std::uint64_t>(raw));
        } else {
            lwe_key[i] = lbcrypto::NativeInteger(static_cast<std::uint64_t>(raw));
        }
    }
    return std::make_shared<lbcrypto::LWEPrivateKeyImpl>(std::move(lwe_key));
}

template <class LweP>
std::vector<lbcrypto::LWECiphertext> TFHEppLWEsToOpenFHE(
    const std::vector<TFHEpp::TLWE<LweP>>& lwes,
    std::uint64_t q_target) {
    const std::size_t dimension = LweP::k * LweP::n;
    const long double q_source =
        std::ldexp(1.0L, std::numeric_limits<typename LweP::T>::digits);
    const long double ratio = static_cast<long double>(q_target) / q_source;

    auto switch_coeff = [ratio, q_target](typename LweP::T coeff) {
        const long double scaled = static_cast<long double>(coeff) * ratio;
        auto rounded = static_cast<std::int64_t>(std::llround(scaled));
        rounded %= static_cast<std::int64_t>(q_target);
        if (rounded < 0) rounded += static_cast<std::int64_t>(q_target);
        return lbcrypto::NativeInteger(static_cast<std::uint64_t>(rounded));
    };

    std::vector<lbcrypto::LWECiphertext> out;
    out.reserve(lwes.size());
    for (const auto& lwe : lwes) {
        lbcrypto::NativeVector a(dimension, lbcrypto::NativeInteger(q_target));
        for (std::size_t i = 0; i < dimension; i++)
            a[i] = switch_coeff(lwe[i]);
        out.push_back(std::make_shared<lbcrypto::LWECiphertextImpl>(
            std::move(a), switch_coeff(lwe[dimension])));
    }
    return out;
}

template <class LweP>
struct RepackContext {
    CkksContext cc;
    CkksKeyPair keys;
    std::uint32_t slots = 0;
    std::uint32_t log_q_lwe = 28;
    std::uint64_t q_target = 1ULL << 28;
    bool ready = false;
};

template <class LweP>
RepackContext<LweP> RepackSetup(
    CkksContext& cc,
    const CkksKeyPair& keys,
    const TFHEpp::Key<LweP>& tfhe_secret_key,
    std::uint32_t slots,
    std::uint32_t log_q_lwe = 28) {
    RepackContext<LweP> context;
    context.cc = cc;
    context.keys = keys;
    context.slots = slots;
    context.log_q_lwe = log_q_lwe;
    context.q_target = 1ULL << log_q_lwe;

    auto lwe_context = std::make_shared<lbcrypto::BinFHEContext>();
    lwe_context->GenerateBinFHEContext(
        lbcrypto::STD128, false, log_q_lwe, 0, lbcrypto::GINX, false);

    auto lwe_secret_key =
        TFHEppKeyToOpenFHE<LweP>(tfhe_secret_key, context.q_target);

    cc->EvalFHEWtoCKKSSetup(lwe_context, slots, log_q_lwe);
    cc->SetBinCCForSchemeSwitch(lwe_context);
    cc->EvalFHEWtoCKKSKeyGen(keys, lwe_secret_key);

    context.ready = true;
    return context;
}

template <class LweP>
CkksCiphertext RepackExecute(
    RepackContext<LweP>& context,
    const std::vector<TFHEpp::TLWE<LweP>>& lwes) {
    if (!context.ready)
        throw std::logic_error("RepackExecute called before RepackSetup");
    auto openfhe_lwes = TFHEppLWEsToOpenFHE<LweP>(lwes, context.q_target);
    auto result =
        context.cc->EvalFHEWtoCKKS(openfhe_lwes, context.slots, context.slots);
    return context.cc->EvalMult(result, std::sqrt(2.0));
}

}  // namespace PaperReview
