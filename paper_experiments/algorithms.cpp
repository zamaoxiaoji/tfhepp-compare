#include "algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace PaperReview {

namespace {

double IntPow(double base, int exp) {
    if (exp < 0) throw std::invalid_argument("negative exponent");
    double out = 1.0;
    for (int i = 0; i < exp; i++) out *= base;
    return out;
}

std::vector<std::vector<double>> InvertSquareMatrix(
    const std::vector<std::vector<double>>& matrix) {
    const std::size_t n = matrix.size();
    if (n == 0) throw std::invalid_argument("cannot invert empty matrix");
    for (const auto& row : matrix)
        if (row.size() != n)
            throw std::invalid_argument("basis matrix must be square");

    std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n, 0.0));
    for (std::size_t r = 0; r < n; r++) {
        for (std::size_t c = 0; c < n; c++) aug[r][c] = matrix[r][c];
        aug[r][n + r] = 1.0;
    }

    for (std::size_t col = 0; col < n; col++) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < n; r++)
            if (std::abs(aug[r][col]) > std::abs(aug[pivot][col]))
                pivot = r;
        if (std::abs(aug[pivot][col]) < 1e-12)
            throw std::invalid_argument("basis matrix is singular");
        if (pivot != col) std::swap(aug[pivot], aug[col]);

        const double pivot_value = aug[col][col];
        for (double& value : aug[col]) value /= pivot_value;

        for (std::size_t r = 0; r < n; r++) {
            if (r == col) continue;
            const double factor = aug[r][col];
            if (std::abs(factor) < 1e-18) continue;
            for (std::size_t c = 0; c < 2 * n; c++)
                aug[r][c] -= factor * aug[col][c];
        }
    }

    std::vector<std::vector<double>> inverse(n, std::vector<double>(n));
    for (std::size_t r = 0; r < n; r++)
        for (std::size_t c = 0; c < n; c++)
            inverse[r][c] = aug[r][n + c];
    return inverse;
}

std::vector<std::vector<double>> BuildLagrangeCoefficients(int domain_size) {
    if (domain_size <= 0)
        throw std::invalid_argument("domain_size must be positive");
    std::vector<std::vector<double>> coeffs(
        static_cast<std::size_t>(domain_size),
        std::vector<double>(static_cast<std::size_t>(domain_size), 0.0));

    for (int target = 0; target < domain_size; target++) {
        std::vector<double> poly{1.0};
        double denominator = 1.0;
        for (int x = 0; x < domain_size; x++) {
            if (x == target) continue;
            std::vector<double> next(poly.size() + 1, 0.0);
            for (std::size_t d = 0; d < poly.size(); d++) {
                next[d] -= static_cast<double>(x) * poly[d];
                next[d + 1] += poly[d];
            }
            poly.swap(next);
            denominator *= static_cast<double>(target - x);
        }
        for (std::size_t d = 0; d < poly.size(); d++)
            coeffs[static_cast<std::size_t>(target)][d] =
                poly[d] / denominator;
    }
    return coeffs;
}

std::vector<CkksCiphertext> BuildPowers(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const CkksCiphertext& base,
    int max_power,
    std::size_t slots) {
    if (max_power < 0) throw std::invalid_argument("max_power is negative");
    std::vector<CkksCiphertext> powers(static_cast<std::size_t>(max_power) + 1);
    std::vector<bool> ready(static_cast<std::size_t>(max_power) + 1, false);
    powers[0] = EncryptConstant(cc, public_key, 1.0, slots);
    ready[0] = true;
    if (max_power >= 1) {
        powers[1] = base;
        ready[1] = true;
    }

    std::function<const CkksCiphertext&(int)> compute_power =
        [&](int degree) -> const CkksCiphertext& {
        const auto idx = static_cast<std::size_t>(degree);
        if (ready[idx]) return powers[idx];
        const int left = degree / 2;
        const int right = degree - left;
        powers[idx] = cc->EvalMult(compute_power(left), compute_power(right));
        ready[idx] = true;
        return powers[idx];
    };

    for (int d = 2; d <= max_power; d++) compute_power(d);
    return powers;
}

}  // namespace

MetaPBS2::Algorithm1Config Chapter3MetaPBSConfig() {
    return MetaPBS2::PaperRowT2NConfig();
}

MetaPBS2::BitExtractOptions ToBitExtractOptions(const Chapter3Params& params) {
    return MetaPBS2::BitExtractOptions{
        .p = params.p,
        .k = params.k,
        .enable_periodic_pruning = params.enable_periodic_pruning,
        .period = 0,
    };
}

MetaPBS2::GapMSBOptions ToGapMSBOptions(const Chapter3Params& params) {
    MetaPBS2::GapMSBOptions options;
    options.p = params.p;
    options.k = params.k;
    options.enable_periodic_pruning = params.enable_periodic_pruning;
    return options;
}

std::string Chapter3ParameterSummary(const Chapter3Params& params) {
    return "p=" + std::to_string(params.p) +
           ", k=" + std::to_string(params.k) +
           ", kappa=" + std::to_string(params.kappa) +
           ", R=" + std::to_string(params.rounds) +
           ", periodic_pruning=" +
           (params.enable_periodic_pruning ? "on" : "off");
}

std::vector<std::string> Chapter3AlgorithmTrace() {
    return {
        "BitExtract: build 0/Q/2 bit LUT, run Meta-PBS Algorithm 1, apply M_k=2^(p-k) pruning only in the first blind rotation.",
        "GapMSB: extract bit_k, convert it to arithmetic weight w_k*Delta, subtract it from ct_I, add offset'=(w_k+1)Delta/2, then run sign PBS.",
        "HomCompare: subtract operands, run recursive HomMSB on the signed difference; each recursive clear adds offset2=w_k*Delta/2 immediately, and the base sign PBS uses offset1=Delta/2.",
    };
}

std::string SecuritySummary() {
    return "TFHEpp default params/128bit.hpp and OpenFHE HEStd_128_classic/STD128 are used; no 80-bit parameter option is enabled by paper_experiments.";
}

std::ostream* OpenOptionalOutputFile(
    int& index,
    int argc,
    char** argv,
    std::unique_ptr<std::ofstream>& file) {
    const std::string arg = argv[index];
    if (arg != "--output" && arg != "--output-file") return nullptr;
    if (index + 1 >= argc)
        throw std::invalid_argument("--output requires a file path");
    file = std::make_unique<std::ofstream>(argv[++index]);
    if (!*file)
        throw std::runtime_error("failed to open output file: " + std::string(argv[index]));
    return file.get();
}

void WriteOutputLine(
    std::ostream& primary,
    std::ostream* secondary,
    const std::string& line) {
    primary << line << "\n";
    if (secondary) *secondary << line << "\n";
}

std::size_t TotalDomainSize(const std::vector<AttributeSpec>& attrs) {
    if (attrs.empty()) throw std::invalid_argument("at least one attribute is required");
    std::size_t total = 1;
    for (const auto& attr : attrs) {
        if (attr.num_categories == 0)
            throw std::invalid_argument("attribute domain must be nonempty");
        if (total > std::numeric_limits<std::size_t>::max() / attr.num_categories)
            throw std::overflow_error("domain size overflow");
        total *= attr.num_categories;
    }
    return total;
}

std::vector<DomainPoint> EnumerateDomainTuples(
    const std::vector<AttributeSpec>& attrs) {
    const std::size_t total = TotalDomainSize(attrs);
    std::vector<DomainPoint> points(total);
    for (std::size_t index = 0; index < total; index++) {
        std::size_t x = index;
        points[index].coordinates.resize(attrs.size());
        for (std::size_t i = 0; i < attrs.size(); i++) {
            points[index].coordinates[i] =
                static_cast<int>(x % attrs[i].num_categories);
            x /= attrs[i].num_categories;
        }
    }
    return points;
}

std::vector<BasisFunction> EnumerateMonomialBasis(
    const std::vector<AttributeSpec>& attrs) {
    const auto points = EnumerateDomainTuples(attrs);
    std::vector<BasisFunction> basis(points.size());
    for (std::size_t i = 0; i < points.size(); i++)
        basis[i].exponents = points[i].coordinates;
    return basis;
}

double EvaluateBasisAtPoint(
    const BasisFunction& basis,
    const DomainPoint& point) {
    if (basis.exponents.size() != point.coordinates.size())
        throw std::invalid_argument("basis dimension mismatch");
    double product = 1.0;
    for (std::size_t i = 0; i < basis.exponents.size(); i++)
        product *= IntPow(static_cast<double>(point.coordinates[i]), basis.exponents[i]);
    return product;
}

BasisMatrix BuildBasisMatrix(
    const std::vector<DomainPoint>& domain,
    const std::vector<BasisFunction>& basis) {
    if (domain.empty() || basis.empty())
        throw std::invalid_argument("domain and basis must be nonempty");
    BasisMatrix matrix;
    matrix.values.assign(domain.size(), std::vector<double>(basis.size(), 0.0));
    for (std::size_t row = 0; row < domain.size(); row++)
        for (std::size_t col = 0; col < basis.size(); col++)
            matrix.values[row][col] = EvaluateBasisAtPoint(basis[col], domain[row]);
    return matrix;
}

void BuildBasis(
    const std::vector<std::vector<int>>& domains,
    const std::vector<BasisFunction>& basis,
    BasisMatrix& matrix,
    std::vector<DomainPoint>& points) {
    std::vector<AttributeSpec> attrs;
    attrs.reserve(domains.size());
    for (const auto& domain : domains) {
        if (domain.empty()) throw std::invalid_argument("domain is empty");
        attrs.push_back(AttributeSpec{domain.size()});
    }

    const auto enumerated = EnumerateDomainTuples(attrs);
    points.resize(enumerated.size());
    for (std::size_t i = 0; i < enumerated.size(); i++) {
        points[i].coordinates.resize(domains.size());
        for (std::size_t d = 0; d < domains.size(); d++)
            points[i].coordinates[d] =
                domains[d][static_cast<std::size_t>(enumerated[i].coordinates[d])];
    }
    matrix = BuildBasisMatrix(points, basis);
}

AlphaTable SolveCoeffTable(
    const BasisMatrix& matrix,
    const std::vector<std::size_t>& target_indices) {
    const auto inverse = InvertSquareMatrix(matrix.values);
    const std::size_t n = inverse.size();
    AlphaTable table;
    table.target_indices = target_indices;
    table.coefficients.assign(target_indices.size(), std::vector<double>(n, 0.0));
    for (std::size_t t = 0; t < target_indices.size(); t++) {
        if (target_indices[t] >= n)
            throw std::out_of_range("target index outside basis matrix");
        for (std::size_t q = 0; q < n; q++)
            table.coefficients[t][q] = inverse[q][target_indices[t]];
    }
    return table;
}

AlphaTable SolveCoeffTableForAllTargets(const BasisMatrix& matrix) {
    std::vector<std::size_t> targets(matrix.values.size());
    std::iota(targets.begin(), targets.end(), std::size_t{0});
    return SolveCoeffTable(matrix, targets);
}

bool VerifyAlphaTable(
    const BasisMatrix& matrix,
    const AlphaTable& alpha_table,
    double tolerance) {
    const std::size_t rows = matrix.values.size();
    if (rows == 0) return false;
    for (std::size_t t = 0; t < alpha_table.target_indices.size(); t++) {
        const auto target = alpha_table.target_indices[t];
        for (std::size_t row = 0; row < rows; row++) {
            double value = 0.0;
            for (std::size_t q = 0; q < matrix.values[row].size(); q++)
                value += matrix.values[row][q] * alpha_table.coefficients[t][q];
            const double expected = row == target ? 1.0 : 0.0;
            if (std::abs(value - expected) > tolerance) return false;
        }
    }
    return true;
}

CkksRuntime MakeCkksRuntime(
    std::size_t slots,
    std::uint32_t multiplicative_depth,
    bool enable_scheme_switch,
    std::uint32_t scaling_mod_size) {
    if (slots == 0) throw std::invalid_argument("slots must be nonzero");
    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
    params.SetMultiplicativeDepth(multiplicative_depth);
    params.SetScalingModSize(scaling_mod_size);
    params.SetScalingTechnique(lbcrypto::FIXEDAUTO);
    params.SetSecurityLevel(lbcrypto::HEStd_128_classic);
    params.SetBatchSize(slots);

    CkksRuntime runtime;
    runtime.cc = lbcrypto::GenCryptoContext(params);
    runtime.slots = slots;
    runtime.cc->Enable(lbcrypto::PKE);
    runtime.cc->Enable(lbcrypto::KEYSWITCH);
    runtime.cc->Enable(lbcrypto::LEVELEDSHE);
    runtime.cc->Enable(lbcrypto::ADVANCEDSHE);
    if (enable_scheme_switch) runtime.cc->Enable(lbcrypto::SCHEMESWITCH);

    runtime.keys = runtime.cc->KeyGen();
    runtime.cc->EvalMultKeyGen(runtime.keys.secretKey);
    runtime.cc->EvalRotateKeyGen(
        runtime.keys.secretKey, RotationIndicesForPowerOfTwoSum(slots));
    return runtime;
}

std::vector<int32_t> RotationIndicesForPowerOfTwoSum(std::size_t slots) {
    std::vector<int32_t> rotations;
    for (std::size_t step = 1; step < slots; step <<= 1)
        rotations.push_back(static_cast<int32_t>(step));
    return rotations;
}

CkksCiphertext EncryptVector(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<double>& values,
    std::size_t slots) {
    std::vector<double> padded(slots, 0.0);
    if (values.size() > slots)
        throw std::invalid_argument("values exceed CKKS slot count");
    std::copy(values.begin(), values.end(), padded.begin());
    return cc->Encrypt(public_key, cc->MakeCKKSPackedPlaintext(padded));
}

CkksCiphertext EncryptConstant(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    double value,
    std::size_t slots) {
    return EncryptVector(cc, public_key, std::vector<double>(slots, value), slots);
}

CkksCiphertext EvalRotateAndSum(
    const CkksContext& cc,
    CkksCiphertext ct,
    std::size_t slots) {
    for (std::size_t step = 1; step < slots; step <<= 1)
        ct = cc->EvalAdd(ct, cc->EvalRotate(ct, static_cast<int32_t>(step)));
    return ct;
}

std::vector<std::vector<CkksCiphertext>> PrecomputeAttributePowers(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<AttributeSpec>& attrs,
    const std::vector<CkksCiphertext>& encrypted_attributes,
    std::size_t slots) {
    if (attrs.size() != encrypted_attributes.size())
        throw std::invalid_argument("attribute/ciphertext count mismatch");
    std::vector<std::vector<CkksCiphertext>> powers(attrs.size());
    for (std::size_t i = 0; i < attrs.size(); i++) {
        const int max_exp = static_cast<int>(attrs[i].num_categories) - 1;
        powers[i] = BuildPowers(
            cc, public_key, encrypted_attributes[i], max_exp, slots);
    }
    return powers;
}

CkksCiphertext MultiplyBalanced(
    const CkksContext& cc,
    const std::vector<CkksCiphertext>& factors) {
    if (factors.empty()) throw std::invalid_argument("no factors to multiply");
    std::vector<CkksCiphertext> layer = factors;
    while (layer.size() > 1) {
        std::vector<CkksCiphertext> next;
        next.reserve((layer.size() + 1) / 2);
        for (std::size_t i = 0; i < layer.size(); i += 2) {
            if (i + 1 < layer.size())
                next.push_back(cc->EvalMult(layer[i], layer[i + 1]));
            else
                next.push_back(layer[i]);
        }
        layer.swap(next);
    }
    return layer.front();
}

std::vector<CkksCiphertext> BuildEncryptedBasis(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<AttributeSpec>& attrs,
    const std::vector<CkksCiphertext>& encrypted_attributes,
    const std::vector<BasisFunction>& basis,
    std::size_t slots) {
    const auto powers =
        PrecomputeAttributePowers(cc, public_key, attrs, encrypted_attributes, slots);
    const auto one = EncryptConstant(cc, public_key, 1.0, slots);
    std::vector<CkksCiphertext> encrypted_basis;
    encrypted_basis.reserve(basis.size());

    for (const auto& phi : basis) {
        if (phi.exponents.size() != attrs.size())
            throw std::invalid_argument("basis dimension mismatch");
        std::vector<CkksCiphertext> factors;
        for (std::size_t attr = 0; attr < attrs.size(); attr++) {
            const int exp = phi.exponents[attr];
            if (exp < 0 || exp >= static_cast<int>(powers[attr].size()))
                throw std::invalid_argument("basis exponent out of range");
            if (exp != 0)
                factors.push_back(powers[attr][static_cast<std::size_t>(exp)]);
        }
        encrypted_basis.push_back(factors.empty() ? one : MultiplyBalanced(cc, factors));
    }
    return encrypted_basis;
}

CkksCiphertext LinearCombination(
    const CkksContext& cc,
    const std::vector<CkksCiphertext>& terms,
    const std::vector<double>& coefficients) {
    if (terms.empty()) throw std::invalid_argument("empty linear-combination terms");
    if (terms.size() != coefficients.size())
        throw std::invalid_argument("linear-combination size mismatch");

    auto result = cc->EvalMult(terms.front(), 0.0);
    for (std::size_t i = 0; i < terms.size(); i++) {
        if (std::abs(coefficients[i]) < 1e-12) continue;
        result = cc->EvalAdd(result, cc->EvalMult(terms[i], coefficients[i]));
    }
    return result;
}

std::vector<CkksCiphertext> MatrixGroupBySum(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<AttributeSpec>& group_attributes,
    const std::vector<CkksCiphertext>& encrypted_group_columns,
    const CkksCiphertext& encrypted_value_column,
    const std::vector<BasisFunction>& basis,
    const AlphaTable& alpha_table,
    std::size_t slots) {
    const auto encrypted_basis = BuildEncryptedBasis(
        cc, public_key, group_attributes, encrypted_group_columns, basis, slots);

    std::vector<CkksCiphertext> weighted_basis(encrypted_basis.size());
    for (std::size_t q = 0; q < encrypted_basis.size(); q++)
        weighted_basis[q] = cc->EvalMult(encrypted_basis[q], encrypted_value_column);

    std::vector<CkksCiphertext> sums;
    sums.reserve(alpha_table.coefficients.size());
    for (const auto& alpha : alpha_table.coefficients) {
        auto masked_values = LinearCombination(cc, weighted_basis, alpha);
        sums.push_back(EvalRotateAndSum(cc, masked_values, slots));
    }
    return sums;
}

std::vector<CkksCiphertext> LookupJoin(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const CkksCiphertext& left_key,
    const std::vector<CkksCiphertext>& left_payload_columns,
    const CkksCiphertext& right_key,
    const std::vector<CkksCiphertext>& right_payload_columns,
    const std::vector<BasisFunction>& basis,
    const AlphaTable& alpha_table,
    std::size_t slots) {
    (void)left_payload_columns;
    if (right_payload_columns.empty())
        throw std::invalid_argument("JOIN requires at least one right payload column");

    const std::vector<AttributeSpec> key_attrs{{basis.size()}};
    const auto basis_left =
        BuildEncryptedBasis(cc, public_key, key_attrs, {left_key}, basis, slots);
    const auto basis_right =
        BuildEncryptedBasis(cc, public_key, key_attrs, {right_key}, basis, slots);

    std::vector<std::vector<CkksCiphertext>> right_bases(right_payload_columns.size());
    for (std::size_t v = 0; v < right_payload_columns.size(); v++) {
        right_bases[v].resize(basis.size());
        for (std::size_t q = 0; q < basis.size(); q++) {
            auto weighted = cc->EvalMult(basis_right[q], right_payload_columns[v]);
            right_bases[v][q] = EvalRotateAndSum(cc, weighted, slots);
        }
    }

    std::vector<CkksCiphertext> joined(right_payload_columns.size());
    for (std::size_t v = 0; v < right_payload_columns.size(); v++)
        joined[v] = cc->EvalMult(right_payload_columns[v], 0.0);

    for (const auto& alpha : alpha_table.coefficients) {
        auto left_mask = LinearCombination(cc, basis_left, alpha);
        for (std::size_t v = 0; v < right_payload_columns.size(); v++) {
            auto payload = LinearCombination(cc, right_bases[v], alpha);
            joined[v] = cc->EvalAdd(joined[v], cc->EvalMult(left_mask, payload));
        }
    }
    return joined;
}

std::vector<int> BaseDigits(int value, int base, int digit_count) {
    if (value < 0 || base <= 1 || digit_count <= 0)
        throw std::invalid_argument("invalid digit decomposition parameters");
    std::vector<int> digits(static_cast<std::size_t>(digit_count), 0);
    for (int i = 0; i < digit_count; i++) {
        digits[static_cast<std::size_t>(i)] = value % base;
        value /= base;
    }
    if (value != 0) throw std::invalid_argument("value exceeds digit capacity");
    return digits;
}

std::vector<std::vector<double>> PlainDigitColumns(
    const std::vector<int>& values,
    int base,
    int digit_count,
    std::size_t slots) {
    if (values.size() > slots)
        throw std::invalid_argument("values exceed CKKS slot count");
    std::vector<std::vector<double>> columns(
        static_cast<std::size_t>(digit_count), std::vector<double>(slots, 0.0));
    for (std::size_t row = 0; row < values.size(); row++) {
        const auto digits = BaseDigits(values[row], base, digit_count);
        for (int d = 0; d < digit_count; d++)
            columns[static_cast<std::size_t>(d)][row] =
                static_cast<double>(digits[static_cast<std::size_t>(d)]);
    }
    return columns;
}

std::vector<CkksCiphertext> EncryptDigitColumns(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<int>& values,
    int base,
    int digit_count,
    std::size_t slots) {
    const auto plain_columns = PlainDigitColumns(values, base, digit_count, slots);
    std::vector<CkksCiphertext> encrypted;
    encrypted.reserve(plain_columns.size());
    for (const auto& col : plain_columns)
        encrypted.push_back(EncryptVector(cc, public_key, col, slots));
    return encrypted;
}

CkksCiphertext DigitMaskFromDigits(
    const CkksContext& cc,
    const CkksPublicKey& public_key,
    const std::vector<CkksCiphertext>& digit_ciphertexts,
    int target_value,
    int base,
    int digit_count,
    std::size_t slots) {
    if (static_cast<int>(digit_ciphertexts.size()) != digit_count)
        throw std::invalid_argument("digit ciphertext count mismatch");
    const auto target_digits = BaseDigits(target_value, base, digit_count);
    const auto coeffs = BuildLagrangeCoefficients(base);

    std::vector<CkksCiphertext> digit_masks;
    digit_masks.reserve(digit_ciphertexts.size());
    for (int digit = 0; digit < digit_count; digit++) {
        auto powers = BuildPowers(
            cc, public_key, digit_ciphertexts[static_cast<std::size_t>(digit)],
            base - 1, slots);
        digit_masks.push_back(LinearCombination(
            cc, powers,
            coeffs[static_cast<std::size_t>(
                target_digits[static_cast<std::size_t>(digit)])]));
    }
    return MultiplyBalanced(cc, digit_masks);
}

}  // namespace PaperReview
