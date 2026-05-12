#include "ckks_relational.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace tfhepp_ckks
{
    namespace
    {
        void RequireActiveSlots(std::size_t active_slots, const char *name)
        {
            if (active_slots == 0)
                throw std::invalid_argument(std::string(name) +
                                            ": active_slots must be positive");
        }

        void ModSwitchToCommonLevel(seal::Ciphertext &lhs,
                                    seal::Ciphertext &rhs,
                                    seal::Evaluator &evaluator)
        {
            if (lhs.parms_id() == rhs.parms_id()) return;
            if (lhs.coeff_modulus_size() > rhs.coeff_modulus_size()) {
                evaluator.mod_switch_to_inplace(lhs, rhs.parms_id());
            }
            else if (rhs.coeff_modulus_size() > lhs.coeff_modulus_size()) {
                evaluator.mod_switch_to_inplace(rhs, lhs.parms_id());
            }
            else {
                throw std::invalid_argument(
                    "ciphertexts have different parameter ids at the same level");
            }
        }

        seal::Ciphertext MultiplyAndRescale(const seal::Ciphertext &lhs_in,
                                            const seal::Ciphertext &rhs_in,
                                            const seal::RelinKeys &relin_keys,
                                            seal::Evaluator &evaluator)
        {
            seal::Ciphertext lhs = lhs_in;
            seal::Ciphertext rhs = rhs_in;
            ModSwitchToCommonLevel(lhs, rhs, evaluator);

            const double target_scale = std::min(lhs.scale(), rhs.scale());
            seal::Ciphertext result;
            evaluator.multiply(lhs, rhs, result);
            evaluator.relinearize_inplace(result, relin_keys);
            evaluator.rescale_to_next_inplace(result);
            result.scale() = target_scale;
            return result;
        }

        seal::Ciphertext MultiplyPlainScalar(
            const seal::Ciphertext &cipher, double scalar,
            seal::CKKSEncoder &encoder, seal::Evaluator &evaluator)
        {
            if (scalar == 1.0) return cipher;

            const double target_scale = cipher.scale();
            seal::Plaintext plain;
            encoder.encode(scalar, cipher.parms_id(), cipher.scale(), plain);
            seal::Ciphertext result = cipher;
            evaluator.multiply_plain_inplace(result, plain);
            evaluator.rescale_to_next_inplace(result);
            result.scale() = target_scale;
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

        void AddPlainScalarInPlace(seal::Ciphertext &cipher, double scalar,
                                   seal::CKKSEncoder &encoder,
                                   seal::Evaluator &evaluator)
        {
            if (scalar == 0.0) return;

            seal::Plaintext plain;
            encoder.encode(scalar, cipher.parms_id(), cipher.scale(), plain);
            evaluator.add_plain_inplace(cipher, plain);
        }

        seal::Ciphertext ApplyActiveSlotMask(
            const seal::Ciphertext &cipher, std::size_t active_slots,
            seal::CKKSEncoder &encoder, seal::Evaluator &evaluator)
        {
            if (active_slots > encoder.slot_count())
                throw std::invalid_argument(
                    "ApplyActiveSlotMask: active_slots exceeds slot count");
            if (active_slots == encoder.slot_count()) return cipher;

            std::vector<double> slots(encoder.slot_count(), 0.0);
            std::fill(slots.begin(), slots.begin() + active_slots, 1.0);

            const double target_scale = cipher.scale();
            seal::Plaintext plain;
            encoder.encode(slots, cipher.parms_id(), cipher.scale(), plain);

            seal::Ciphertext result = cipher;
            evaluator.multiply_plain_inplace(result, plain);
            evaluator.rescale_to_next_inplace(result);
            result.scale() = target_scale;
            return result;
        }

        std::vector<std::size_t> DomainSizes(
            const std::vector<std::vector<double>> &domains,
            const char *name)
        {
            if (domains.empty())
                throw std::invalid_argument(std::string(name) +
                                            ": empty domain list");

            std::vector<std::size_t> sizes;
            sizes.reserve(domains.size());
            for (const auto &domain : domains) {
                if (domain.empty())
                    throw std::invalid_argument(std::string(name) +
                                                ": empty attribute domain");
                sizes.push_back(domain.size());
            }
            return sizes;
        }

        std::size_t ProductSize(const std::vector<std::size_t> &sizes)
        {
            std::size_t product = 1;
            for (const auto size : sizes) product *= size;
            return product;
        }

        std::vector<std::size_t> DecodeMixedRadix(
            std::size_t index, const std::vector<std::size_t> &radices)
        {
            std::vector<std::size_t> digits(radices.size(), 0);
            for (std::size_t rev = 0; rev < radices.size(); ++rev) {
                const std::size_t pos = radices.size() - 1 - rev;
                digits[pos] = index % radices[pos];
                index /= radices[pos];
            }
            return digits;
        }

        void ValidateSquareMatrix(const DenseMatrix &matrix, const char *name)
        {
            if (matrix.empty())
                throw std::invalid_argument(std::string(name) +
                                            ": empty matrix");
            const std::size_t n = matrix.size();
            for (const auto &row : matrix) {
                if (row.size() != n)
                    throw std::invalid_argument(std::string(name) +
                                                ": matrix must be square");
            }
        }

        DenseMatrix InvertMatrix(const DenseMatrix &matrix)
        {
            ValidateSquareMatrix(matrix, "InvertMatrix");

            const std::size_t n = matrix.size();
            DenseMatrix augmented(n, std::vector<double>(2 * n, 0.0));
            for (std::size_t row = 0; row < n; ++row) {
                for (std::size_t col = 0; col < n; ++col)
                    augmented[row][col] = matrix[row][col];
                augmented[row][n + row] = 1.0;
            }

            for (std::size_t col = 0; col < n; ++col) {
                std::size_t pivot = col;
                double pivot_abs = std::abs(augmented[col][col]);
                for (std::size_t row = col + 1; row < n; ++row) {
                    const double candidate = std::abs(augmented[row][col]);
                    if (candidate > pivot_abs) {
                        pivot = row;
                        pivot_abs = candidate;
                    }
                }
                if (pivot_abs <= std::numeric_limits<double>::epsilon())
                    throw std::invalid_argument(
                        "InvertMatrix: singular basis matrix");

                if (pivot != col) std::swap(augmented[pivot], augmented[col]);

                const double pivot_value = augmented[col][col];
                for (double &value : augmented[col]) value /= pivot_value;

                for (std::size_t row = 0; row < n; ++row) {
                    if (row == col) continue;
                    const double factor = augmented[row][col];
                    if (factor == 0.0) continue;
                    for (std::size_t k = 0; k < 2 * n; ++k)
                        augmented[row][k] -= factor * augmented[col][k];
                }
            }

            DenseMatrix inverse(n, std::vector<double>(n, 0.0));
            for (std::size_t row = 0; row < n; ++row)
                for (std::size_t col = 0; col < n; ++col)
                    inverse[row][col] = augmented[row][n + col];
            return inverse;
        }
    } // namespace

    ExponentTable EnumerateExponents(
        const std::vector<std::vector<double>> &domains)
    {
        const auto sizes = DomainSizes(domains, "EnumerateExponents");
        const std::size_t total_size = ProductSize(sizes);

        ExponentTable exponents;
        exponents.reserve(total_size);
        for (std::size_t index = 0; index < total_size; ++index)
            exponents.push_back(DecodeMixedRadix(index, sizes));
        return exponents;
    }

    DenseMatrix EnumerateDomainPoints(
        const std::vector<std::vector<double>> &domains)
    {
        const auto sizes = DomainSizes(domains, "EnumerateDomainPoints");
        const std::size_t total_size = ProductSize(sizes);

        DenseMatrix points;
        points.reserve(total_size);
        for (std::size_t index = 0; index < total_size; ++index) {
            const auto digits = DecodeMixedRadix(index, sizes);
            std::vector<double> point(domains.size(), 0.0);
            for (std::size_t attr = 0; attr < domains.size(); ++attr)
                point[attr] = domains[attr][digits[attr]];
            points.push_back(std::move(point));
        }
        return points;
    }

    DenseMatrix BuildBasisMatrix(
        const std::vector<std::vector<double>> &domains)
    {
        const auto points = EnumerateDomainPoints(domains);
        const auto exponents = EnumerateExponents(domains);

        DenseMatrix basis(points.size(),
                          std::vector<double>(exponents.size(), 1.0));
        for (std::size_t row = 0; row < points.size(); ++row) {
            for (std::size_t col = 0; col < exponents.size(); ++col) {
                double value = 1.0;
                for (std::size_t attr = 0; attr < domains.size(); ++attr) {
                    const std::size_t exponent = exponents[col][attr];
                    if (exponent != 0)
                        value *= std::pow(points[row][attr],
                                          static_cast<double>(exponent));
                }
                basis[row][col] = value;
            }
        }
        return basis;
    }

    std::vector<std::vector<double>> SolveAlphaTable(
        const DenseMatrix &basis_matrix, double verification_tolerance)
    {
        ValidateSquareMatrix(basis_matrix, "SolveAlphaTable");
        if (verification_tolerance <= 0.0)
            throw std::invalid_argument(
                "SolveAlphaTable: verification tolerance must be positive");

        const auto inverse = InvertMatrix(basis_matrix);
        const std::size_t n = basis_matrix.size();

        std::vector<std::vector<double>> alpha_table(
            n, std::vector<double>(n, 0.0));
        for (std::size_t theta = 0; theta < n; ++theta)
            for (std::size_t q = 0; q < n; ++q)
                alpha_table[theta][q] = inverse[q][theta];

        double max_abs_error = 0.0;
        for (std::size_t theta = 0; theta < n; ++theta) {
            for (std::size_t row = 0; row < n; ++row) {
                double got = 0.0;
                for (std::size_t q = 0; q < n; ++q)
                    got += basis_matrix[row][q] * alpha_table[theta][q];
                const double expected = (row == theta) ? 1.0 : 0.0;
                max_abs_error =
                    std::max(max_abs_error, std::abs(got - expected));
            }
        }
        if (max_abs_error > verification_tolerance)
            throw std::invalid_argument(
                "SolveAlphaTable: basis verification failed");

        return alpha_table;
    }

    std::vector<std::vector<double>> SolveLagrangeAlphaTable(
        const std::vector<double> &domain)
    {
        if (domain.empty())
            throw std::invalid_argument(
                "SolveLagrangeAlphaTable: empty domain");

        const std::size_t degree = domain.size() - 1;
        std::vector<std::vector<double>> alpha_table(
            domain.size(), std::vector<double>(degree + 1, 0.0));

        for (std::size_t target = 0; target < domain.size(); ++target) {
            std::vector<double> coeff{1.0};
            double denominator = 1.0;
            for (std::size_t k = 0; k < domain.size(); ++k) {
                if (k == target) continue;
                denominator *= domain[target] - domain[k];

                std::vector<double> next(coeff.size() + 1, 0.0);
                for (std::size_t i = 0; i < coeff.size(); ++i) {
                    next[i] -= domain[k] * coeff[i];
                    next[i + 1] += coeff[i];
                }
                coeff = std::move(next);
            }

            for (double &value : coeff) value /= denominator;
            alpha_table[target] = std::move(coeff);
        }
        return alpha_table;
    }

    std::vector<std::vector<double>> SolveTensorLagrangeAlphaTable(
        const std::vector<std::vector<double>> &domains)
    {
        return SolveAlphaTable(BuildBasisMatrix(domains));
    }

    std::vector<seal::Ciphertext> BuildMonomialBasis(
        const seal::Ciphertext &column, std::size_t max_degree,
        const seal::RelinKeys &relin_keys, seal::Evaluator &evaluator)
    {
        std::vector<seal::Ciphertext> basis;
        if (max_degree == 0) return basis;

        basis.reserve(max_degree);
        basis.push_back(column);
        for (std::size_t degree = 2; degree <= max_degree; ++degree)
            basis.push_back(MultiplyAndRescale(basis.back(), column,
                                               relin_keys, evaluator));
        return basis;
    }

    std::vector<seal::Ciphertext> BuildTensorMonomialBasis(
        const std::vector<seal::Ciphertext> &columns,
        const std::vector<std::vector<double>> &domains,
        const seal::RelinKeys &relin_keys, seal::Evaluator &evaluator)
    {
        const auto sizes = DomainSizes(domains, "BuildTensorMonomialBasis");
        if (columns.size() != domains.size())
            throw std::invalid_argument(
                "BuildTensorMonomialBasis: column/domain count mismatch");

        const std::size_t total_size = ProductSize(sizes);
        if (total_size <= 1)
            throw std::invalid_argument(
                "BuildTensorMonomialBasis: tensor basis is constant-only");

        std::vector<std::vector<seal::Ciphertext>> powers;
        powers.reserve(columns.size());
        for (std::size_t attr = 0; attr < columns.size(); ++attr) {
            powers.push_back(BuildMonomialBasis(columns[attr],
                                                domains[attr].size() - 1,
                                                relin_keys, evaluator));
        }

        std::vector<seal::Ciphertext> basis;
        basis.reserve(total_size - 1);
        for (std::size_t basis_index = 1; basis_index < total_size;
             ++basis_index) {
            const auto exponents = DecodeMixedRadix(basis_index, sizes);

            seal::Ciphertext term;
            bool initialized = false;
            for (std::size_t attr = 0; attr < exponents.size(); ++attr) {
                const std::size_t exponent = exponents[attr];
                if (exponent == 0) continue;

                if (!initialized) {
                    term = powers[attr][exponent - 1];
                    initialized = true;
                }
                else {
                    term = MultiplyAndRescale(term,
                                              powers[attr][exponent - 1],
                                              relin_keys, evaluator);
                }
            }
            if (!initialized)
                throw std::logic_error(
                    "BuildTensorMonomialBasis: unexpected constant term");
            basis.push_back(std::move(term));
        }
        return basis;
    }

    void RotateAndSumInPlace(seal::Ciphertext &cipher,
                             std::size_t active_slots,
                             const seal::GaloisKeys &galois_keys,
                             seal::Evaluator &evaluator)
    {
        RequireActiveSlots(active_slots, "RotateAndSumInPlace");

        for (std::size_t step = 1; step < active_slots; step <<= 1) {
            seal::Ciphertext rotated;
            evaluator.rotate_vector(cipher, static_cast<int>(step),
                                    galois_keys, rotated);
            evaluator.add_inplace(cipher, rotated);
        }
    }

    seal::Ciphertext RotateAndSum(const seal::Ciphertext &cipher,
                                  std::size_t active_slots,
                                  const seal::GaloisKeys &galois_keys,
                                  seal::Evaluator &evaluator)
    {
        seal::Ciphertext result = cipher;
        RotateAndSumInPlace(result, active_slots, galois_keys, evaluator);
        return result;
    }

    std::vector<seal::Ciphertext> ReconstructMasksFromBasis(
        const std::vector<seal::Ciphertext> &basis,
        const std::vector<std::vector<double>> &alpha_table,
        seal::CKKSEncoder &encoder, seal::Evaluator &evaluator)
    {
        if (basis.empty())
            throw std::invalid_argument("ReconstructMasksFromBasis: empty basis");

        std::vector<seal::Ciphertext> masks;
        masks.reserve(alpha_table.size());
        for (const auto &row : alpha_table) {
            if (row.size() != basis.size())
                throw std::invalid_argument(
                    "ReconstructMasksFromBasis: alpha row size mismatch");

            seal::Ciphertext acc;
            bool initialized = false;
            for (std::size_t q = 0; q < basis.size(); ++q) {
                if (std::abs(row[q]) == 0.0) continue;
                seal::Ciphertext term =
                    MultiplyPlainScalar(basis[q], row[q], encoder, evaluator);
                if (!initialized) {
                    acc = std::move(term);
                    initialized = true;
                }
                else {
                    AddAlignedInPlace(acc, term, evaluator);
                }
            }
            if (!initialized)
                throw std::invalid_argument(
                    "ReconstructMasksFromBasis: zero alpha row would produce "
                    "a transparent ciphertext");
            masks.push_back(std::move(acc));
        }
        return masks;
    }

    std::vector<seal::Ciphertext> ReconstructMasksFromMonomialBasis(
        const std::vector<seal::Ciphertext> &monomial_basis,
        const std::vector<std::vector<double>> &alpha_table,
        seal::CKKSEncoder &encoder, seal::Evaluator &evaluator)
    {
        std::vector<seal::Ciphertext> masks;
        masks.reserve(alpha_table.size());
        for (const auto &row : alpha_table) {
            if (row.size() != monomial_basis.size() + 1)
                throw std::invalid_argument(
                    "ReconstructMasksFromMonomialBasis: alpha row size "
                    "mismatch");

            seal::Ciphertext acc;
            bool initialized = false;
            for (std::size_t degree = 1; degree < row.size(); ++degree) {
                if (std::abs(row[degree]) == 0.0) continue;
                seal::Ciphertext term = MultiplyPlainScalar(
                    monomial_basis[degree - 1], row[degree], encoder,
                    evaluator);
                if (!initialized) {
                    acc = std::move(term);
                    initialized = true;
                }
                else {
                    AddAlignedInPlace(acc, term, evaluator);
                }
            }
            if (!initialized)
                throw std::invalid_argument(
                    "ReconstructMasksFromMonomialBasis: constant-only "
                    "polynomial needs an encrypted-one anchor");

            AddPlainScalarInPlace(acc, row[0], encoder, evaluator);
            masks.push_back(std::move(acc));
        }
        return masks;
    }

    std::vector<seal::Ciphertext> BuildLagrangeMasks(
        const seal::Ciphertext &column, const std::vector<double> &domain,
        const seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder,
        seal::Evaluator &evaluator)
    {
        const auto alpha_table = SolveLagrangeAlphaTable(domain);
        auto basis =
            BuildMonomialBasis(column, domain.size() - 1, relin_keys,
                               evaluator);
        return ReconstructMasksFromMonomialBasis(basis, alpha_table, encoder,
                                                 evaluator);
    }

    std::vector<seal::Ciphertext> BuildTensorLagrangeMasks(
        const std::vector<seal::Ciphertext> &columns,
        const std::vector<std::vector<double>> &domains,
        const seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder,
        seal::Evaluator &evaluator)
    {
        const auto alpha_table = SolveTensorLagrangeAlphaTable(domains);
        auto basis =
            BuildTensorMonomialBasis(columns, domains, relin_keys, evaluator);
        return ReconstructMasksFromMonomialBasis(basis, alpha_table, encoder,
                                                 evaluator);
    }

    std::vector<seal::Ciphertext> GroupByCountFromEncryptedMasks(
        const std::vector<seal::Ciphertext> &group_masks,
        std::size_t active_slots, seal::CKKSEncoder &encoder,
        const seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        if (group_masks.empty())
            throw std::invalid_argument(
                "GroupByCountFromEncryptedMasks: empty group mask set");
        RequireActiveSlots(active_slots, "GroupByCountFromEncryptedMasks");

        std::vector<seal::Ciphertext> counts;
        counts.reserve(group_masks.size());
        for (const auto &mask : group_masks) {
            seal::Ciphertext active_mask =
                ApplyActiveSlotMask(mask, active_slots, encoder, evaluator);
            RotateAndSumInPlace(active_mask, active_slots, galois_keys,
                                evaluator);
            counts.push_back(std::move(active_mask));
        }
        return counts;
    }

    std::vector<seal::Ciphertext> GroupBySumFromEncryptedMasks(
        const seal::Ciphertext &values,
        const std::vector<seal::Ciphertext> &group_masks,
        std::size_t slot_count, const seal::RelinKeys &relin_keys,
        const seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        if (group_masks.empty())
            throw std::invalid_argument(
                "GroupBySumFromEncryptedMasks: empty group mask set");

        std::vector<seal::Ciphertext> sums;
        sums.reserve(group_masks.size());
        for (const auto &mask : group_masks) {
            seal::Ciphertext selected =
                MultiplyAndRescale(values, mask, relin_keys, evaluator);
            RotateAndSumInPlace(selected, slot_count, galois_keys, evaluator);
            sums.push_back(std::move(selected));
        }
        return sums;
    }

    seal::Ciphertext LookupJoinFromEncryptedMasks(
        const std::vector<seal::Ciphertext> &left_key_masks,
        const std::vector<seal::Ciphertext> &right_key_masks,
        const seal::Ciphertext &right_payload, std::size_t slot_count,
        const seal::RelinKeys &relin_keys,
        const seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator)
    {
        if (left_key_masks.empty())
            throw std::invalid_argument(
                "LookupJoinFromEncryptedMasks: empty key domain");
        if (left_key_masks.size() != right_key_masks.size())
            throw std::invalid_argument(
                "LookupJoinFromEncryptedMasks: mask domain size mismatch");

        seal::Ciphertext joined;
        bool initialized = false;
        for (std::size_t j = 0; j < left_key_masks.size(); ++j) {
            seal::Ciphertext payload =
                MultiplyAndRescale(right_payload, right_key_masks[j],
                                   relin_keys, evaluator);
            RotateAndSumInPlace(payload, slot_count, galois_keys, evaluator);

            seal::Ciphertext contribution =
                MultiplyAndRescale(left_key_masks[j], payload, relin_keys,
                                   evaluator);
            if (!initialized) {
                joined = std::move(contribution);
                initialized = true;
            }
            else {
                AddAlignedInPlace(joined, contribution, evaluator);
            }
        }
        return joined;
    }

} // namespace tfhepp_ckks
