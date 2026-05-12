#pragma once

#include <cstddef>
#include <vector>

#include "seal/seal.h"

namespace tfhepp_ckks
{
    using DenseMatrix = std::vector<std::vector<double>>;
    using ExponentTable = std::vector<std::vector<std::size_t>>;

    ExponentTable EnumerateExponents(
        const std::vector<std::vector<double>> &domains);

    DenseMatrix EnumerateDomainPoints(
        const std::vector<std::vector<double>> &domains);

    DenseMatrix BuildBasisMatrix(
        const std::vector<std::vector<double>> &domains);

    std::vector<std::vector<double>> SolveAlphaTable(
        const DenseMatrix &basis_matrix,
        double verification_tolerance = 1e-8);

    std::vector<std::vector<double>> SolveLagrangeAlphaTable(
        const std::vector<double> &domain);

    std::vector<std::vector<double>> SolveTensorLagrangeAlphaTable(
        const std::vector<std::vector<double>> &domains);

    std::vector<seal::Ciphertext> BuildMonomialBasis(
        const seal::Ciphertext &column, std::size_t max_degree,
        const seal::RelinKeys &relin_keys, seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> BuildTensorMonomialBasis(
        const std::vector<seal::Ciphertext> &columns,
        const std::vector<std::vector<double>> &domains,
        const seal::RelinKeys &relin_keys, seal::Evaluator &evaluator);

    // active_slots is public metadata. For non-power-of-two active regions,
    // slot 0 holds the exact aggregate when slots after active_slots are zeroed.
    void RotateAndSumInPlace(seal::Ciphertext &cipher,
                             std::size_t active_slots,
                             const seal::GaloisKeys &galois_keys,
                             seal::Evaluator &evaluator);

    seal::Ciphertext RotateAndSum(const seal::Ciphertext &cipher,
                                  std::size_t active_slots,
                                  const seal::GaloisKeys &galois_keys,
                                  seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> ReconstructMasksFromBasis(
        const std::vector<seal::Ciphertext> &basis,
        const std::vector<std::vector<double>> &alpha_table,
        seal::CKKSEncoder &encoder, seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> ReconstructMasksFromMonomialBasis(
        const std::vector<seal::Ciphertext> &monomial_basis,
        const std::vector<std::vector<double>> &alpha_table,
        seal::CKKSEncoder &encoder, seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> BuildLagrangeMasks(
        const seal::Ciphertext &column, const std::vector<double> &domain,
        const seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder,
        seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> BuildTensorLagrangeMasks(
        const std::vector<seal::Ciphertext> &columns,
        const std::vector<std::vector<double>> &domains,
        const seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder,
        seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> GroupByCountFromEncryptedMasks(
        const std::vector<seal::Ciphertext> &group_masks,
        std::size_t active_slots, seal::CKKSEncoder &encoder,
        const seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> GroupBySumFromEncryptedMasks(
        const seal::Ciphertext &values,
        const std::vector<seal::Ciphertext> &group_masks,
        std::size_t slot_count, const seal::RelinKeys &relin_keys,
        const seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);

    seal::Ciphertext LookupJoinFromEncryptedMasks(
        const std::vector<seal::Ciphertext> &left_key_masks,
        const std::vector<seal::Ciphertext> &right_key_masks,
        const seal::Ciphertext &right_payload, std::size_t slot_count,
        const seal::RelinKeys &relin_keys,
        const seal::GaloisKeys &galois_keys, seal::Evaluator &evaluator);

} // namespace tfhepp_ckks
