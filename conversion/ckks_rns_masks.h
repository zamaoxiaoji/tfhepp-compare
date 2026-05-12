#pragma once

#include <cstddef>
#include <vector>

#include "seal/seal.h"

namespace tfhepp_ckks
{
    struct RnsMaskConfig {
        std::size_t radix = 0;
        std::size_t digit_count = 0;
    };

    std::size_t RnsDigitCount(std::size_t domain_size, std::size_t radix);

    std::vector<std::size_t> RnsDigits(std::size_t value,
                                       const RnsMaskConfig &config);

    seal::Ciphertext RnsDigitEqualityMask(
        const seal::Ciphertext &digit_column, std::size_t target_digit,
        std::size_t radix, seal::CKKSEncoder &encoder,
        const seal::RelinKeys &relin_keys, seal::Evaluator &evaluator);

    seal::Ciphertext RnsEqualityMask(
        const std::vector<seal::Ciphertext> &digit_columns,
        std::size_t target_value, const RnsMaskConfig &config,
        seal::CKKSEncoder &encoder, const seal::RelinKeys &relin_keys,
        seal::Evaluator &evaluator);

    std::vector<seal::Ciphertext> RnsEqualityMasks(
        const std::vector<seal::Ciphertext> &digit_columns,
        std::size_t domain_size, const RnsMaskConfig &config,
        seal::CKKSEncoder &encoder, const seal::RelinKeys &relin_keys,
        seal::Evaluator &evaluator);

} // namespace tfhepp_ckks
