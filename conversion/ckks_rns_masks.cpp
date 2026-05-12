#include "ckks_rns_masks.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace tfhepp_ckks
{
    namespace
    {
        void ValidateRadix(std::size_t radix, const char *name)
        {
            if (radix < 2)
                throw std::invalid_argument(std::string(name) +
                                            ": radix must be at least 2");
        }

        void ValidateConfig(const RnsMaskConfig &config, const char *name)
        {
            ValidateRadix(config.radix, name);
            if (config.digit_count == 0)
                throw std::invalid_argument(std::string(name) +
                                            ": digit_count must be positive");
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

        void SubPlainScalarInPlace(seal::Ciphertext &cipher, double scalar,
                                   seal::CKKSEncoder &encoder,
                                   seal::Evaluator &evaluator)
        {
            seal::Plaintext plain;
            encoder.encode(scalar, cipher.parms_id(), cipher.scale(), plain);
            evaluator.sub_plain_inplace(cipher, plain);
        }

        void MultiplyPlainScalarInPlace(seal::Ciphertext &cipher,
                                        double scalar,
                                        seal::CKKSEncoder &encoder,
                                        seal::Evaluator &evaluator)
        {
            if (scalar == 1.0) return;

            const double target_scale = cipher.scale();
            seal::Plaintext plain;
            encoder.encode(scalar, cipher.parms_id(), cipher.scale(), plain);
            evaluator.multiply_plain_inplace(cipher, plain);
            evaluator.rescale_to_next_inplace(cipher);
            cipher.scale() = target_scale;
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

        seal::Ciphertext DigitMinusConstant(
            const seal::Ciphertext &digit_column, std::size_t constant,
            seal::CKKSEncoder &encoder, seal::Evaluator &evaluator)
        {
            seal::Ciphertext result = digit_column;
            SubPlainScalarInPlace(result, static_cast<double>(constant),
                                  encoder, evaluator);
            return result;
        }
    } // namespace

    std::size_t RnsDigitCount(std::size_t domain_size, std::size_t radix)
    {
        ValidateRadix(radix, "RnsDigitCount");
        if (domain_size == 0)
            throw std::invalid_argument(
                "RnsDigitCount: domain_size must be positive");

        std::size_t digits = 0;
        std::size_t capacity = 1;
        while (capacity < domain_size) {
            if (capacity > static_cast<std::size_t>(-1) / radix)
                throw std::overflow_error("RnsDigitCount: radix power overflow");
            capacity *= radix;
            ++digits;
        }
        return std::max<std::size_t>(digits, 1);
    }

    std::vector<std::size_t> RnsDigits(std::size_t value,
                                       const RnsMaskConfig &config)
    {
        ValidateConfig(config, "RnsDigits");

        std::vector<std::size_t> digits(config.digit_count, 0);
        for (std::size_t i = 0; i < config.digit_count; ++i) {
            digits[i] = value % config.radix;
            value /= config.radix;
        }
        if (value != 0)
            throw std::invalid_argument(
                "RnsDigits: value does not fit in digit_count");
        return digits;
    }

    seal::Ciphertext RnsDigitEqualityMask(
        const seal::Ciphertext &digit_column, std::size_t target_digit,
        std::size_t radix, seal::CKKSEncoder &encoder,
        const seal::RelinKeys &relin_keys, seal::Evaluator &evaluator)
    {
        ValidateRadix(radix, "RnsDigitEqualityMask");
        if (target_digit >= radix)
            throw std::invalid_argument(
                "RnsDigitEqualityMask: target digit is outside radix");

        seal::Ciphertext mask;
        bool initialized = false;
        double denominator = 1.0;
        for (std::size_t k = 0; k < radix; ++k) {
            if (k == target_digit) continue;
            denominator *= static_cast<double>(
                static_cast<long long>(target_digit) -
                static_cast<long long>(k));

            seal::Ciphertext term =
                DigitMinusConstant(digit_column, k, encoder, evaluator);
            if (!initialized) {
                mask = std::move(term);
                initialized = true;
            }
            else {
                mask = MultiplyAndRescale(mask, term, relin_keys, evaluator);
            }
        }

        MultiplyPlainScalarInPlace(mask, 1.0 / denominator, encoder,
                                   evaluator);
        return mask;
    }

    seal::Ciphertext RnsEqualityMask(
        const std::vector<seal::Ciphertext> &digit_columns,
        std::size_t target_value, const RnsMaskConfig &config,
        seal::CKKSEncoder &encoder, const seal::RelinKeys &relin_keys,
        seal::Evaluator &evaluator)
    {
        ValidateConfig(config, "RnsEqualityMask");
        if (digit_columns.size() != config.digit_count)
            throw std::invalid_argument(
                "RnsEqualityMask: digit column count mismatch");

        const std::vector<std::size_t> target_digits =
            RnsDigits(target_value, config);

        seal::Ciphertext mask;
        bool initialized = false;
        for (std::size_t i = 0; i < config.digit_count; ++i) {
            seal::Ciphertext digit_mask = RnsDigitEqualityMask(
                digit_columns[i], target_digits[i], config.radix, encoder,
                relin_keys, evaluator);
            if (!initialized) {
                mask = std::move(digit_mask);
                initialized = true;
            }
            else {
                mask = MultiplyAndRescale(mask, digit_mask, relin_keys,
                                          evaluator);
            }
        }
        return mask;
    }

    std::vector<seal::Ciphertext> RnsEqualityMasks(
        const std::vector<seal::Ciphertext> &digit_columns,
        std::size_t domain_size, const RnsMaskConfig &config,
        seal::CKKSEncoder &encoder, const seal::RelinKeys &relin_keys,
        seal::Evaluator &evaluator)
    {
        ValidateConfig(config, "RnsEqualityMasks");
        if (domain_size == 0)
            throw std::invalid_argument(
                "RnsEqualityMasks: domain_size must be positive");

        std::vector<seal::Ciphertext> masks;
        masks.reserve(domain_size);
        for (std::size_t target = 0; target < domain_size; ++target)
            masks.push_back(RnsEqualityMask(digit_columns, target, config,
                                            encoder, relin_keys, evaluator));
        return masks;
    }

} // namespace tfhepp_ckks
