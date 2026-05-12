#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "key.hpp"
#include "tlwe.hpp"
#include "seal/seal.h"

namespace tfhepp_ckks
{
    struct RepackEvaluationKey {
        seal::Ciphertext encrypted_secret_key;
        std::vector<seal::Ciphertext> rotated_secret_keys;
        std::size_t tfhe_dimension = 0;
        double key_scale = 0.0;
    };

    struct RepackConfig {
        double message_scale = 0.0;
        double lwe_modulus = 0.0;
        double linear_rescale = 1.0;
        double key_scale = 0.0;
    };

    inline double Pow2(int exponent)
    {
        return std::ldexp(1.0, exponent);
    }

    template <class P>
    RepackConfig DefaultRepackConfig(uint32_t message_scale_bits = 29,
                                     uint32_t key_scale_bits = 45)
    {
        const uint32_t torus_bits =
            std::numeric_limits<typename P::T>::digits;
        return RepackConfig{
            Pow2(static_cast<int>(message_scale_bits)),
            Pow2(static_cast<int>(torus_bits)),
            Pow2(static_cast<int>(key_scale_bits) -
                 static_cast<int>(torus_bits)),
            Pow2(static_cast<int>(key_scale_bits)),
        };
    }

    std::vector<int> DefaultCoeffModulusBits();

    seal::EncryptionParameters MakeDefaultCKKSParameters(
        std::size_t poly_modulus_degree = 65536,
        const std::vector<int> &coeff_modulus_bits = DefaultCoeffModulusBits());

    void EncodeSlots(const std::vector<double> &input, double scale,
                     seal::Plaintext &plain, seal::CKKSEncoder &encoder);

    void EncodeSlots(const std::vector<double> &input,
                     seal::parms_id_type parms_id, double scale,
                     seal::Plaintext &plain, seal::CKKSEncoder &encoder);

    void DecodeSlots(std::vector<double> &result, const seal::Plaintext &plain,
                     seal::CKKSEncoder &encoder);

    void GenerateRepackKeyFromSecret(RepackEvaluationKey &eval_key,
                                     const std::vector<double> &secret_key,
                                     double key_scale,
                                     seal::CKKSEncoder &encoder,
                                     const seal::Encryptor &encryptor,
                                     const seal::SEALContext &context);

    void LinearTransform(seal::Ciphertext &result,
                         const std::vector<std::vector<double>> &matrix,
                         double scale, const RepackEvaluationKey &eval_key,
                         seal::CKKSEncoder &encoder,
                         const seal::GaloisKeys &galois_keys,
                         seal::Evaluator &evaluator);

    void PackLWECoefficientsToCKKS(
        seal::Ciphertext &result,
        const std::vector<std::vector<double>> &a_matrix,
        const std::vector<double> &b_vector, const RepackEvaluationKey &eval_key,
        const RepackConfig &config, seal::CKKSEncoder &encoder,
        const seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys,
        seal::Evaluator &evaluator, seal::SEALContext &context);

    void HomomorphicMod(seal::Ciphertext &cipher, double scale, double q0,
                        seal::CKKSEncoder &encoder,
                        seal::Evaluator &evaluator,
                        seal::RelinKeys &relin_keys,
                        seal::SEALContext &context);

    void HomomorphicRound(seal::Ciphertext &cipher, double scale,
                          seal::CKKSEncoder &encoder,
                          seal::RelinKeys &relin_keys,
                          seal::Evaluator &evaluator,
                          seal::SEALContext &context);

    namespace detail
    {
        template <class P>
        double SecretKeySlot(typename P::T value)
        {
            if constexpr (std::is_unsigned_v<typename P::T>) {
                return value > typename P::T(1) ? -1.0
                                                : static_cast<double>(value);
            }
            else {
                return static_cast<double>(value);
            }
        }

        template <class P>
        std::vector<double> SecretKeySlots(const TFHEpp::SecretKey &secret_key)
        {
            const auto key = secret_key.key.template get<P>();
            std::vector<double> slots(P::k * P::n);
            for (std::size_t i = 0; i < slots.size(); ++i)
                slots[i] = SecretKeySlot<P>(key[i]);
            return slots;
        }

        template <class P>
        void BuildLWECoefficients(
            const std::vector<TFHEpp::TLWE<P>> &lwe_ciphers,
            double linear_rescale,
            std::vector<std::vector<double>> &a_matrix,
            std::vector<double> &b_vector)
        {
            if (lwe_ciphers.empty())
                throw std::invalid_argument("PackLWEsToCKKS: empty input");

            constexpr double kHomModBound = 41.0;
            using SignedT = std::make_signed_t<typename P::T>;
            const std::size_t dimension = P::k * P::n;

            a_matrix.assign(lwe_ciphers.size(),
                            std::vector<double>(dimension, 0.0));
            b_vector.assign(lwe_ciphers.size(), 0.0);
            for (std::size_t row = 0; row < lwe_ciphers.size(); ++row) {
                const auto &lwe = lwe_ciphers[row];
                for (std::size_t col = 0; col < dimension; ++col) {
                    const auto neg_a =
                        static_cast<SignedT>(typename P::T(0) - lwe[col]);
                    a_matrix[row][col] =
                        static_cast<double>(neg_a) * linear_rescale /
                        kHomModBound;
                }
                b_vector[row] =
                    static_cast<double>(
                        static_cast<SignedT>(lwe[dimension])) *
                    linear_rescale / kHomModBound;
            }
        }
    } // namespace detail

    template <class P>
    void GenerateRepackKey(RepackEvaluationKey &eval_key,
                           const TFHEpp::SecretKey &tfhe_secret_key,
                           double key_scale, seal::CKKSEncoder &encoder,
                           const seal::Encryptor &encryptor,
                           const seal::SEALContext &context)
    {
        GenerateRepackKeyFromSecret(eval_key,
                                    detail::SecretKeySlots<P>(tfhe_secret_key),
                                    key_scale, encoder, encryptor, context);
    }

    template <class P>
    void PackLWEsToCKKS(seal::Ciphertext &result,
                        const std::vector<TFHEpp::TLWE<P>> &lwe_ciphers,
                        const RepackEvaluationKey &eval_key,
                        const RepackConfig &config,
                        seal::CKKSEncoder &encoder,
                        const seal::GaloisKeys &galois_keys,
                        seal::RelinKeys &relin_keys,
                        seal::Evaluator &evaluator,
                        seal::SEALContext &context)
    {
        std::vector<std::vector<double>> a_matrix;
        std::vector<double> b_vector;
        detail::BuildLWECoefficients<P>(lwe_ciphers, config.linear_rescale,
                                        a_matrix, b_vector);
        PackLWECoefficientsToCKKS(result, a_matrix, b_vector, eval_key, config,
                                  encoder, galois_keys, relin_keys, evaluator,
                                  context);
    }

    template <class P>
    void PackBinaryLWEsToCKKS(seal::Ciphertext &result,
                              const std::vector<TFHEpp::TLWE<P>> &lwe_ciphers,
                              const RepackEvaluationKey &eval_key,
                              const RepackConfig &config,
                              seal::CKKSEncoder &encoder,
                              const seal::GaloisKeys &galois_keys,
                              seal::RelinKeys &relin_keys,
                              seal::Evaluator &evaluator,
                              seal::SEALContext &context)
    {
        PackLWEsToCKKS(result, lwe_ciphers, eval_key, config, encoder,
                       galois_keys, relin_keys, evaluator, context);
        HomomorphicRound(result, result.scale(), encoder, relin_keys, evaluator,
                         context);
    }

} // namespace tfhepp_ckks
