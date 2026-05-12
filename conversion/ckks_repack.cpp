#include "ckks_repack.h"

#include <numbers>

namespace tfhepp_ckks
{
    namespace
    {
        constexpr double kHomModBound = 41.0;

        template <typename T>
        T CeilSqrt(T value)
        {
            return static_cast<T>(std::ceil(std::sqrt(static_cast<double>(value))));
        }

        template <typename T>
        T CeilDiv(T a, T b)
        {
            return (a + b - 1) / b;
        }

        template <typename T>
        uint64_t CeilLog2(T value)
        {
            if (value <= 1) return 0;
            return static_cast<uint64_t>(
                std::ceil(std::log2(static_cast<double>(value))));
        }

        bool IsPowerOfTwo(std::size_t value)
        {
            return value != 0 && (value & (value - 1)) == 0;
        }

        void AddScalar(seal::Ciphertext &result, double scalar,
                       seal::CKKSEncoder &encoder,
                       seal::Evaluator &evaluator)
        {
            seal::Plaintext plain;
            encoder.encode(scalar, result.parms_id(), result.scale(), plain);
            evaluator.add_plain_inplace(result, plain);
        }

        void MultiplyScalar(seal::Ciphertext &result, double scalar,
                            double scale, seal::CKKSEncoder &encoder,
                            seal::Evaluator &evaluator)
        {
            seal::Plaintext plain;
            encoder.encode(scalar, result.parms_id(), scale, plain);
            evaluator.multiply_plain_inplace(result, plain);
        }

        void MultiplyAndRelinearize(const seal::Ciphertext &lhs_in,
                                    const seal::Ciphertext &rhs_in,
                                    seal::Ciphertext &result,
                                    seal::Evaluator &evaluator,
                                    seal::RelinKeys &relin_keys)
        {
            seal::Ciphertext lhs = lhs_in;
            seal::Ciphertext rhs = rhs_in;
            if (lhs.coeff_modulus_size() > rhs.coeff_modulus_size())
                evaluator.mod_switch_to_inplace(lhs, rhs.parms_id());
            else
                evaluator.mod_switch_to_inplace(rhs, lhs.parms_id());
            evaluator.multiply(lhs, rhs, result);
            evaluator.relinearize_inplace(result, relin_keys);
        }

        struct ChebyshevPoly {
            std::vector<double> coefficients;

            ChebyshevPoly() = default;
            explicit ChebyshevPoly(std::size_t degree)
                : coefficients(degree + 1, 0.0)
            {}

            bool empty() const { return coefficients.empty(); }
            uint64_t degree() const { return coefficients.size() - 1; }
        };

        const std::vector<double> &SinApproxCoefficients()
        {
            static const std::vector<double> coeff{
                0.07043813426689680968, 0.00000000000000021373,
                0.13648426242848080148, -0.00000000000000065511,
                0.12251850008297035521, 0.00000000000000138338,
                0.09687911813122730698, -0.00000000000000086808,
                0.05705995636145275163, -0.00000000000000031285,
                0.00190187312696625684, 0.00000000000000044354,
                -0.06571531807303646056, 0.00000000000000006187,
                -0.13574001749816311246, 0.00000000000000010202,
                -0.18904743240584598318, 0.00000000000000071687,
                -0.19987271926250091414, -0.00000000000000041061,
                -0.14604929968900889903, -0.00000000000000043404,
                -0.02740391745089436684, 0.00000000000000064004,
                0.11591361998332121164, -0.00000000000000074717,
                0.20462201735538801395, -0.00000000000000080824,
                0.16189767223136858343, -0.00000000000000026978,
                -0.01077097508287188157, -0.00000000000000003535,
                -0.18568750764013738919, 0.00000000000000053534,
                -0.18278566848347499452, -0.00000000000000019242,
                0.03006119366883018562, -0.00000000000000096166,
                0.21645511738551495573, 0.00000000000000012574,
                0.10356123542751023703, 0.00000000000000037192,
                -0.17891441905512733834, 0.00000000000000035345,
                -0.16355634150343503763, -0.00000000000000053488,
                0.16482549804489621259, -0.00000000000000118865,
                0.16413921800264116846, 0.00000000000000026596,
                -0.20888412226960773044, 0.00000000000000142608,
                -0.08344642688900480443, -0.00000000000000024583,
                0.26869912801153833515, 0.00000000000000000804,
                -0.13548482981129844616, 0.00000000000000109860,
                -0.13729914741849616466, 0.00000000000000277430,
                0.31392945514509801308, 0.00000000000000319002,
                -0.32761217572456985403, 0.00000000000000107757,
                0.24389625486353291861, 0.00000000000000062732,
                -0.14493001079259371089, 0.00000000000000155720,
                0.07233834282998567733, -0.00000000000000110679,
                -0.03123261034259527599, -0.00000000000000262276,
                0.01189279261890215621, -0.00000000000000323379,
                -0.00404956978046461594, -0.00000000000000441352,
                0.00124605558777123146, -0.00000000000000212377,
                -0.00034935166504105009, -0.00000000000000081908,
                0.00008984823446769276, -0.00000000000000558015,
                -0.00002131648317308570, 0.00000000000000019052,
                0.00000468762132811721, 0.00000000000000103557,
                -0.00000095941952191024, 0.00000000000000107394,
                0.00000018342015938494, 0.00000000000000161096,
                -0.00000003285866440681, 0.00000000000000201318,
                0.00000000553158897361, 0.00000000000000262726,
                -0.00000000087731624365, 0.00000000000000194294,
                0.00000000013139619346, 0.00000000000000145072,
                -0.00000000001866069215, -0.00000000000000088172,
                0.00000000000282372220};
            return coeff;
        }

        ChebyshevPoly GenerateModPoly(uint32_t r)
        {
            ChebyshevPoly poly;
            poly.coefficients = SinApproxCoefficients();
            const double cnst_scale =
                std::pow(0.5 / std::numbers::pi_v<double>,
                         1.0 / static_cast<double>(uint64_t{1} << r));
            for (double &c : poly.coefficients) c *= cnst_scale;
            return poly;
        }

        void EvalChebyshevBasic(std::vector<seal::Ciphertext> &t1,
                                std::vector<seal::Ciphertext> &t2,
                                seal::Ciphertext &x, uint64_t degree,
                                seal::CKKSEncoder &encoder,
                                seal::Evaluator &evaluator,
                                seal::RelinKeys &relin_keys)
        {
            if (degree < 1)
                throw std::invalid_argument(
                    "EvalChebyshevBasic: degree is less than 1");

            const uint64_t m = CeilLog2(degree + 1);
            const uint64_t l = m / 2;
            t1.resize((uint64_t{1} << l) - 1);

            t1[0] = x;
            for (std::size_t i = 2; i < (uint64_t{1} << l); ++i) {
                const uint64_t a = (i + 1) / 2;
                const uint64_t b = i - a;
                const uint64_t c = a > b ? a - b : b - a;
                MultiplyAndRelinearize(t1[a - 1], t1[b - 1], t1[i - 1],
                                        evaluator, relin_keys);
                evaluator.add_inplace(t1[i - 1], t1[i - 1]);

                if (c == 0) {
                    AddScalar(t1[i - 1], -1.0, encoder, evaluator);
                    evaluator.rescale_to_next_inplace(t1[i - 1]);
                }
                else {
                    seal::Ciphertext temp = t1[0];
                    MultiplyScalar(temp, 1.0,
                                   t1[i - 1].scale() / t1[0].scale(), encoder,
                                   evaluator);
                    evaluator.mod_switch_to_inplace(temp,
                                                    t1[i - 1].parms_id());
                    evaluator.sub_inplace(t1[i - 1], temp);
                    evaluator.rescale_to_next_inplace(t1[i - 1]);
                }
            }

            t2.resize(m - l);
            seal::Ciphertext temp = t1[(uint64_t{1} << (l - 1)) - 1];
            evaluator.mod_switch_to_inplace(temp,
                                            t1[(uint64_t{1} << l) - 2]
                                                .parms_id());
            for (std::size_t i = 0; i < m - l; ++i) {
                MultiplyAndRelinearize(temp, temp, t2[i], evaluator,
                                        relin_keys);
                evaluator.add_inplace(t2[i], t2[i]);
                AddScalar(t2[i], -1.0, encoder, evaluator);
                evaluator.rescale_to_next_inplace(t2[i]);
                temp = t2[i];
            }
        }

        void DivisionChebyshevLazy(const ChebyshevPoly &f,
                                   uint64_t chebyshev_degree,
                                   ChebyshevPoly &quotient,
                                   ChebyshevPoly &remainder)
        {
            if (f.degree() < chebyshev_degree) {
                quotient = ChebyshevPoly();
                quotient.coefficients.resize(1, 0.0);
                remainder = f;
                return;
            }

            remainder = ChebyshevPoly(chebyshev_degree - 1);
            std::copy_n(f.coefficients.data(), chebyshev_degree,
                        remainder.coefficients.data());
            quotient = ChebyshevPoly(f.degree() - chebyshev_degree);
            quotient.coefficients[0] = f.coefficients[chebyshev_degree];

            for (std::size_t i = chebyshev_degree + 1, j = 1;
                 i <= f.degree(); ++i, ++j) {
                quotient.coefficients[i - chebyshev_degree] =
                    2.0 * f.coefficients[i];
                remainder.coefficients[chebyshev_degree - j] -=
                    f.coefficients[i];
            }
        }

        void EvalRecurseLazy(double target_scale, uint64_t m, uint64_t l,
                             const ChebyshevPoly &poly,
                             std::vector<seal::Ciphertext> &t1,
                             std::vector<seal::Ciphertext> &t2,
                             seal::Ciphertext &result,
                             seal::SEALContext &context,
                             seal::CKKSEncoder &encoder,
                             seal::Evaluator &evaluator,
                             seal::RelinKeys &relin_keys)
        {
            if (poly.empty()) {
                result.release();
                return;
            }

            const uint64_t d = poly.degree();
            if (d < (uint64_t{1} << l)) {
                if (d == 0) {
                    if (std::abs(std::round(poly.coefficients[0] *
                                            target_scale)) > 1.0) {
                        result.release();
                        result.resize(context, t1[0].parms_id(), 2);
                        result.scale() = target_scale;
                        result.is_ntt_form() = true;
                        AddScalar(result, poly.coefficients[0], encoder,
                                  evaluator);
                    }
                    else {
                        result.release();
                    }
                    return;
                }

                auto context_data = context.get_context_data(t1[d - 1].parms_id());
                if (!context_data)
                    throw std::invalid_argument(
                        "EvalRecurseLazy: invalid ciphertext parameter id");
                const auto &parms = context_data->parms();
                const double q_t_d =
                    parms.coeff_modulus()[t1[d - 1].coeff_modulus_size() - 1]
                        .value();
                const double cipher_scale = target_scale * q_t_d;
                seal::Ciphertext temp;
                bool all_zero = true;
                result.release();
                result.resize(context, t1[d - 1].parms_id(), 2);
                result.scale() = cipher_scale;
                result.is_ntt_form() = true;

                for (std::size_t i = d; i > 0; --i) {
                    if (std::abs(std::round(poly.coefficients[i] *
                                            target_scale)) > 1.0) {
                        temp = t1[i - 1];
                        evaluator.mod_switch_to_inplace(temp,
                                                        t1[d - 1].parms_id());
                        MultiplyScalar(temp, poly.coefficients[i],
                                       cipher_scale / t1[i - 1].scale(),
                                       encoder, evaluator);
                        evaluator.add_inplace(result, temp);
                        all_zero = false;
                    }
                }

                if (std::abs(std::round(poly.coefficients[0] *
                                        target_scale)) > 1.0) {
                    AddScalar(result, poly.coefficients[0], encoder, evaluator);
                    all_zero = false;
                }

                if (all_zero) {
                    result.release();
                }
                else {
                    evaluator.rescale_to_next_inplace(result);
                }
                return;
            }

            const uint64_t chebyshev_degree = uint64_t{1} << (m - 1);
            const uint64_t chebyshev_index = m - 1 - l;
            ChebyshevPoly quotient, remainder;
            DivisionChebyshevLazy(poly, chebyshev_degree, quotient, remainder);

            const std::size_t level = t2[chebyshev_index].coeff_modulus_size() - 1;
            auto context_data = context.get_context_data(context.first_parms_id());
            if (!context_data)
                throw std::invalid_argument(
                    "EvalRecurseLazy: invalid first parameter id");
            const auto &parms = context_data->parms();
            if (level >= parms.coeff_modulus().size())
                throw std::invalid_argument("EvalRecurseLazy: level overflow");
            const double q = parms.coeff_modulus()[level].value();
            const double quotient_scale =
                target_scale * q / t2[chebyshev_index].scale();

            seal::Ciphertext cipher_quotient, cipher_remainder;
            EvalRecurseLazy(target_scale, m - 1, l, remainder, t1, t2,
                            cipher_remainder, context, encoder, evaluator,
                            relin_keys);
            EvalRecurseLazy(quotient_scale, m - 1, l, quotient, t1, t2,
                            cipher_quotient, context, encoder, evaluator,
                            relin_keys);

            if (cipher_quotient.size() && cipher_remainder.size()) {
                MultiplyAndRelinearize(cipher_quotient, t2[chebyshev_index],
                                        result, evaluator, relin_keys);
                if (result.coeff_modulus_size() <=
                    cipher_remainder.coeff_modulus_size()) {
                    MultiplyScalar(cipher_remainder, 1.0,
                                   result.scale() /
                                       cipher_remainder.scale(),
                                   encoder, evaluator);
                    evaluator.mod_switch_to_inplace(cipher_remainder,
                                                    result.parms_id());
                    evaluator.add_inplace(result, cipher_remainder);
                    evaluator.rescale_to_next_inplace(result);
                }
                else {
                    evaluator.rescale_to_next_inplace(result);
                    if (!seal::util::are_close(result.scale(),
                                               cipher_remainder.scale()))
                        result.scale() = cipher_remainder.scale();
                    evaluator.mod_switch_to_inplace(result,
                                                    cipher_remainder.parms_id());
                    evaluator.add_inplace(result, cipher_remainder);
                }
            }
            else if (cipher_quotient.size()) {
                MultiplyAndRelinearize(cipher_quotient, t2[chebyshev_index],
                                        result, evaluator, relin_keys);
                evaluator.rescale_to_next_inplace(result);
            }
            else if (cipher_remainder.size()) {
                result = cipher_remainder;
            }
            else {
                result.release();
            }
        }

        void PolyEvaluateBSGSLazy(seal::Ciphertext &result,
                                  double target_scale, seal::Ciphertext &x,
                                  const ChebyshevPoly &poly,
                                  seal::SEALContext &context,
                                  seal::CKKSEncoder &encoder,
                                  seal::Evaluator &evaluator,
                                  seal::RelinKeys &relin_keys)
        {
            const uint64_t degree = poly.degree();
            const uint64_t m = CeilLog2(degree + 1);
            const uint64_t l = m / 2;
            std::vector<seal::Ciphertext> t1, t2;
            EvalChebyshevBasic(t1, t2, x, degree, encoder, evaluator,
                               relin_keys);
            EvalRecurseLazy(target_scale, m, l, poly, t1, t2, result, context,
                            encoder, evaluator, relin_keys);
        }

        void EvalPower(double target_scale,
                       const std::vector<double> &coefficients,
                       std::vector<seal::Ciphertext> &power_basis,
                       seal::Ciphertext &result, seal::SEALContext &context,
                       seal::CKKSEncoder &encoder,
                       seal::Evaluator &evaluator,
                       seal::RelinKeys &relin_keys)
        {
            if (coefficients.size() == 1) {
                if (std::abs(std::round(coefficients[0] * target_scale)) >
                    1.0) {
                    auto context_data =
                        context.get_context_data(power_basis[0].parms_id());
                    if (!context_data)
                        throw std::invalid_argument(
                            "EvalPower: invalid ciphertext parameter id");
                    const auto &parms = context_data->parms();
                    const double q_t_d =
                        parms
                            .coeff_modulus()[power_basis[0]
                                                 .coeff_modulus_size() -
                                             1]
                            .value();
                    const double cipher_scale = target_scale * q_t_d;
                    result = power_basis[0];
                    MultiplyScalar(result, coefficients[0],
                                   cipher_scale / power_basis[0].scale(),
                                   encoder, evaluator);
                    evaluator.rescale_to_next_inplace(result);
                    return;
                }
                result.release();
                return;
            }

            const uint64_t degree = coefficients.size() * 2 - 1;
            const uint64_t m = CeilLog2(degree + 1);
            std::vector<double> remainder((uint64_t{1} << (m - 1)) / 2);
            std::vector<double> quotient(coefficients.size() -
                                         remainder.size());
            for (std::size_t i = 0; i < remainder.size(); ++i)
                remainder[i] = coefficients[i];
            for (std::size_t i = 0; i < quotient.size(); ++i)
                quotient[i] = coefficients[i + remainder.size()];

            auto context_data =
                context.get_context_data(power_basis[m - 1].parms_id());
            if (!context_data)
                throw std::invalid_argument(
                    "EvalPower: invalid basis parameter id");
            const auto &parms = context_data->parms();
            const double q_t_d =
                parms
                    .coeff_modulus()[power_basis[m - 1].coeff_modulus_size() -
                                     1]
                    .value();
            seal::Ciphertext cipher_quotient, cipher_remainder;
            const double quotient_scale =
                target_scale * q_t_d / power_basis[m - 1].scale();
            EvalPower(quotient_scale, quotient, power_basis, cipher_quotient,
                      context, encoder, evaluator, relin_keys);
            EvalPower(target_scale, remainder, power_basis, cipher_remainder,
                      context, encoder, evaluator, relin_keys);
            MultiplyAndRelinearize(cipher_quotient, power_basis[m - 1], result,
                                    evaluator, relin_keys);
            evaluator.rescale_to_next_inplace(result);
            evaluator.mod_switch_to_inplace(cipher_remainder,
                                            result.parms_id());
            cipher_remainder.scale() = result.scale();
            evaluator.add_inplace(result, cipher_remainder);
        }

        void PolyEvaluatePower(seal::Ciphertext &result, double target_scale,
                               seal::Ciphertext &x,
                               const std::vector<double> &coefficients,
                               seal::SEALContext &context,
                               seal::CKKSEncoder &encoder,
                               seal::Evaluator &evaluator,
                               seal::RelinKeys &relin_keys)
        {
            const uint64_t degree = coefficients.size() * 2 - 1;
            const uint64_t m = CeilLog2(degree + 1);
            std::vector<seal::Ciphertext> power_basis(m);
            power_basis[0] = x;
            for (std::size_t i = 1; i < m; ++i) {
                MultiplyAndRelinearize(power_basis[i - 1], power_basis[i - 1],
                                        power_basis[i], evaluator, relin_keys);
                evaluator.rescale_to_next_inplace(power_basis[i]);
            }
            EvalPower(target_scale, coefficients, power_basis, result, context,
                      encoder, evaluator, relin_keys);
        }
    } // namespace

    std::vector<int> DefaultCoeffModulusBits()
    {
        return {45, 42, 42, 42, 42, 42, 42, 45, 45, 45,
                45, 45, 45, 45, 45, 45, 45, 45, 59};
    }

    seal::EncryptionParameters MakeDefaultCKKSParameters(
        std::size_t poly_modulus_degree,
        const std::vector<int> &coeff_modulus_bits)
    {
        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(
            seal::CoeffModulus::Create(poly_modulus_degree,
                                       coeff_modulus_bits));
        return parms;
    }

    void EncodeSlots(const std::vector<double> &input, double scale,
                     seal::Plaintext &plain, seal::CKKSEncoder &encoder)
    {
        const std::size_t slot_count = encoder.slot_count();
        if (input.empty() || input.size() > slot_count)
            throw std::invalid_argument("EncodeSlots: invalid input size");

        std::vector<double> plain_input(slot_count, 0.0);
        for (std::size_t i = 0; i < slot_count; ++i)
            plain_input[i] = input[i % input.size()];
        encoder.encode(plain_input, scale, plain);
    }

    void EncodeSlots(const std::vector<double> &input,
                     seal::parms_id_type parms_id, double scale,
                     seal::Plaintext &plain, seal::CKKSEncoder &encoder)
    {
        const std::size_t slot_count = encoder.slot_count();
        if (input.empty() || input.size() > slot_count)
            throw std::invalid_argument("EncodeSlots: invalid input size");

        std::vector<double> plain_input(slot_count, 0.0);
        for (std::size_t i = 0; i < slot_count; ++i)
            plain_input[i] = input[i % input.size()];
        encoder.encode(plain_input, parms_id, scale, plain);
    }

    void DecodeSlots(std::vector<double> &result, const seal::Plaintext &plain,
                     seal::CKKSEncoder &encoder)
    {
        const std::size_t slot_count = encoder.slot_count();
        if (result.size() > slot_count)
            throw std::invalid_argument("DecodeSlots: output is too large");

        std::vector<double> plain_output(slot_count, 0.0);
        encoder.decode(plain, plain_output);
        for (std::size_t i = 0; i < result.size(); ++i)
            result[i] = plain_output[i];
    }

    void GenerateRepackKeyFromSecret(RepackEvaluationKey &eval_key,
                                     const std::vector<double> &secret_key,
                                     double key_scale,
                                     seal::CKKSEncoder &encoder,
                                     const seal::Encryptor &encryptor,
                                     const seal::SEALContext &context)
    {
        const std::size_t slot_count = encoder.slot_count();
        if (secret_key.empty() || secret_key.size() > slot_count)
            throw std::invalid_argument(
                "GenerateRepackKeyFromSecret: invalid secret-key size");

        std::vector<double> slots = secret_key;
        seal::Plaintext plain_sk;
        plain_sk.parms_id() = context.first_parms_id();
        EncodeSlots(slots, key_scale, plain_sk, encoder);
        encryptor.encrypt_symmetric(plain_sk, eval_key.encrypted_secret_key);

        const std::size_t baby_steps = CeilSqrt(secret_key.size());
        eval_key.rotated_secret_keys.resize(baby_steps);
        eval_key.rotated_secret_keys[0] = eval_key.encrypted_secret_key;
        for (std::size_t j = 1; j < baby_steps; ++j) {
            std::rotate(slots.begin(), slots.begin() + 1, slots.end());
            EncodeSlots(slots, key_scale, plain_sk, encoder);
            encryptor.encrypt_symmetric(plain_sk,
                                        eval_key.rotated_secret_keys[j]);
        }

        seal::util::seal_memzero(slots.data(),
                                 slots.size() * sizeof(slots[0]));
        eval_key.tfhe_dimension = secret_key.size();
        eval_key.key_scale = key_scale;
    }

    void LinearTransform(seal::Ciphertext &result,
                         const std::vector<std::vector<double>> &matrix,
                         double scale, const RepackEvaluationKey &eval_key,
                         seal::CKKSEncoder &encoder,
                         const seal::GaloisKeys &galois_keys,
                         seal::Evaluator &evaluator)
    {
        if (matrix.empty() || matrix.front().empty())
            throw std::invalid_argument("LinearTransform: empty matrix");

        const std::size_t rows = matrix.size();
        const std::size_t columns = matrix.front().size();
        for (const auto &row : matrix)
            if (row.size() != columns)
                throw std::invalid_argument("LinearTransform: ragged matrix");
        if (columns != eval_key.tfhe_dimension)
            throw std::invalid_argument(
                "LinearTransform: matrix/key dimension mismatch");

        const std::size_t slot_count = encoder.slot_count();
        const std::size_t max_len = std::max(rows, columns);
        const std::size_t min_len = std::min(rows, columns);
        if (max_len > slot_count)
            throw std::invalid_argument(
                "LinearTransform: matrix exceeds CKKS slot count");

        const std::size_t g_tilde = CeilSqrt(min_len);
        const std::size_t b_tilde = CeilDiv(min_len, g_tilde);
        if (eval_key.rotated_secret_keys.size() < g_tilde)
            throw std::invalid_argument(
                "LinearTransform: repack key has too few baby-step rotations");

        std::vector<double> diag(max_len, 0.0);
        seal::Plaintext plain;
        plain.parms_id() = eval_key.rotated_secret_keys[0].parms_id();
        for (std::size_t b = 0; b < b_tilde && g_tilde * b < min_len; ++b) {
            seal::Ciphertext sum;
            for (std::size_t g = 0;
                 g < g_tilde && b * g_tilde + g < min_len; ++g) {
                const std::size_t j = b * g_tilde + g;
                for (std::size_t r = 0; r < max_len; ++r)
                    diag[r] = matrix[r % rows][(r + j) % columns];
                std::rotate(diag.rbegin(), diag.rbegin() + b * g_tilde,
                            diag.rend());
                EncodeSlots(diag, scale, plain, encoder);
                if (g == 0) {
                    evaluator.multiply_plain(eval_key.rotated_secret_keys[g],
                                             plain, sum);
                }
                else {
                    seal::Ciphertext temp;
                    evaluator.multiply_plain(eval_key.rotated_secret_keys[g],
                                             plain, temp);
                    evaluator.add_inplace(sum, temp);
                }
            }
            if (b == 0) {
                result = sum;
            }
            else {
                evaluator.rotate_vector_inplace(sum,
                                                static_cast<int>(b * g_tilde),
                                                galois_keys);
                evaluator.add_inplace(result, sum);
            }
        }

        if (rows < columns) {
            if (columns % rows != 0 || !IsPowerOfTwo(columns / rows))
                throw std::invalid_argument(
                    "LinearTransform: rows must divide columns by a power of 2");
            const std::size_t rounds =
                static_cast<std::size_t>(std::log2(columns / rows));
            for (std::size_t j = 0; j < rounds; ++j) {
                seal::Ciphertext temp = result;
                evaluator.rotate_vector_inplace(
                    temp, static_cast<int>((std::size_t{1} << j) * rows),
                    galois_keys);
                evaluator.add_inplace(result, temp);
            }
        }
    }

    void PackLWECoefficientsToCKKS(
        seal::Ciphertext &result,
        const std::vector<std::vector<double>> &a_matrix,
        const std::vector<double> &b_vector, const RepackEvaluationKey &eval_key,
        const RepackConfig &config, seal::CKKSEncoder &encoder,
        const seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys,
        seal::Evaluator &evaluator, seal::SEALContext &context)
    {
        if (a_matrix.empty() || b_vector.size() != a_matrix.size())
            throw std::invalid_argument(
                "PackLWECoefficientsToCKKS: invalid LWE coefficient sizes");
        if (config.message_scale <= 0.0 || config.lwe_modulus <= 0.0 ||
            config.linear_rescale <= 0.0)
            throw std::invalid_argument(
                "PackLWECoefficientsToCKKS: invalid repack config");

        LinearTransform(result, a_matrix, 1.0, eval_key, encoder, galois_keys,
                        evaluator);
        evaluator.rescale_to_next_inplace(result);
        result.scale() = 1.0;

        seal::Plaintext plain;
        plain.parms_id() = result.parms_id();
        EncodeSlots(b_vector, result.parms_id(), 1.0, plain, encoder);
        evaluator.add_plain_inplace(result, plain);

        result.scale() = config.message_scale * config.linear_rescale;
        HomomorphicMod(result, config.message_scale * config.linear_rescale,
                       config.lwe_modulus * config.linear_rescale, encoder,
                       evaluator, relin_keys, context);
    }

    void HomomorphicMod(seal::Ciphertext &cipher, double scale, double q0,
                        seal::CKKSEncoder &encoder,
                        seal::Evaluator &evaluator,
                        seal::RelinKeys &relin_keys,
                        seal::SEALContext &context)
    {
        constexpr uint32_t r = 2;
        ChebyshevPoly poly = GenerateModPoly(r);
        const uint32_t depth =
            static_cast<uint32_t>(CeilLog2(poly.degree())) + 1 + r;
        if (cipher.coeff_modulus_size() <= depth + 1)
            throw std::invalid_argument("HomomorphicMod: level is too small");

        auto context_data = context.get_context_data(cipher.parms_id());
        if (!context_data)
            throw std::invalid_argument(
                "HomomorphicMod: invalid ciphertext parameter id");
        const auto &parms = context_data->parms();

        double target_scale = cipher.scale() * std::round(q0 / scale);
        cipher.scale() = target_scale;
        const std::size_t output_level =
            cipher.coeff_modulus_size() - 1 - depth;

        for (std::size_t i = 1; i <= r; ++i) {
            const uint64_t qi =
                parms.coeff_modulus()[output_level + i].value();
            target_scale *= static_cast<double>(qi);
            target_scale = std::sqrt(target_scale);
        }
        cipher.scale() = target_scale;

        AddScalar(cipher, -0.25 / kHomModBound, encoder, evaluator);
        seal::Ciphertext cipher_result;
        PolyEvaluateBSGSLazy(cipher_result, target_scale, cipher, poly, context,
                             encoder, evaluator, relin_keys);
        if (!seal::util::are_close(target_scale, cipher_result.scale()))
            throw std::invalid_argument("HomomorphicMod: scale mismatch");

        double theta =
            std::pow(0.5 / std::numbers::pi_v<double>,
                     1.0 / static_cast<double>(uint64_t{1} << r));
        cipher = cipher_result;
        for (std::size_t i = 0; i < r; ++i) {
            theta *= theta;
            MultiplyAndRelinearize(cipher, cipher, cipher, evaluator,
                                    relin_keys);
            evaluator.add_inplace(cipher, cipher);
            AddScalar(cipher, -theta, encoder, evaluator);
            evaluator.rescale_to_next_inplace(cipher);
        }

        cipher.scale() /= std::round(q0 / scale);
    }

    void HomomorphicRound(seal::Ciphertext &cipher, double scale,
                          seal::CKKSEncoder &encoder,
                          seal::RelinKeys &relin_keys,
                          seal::Evaluator &evaluator,
                          seal::SEALContext &context)
    {
        static const std::vector<double> coeff1{1.5, -0.5};
        static const std::vector<double> coeff5{
            2.4609375 / 2.0, -3.28125 / 2.0, 2.953125 / 2.0,
            -1.40625 / 2.0, 0.2734375 / 2.0};

        MultiplyScalar(cipher, 2.0, 1.0, encoder, evaluator);
        AddScalar(cipher, -1.0, encoder, evaluator);
        PolyEvaluatePower(cipher, scale, cipher, coeff1, context, encoder,
                          evaluator, relin_keys);
        PolyEvaluatePower(cipher, scale, cipher, coeff5, context, encoder,
                          evaluator, relin_keys);
        AddScalar(cipher, 0.5, encoder, evaluator);
    }

} // namespace tfhepp_ckks
