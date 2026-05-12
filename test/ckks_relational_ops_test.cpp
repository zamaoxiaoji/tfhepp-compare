#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "ckks_relational.h"
#include "ckks_repack.h"

namespace
{
    std::vector<int> RotationSteps(std::size_t slot_count)
    {
        std::vector<int> steps;
        for (std::size_t step = 1; step < slot_count; step <<= 1)
            steps.push_back(static_cast<int>(step));
        return steps;
    }

    seal::Ciphertext EncryptSlots(const std::vector<double> &slots,
                                  double scale, seal::CKKSEncoder &encoder,
                                  seal::Encryptor &encryptor)
    {
        seal::Plaintext plain;
        encoder.encode(slots, scale, plain);
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

    void ExpectNear(double got, double expected, double tolerance,
                    const char *label)
    {
        if (std::abs(got - expected) > tolerance) {
            std::cerr << label << ": got " << got << ", expected "
                      << expected << std::endl;
            throw std::runtime_error(label);
        }
    }

    double EvalPowerPolynomial(const std::vector<double> &coefficients,
                               double x)
    {
        double result = 0.0;
        double power = 1.0;
        for (const double coefficient : coefficients) {
            result += coefficient * power;
            power *= x;
        }
        return result;
    }

    void VerifySingleAttributeLagrange()
    {
        const std::vector<double> domain{0.0, 1.0, 2.0};
        const auto alpha = tfhepp_ckks::SolveLagrangeAlphaTable(domain);
        for (std::size_t target = 0; target < domain.size(); ++target) {
            for (std::size_t point = 0; point < domain.size(); ++point) {
                const double expected = (point == target) ? 1.0 : 0.0;
                ExpectNear(EvalPowerPolynomial(alpha[target], domain[point]),
                           expected, 1e-10, "single-attribute lagrange");
            }
        }
    }

    void VerifyTensorBasisCompiler()
    {
        const std::vector<std::vector<double>> domains{{0.0, 1.0},
                                                       {0.0, 1.0, 2.0}};
        const tfhepp_ckks::ExponentTable expected_exponents{
            {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}};
        const tfhepp_ckks::DenseMatrix expected_points{
            {0.0, 0.0}, {0.0, 1.0}, {0.0, 2.0},
            {1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0}};
        const tfhepp_ckks::DenseMatrix expected_basis{
            {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0, 0.0, 0.0, 0.0},
            {1.0, 2.0, 4.0, 0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0, 1.0, 0.0, 0.0},
            {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
            {1.0, 2.0, 4.0, 1.0, 2.0, 4.0}};

        if (tfhepp_ckks::EnumerateExponents(domains) != expected_exponents)
            throw std::runtime_error("tensor exponent order");
        if (tfhepp_ckks::EnumerateDomainPoints(domains) != expected_points)
            throw std::runtime_error("tensor point order");

        const auto basis = tfhepp_ckks::BuildBasisMatrix(domains);
        if (basis != expected_basis)
            throw std::runtime_error("tensor basis matrix");

        const auto alpha = tfhepp_ckks::SolveAlphaTable(basis);
        for (std::size_t theta = 0; theta < basis.size(); ++theta) {
            for (std::size_t row = 0; row < basis.size(); ++row) {
                double got = 0.0;
                for (std::size_t q = 0; q < basis.size(); ++q)
                    got += basis[row][q] * alpha[theta][q];
                const double expected = (row == theta) ? 1.0 : 0.0;
                ExpectNear(got, expected, 1e-10, "tensor alpha residual");
            }
        }

        const auto tensor_alpha =
            tfhepp_ckks::SolveTensorLagrangeAlphaTable(domains);
        for (std::size_t row = 0; row < alpha.size(); ++row)
            for (std::size_t col = 0; col < alpha[row].size(); ++col)
                ExpectNear(tensor_alpha[row][col], alpha[row][col], 1e-10,
                           "tensor alpha table");
    }

    std::vector<double> PlainGroupCounts(
        const std::vector<std::size_t> &keys, std::size_t domain_size)
    {
        std::vector<double> counts(domain_size, 0.0);
        for (const auto key : keys) counts[key] += 1.0;
        return counts;
    }

    std::vector<double> PlainGroupSums(
        const std::vector<std::size_t> &keys,
        const std::vector<double> &values, std::size_t domain_size)
    {
        std::vector<double> sums(domain_size, 0.0);
        for (std::size_t row = 0; row < keys.size(); ++row)
            sums[keys[row]] += values[row];
        return sums;
    }

    std::vector<double> PlainLookupJoin(
        const std::vector<std::size_t> &left_keys,
        const std::vector<std::size_t> &right_keys,
        const std::vector<double> &right_payload,
        std::size_t domain_size)
    {
        std::vector<double> payload_by_key(domain_size, 0.0);
        std::vector<bool> seen(domain_size, false);
        for (std::size_t row = 0; row < right_keys.size(); ++row) {
            const auto key = right_keys[row];
            if (key >= domain_size || seen[key])
                throw std::runtime_error("right key must be unique in domain");
            payload_by_key[key] = right_payload[row];
            seen[key] = true;
        }

        std::vector<double> result(left_keys.size(), 0.0);
        for (std::size_t row = 0; row < left_keys.size(); ++row)
            result[row] = payload_by_key[left_keys[row]];
        return result;
    }

    std::size_t FlattenCompositeKey(
        const std::vector<std::vector<std::size_t>> &columns,
        const std::vector<std::size_t> &domain_sizes, std::size_t row)
    {
        std::size_t index = 0;
        for (std::size_t attr = 0; attr < columns.size(); ++attr) {
            const auto key = columns[attr][row];
            if (key >= domain_sizes[attr])
                throw std::runtime_error("composite key is outside domain");
            index = index * domain_sizes[attr] + key;
        }
        return index;
    }

    std::vector<double> PlainCompositeGroupCounts(
        const std::vector<std::vector<std::size_t>> &columns,
        const std::vector<std::size_t> &domain_sizes)
    {
        std::size_t group_count = 1;
        for (const auto size : domain_sizes) group_count *= size;

        std::vector<double> counts(group_count, 0.0);
        for (std::size_t row = 0; row < columns.front().size(); ++row)
            counts[FlattenCompositeKey(columns, domain_sizes, row)] += 1.0;
        return counts;
    }

    std::vector<double> PlainCompositeGroupSums(
        const std::vector<std::vector<std::size_t>> &columns,
        const std::vector<double> &values,
        const std::vector<std::size_t> &domain_sizes)
    {
        std::size_t group_count = 1;
        for (const auto size : domain_sizes) group_count *= size;

        std::vector<double> sums(group_count, 0.0);
        for (std::size_t row = 0; row < values.size(); ++row)
            sums[FlattenCompositeKey(columns, domain_sizes, row)] +=
                values[row];
        return sums;
    }

    std::vector<double> PlainCompositeLookupJoin(
        const std::vector<std::vector<std::size_t>> &left_columns,
        const std::vector<std::vector<std::size_t>> &right_columns,
        const std::vector<double> &right_payload,
        const std::vector<std::size_t> &domain_sizes)
    {
        std::size_t group_count = 1;
        for (const auto size : domain_sizes) group_count *= size;

        std::vector<double> payload_by_key(group_count, 0.0);
        std::vector<bool> seen(group_count, false);
        for (std::size_t row = 0; row < right_payload.size(); ++row) {
            const auto key =
                FlattenCompositeKey(right_columns, domain_sizes, row);
            if (seen[key])
                throw std::runtime_error(
                    "right composite key must be unique");
            payload_by_key[key] = right_payload[row];
            seen[key] = true;
        }

        std::vector<double> result(left_columns.front().size(), 0.0);
        for (std::size_t row = 0; row < result.size(); ++row)
            result[row] =
                payload_by_key[FlattenCompositeKey(left_columns, domain_sizes,
                                                   row)];
        return result;
    }

    std::vector<double> SlotVector(const std::vector<double> &active_values,
                                   std::size_t slot_count)
    {
        std::vector<double> slots(slot_count, 0.0);
        for (std::size_t row = 0; row < active_values.size(); ++row)
            slots[row] = active_values[row];
        return slots;
    }

    std::vector<double> KeySlotVector(
        const std::vector<std::size_t> &active_keys, std::size_t slot_count)
    {
        std::vector<double> slots(slot_count, 0.0);
        for (std::size_t row = 0; row < active_keys.size(); ++row)
            slots[row] = static_cast<double>(active_keys[row]);
        return slots;
    }
} // namespace

int main()
{
    VerifySingleAttributeLagrange();
    VerifyTensorBasisCompiler();

    constexpr double scale = static_cast<double>(uint64_t{1} << 40);
    constexpr double tolerance = 0.05;
    constexpr std::size_t row_count = 4;
    constexpr std::size_t group_domain_size = 2;
    constexpr std::size_t join_domain_size = 3;

    auto parms = tfhepp_ckks::MakeDefaultCKKSParameters(
        16384, {60, 40, 40, 40, 40, 40, 40, 60});
    seal::SEALContext context(parms);
    seal::KeyGenerator keygen(context);
    const auto secret_key = keygen.secret_key();
    seal::PublicKey public_key;
    keygen.create_public_key(public_key);
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    seal::CKKSEncoder encoder(context);
    const std::size_t slot_count = encoder.slot_count();
    seal::GaloisKeys galois_keys;
    keygen.create_galois_keys(RotationSteps(slot_count), galois_keys);

    seal::Encryptor encryptor(context, public_key);
    seal::Decryptor decryptor(context, secret_key);
    seal::Evaluator evaluator(context);

    const std::vector<std::size_t> group_keys{0, 1, 0, 1};
    const std::vector<double> group_values_plain{10.0, 20.0, 30.0, 40.0};
    const std::vector<double> group_domain{0.0, 1.0};
    const auto expected_counts =
        PlainGroupCounts(group_keys, group_domain_size);
    const auto expected_sums =
        PlainGroupSums(group_keys, group_values_plain, group_domain_size);

    auto group_key_ct = EncryptSlots(KeySlotVector(group_keys, slot_count),
                                     scale, encoder, encryptor);
    auto value_ct = EncryptSlots(SlotVector(group_values_plain, slot_count),
                                 scale, encoder, encryptor);
    auto group_masks = tfhepp_ckks::BuildLagrangeMasks(
        group_key_ct, group_domain, relin_keys, encoder, evaluator);

    auto group_counts = tfhepp_ckks::GroupByCountFromEncryptedMasks(
        group_masks, row_count, encoder, galois_keys, evaluator);
    auto count0 = DecryptSlots(group_counts[0], decryptor, encoder);
    auto count1 = DecryptSlots(group_counts[1], decryptor, encoder);
    ExpectNear(count0[0], expected_counts[0], tolerance, "group0 count");
    ExpectNear(count1[0], expected_counts[1], tolerance, "group1 count");

    auto group_sums = tfhepp_ckks::GroupBySumFromEncryptedMasks(
        value_ct, group_masks, row_count, relin_keys, galois_keys, evaluator);
    auto sum0 = DecryptSlots(group_sums[0], decryptor, encoder);
    auto sum1 = DecryptSlots(group_sums[1], decryptor, encoder);
    ExpectNear(sum0[0], expected_sums[0], tolerance, "group0 sum");
    ExpectNear(sum1[0], expected_sums[1], tolerance, "group1 sum");

    const std::vector<std::size_t> left_keys{0, 1, 0, 2};
    const std::vector<std::size_t> right_keys{0, 1, 2};
    const std::vector<double> right_payload_plain{100.0, 200.0, 300.0};
    const std::vector<double> join_domain{0.0, 1.0, 2.0};
    const auto expected_join =
        PlainLookupJoin(left_keys, right_keys, right_payload_plain,
                        join_domain_size);
    auto left_key_ct = EncryptSlots(KeySlotVector(left_keys, slot_count),
                                    scale, encoder, encryptor);
    auto right_key_ct = EncryptSlots(KeySlotVector(right_keys, slot_count),
                                     scale, encoder, encryptor);
    auto left_masks = tfhepp_ckks::BuildLagrangeMasks(
        left_key_ct, join_domain, relin_keys, encoder, evaluator);
    auto right_masks = tfhepp_ckks::BuildLagrangeMasks(
        right_key_ct, join_domain, relin_keys, encoder, evaluator);

    auto joined = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        left_masks, right_masks,
        EncryptSlots(SlotVector(right_payload_plain, slot_count), scale,
                     encoder, encryptor),
        slot_count,
        relin_keys, galois_keys, evaluator);
    auto decoded_join = DecryptSlots(joined, decryptor, encoder);
    for (std::size_t row = 0; row < expected_join.size(); ++row)
        ExpectNear(decoded_join[row], expected_join[row], tolerance,
                   "join slot");

    const std::vector<std::size_t> gender{0, 0, 1, 1, 0, 1, 0};
    const std::vector<std::size_t> dept{0, 1, 0, 2, 2, 1, 1};
    const std::vector<double> salaries{5.0, 10.0, 20.0, 30.0,
                                       40.0, 50.0, 7.0};
    const std::vector<std::vector<std::size_t>> group_columns{gender, dept};
    const std::vector<std::size_t> multi_domain_sizes{2, 3};
    const std::vector<std::vector<double>> multi_domains{{0.0, 1.0},
                                                         {0.0, 1.0, 2.0}};
    const auto expected_multi_sums =
        PlainCompositeGroupSums(group_columns, salaries, multi_domain_sizes);
    const auto expected_multi_counts =
        PlainCompositeGroupCounts(group_columns, multi_domain_sizes);

    std::vector<seal::Ciphertext> group_column_cts{
        EncryptSlots(KeySlotVector(gender, slot_count), scale, encoder,
                     encryptor),
        EncryptSlots(KeySlotVector(dept, slot_count), scale, encoder,
                     encryptor)};
    auto multi_masks = tfhepp_ckks::BuildTensorLagrangeMasks(
        group_column_cts, multi_domains, relin_keys, encoder, evaluator);
    auto multi_counts = tfhepp_ckks::GroupByCountFromEncryptedMasks(
        multi_masks, salaries.size(), encoder, galois_keys, evaluator);
    for (std::size_t group = 0; group < expected_multi_counts.size(); ++group) {
        const auto decoded = DecryptSlots(multi_counts[group], decryptor,
                                          encoder);
        ExpectNear(decoded[0], expected_multi_counts[group], tolerance,
                   "multi group count");
    }

    auto multi_sums = tfhepp_ckks::GroupBySumFromEncryptedMasks(
        EncryptSlots(SlotVector(salaries, slot_count), scale, encoder,
                     encryptor),
        multi_masks, salaries.size(), relin_keys, galois_keys, evaluator);
    for (std::size_t group = 0; group < expected_multi_sums.size(); ++group) {
        const auto decoded = DecryptSlots(multi_sums[group], decryptor,
                                          encoder);
        ExpectNear(decoded[0], expected_multi_sums[group], tolerance,
                   "multi group sum");
    }

    const std::vector<std::size_t> left_a{0, 0, 1, 1};
    const std::vector<std::size_t> left_b{0, 2, 1, 2};
    const std::vector<std::size_t> right_a{0, 0, 1, 1};
    const std::vector<std::size_t> right_b{0, 2, 1, 2};
    const std::vector<double> composite_payload{100.0, 200.0, 300.0, 400.0};
    const std::vector<std::vector<std::size_t>> left_composite{left_a,
                                                               left_b};
    const std::vector<std::vector<std::size_t>> right_composite{right_a,
                                                                right_b};
    const auto expected_composite_join = PlainCompositeLookupJoin(
        left_composite, right_composite, composite_payload,
        multi_domain_sizes);

    std::vector<seal::Ciphertext> left_composite_cts{
        EncryptSlots(KeySlotVector(left_a, slot_count), scale, encoder,
                     encryptor),
        EncryptSlots(KeySlotVector(left_b, slot_count), scale, encoder,
                     encryptor)};
    std::vector<seal::Ciphertext> right_composite_cts{
        EncryptSlots(KeySlotVector(right_a, slot_count), scale, encoder,
                     encryptor),
        EncryptSlots(KeySlotVector(right_b, slot_count), scale, encoder,
                     encryptor)};
    auto left_composite_masks = tfhepp_ckks::BuildTensorLagrangeMasks(
        left_composite_cts, multi_domains, relin_keys, encoder, evaluator);
    auto right_composite_masks = tfhepp_ckks::BuildTensorLagrangeMasks(
        right_composite_cts, multi_domains, relin_keys, encoder, evaluator);
    auto composite_join = tfhepp_ckks::LookupJoinFromEncryptedMasks(
        left_composite_masks, right_composite_masks,
        EncryptSlots(SlotVector(composite_payload, slot_count), scale,
                     encoder, encryptor),
        slot_count, relin_keys, galois_keys, evaluator);
    const auto decoded_composite_join =
        DecryptSlots(composite_join, decryptor, encoder);
    for (std::size_t row = 0; row < expected_composite_join.size(); ++row)
        ExpectNear(decoded_composite_join[row], expected_composite_join[row],
                   tolerance, "composite join slot");

    std::cout << "ckks_relational_ops_test ok" << std::endl;
    return 0;
}
