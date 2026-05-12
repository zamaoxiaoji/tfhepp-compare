#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "ckks_repack.h"
#include "ckks_rns_masks.h"

namespace
{
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
} // namespace

int main()
{
    constexpr double scale = static_cast<double>(uint64_t{1} << 40);
    constexpr double tolerance = 0.08;

    auto parms = tfhepp_ckks::MakeDefaultCKKSParameters(
        16384, {60, 40, 40, 40, 40, 60});
    seal::SEALContext context(parms);
    seal::KeyGenerator keygen(context);
    const auto secret_key = keygen.secret_key();
    seal::PublicKey public_key;
    keygen.create_public_key(public_key);
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    seal::CKKSEncoder encoder(context);
    seal::Encryptor encryptor(context, public_key);
    seal::Decryptor decryptor(context, secret_key);
    seal::Evaluator evaluator(context);

    const std::size_t slot_count = encoder.slot_count();
    const tfhepp_ckks::RnsMaskConfig config{
        4, tfhepp_ckks::RnsDigitCount(16, 4)};
    if (config.digit_count != 2) return 1;

    const std::vector<std::size_t> keys{0, 6, 10, 6, 15, 1, 7, 8};
    std::vector<std::vector<double>> digit_plain(
        config.digit_count, std::vector<double>(slot_count, 0.0));
    for (std::size_t slot = 0; slot < keys.size(); ++slot) {
        const auto digits = tfhepp_ckks::RnsDigits(keys[slot], config);
        for (std::size_t i = 0; i < config.digit_count; ++i)
            digit_plain[i][slot] = static_cast<double>(digits[i]);
    }

    std::vector<seal::Ciphertext> digit_columns;
    for (const auto &digit : digit_plain)
        digit_columns.push_back(
            EncryptSlots(digit, scale, encoder, encryptor));

    auto low_digit_is_two = tfhepp_ckks::RnsDigitEqualityMask(
        digit_columns[0], 2, config.radix, encoder, relin_keys, evaluator);
    const auto decoded_digit =
        DecryptSlots(low_digit_is_two, decryptor, encoder);
    ExpectNear(decoded_digit[0], 0.0, tolerance, "digit slot0");
    ExpectNear(decoded_digit[1], 1.0, tolerance, "digit slot1");
    ExpectNear(decoded_digit[2], 1.0, tolerance, "digit slot2");
    ExpectNear(decoded_digit[3], 1.0, tolerance, "digit slot3");
    ExpectNear(decoded_digit[4], 0.0, tolerance, "digit slot4");

    auto target_six = tfhepp_ckks::RnsEqualityMask(
        digit_columns, 6, config, encoder, relin_keys, evaluator);
    const auto decoded_six = DecryptSlots(target_six, decryptor, encoder);
    ExpectNear(decoded_six[0], 0.0, tolerance, "target6 slot0");
    ExpectNear(decoded_six[1], 1.0, tolerance, "target6 slot1");
    ExpectNear(decoded_six[2], 0.0, tolerance, "target6 slot2");
    ExpectNear(decoded_six[3], 1.0, tolerance, "target6 slot3");
    ExpectNear(decoded_six[4], 0.0, tolerance, "target6 slot4");

    auto all_masks = tfhepp_ckks::RnsEqualityMasks(
        digit_columns, 16, config, encoder, relin_keys, evaluator);
    const auto decoded_fifteen =
        DecryptSlots(all_masks[15], decryptor, encoder);
    ExpectNear(decoded_fifteen[0], 0.0, tolerance, "target15 slot0");
    ExpectNear(decoded_fifteen[4], 1.0, tolerance, "target15 slot4");
    ExpectNear(decoded_fifteen[7], 0.0, tolerance, "target15 slot7");

    std::cout << "ckks_rns_masks_test ok" << std::endl;
    return 0;
}
