#include "algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::size_t NextPowerOfTwo(std::size_t x) {
    std::size_t out = 1;
    while (out < x) out <<= 1;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    std::size_t rows = 4;
    std::size_t key_domain = 4;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--rows" && i + 1 < argc) {
                rows = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--domain" && i + 1 < argc) {
                key_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap4_join [--rows N] [--domain M]");
            }
        }
        if (rows == 0 || key_domain == 0)
            throw std::invalid_argument("rows and domain must be nonzero");

        const std::size_t slots = NextPowerOfTwo(std::max(rows, key_domain));
        std::vector<double> left_key(rows);
        std::vector<double> expected(rows);
        std::vector<double> right_key(slots, 0.0);
        std::vector<double> right_payload(slots, 0.0);

        for (std::size_t j = 0; j < key_domain; j++) {
            right_key[j] = static_cast<double>(j);
            right_payload[j] = 100.0 + static_cast<double>(10 * j);
        }
        for (std::size_t i = 0; i < rows; i++) {
            const std::size_t key = (rows - 1 - i) % key_domain;
            left_key[i] = static_cast<double>(key);
            expected[i] = right_payload[key];
        }

        auto runtime = PaperReview::MakeCkksRuntime(slots, 12);
        const std::vector<PaperReview::AttributeSpec> attrs{{key_domain}};
        const auto domain = PaperReview::EnumerateDomainTuples(attrs);
        const auto basis = PaperReview::EnumerateMonomialBasis(attrs);
        const auto matrix = PaperReview::BuildBasisMatrix(domain, basis);
        const auto alpha = PaperReview::SolveCoeffTableForAllTargets(matrix);
        if (!PaperReview::VerifyAlphaTable(matrix, alpha)) {
            std::cerr << "AlphaTable verification failed\n";
            return 1;
        }

        const auto ct_left_key = PaperReview::EncryptVector(
            runtime.cc, runtime.keys.publicKey, left_key, slots);
        const auto ct_right_key = PaperReview::EncryptVector(
            runtime.cc, runtime.keys.publicKey, right_key, slots);
        const auto ct_payload = PaperReview::EncryptVector(
            runtime.cc, runtime.keys.publicKey, right_payload, slots);
        const auto joined = PaperReview::LookupJoin(
            runtime.cc, runtime.keys.publicKey, ct_left_key, {},
            ct_right_key, {ct_payload}, basis, alpha, slots);

        lbcrypto::Plaintext plain;
        runtime.cc->Decrypt(runtime.keys.secretKey, joined.front(), &plain);
        plain->SetLength(slots);
        const auto got = plain->GetRealPackedValue();

        double max_err = 0.0;
        std::cout << "rows=" << rows << " domain=" << key_domain
                  << " slots=" << slots << "\n";
        std::cout << "slot,left_key,expected_payload,encrypted_payload,abs_error\n";
        for (std::size_t i = 0; i < rows; i++) {
            const double err = std::abs(expected[i] - got[i]);
            max_err = std::max(max_err, err);
            std::cout << i << "," << left_key[i] << "," << expected[i]
                      << "," << got[i] << "," << err << "\n";
        }
        std::cout << "max_abs_error=" << max_err << "\n";
        return max_err < 0.5 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap4_join failed: " << e.what() << "\n";
        return 1;
    }
}
