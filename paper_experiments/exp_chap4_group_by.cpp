#include "algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::size_t NextPowerOfTwo(std::size_t x) {
    std::size_t out = 1;
    while (out < x) out <<= 1;
    return out;
}

std::vector<double> DecryptSlot0(
    const PaperReview::CkksRuntime& runtime,
    const std::vector<PaperReview::CkksCiphertext>& cts) {
    std::vector<double> out;
    out.reserve(cts.size());
    for (const auto& ct : cts) {
        lbcrypto::Plaintext plain;
        runtime.cc->Decrypt(runtime.keys.secretKey, ct, &plain);
        plain->SetLength(runtime.slots);
        out.push_back(plain->GetRealPackedValue().front());
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    std::size_t rows = 8;
    std::size_t groups = 4;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--rows" && i + 1 < argc) {
                rows = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--groups" && i + 1 < argc) {
                groups = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap4_group_by [--rows N] [--groups M]");
            }
        }
        if (rows == 0 || groups == 0)
            throw std::invalid_argument("rows and groups must be nonzero");

        const std::size_t slots = NextPowerOfTwo(rows);
        std::vector<double> group_key(rows), value(rows);
        std::vector<double> plain(groups, 0.0);
        for (std::size_t i = 0; i < rows; i++) {
            group_key[i] = static_cast<double>(i % groups);
            value[i] = static_cast<double>(i + 1);
            plain[static_cast<std::size_t>(group_key[i])] += value[i];
        }

        auto runtime = PaperReview::MakeCkksRuntime(slots, 12);
        const std::vector<PaperReview::AttributeSpec> attrs{{groups}};
        const auto domain = PaperReview::EnumerateDomainTuples(attrs);
        const auto basis = PaperReview::EnumerateMonomialBasis(attrs);
        const auto matrix = PaperReview::BuildBasisMatrix(domain, basis);
        const auto alpha = PaperReview::SolveCoeffTableForAllTargets(matrix);
        if (!PaperReview::VerifyAlphaTable(matrix, alpha)) {
            std::cerr << "AlphaTable verification failed\n";
            return 1;
        }

        const auto ct_key = PaperReview::EncryptVector(
            runtime.cc, runtime.keys.publicKey, group_key, slots);
        const auto ct_value = PaperReview::EncryptVector(
            runtime.cc, runtime.keys.publicKey, value, slots);
        const auto encrypted_sums = PaperReview::MatrixGroupBySum(
            runtime.cc, runtime.keys.publicKey, attrs, {ct_key},
            ct_value, basis, alpha, slots);
        const auto got = DecryptSlot0(runtime, encrypted_sums);

        double max_err = 0.0;
        std::cout << "rows=" << rows << " groups=" << groups
                  << " slots=" << slots << "\n";
        std::cout << "group,plain,encrypted,abs_error\n";
        for (std::size_t g = 0; g < groups; g++) {
            const double err = std::abs(plain[g] - got[g]);
            max_err = std::max(max_err, err);
            std::cout << g << "," << plain[g] << "," << got[g] << ","
                      << err << "\n";
        }
        std::cout << "max_abs_error=" << max_err << "\n";
        return max_err < 0.5 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap4_group_by failed: " << e.what() << "\n";
        return 1;
    }
}
