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

int PowInt(int base, int exp) {
    int out = 1;
    for (int i = 0; i < exp; i++) out *= base;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    std::size_t rows = 8;
    int base = 4;
    int digit_count = 3;
    int target = 15;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--rows" && i + 1 < argc) {
                rows = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--base" && i + 1 < argc) {
                base = std::stoi(argv[++i]);
            } else if (arg == "--digits" && i + 1 < argc) {
                digit_count = std::stoi(argv[++i]);
            } else if (arg == "--target" && i + 1 < argc) {
                target = std::stoi(argv[++i]);
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap4_rns [--rows N] [--base P] [--digits L] [--target V]");
            }
        }
        if (rows == 0 || base <= 1 || digit_count <= 0)
            throw std::invalid_argument("rows, base, and digits must be valid");
        const int capacity = PowInt(base, digit_count);
        if (target < 0 || target >= capacity)
            throw std::invalid_argument("target is outside the digit capacity");

        const std::size_t slots = NextPowerOfTwo(rows);
        std::vector<int> values(rows);
        std::vector<double> expected(slots, 0.0);
        for (std::size_t i = 0; i < rows; i++) {
            values[i] = static_cast<int>(i % static_cast<std::size_t>(capacity));
            if (i == 0) values[i] = target;
            expected[i] = values[i] == target ? 1.0 : 0.0;
        }

        auto runtime = PaperReview::MakeCkksRuntime(slots, 12);
        const auto digit_cts = PaperReview::EncryptDigitColumns(
            runtime.cc, runtime.keys.publicKey, values, base, digit_count, slots);
        const auto mask = PaperReview::DigitMaskFromDigits(
            runtime.cc, runtime.keys.publicKey, digit_cts,
            target, base, digit_count, slots);

        lbcrypto::Plaintext plain;
        runtime.cc->Decrypt(runtime.keys.secretKey, mask, &plain);
        plain->SetLength(slots);
        const auto got = plain->GetRealPackedValue();

        const auto target_digits = PaperReview::BaseDigits(target, base, digit_count);
        std::cout << "rows=" << rows << " base=" << base
                  << " digits=" << digit_count << " target=" << target
                  << " slots=" << slots << "\n";
        std::cout << "target_digits";
        for (const auto d : target_digits) std::cout << "," << d;
        std::cout << "\n";
        std::cout << "slot,value,expected_mask,encrypted_mask,abs_error\n";
        double max_err = 0.0;
        for (std::size_t i = 0; i < rows; i++) {
            const double err = std::abs(expected[i] - got[i]);
            max_err = std::max(max_err, err);
            std::cout << i << "," << values[i] << "," << expected[i]
                      << "," << got[i] << "," << err << "\n";
        }
        std::cout << "max_abs_error=" << max_err << "\n";
        return max_err < 0.5 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap4_rns failed: " << e.what() << "\n";
        return 1;
    }
}
