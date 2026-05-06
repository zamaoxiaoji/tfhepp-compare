#include "algorithms.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kCkksTolerance = 1e-6;

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
    int domain_size = 64;
    int target = 15;
    std::uint32_t ckks_depth = 12;
    bool digits_set = false;
    std::unique_ptr<std::ofstream> output_file;
    std::ostream* output = nullptr;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--rows" && i + 1 < argc) {
                rows = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--base" && i + 1 < argc) {
                base = std::stoi(argv[++i]);
            } else if (arg == "--digits" && i + 1 < argc) {
                digit_count = std::stoi(argv[++i]);
                digits_set = true;
            } else if (arg == "--domain-size" && i + 1 < argc) {
                domain_size = std::stoi(argv[++i]);
            } else if (arg == "--target" && i + 1 < argc) {
                target = std::stoi(argv[++i]);
            } else if (arg == "--ckks-depth" && i + 1 < argc) {
                ckks_depth = static_cast<std::uint32_t>(std::stoul(argv[++i]));
            } else if (arg == "--output" || arg == "--output-file") {
                output = PaperReview::OpenOptionalOutputFile(
                    i, argc, argv, output_file);
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap4_rns [--rows N] [--domain-size M] "
                    "[--base P] [--digits L] [--target V] [--ckks-depth D] "
                    "[--output PATH]");
            }
        }
        if (rows == 0 || base <= 1 || digit_count <= 0 || domain_size <= 0)
            throw std::invalid_argument("rows, base, and digits must be valid");
        if (!digits_set) {
            digit_count = 1;
            while (PowInt(base, digit_count) < domain_size) digit_count++;
        }
        const int capacity = PowInt(base, digit_count);
        if (target < 0 || target >= capacity)
            throw std::invalid_argument("target is outside the digit capacity");
        if (domain_size > capacity)
            throw std::invalid_argument("domain size exceeds digit capacity");

        const std::size_t slots = NextPowerOfTwo(rows);
        std::vector<int> values(rows);
        std::vector<double> expected(slots, 0.0);
        for (std::size_t i = 0; i < rows; i++) {
            values[i] = static_cast<int>(i % static_cast<std::size_t>(domain_size));
            if (i == 0) values[i] = target;
            expected[i] = values[i] == target ? 1.0 : 0.0;
        }

        auto runtime = PaperReview::MakeCkksRuntime(
            slots, ckks_depth, false, 55);
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
        {
            std::ostringstream line;
            line << "rows=" << rows << " domain_size=" << domain_size
                 << " direct_degree=" << (domain_size - 1)
                 << " base=" << base
                 << " digits=" << digit_count << " target=" << target
                 << " slots=" << slots << " ckks_depth=" << ckks_depth
                 << " tolerance=" << kCkksTolerance
                 << " scaling_mod_size=55";
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        PaperReview::WriteOutputLine(std::cout, output, "metric,value");
        PaperReview::WriteOutputLine(std::cout, output, "operator,rns_digit_mask");
        PaperReview::WriteOutputLine(std::cout, output, "rows," + std::to_string(rows));
        PaperReview::WriteOutputLine(std::cout, output, "domain_size," + std::to_string(domain_size));
        PaperReview::WriteOutputLine(std::cout, output, "direct_degree," + std::to_string(domain_size - 1));
        PaperReview::WriteOutputLine(std::cout, output, "base," + std::to_string(base));
        PaperReview::WriteOutputLine(std::cout, output, "digits," + std::to_string(digit_count));
        PaperReview::WriteOutputLine(std::cout, output, "rns_single_degree," + std::to_string(base - 1));
        PaperReview::WriteOutputLine(std::cout, output, "target," + std::to_string(target));
        PaperReview::WriteOutputLine(std::cout, output, "slots," + std::to_string(slots));
        PaperReview::WriteOutputLine(std::cout, output, "ckks_depth," + std::to_string(ckks_depth));
        PaperReview::WriteOutputLine(std::cout, output, "tolerance,1e-6");
        {
            std::ostringstream line;
            line << "target_digits";
            for (const auto d : target_digits) line << "," << d;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        PaperReview::WriteOutputLine(std::cout, output, "slot,value,expected_mask,encrypted_mask,abs_error");
        double max_err = 0.0;
        std::size_t correct = 0;
        for (std::size_t i = 0; i < rows; i++) {
            const double err = std::abs(expected[i] - got[i]);
            max_err = std::max(max_err, err);
            if (err <= kCkksTolerance) correct++;
            std::ostringstream line;
            line << i << "," << values[i] << "," << expected[i]
                 << "," << got[i] << "," << err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        const double accuracy =
            rows == 0 ? 0.0 : static_cast<double>(correct) / rows;
        PaperReview::WriteOutputLine(
            std::cout, output,
            "correct_slots=" + std::to_string(correct) + "/" + std::to_string(rows) +
                " accuracy=" + std::to_string(accuracy));
        PaperReview::WriteOutputLine(
            std::cout, output, "correct_slots," + std::to_string(correct));
        PaperReview::WriteOutputLine(
            std::cout, output, "accuracy," + std::to_string(accuracy));
        {
            std::ostringstream line;
            line << "max_abs_error," << std::scientific << max_err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        return correct == rows ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap4_rns failed: " << e.what() << "\n";
        return 1;
    }
}
