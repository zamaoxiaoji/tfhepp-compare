#include "algorithms.hpp"
#include "native_baseline_experiments.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
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

}  // namespace

int main(int argc, char** argv) {
    std::size_t rows = 4;
    std::size_t key_domain = 4;
    std::uint64_t seed = 42;
    std::uint32_t ckks_depth = 12;
    auto system = PaperReview::NativeBaselineSystem::Ours;
    std::unique_ptr<std::ofstream> output_file;
    std::ostream* output = nullptr;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--rows" && i + 1 < argc) {
                rows = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--domain" && i + 1 < argc) {
                key_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--key-domain" && i + 1 < argc) {
                key_domain = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--seed" && i + 1 < argc) {
                seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            } else if (arg == "--ckks-depth" && i + 1 < argc) {
                ckks_depth = static_cast<std::uint32_t>(std::stoul(argv[++i]));
            } else if (arg == "--system" && i + 1 < argc) {
                system = PaperReview::ParseNativeBaselineSystem(argv[++i]);
            } else if (arg == "--output" || arg == "--output-file") {
                output = PaperReview::OpenOptionalOutputFile(
                    i, argc, argv, output_file);
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap4_join [--rows N] [--domain M|--key-domain M] [--seed S] "
                    "[--system ours|he3db|arcedb|all] [--ckks-depth D] [--output PATH]");
            }
        }
        if (rows == 0 || key_domain == 0)
            throw std::invalid_argument("rows and domain must be nonzero");

        if (system == PaperReview::NativeBaselineSystem::HE3DB ||
            system == PaperReview::NativeBaselineSystem::ArcEDB ||
            system == PaperReview::NativeBaselineSystem::All) {
            PaperReview::NativeBaselineOptions native_options;
            native_options.system = system;
            native_options.rows = rows;
            native_options.key_domain = key_domain;
            native_options.seed = seed;
            for (const auto& result :
                 PaperReview::RunChap4JoinBaselines(native_options)) {
                PaperReview::PrintNativeBaselineResult(result, std::cout);
                if (output) PaperReview::PrintNativeBaselineResult(result, *output);
            }
            if (system != PaperReview::NativeBaselineSystem::All) return 0;
        }

        if (system == PaperReview::NativeBaselineSystem::HE3DB ||
            system == PaperReview::NativeBaselineSystem::ArcEDB)
            return 0;

        const std::size_t slots = NextPowerOfTwo(std::max(rows, key_domain));
        std::vector<double> left_key(rows);
        std::vector<double> expected(rows);
        std::vector<double> right_key(slots, 0.0);
        std::vector<double> right_payload(slots, 0.0);
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<int> key_dist(
            0, static_cast<int>(key_domain - 1));
        std::uniform_int_distribution<int> payload_dist(1, 1000);

        for (std::size_t j = 0; j < key_domain; j++) {
            right_key[j] = static_cast<double>(j);
            right_payload[j] = static_cast<double>(payload_dist(rng));
        }
        for (std::size_t i = 0; i < rows; i++) {
            const std::size_t key = static_cast<std::size_t>(key_dist(rng));
            left_key[i] = static_cast<double>(key);
            expected[i] = right_payload[key];
        }

        auto runtime = PaperReview::MakeCkksRuntime(
            slots, ckks_depth, false, 55);
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
        double sum_abs_err = 0.0;
        std::size_t correct = 0;
        {
            std::ostringstream line;
            line << "rows=" << rows << " domain=" << key_domain
                 << " slots=" << slots << " seed=" << seed
                 << " ckks_depth=" << ckks_depth
                 << " tolerance=" << kCkksTolerance
                 << " scaling_mod_size=55";
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        PaperReview::WriteOutputLine(std::cout, output, "metric,value");
        PaperReview::WriteOutputLine(std::cout, output, "operator,join");
        PaperReview::WriteOutputLine(std::cout, output, "system,ours");
        PaperReview::WriteOutputLine(std::cout, output, "rows," + std::to_string(rows));
        PaperReview::WriteOutputLine(std::cout, output, "key_domain," + std::to_string(key_domain));
        PaperReview::WriteOutputLine(std::cout, output, "slots," + std::to_string(slots));
        PaperReview::WriteOutputLine(std::cout, output, "seed," + std::to_string(seed));
        PaperReview::WriteOutputLine(std::cout, output, "ckks_depth," + std::to_string(ckks_depth));
        PaperReview::WriteOutputLine(std::cout, output, "tolerance,1e-6");
        PaperReview::WriteOutputLine(std::cout, output, "slot,left_key,expected_payload,encrypted_payload,abs_error");
        for (std::size_t i = 0; i < rows; i++) {
            const double err = std::abs(expected[i] - got[i]);
            max_err = std::max(max_err, err);
            sum_abs_err += err;
            if (err <= kCkksTolerance) correct++;
            std::ostringstream line;
            line << i << "," << left_key[i] << "," << expected[i]
                 << "," << got[i] << "," << err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        const double accuracy =
            rows == 0 ? 0.0 : static_cast<double>(correct) / rows;
        const double avg_err = rows == 0 ? 0.0 : sum_abs_err / rows;
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
            line << "avg_abs_error," << std::scientific << avg_err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        {
            std::ostringstream line;
            line << "max_abs_error," << std::scientific << max_err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        return correct == rows ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap4_join failed: " << e.what() << "\n";
        return 1;
    }
}
