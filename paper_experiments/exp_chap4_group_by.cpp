#include "algorithms.hpp"
#include "native_baseline_experiments.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
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
    std::uint32_t ckks_depth = 12;
    auto system = PaperReview::NativeBaselineSystem::Ours;
    std::unique_ptr<std::ofstream> output_file;
    std::ostream* output = nullptr;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--rows" && i + 1 < argc) {
                rows = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--groups" && i + 1 < argc) {
                groups = static_cast<std::size_t>(std::stoull(argv[++i]));
            } else if (arg == "--group-bits" && i + 1 < argc) {
                const auto bits = std::stoul(argv[++i]);
                groups = std::size_t{1} << bits;
            } else if (arg == "--ckks-depth" && i + 1 < argc) {
                ckks_depth = static_cast<std::uint32_t>(std::stoul(argv[++i]));
            } else if (arg == "--system" && i + 1 < argc) {
                system = PaperReview::ParseNativeBaselineSystem(argv[++i]);
            } else if (arg == "--output" || arg == "--output-file") {
                output = PaperReview::OpenOptionalOutputFile(
                    i, argc, argv, output_file);
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap4_group_by [--rows N] [--groups M|--group-bits B] "
                    "[--system ours|he3db|arcedb|all] [--ckks-depth D] [--output PATH]");
            }
        }
        if (rows == 0 || groups == 0)
            throw std::invalid_argument("rows and groups must be nonzero");

        if (system == PaperReview::NativeBaselineSystem::HE3DB ||
            system == PaperReview::NativeBaselineSystem::ArcEDB ||
            system == PaperReview::NativeBaselineSystem::All) {
            PaperReview::NativeBaselineOptions native_options;
            native_options.system = system;
            native_options.rows = rows;
            native_options.groups = groups;
            for (const auto& result :
                 PaperReview::RunChap4GroupByBaselines(native_options)) {
                PaperReview::PrintNativeBaselineResult(result, std::cout);
                if (output) PaperReview::PrintNativeBaselineResult(result, *output);
            }
            if (system != PaperReview::NativeBaselineSystem::All) return 0;
        }

        if (system == PaperReview::NativeBaselineSystem::HE3DB ||
            system == PaperReview::NativeBaselineSystem::ArcEDB)
            return 0;

        const std::size_t slots = NextPowerOfTwo(rows);
        std::vector<double> group_key(rows), value(rows);
        std::vector<double> plain(groups, 0.0);
        for (std::size_t i = 0; i < rows; i++) {
            group_key[i] = static_cast<double>(i % groups);
            value[i] = static_cast<double>(i + 1);
            plain[static_cast<std::size_t>(group_key[i])] += value[i];
        }

        auto runtime = PaperReview::MakeCkksRuntime(
            slots, ckks_depth, false, 55);
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
        std::size_t correct = 0;
        {
            std::ostringstream line;
            line << "rows=" << rows << " groups=" << groups
                 << " slots=" << slots << " ckks_depth=" << ckks_depth
                 << " tolerance=" << kCkksTolerance
                 << " scaling_mod_size=55";
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        PaperReview::WriteOutputLine(std::cout, output, "metric,value");
        PaperReview::WriteOutputLine(std::cout, output, "operator,group_by");
        PaperReview::WriteOutputLine(std::cout, output, "system,ours");
        PaperReview::WriteOutputLine(std::cout, output, "rows," + std::to_string(rows));
        PaperReview::WriteOutputLine(std::cout, output, "groups," + std::to_string(groups));
        PaperReview::WriteOutputLine(std::cout, output, "slots," + std::to_string(slots));
        PaperReview::WriteOutputLine(std::cout, output, "ckks_depth," + std::to_string(ckks_depth));
        PaperReview::WriteOutputLine(std::cout, output, "tolerance,1e-6");
        PaperReview::WriteOutputLine(std::cout, output, "group,plain,encrypted,abs_error");
        for (std::size_t g = 0; g < groups; g++) {
            const double err = std::abs(plain[g] - got[g]);
            max_err = std::max(max_err, err);
            if (err <= kCkksTolerance) correct++;
            std::ostringstream line;
            line << g << "," << plain[g] << "," << got[g] << "," << err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        const double accuracy =
            groups == 0 ? 0.0 : static_cast<double>(correct) / groups;
        PaperReview::WriteOutputLine(
            std::cout, output,
            "correct_groups=" + std::to_string(correct) + "/" + std::to_string(groups) +
                " accuracy=" + std::to_string(accuracy));
        PaperReview::WriteOutputLine(
            std::cout, output, "correct_groups," + std::to_string(correct));
        PaperReview::WriteOutputLine(
            std::cout, output, "accuracy," + std::to_string(accuracy));
        {
            std::ostringstream line;
            line << "max_abs_error," << std::scientific << max_err;
            PaperReview::WriteOutputLine(std::cout, output, line.str());
        }
        return correct == groups ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap4_group_by failed: " << e.what() << "\n";
        return 1;
    }
}
