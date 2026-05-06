#include "chap3_compare_experiment.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<PaperReview::Chapter3CompareCase> BoundaryCases(
    int bits,
    int radius,
    std::uint64_t threshold) {
    const auto max = (std::uint64_t{1} << bits) - 1;
    std::vector<PaperReview::Chapter3CompareCase> cases;
    for (int offset = -radius; offset <= radius; offset++) {
        const auto rhs_signed =
            static_cast<long long>(threshold) + static_cast<long long>(offset);
        if (rhs_signed < 0 || static_cast<std::uint64_t>(rhs_signed) > max)
            continue;
        cases.push_back({
            bits,
            threshold,
            static_cast<std::uint64_t>(rhs_signed),
        });
    }
    return cases;
}

}  // namespace

int main(int argc, char** argv) {
    int bits = 4;
    int radius = 3;
    std::uint64_t threshold = 0;
    PaperReview::Chapter3CompareExperimentOptions options;
    options.label = "Chapter 3 comparison boundary encrypted test";
    options.kappas = {5};

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--full") {
                bits = 16;
                radius = 4;
            } else if (arg == "--bits" && i + 1 < argc) {
                bits = std::stoi(argv[++i]);
            } else if (arg == "--radius" && i + 1 < argc) {
                radius = std::stoi(argv[++i]);
            } else if (arg == "--threshold" && i + 1 < argc) {
                threshold = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            } else if (arg == "--kappa" && i + 1 < argc) {
                options.kappas = {std::stoi(argv[++i])};
            } else if ((arg == "--trials" || arg == "--random") && i + 1 < argc) {
                options.random_cases_per_width = std::stoi(argv[++i]);
            } else if (arg == "--seed" && i + 1 < argc) {
                options.seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            } else if ((arg == "--output" || arg == "--output-file") && i + 1 < argc) {
                options.output_path = argv[++i];
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap3_boundary [--full|--bits N] [--radius N] "
                    "[--threshold T] [--kappa N] [--trials N] [--seed S] "
                    "[--output PATH]");
            }
        }
        options.bit_widths = {bits};
        if (threshold == 0) threshold = std::uint64_t{1} << (bits - 1);
        options.cases = BoundaryCases(bits, radius, threshold);
        return PaperReview::RunChapter3CompareExperiment(options, std::cout);
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_boundary failed: " << e.what() << "\n";
        return 1;
    }
}
