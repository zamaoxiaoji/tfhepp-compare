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
    int radius) {
    const auto max = (std::uint64_t{1} << (bits - 1)) - 1;
    const auto rhs = max / 2;
    std::vector<PaperReview::Chapter3CompareCase> cases;
    for (int offset = -radius; offset <= radius; offset++) {
        if (offset == 0) {
            cases.push_back({bits, rhs, rhs});
            continue;
        }
        const auto lhs_signed =
            static_cast<long long>(rhs) + static_cast<long long>(offset);
        if (lhs_signed < 0 || static_cast<std::uint64_t>(lhs_signed) > max)
            continue;
        cases.push_back({
            bits,
            static_cast<std::uint64_t>(lhs_signed),
            rhs,
        });
    }
    return cases;
}

}  // namespace

int main(int argc, char** argv) {
    int bits = 4;
    int radius = 3;
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
            } else if (arg == "--kappa" && i + 1 < argc) {
                options.kappas = {std::stoi(argv[++i])};
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap3_boundary [--full|--bits N] [--radius N] [--kappa N]");
            }
        }
        options.bit_widths = {bits};
        options.cases = BoundaryCases(bits, radius);
        return PaperReview::RunChapter3CompareExperiment(options, std::cout);
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_boundary failed: " << e.what() << "\n";
        return 1;
    }
}
