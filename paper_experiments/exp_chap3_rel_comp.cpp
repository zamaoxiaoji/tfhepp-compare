#include "chap3_compare_experiment.hpp"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    PaperReview::Chapter3CompareExperimentOptions options;
    options.label = "Chapter 3 relational comparison encrypted benchmark";
    options.bit_widths = {4};
    options.kappas = {5};

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--full") {
                options.bit_widths = {4, 8, 16};
                options.random_cases_per_width = 4;
            } else if (arg == "--bits" && i + 1 < argc) {
                options.bit_widths = {std::stoi(argv[++i])};
            } else if (arg == "--random" && i + 1 < argc) {
                options.random_cases_per_width = std::stoi(argv[++i]);
            } else if (arg == "--kappa" && i + 1 < argc) {
                options.kappas = {std::stoi(argv[++i])};
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap3_rel_comp [--full|--bits N] [--random N] [--kappa N]");
            }
        }
        return PaperReview::RunChapter3CompareExperiment(options, std::cout);
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_rel_comp failed: " << e.what() << "\n";
        return 1;
    }
}
