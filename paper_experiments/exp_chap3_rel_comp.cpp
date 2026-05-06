#include "chap3_compare_experiment.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    PaperReview::Chapter3CompareExperimentOptions options;
    options.label = "Chapter 3 relational comparison encrypted benchmark";
    options.bit_widths = {4};
    options.kappas = {5};

    try {
        PaperReview::ParseChapter3CompareArgs(options, argc, argv);
        return PaperReview::RunChapter3CompareExperiment(options, std::cout);
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_rel_comp failed: " << e.what() << "\n";
        return 1;
    }
}
