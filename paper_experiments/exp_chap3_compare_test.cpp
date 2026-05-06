#include "chap3_compare_experiment.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    try {
        PaperReview::Chapter3CompareExperimentOptions options;
        options.label = "Chapter 3 HomCompare encrypted test";
        PaperReview::ParseChapter3CompareArgs(options, argc, argv);
        return PaperReview::RunChapter3CompareExperiment(options, std::cout);
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_compare_test failed: " << e.what() << "\n";
        return 1;
    }
}
