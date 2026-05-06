#include "tpch_experiments.hpp"

#include <iostream>

int main(int argc, char** argv) {
    try {
        const auto options = PaperReview::ParseExperimentOptions(argc, argv);
        const auto result = PaperReview::RunTpchQ14Experiment(options);
        PaperReview::PrintExperimentResult(result, std::cout);
        PaperReview::WriteExperimentResultIfRequested(result, options);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query14 failed: " << e.what() << "\n";
        return 1;
    }
}
