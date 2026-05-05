#include "tpch_experiments.hpp"

#include <iostream>

int main(int argc, char** argv) {
    try {
        const auto options = PaperReview::ParseExperimentOptions(argc, argv);
        const auto result = PaperReview::RunTpchQ6Experiment(options);
        PaperReview::PrintExperimentResult(result, std::cout);
        return result.abs_error < 1.0 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query6 failed: " << e.what() << "\n";
        return 1;
    }
}
