#include "tpch_experiments.hpp"

#include <iostream>

int main(int argc, char** argv) {
    try {
        auto options = PaperReview::ParseExperimentOptions(argc, argv);
        const auto result = PaperReview::RunTpchQ3Experiment(options);
        PaperReview::PrintExperimentResult(result, std::cout);
        return result.abs_error < 2.0 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query3 failed: " << e.what() << "\n";
        return 1;
    }
}
