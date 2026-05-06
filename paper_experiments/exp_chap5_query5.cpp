#include "tpch_experiments.hpp"
#include "native_baseline_experiments.hpp"

#include <iostream>

int main(int argc, char** argv) {
    try {
        auto options = PaperReview::ParseExperimentOptions(argc, argv);
        const auto system = PaperReview::ParseNativeBaselineSystem(options.system);
        if (system == PaperReview::NativeBaselineSystem::HE3DB ||
            system == PaperReview::NativeBaselineSystem::ArcEDB ||
            system == PaperReview::NativeBaselineSystem::All) {
            PaperReview::NativeBaselineOptions native_options;
            native_options.system = system;
            native_options.rows = options.rows;
            native_options.key_domain = options.q5_key_domain;
            native_options.nation_domain = options.q5_nation_domain;
            native_options.seed = options.seed;
            for (const auto& result :
                 PaperReview::RunChap5Q5Baselines(native_options)) {
                PaperReview::PrintNativeBaselineResult(result, std::cout);
                PaperReview::WriteNativeBaselineResultIfRequested(
                    result, options.output_path);
            }
            if (system != PaperReview::NativeBaselineSystem::All) return 0;
        }
        const auto result = PaperReview::RunTpchQ5Experiment(options);
        PaperReview::PrintExperimentResult(result, std::cout);
        PaperReview::WriteExperimentResultIfRequested(result, options);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query5 failed: " << e.what() << "\n";
        return 1;
    }
}
