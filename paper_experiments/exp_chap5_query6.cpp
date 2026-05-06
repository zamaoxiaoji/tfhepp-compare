#include "tpch_reference_results.hpp"

#include <iostream>

int main(int argc, char** argv) {
    try {
        constexpr auto query = PaperReview::TpchReferenceQuery::Q6;
        const auto options =
            PaperReview::ParseTpchReferenceOptions(argc, argv, "exp_chap5_query6", query);
        PaperReview::PrintTpchReferenceResult(query, options, std::cout);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query6 failed: " << e.what() << "\n";
        return 1;
    }
}
