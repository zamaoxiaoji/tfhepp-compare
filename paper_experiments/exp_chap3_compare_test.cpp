#include "chap3_compare_experiment.hpp"

#include <iostream>
#include <string>

int main(int argc, char** argv) {
    bool full = false;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--full") full = true;
    }
    return PaperReview::RunChapter3CompareExperiment(full, std::cout);
}
