#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace PaperReview {

enum class Chapter3PredicateKind {
    GreaterThan,
    Equal,
    NotEqual,
};

struct Chapter3CompareCase {
    int bits = 0;
    std::uint64_t lhs = 0;
    std::uint64_t rhs = 0;
};

struct Chapter3CompareExperimentOptions {
    std::string label = "Chapter 3 HomCompare encrypted test";
    std::vector<int> bit_widths{4};
    std::vector<int> kappas{5};
    std::vector<Chapter3CompareCase> cases;
    int random_cases_per_width = 0;
    std::uint64_t seed = 0xC0FFEE;
    Chapter3PredicateKind predicate = Chapter3PredicateKind::GreaterThan;
    std::string output_path;
};

void ParseChapter3CompareArgs(
    Chapter3CompareExperimentOptions& options,
    int argc,
    char** argv);

int RunChapter3CompareExperiment(bool full, std::ostream& os);
int RunChapter3CompareExperiment(
    const Chapter3CompareExperimentOptions& options,
    std::ostream& os);

}  // namespace PaperReview
