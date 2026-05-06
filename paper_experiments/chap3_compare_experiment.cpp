#include "chap3_compare_experiment.hpp"

#include "algorithms.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <gate.hpp>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace PaperReview {
namespace {

struct br_lvl22param {
    using domainP = TFHEpp::lvl2param;
    using targetP = TFHEpp::lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

using brP_meta = br_lvl22param;
using brP_logari = TFHEpp::lvl02param;
using brP_base = TFHEpp::lvl01param;
using iksP_t = TFHEpp::lvl20param;
using P_in = TFHEpp::lvl2param;
using P_out = TFHEpp::lvl1param;

template <class P>
typename P::T EncodeMessage(std::uint64_t message, int precision_bits) {
    return static_cast<typename P::T>(message)
           << (std::numeric_limits<typename P::T>::digits - precision_bits);
}

template <class P>
int DecodeSign(typename P::T phase) {
    return static_cast<std::make_signed_t<typename P::T>>(phase) < 0 ? 1 : 0;
}

double MsSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

std::uint64_t MaxInputValue(int bits) {
    if (bits < 2 || bits >= 63)
        throw std::invalid_argument("comparison bits must be in [2, 62]");
    return (std::uint64_t{1} << bits) - 1;
}

std::vector<Chapter3CompareCase> StandardCases(int bits) {
    const std::uint64_t max = MaxInputValue(bits);
    return {
        Chapter3CompareCase{bits, max / 2 + 1, max / 4},
        Chapter3CompareCase{bits, max / 4, max / 2 + 1},
        Chapter3CompareCase{bits, max, max / 2},
        Chapter3CompareCase{bits, max / 2, max},
        Chapter3CompareCase{bits, max / 2, max / 2},
    };
}

std::vector<int> ParseIntList(const std::string& text) {
    std::vector<int> values;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) throw std::invalid_argument("empty integer list entry");
        values.push_back(std::stoi(item));
    }
    if (values.empty()) throw std::invalid_argument("empty integer list");
    return values;
}

Chapter3PredicateKind ParsePredicate(const std::string& text) {
    if (text == "gt" || text == "greater_than" || text == "relation" ||
        text == "rel")
        return Chapter3PredicateKind::GreaterThan;
    if (text == "eq" || text == "equal" || text == "equality")
        return Chapter3PredicateKind::Equal;
    if (text == "ne" || text == "neq" || text == "not_equal" ||
        text == "inequality")
        return Chapter3PredicateKind::NotEqual;
    throw std::invalid_argument(
        "predicate must be relation/gt, equality/eq, or inequality/ne");
}

const char* PredicateName(Chapter3PredicateKind predicate) {
    switch (predicate) {
    case Chapter3PredicateKind::GreaterThan:
        return "greater_than";
    case Chapter3PredicateKind::Equal:
        return "equal";
    case Chapter3PredicateKind::NotEqual:
        return "not_equal";
    }
    return "unknown";
}

std::vector<Chapter3CompareCase> BuildCases(
    const Chapter3CompareExperimentOptions& options) {
    std::vector<Chapter3CompareCase> cases = options.cases;
    if (cases.empty()) {
        for (const int bits : options.bit_widths) {
            MaxInputValue(bits);
            const auto standard = StandardCases(bits);
            cases.insert(cases.end(), standard.begin(), standard.end());
        }
    } else {
        for (const auto& c : cases) {
            const auto max = MaxInputValue(c.bits);
            if (c.lhs > max || c.rhs > max)
                throw std::invalid_argument("comparison case exceeds [0, 2^bits-1]");
        }
    }

    if (options.random_cases_per_width > 0) {
        std::mt19937_64 rng(options.seed);
        for (const int bits : options.bit_widths) {
            const std::uint64_t max = MaxInputValue(bits);
            std::uniform_int_distribution<std::uint64_t> dist(0, max);
            for (int i = 0; i < options.random_cases_per_width; i++)
                cases.push_back(Chapter3CompareCase{
                    bits, dist(rng), dist(rng)});
        }
    }
    return cases;
}

}  // namespace

void ParseChapter3CompareArgs(
    Chapter3CompareExperimentOptions& options,
    int argc,
    char** argv) {
    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        if (arg == "--full") {
            options.bit_widths = {4, 8, 16, 32};
            options.random_cases_per_width = std::max(options.random_cases_per_width, 8);
        } else if (arg == "--bits" && i + 1 < argc) {
            options.bit_widths = {std::stoi(argv[++i])};
        } else if (arg == "--bits-list" && i + 1 < argc) {
            options.bit_widths = ParseIntList(argv[++i]);
        } else if (arg == "--kappa" && i + 1 < argc) {
            options.kappas = {std::stoi(argv[++i])};
        } else if (arg == "--kappas" && i + 1 < argc) {
            options.kappas = ParseIntList(argv[++i]);
        } else if ((arg == "--random" || arg == "--trials") && i + 1 < argc) {
            options.random_cases_per_width = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            options.seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        } else if ((arg == "--output" || arg == "--output-file") && i + 1 < argc) {
            options.output_path = argv[++i];
        } else if (arg == "--predicate" && i + 1 < argc) {
            options.predicate = ParsePredicate(argv[++i]);
        } else if (arg == "--equality") {
            options.predicate = Chapter3PredicateKind::Equal;
        } else if (arg == "--inequality") {
            options.predicate = Chapter3PredicateKind::NotEqual;
        } else if (arg == "--relation") {
            options.predicate = Chapter3PredicateKind::GreaterThan;
        } else {
            throw std::invalid_argument(
                "usage: [--full|--bits N|--bits-list 4,8,16,32] "
                "[--kappa K|--kappas K1,K2] [--trials N] [--seed S] "
                "[--predicate relation|equality|inequality] [--output PATH]");
        }
    }
}

int RunChapter3CompareExperiment(bool full, std::ostream& os) {
    Chapter3CompareExperimentOptions options;
    options.label = "Chapter 3 HomCompare encrypted test";
    options.bit_widths = full ? std::vector<int>{4, 8, 16, 32}
                              : std::vector<int>{4};
    options.random_cases_per_width = full ? 8 : 0;
    return RunChapter3CompareExperiment(options, os);
}

int RunChapter3CompareExperiment(
    const Chapter3CompareExperimentOptions& options,
    std::ostream& os) {
    std::unique_ptr<std::ofstream> output_file;
    std::ostream* output = nullptr;
    if (!options.output_path.empty()) {
        output_file = std::make_unique<std::ofstream>(options.output_path);
        if (!*output_file)
            throw std::runtime_error("failed to open output file: " + options.output_path);
        output = output_file.get();
    }
    auto write_line = [&](const std::string& text) {
        PaperReview::WriteOutputLine(os, output, text);
    };
    if (options.kappas.empty())
        throw std::invalid_argument("at least one kappa value is required");
    const auto cases = BuildCases(options);
    if (cases.empty())
        throw std::invalid_argument("at least one comparison case is required");

    write_line("=== " + options.label + " ===");
    write_line(std::string("predicate=") + PredicateName(options.predicate) +
               " implementation=PaperReview::Chapter3HomCompare");
    std::ostringstream method;
    method << "method=encrypt lhs/rhs as TFHE lvl2 TLWE, run HomCompare, and decrypt only the final predicate ciphertext";
    if (options.predicate == Chapter3PredicateKind::Equal)
        method << "; equal is computed in TFHE as NOT(OR(lhs>rhs,lhs<rhs))";
    if (options.predicate == Chapter3PredicateKind::NotEqual)
        method << "; inequality is computed in TFHE as OR(lhs>rhs,lhs<rhs)";
    write_line(method.str());
    write_line("metric,value");
    write_line("operator,hom_compare");
    write_line(std::string("predicate,") + PredicateName(options.predicate));
    write_line("security=" + SecuritySummary());
    for (const auto& trace : Chapter3AlgorithmTrace())
        write_line("algorithm_trace=" + trace);

    const auto cfg = Chapter3MetaPBSConfig();
    const int metapbs_window_bits =
        MetaPBS2::MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    for (const int kappa : options.kappas) {
        if (kappa <= 0 || kappa >= metapbs_window_bits)
            throw std::invalid_argument("kappa must be in the MetaPBS extraction window");
    }

    TFHEpp::SecretKey sk;

    auto bk_meta = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_meta>>();
    auto bk_logari = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_logari>>();
    auto bk_base = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP_base>>();
    auto iksk = std::make_unique<TFHEpp::KeySwitchingKey<iksP_t>>();
    TFHEpp::EvalKey gate_ek(sk);

    write_line("keygen,start");
    TFHEpp::bkfftgen<brP_meta>(*bk_meta, sk);
    TFHEpp::bkfftgen<brP_logari>(*bk_logari, sk);
    TFHEpp::bkfftgen<brP_base>(*bk_base, sk);
    TFHEpp::ikskgen<iksP_t>(*iksk, sk);
    gate_ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    gate_ek.emplaceiksk<TFHEpp::lvl10param>(sk);

    std::vector<MetaPBS2::TruncRepeatKey<P_in>> trkeys;
    trkeys.push_back(MetaPBS2::GenerateTruncRepeatKey<P_in>(
        sk.key.get<P_in>(), cfg.rounds[0].beta));
    trkeys.push_back(MetaPBS2::GenerateTruncRepeatKey<P_in>(
        sk.key.get<P_in>(), cfg.rounds[1].beta));
    write_line("keygen,done");

    int total = 0;
    int passed = 0;
    write_line("kappa,bits,lhs,rhs,expected,got,ms,prune_skipped,prune_total");
    for (const int kappa : options.kappas) {
        MetaPBS2::HomMSBOptions hom_options;
        hom_options.kappa = kappa;
        for (const auto& c : cases) {
            const int encode_p = c.bits + 1;
            TFHEpp::TLWE<P_in> lhs_ct{}, rhs_ct{};
            TFHEpp::tlweSymEncrypt<P_in>(
                lhs_ct, EncodeMessage<P_in>(c.lhs, encode_p),
                P_in::α, sk.key.get<P_in>());
            TFHEpp::tlweSymEncrypt<P_in>(
                rhs_ct, EncodeMessage<P_in>(c.rhs, encode_p),
                P_in::α, sk.key.get<P_in>());

            MetaPBS2::BlindRotatePruneStats stats{};
            const auto start = std::chrono::steady_clock::now();
            TFHEpp::TLWE<P_out> out{};
            int expected = 0;
            if (options.predicate == Chapter3PredicateKind::GreaterThan) {
                out = Chapter3HomCompare<
                    brP_meta, brP_logari, brP_base, iksP_t>(
                        lhs_ct, rhs_ct, c.bits, ComparePredicate::GreaterThan,
                        *bk_meta, trkeys, cfg, *bk_logari, *bk_base, *iksk,
                        hom_options, &stats);
                expected = c.lhs > c.rhs ? 1 : 0;
            } else {
                const auto gt = Chapter3HomCompare<
                    brP_meta, brP_logari, brP_base, iksP_t>(
                        lhs_ct, rhs_ct, c.bits, ComparePredicate::GreaterThan,
                        *bk_meta, trkeys, cfg, *bk_logari, *bk_base, *iksk,
                        hom_options, &stats);
                const auto lt = Chapter3HomCompare<
                    brP_meta, brP_logari, brP_base, iksP_t>(
                        lhs_ct, rhs_ct, c.bits, ComparePredicate::LessThan,
                        *bk_meta, trkeys, cfg, *bk_logari, *bk_base, *iksk,
                        hom_options, &stats);
                TFHEpp::TLWE<P_out> gt_gate{};
                TFHEpp::TLWE<P_out> lt_gate{};
                TFHEpp::HomNOT<P_out>(gt_gate, gt);
                TFHEpp::HomNOT<P_out>(lt_gate, lt);
                TFHEpp::TLWE<P_out> neq_gate{};
                TFHEpp::HomOR<TFHEpp::lvl10param, TFHEpp::lvl01param, P_out::μ>(
                    neq_gate, gt_gate, lt_gate, gate_ek);
                if (options.predicate == Chapter3PredicateKind::Equal) {
                    out = neq_gate;
                    expected = c.lhs == c.rhs ? 1 : 0;
                } else {
                    TFHEpp::HomNOT<P_out>(out, neq_gate);
                    expected = c.lhs != c.rhs ? 1 : 0;
                }
            }
            const auto phase =
                TFHEpp::tlweSymPhase<P_out>(out, sk.key.get<P_out>());
            const int got = DecodeSign<P_out>(phase);
            const double ms = MsSince(start);

            total++;
            if (got == expected) passed++;
            std::ostringstream row;
            row << kappa << "," << c.bits << "," << c.lhs << "," << c.rhs << ","
                << expected << "," << got << ","
                << std::fixed << std::setprecision(2) << ms << ","
                << stats.skipped << "," << stats.total;
            write_line(row.str());
        }
    }

    const double accuracy =
        total == 0 ? 0.0 : static_cast<double>(passed) / total;
    write_line("summary,passed,total,accuracy");
    {
        std::ostringstream row;
        row << "summary," << passed << "," << total << ","
            << std::setprecision(6) << accuracy;
        write_line(row.str());
    }
    write_line("passed," + std::to_string(passed));
    write_line("total," + std::to_string(total));
    write_line("accuracy," + std::to_string(accuracy));
    return passed == total ? 0 : 1;
}

}  // namespace PaperReview
