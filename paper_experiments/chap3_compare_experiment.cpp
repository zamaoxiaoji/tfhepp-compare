#include "chap3_compare_experiment.hpp"

#include "algorithms.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>
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
    return (std::uint64_t{1} << (bits - 1)) - 1;
}

std::vector<Chapter3CompareCase> StandardCases(int bits) {
    const std::uint64_t max = (std::uint64_t{1} << (bits - 1)) - 1;
    return {
        Chapter3CompareCase{bits, max / 2 + 1, max / 4},
        Chapter3CompareCase{bits, max / 4, max / 2 + 1},
        Chapter3CompareCase{bits, max, max / 2},
        Chapter3CompareCase{bits, max / 2, max},
    };
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
                throw std::invalid_argument("comparison case exceeds signed-safe input range");
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

int RunChapter3CompareExperiment(bool full, std::ostream& os) {
    Chapter3CompareExperimentOptions options;
    options.label = "Chapter 3 HomCompare encrypted test";
    options.bit_widths = full ? std::vector<int>{4, 8, 16}
                              : std::vector<int>{4};
    options.random_cases_per_width = full ? 8 : 0;
    return RunChapter3CompareExperiment(options, os);
}

int RunChapter3CompareExperiment(
    const Chapter3CompareExperimentOptions& options,
    std::ostream& os) {
    if (options.kappas.empty())
        throw std::invalid_argument("at least one kappa value is required");
    const auto cases = BuildCases(options);
    if (cases.empty())
        throw std::invalid_argument("at least one comparison case is required");

    os << "=== " << options.label << " ===\n";
    os << "predicate=greater_than"
       << " implementation=PaperReview::Chapter3HomCompare\n";
    os << "method=encrypt lhs/rhs as TFHE lvl2 TLWE, run HomCompare, decrypt lvl1 sign bit, compare against plaintext predicate\n";
    os << "security=" << SecuritySummary() << "\n";
    for (const auto& trace : Chapter3AlgorithmTrace())
        os << "algorithm_trace=" << trace << "\n";

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

    os << "keygen,start\n";
    TFHEpp::bkfftgen<brP_meta>(*bk_meta, sk);
    TFHEpp::bkfftgen<brP_logari>(*bk_logari, sk);
    TFHEpp::bkfftgen<brP_base>(*bk_base, sk);
    TFHEpp::ikskgen<iksP_t>(*iksk, sk);

    std::vector<MetaPBS2::TruncRepeatKey<P_in>> trkeys;
    trkeys.push_back(MetaPBS2::GenerateTruncRepeatKey<P_in>(
        sk.key.get<P_in>(), cfg.rounds[0].beta));
    trkeys.push_back(MetaPBS2::GenerateTruncRepeatKey<P_in>(
        sk.key.get<P_in>(), cfg.rounds[1].beta));
    os << "keygen,done\n";

    int total = 0;
    int passed = 0;
    os << "kappa,bits,lhs,rhs,expected,got,ms,prune_skipped,prune_total\n";
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
            const auto out = Chapter3HomCompare<
                brP_meta, brP_logari, brP_base, iksP_t>(
                    lhs_ct, rhs_ct, c.bits, ComparePredicate::GreaterThan,
                    *bk_meta, trkeys, cfg, *bk_logari, *bk_base, *iksk,
                    hom_options, &stats);
            const double ms = MsSince(start);

            const auto phase = TFHEpp::tlweSymPhase<P_out>(out, sk.key.get<P_out>());
            const int got = DecodeSign<P_out>(phase);
            const int expected = c.lhs > c.rhs ? 1 : 0;
            total++;
            if (got == expected) passed++;
            os << kappa << "," << c.bits << "," << c.lhs << "," << c.rhs << ","
               << expected << "," << got << ","
               << std::fixed << std::setprecision(2) << ms << ","
               << stats.skipped << "," << stats.total << "\n";
        }
    }

    os << "summary,passed,total,accuracy\n";
    os << "summary," << passed << "," << total << ","
       << std::setprecision(6)
       << (total == 0 ? 0.0 : static_cast<double>(passed) / total)
       << "\n";
    return passed == total ? 0 : 1;
}

}  // namespace PaperReview
