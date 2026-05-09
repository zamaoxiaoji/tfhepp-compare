/**
 * @file my_ethmsb_pruned_fast_test.cpp
 * @brief Strict ETHMSB vs explicit fast/pruned ETHMSB benchmark.
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "ethmsb_compare.h"

using ETHMSB_NS::PruneStats;
using ETHMSB_NS::PrunedETHMSBMode;
using ETHMSB_NS::PrunedETHMSBOptions;
using ETHMSB_NS::tlweSymInt32Encrypt;

namespace {

struct FailureSample {
    uint32_t idx = 0;
    uint64_t p0 = 0;
    uint64_t p1 = 0;
    uint32_t expected = 0;
    uint32_t got = 0;
};

struct RunResult {
    double avg_ms = 0.0;
    uint32_t errors = 0;
    PruneStats stats;
    std::vector<FailureSample> failures;
};

struct Args {
    int cases = 100;
    std::vector<uint32_t> periods = {2, 4, 8, 16, 32};
    bool run_weighted = true;
    bool run_paper = true;
};

template <class P>
struct CipherInput {
    typename P::T p0 = 0;
    typename P::T p1 = 0;
    uint32_t expected = 0;
    TFHEpp::TLWE<P> c0 = {};
    TFHEpp::TLWE<P> c1 = {};
};

void SeedTFHE(uint64_t seed)
{
#ifdef USE_BLAKE3
    TFHEpp::generator = BLAKE3PRNG::BLAKE3PRNG<uint64_t>(seed);
#else
    (void) seed;
#endif
}

std::vector<std::string> SplitComma(const std::string &value)
{
    std::vector<std::string> out;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

Args ParseArgs(int argc, char **argv)
{
    Args args;
    for (int i = 1; i < argc; i++) {
        const std::string key = argv[i];
        if (key == "--cases" && i + 1 < argc) {
            args.cases = std::max(1, std::stoi(argv[++i]));
        }
        else if (key == "--periods" && i + 1 < argc) {
            args.periods.clear();
            for (const std::string &part : SplitComma(argv[++i]))
                args.periods.push_back(static_cast<uint32_t>(std::stoul(part)));
            if (args.periods.empty()) args.periods = {2, 4, 8, 16, 32};
        }
        else if (key == "--modes" && i + 1 < argc) {
            args.run_weighted = false;
            args.run_paper = false;
            for (const std::string &part : SplitComma(argv[++i])) {
                if (part == "weighted") args.run_weighted = true;
                if (part == "paper") args.run_paper = true;
            }
        }
    }
    return args;
}

template <class P>
void SetupEvalKey(TFHEpp::SecretKey &sk, TFHEpp::EvalKey &ek)
{
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
        ek.emplacebkfft<TFHEpp::lvl02param>(sk);
        ek.emplaceiksk<TFHEpp::lvl20param>(sk);
        ek.emplaceiksk<TFHEpp::lvl21param>(sk);
    }
}

template <class P>
std::vector<CipherInput<P>> MakeInputs(uint32_t plain_bits, int cases,
                                       TFHEpp::SecretKey &sk)
{
    std::mt19937_64 engine(0x5052554e45445448ULL + plain_bits);
    const uint32_t scale_bits =
        std::numeric_limits<typename P::T>::digits - plain_bits - 1;
    const uint64_t max_message = (1ULL << (plain_bits - 1)) - 1;
    std::uniform_int_distribution<uint64_t> message(0, max_message);

    std::vector<CipherInput<P>> inputs;
    inputs.reserve(cases);
    for (int i = 0; i < cases; i++) {
        CipherInput<P> in;
        in.p0 = static_cast<typename P::T>(message(engine));
        in.p1 = static_cast<typename P::T>(message(engine));
        in.expected = (in.p0 > in.p1) ? 1 : 0;
        in.c0 = tlweSymInt32Encrypt<P>(
            in.p0, P::α, std::pow(2., scale_bits), sk.key.template get<P>());
        in.c1 = tlweSymInt32Encrypt<P>(
            in.p1, P::α, std::pow(2., scale_bits), sk.key.template get<P>());
        inputs.push_back(in);
    }
    return inputs;
}

template <class P>
RunResult RunStrict(const std::vector<CipherInput<P>> &inputs,
                    uint32_t plain_bits, TFHEpp::SecretKey &sk,
                    TFHEpp::EvalKey &ek)
{
    RunResult result;
    const auto begin = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < inputs.size(); i++) {
        auto c0 = inputs[i].c0;
        auto c1 = inputs[i].c1;
        ETHMSB_NS::TLWELvl1 out;
        ETHMSB_NS::greater_than<P>(c0, c1, out, plain_bits, ek, LOGIC);
        const uint32_t got =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(out, sk.key.lvl1);
        if (got != inputs[i].expected) {
            result.errors++;
            if (result.failures.size() < 5)
                result.failures.push_back(
                    {static_cast<uint32_t>(i),
                     static_cast<uint64_t>(inputs[i].p0),
                     static_cast<uint64_t>(inputs[i].p1), inputs[i].expected,
                     got});
        }
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.avg_ms =
        std::chrono::duration<double, std::milli>(end - begin).count() /
        std::max<size_t>(inputs.size(), 1);
    return result;
}

template <class P>
RunResult RunPruned(const std::vector<CipherInput<P>> &inputs,
                    uint32_t plain_bits, TFHEpp::SecretKey &sk,
                    TFHEpp::EvalKey &ek, PrunedETHMSBMode mode,
                    uint32_t period)
{
    RunResult result;
    PrunedETHMSBOptions options;
    options.mode = mode;
    options.period_idx = period;
    options.window_bits = 5;

    const auto begin = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < inputs.size(); i++) {
        auto c0 = inputs[i].c0;
        auto c1 = inputs[i].c1;
        ETHMSB_NS::TLWELvl1 out;
        ETHMSB_NS::greater_than_pruned_fast<P>(
            c0, c1, out, plain_bits, ek, LOGIC, options, &result.stats);
        const uint32_t got =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(out, sk.key.lvl1);
        if (got != inputs[i].expected) {
            result.errors++;
            if (result.failures.size() < 5)
                result.failures.push_back(
                    {static_cast<uint32_t>(i),
                     static_cast<uint64_t>(inputs[i].p0),
                     static_cast<uint64_t>(inputs[i].p1), inputs[i].expected,
                     got});
        }
    }
    const auto end = std::chrono::high_resolution_clock::now();
    result.avg_ms =
        std::chrono::duration<double, std::milli>(end - begin).count() /
        std::max<size_t>(inputs.size(), 1);
    return result;
}

std::string ModeName(PrunedETHMSBMode mode)
{
    if (mode == PrunedETHMSBMode::WeightedApprox) return "weighted-approx";
    return "paper-qhalf";
}

void PrintFailures(uint32_t bits, const std::string &mode, uint32_t period,
                   const RunResult &result)
{
    for (const FailureSample &f : result.failures) {
        std::cout << "failure,bits=" << bits << ",mode=" << mode
                  << ",period=" << period << ",idx=" << f.idx
                  << ",p0=" << f.p0 << ",p1=" << f.p1
                  << ",expected=" << f.expected << ",got=" << f.got
                  << "\n";
    }
}

void PrintRow(uint32_t bits, const std::string &mode, uint32_t period,
              int cases, const RunResult &result, double strict_avg_ms)
{
    const double error_rate =
        cases == 0 ? 0.0 : static_cast<double>(result.errors) / cases;
    const double speedup =
        result.avg_ms == 0.0 ? 0.0 : strict_avg_ms / result.avg_ms;
    std::cout << std::fixed << std::setprecision(6)
              << bits << "," << mode << "," << period << "," << cases
              << "," << result.avg_ms << "," << result.errors << ","
              << error_rate << "," << speedup << ","
              << result.stats.total_terms << ","
              << result.stats.zero_skipped << ","
              << result.stats.periodic_skipped << ","
              << result.stats.executed_terms << "\n";
    PrintFailures(bits, mode, period, result);
}

template <class P>
void RunBits(uint32_t plain_bits, const Args &args)
{
    SeedTFHE(0x4554484d53424650ULL + plain_bits);
    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    SetupEvalKey<P>(sk, ek);
    auto inputs = MakeInputs<P>(plain_bits, args.cases, sk);

    const RunResult strict = RunStrict<P>(inputs, plain_bits, sk, ek);
    PrintRow(plain_bits, "strict", 0, args.cases, strict, strict.avg_ms);

    if (args.run_weighted) {
        bool have_best = false;
        uint32_t best_period = 0;
        RunResult best;
        for (uint32_t period : args.periods) {
            RunResult result = RunPruned<P>(
                inputs, plain_bits, sk, ek,
                PrunedETHMSBMode::WeightedApprox, period);
            PrintRow(plain_bits, ModeName(PrunedETHMSBMode::WeightedApprox),
                     period, args.cases, result, strict.avg_ms);
            if (!have_best || result.avg_ms < best.avg_ms) {
                have_best = true;
                best_period = period;
                best = result;
            }
        }
        if (have_best)
            PrintRow(plain_bits, "best-weighted-approx", best_period,
                     args.cases, best, strict.avg_ms);
    }

    if (args.run_paper) {
        for (uint32_t period : args.periods) {
            RunResult result = RunPruned<P>(
                inputs, plain_bits, sk, ek, PrunedETHMSBMode::PaperQHalf3PBS,
                period);
            PrintRow(plain_bits, ModeName(PrunedETHMSBMode::PaperQHalf3PBS),
                     period, args.cases, result, strict.avg_ms);
        }
    }
}

}  // namespace

int main(int argc, char **argv)
{
    const Args args = ParseArgs(argc, argv);
    std::cout << "bits,mode,period,cases,avg_ms,errors,error_rate,"
                 "speedup_vs_strict,total_terms,zero_skipped,"
                 "periodic_skipped,executed_terms\n";
    RunBits<TFHEpp::lvl1param>(4, args);
    RunBits<TFHEpp::lvl1param>(8, args);
    RunBits<TFHEpp::lvl2param>(16, args);
    RunBits<TFHEpp::lvl2param>(32, args);
    return 0;
}
