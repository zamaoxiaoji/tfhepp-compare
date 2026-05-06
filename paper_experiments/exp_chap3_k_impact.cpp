#include "algorithms.hpp"

#include <chrono>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

using brP = br_lvl22param;
using P = TFHEpp::lvl2param;

template <class TorusP>
typename TorusP::T EncodeMessage(std::uint64_t message, int precision_bits) {
    return static_cast<typename TorusP::T>(message)
           << (std::numeric_limits<typename TorusP::T>::digits - precision_bits);
}

template <class TorusP>
int DecodeSign(const TFHEpp::TLWE<TorusP>& ct, const TFHEpp::SecretKey& sk) {
    const auto phase = TFHEpp::tlweSymPhase<TorusP>(ct, sk.key.get<TorusP>());
    return static_cast<std::make_signed_t<typename TorusP::T>>(phase) < 0 ? 1 : 0;
}

double MsSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

std::vector<std::uint64_t> BoundarySamples(int p, int radius) {
    const std::uint64_t threshold = std::uint64_t{1} << (p - 1);
    const std::uint64_t max = (std::uint64_t{1} << p) - 1;
    std::vector<std::uint64_t> samples{0, 1, threshold - 1, threshold, max};
    for (int d = -radius; d <= radius; d++) {
        const auto v = static_cast<long long>(threshold) + d;
        if (v >= 0 && static_cast<std::uint64_t>(v) <= max)
            samples.push_back(static_cast<std::uint64_t>(v));
    }
    return samples;
}

void AppendRandomSamples(
    std::vector<std::uint64_t>& samples,
    int p,
    int trials,
    std::uint64_t seed) {
    if (trials <= 0) return;
    const std::uint64_t max = (std::uint64_t{1} << p) - 1;
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::uint64_t> dist(0, max);
    for (int i = 0; i < trials; i++) samples.push_back(dist(rng));
}

std::vector<int> ParseIntList(const std::string& text) {
    std::vector<int> values;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) throw std::invalid_argument("empty k-list entry");
        values.push_back(std::stoi(item));
    }
    if (values.empty()) throw std::invalid_argument("empty k-list");
    return values;
}

}  // namespace

int main(int argc, char** argv) {
    int p = 12;
    int radius = 4;
    int random_trials = 0;
    std::uint64_t seed = 0xC0FFEE;
    std::vector<int> k_values;
    std::unique_ptr<std::ofstream> output_file;
    std::ostream* output = nullptr;

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--k" && i + 1 < argc) {
                k_values = {std::stoi(argv[++i])};
            } else if (arg == "--k-list" && i + 1 < argc) {
                k_values = ParseIntList(argv[++i]);
            } else if (arg == "--sweep") {
                k_values.clear();
            } else if (arg == "--radius" && i + 1 < argc) {
                radius = std::stoi(argv[++i]);
            } else if ((arg == "--trials" || arg == "--random") && i + 1 < argc) {
                random_trials = std::stoi(argv[++i]);
            } else if (arg == "--seed" && i + 1 < argc) {
                seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            } else if (arg == "--output" || arg == "--output-file") {
                output = PaperReview::OpenOptionalOutputFile(
                    i, argc, argv, output_file);
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap3_k_impact [--sweep|--k K|--k-list K1,K2] "
                    "[--radius R] [--trials N] [--seed S] [--output PATH]");
            }
        }

        const auto cfg = PaperReview::Chapter3MetaPBSConfig();
        p = MetaPBS2::MessagePrecisionFromPowerOfTwoModulus(cfg.t);
        if (k_values.empty())
            k_values = {p - 2, p - 3, p - 4, p - 5, 1};
        for (const int k : k_values) {
            if (k <= 0 || k >= p)
                throw std::invalid_argument("k must satisfy 1 <= k < p");
        }

        TFHEpp::SecretKey sk;
        auto bk = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
        TFHEpp::bkfftgen<brP>(*bk, sk);

        std::vector<MetaPBS2::TruncRepeatKey<P>> trkeys;
        for (const auto& round : cfg.rounds)
            trkeys.push_back(MetaPBS2::GenerateTruncRepeatKey<P>(
                sk.key.get<P>(), round.beta));

        auto samples = BoundarySamples(p, radius);
        AppendRandomSamples(samples, p, random_trials, seed);
        auto write_line = [&](const std::string& text) {
            PaperReview::WriteOutputLine(std::cout, output, text);
        };
        write_line("=== Chapter 3 GapMSB zeroed-bit k impact ===");
        write_line("method=encrypt m as TFHE lvl2 TLWE, run GapMSB_k, decrypt MSB sign bit, compare with plaintext bit_0(m)");
        write_line("metric,value");
        write_line("operator,gap_msb_k_impact");
        write_line("security=" + PaperReview::SecuritySummary());
        {
            std::ostringstream line;
            line << "p=" << p << " radius=" << radius
                 << " random_trials=" << random_trials
                 << " seed=" << seed;
            write_line(line.str());
        }
        write_line("p," + std::to_string(p));
        write_line("radius," + std::to_string(radius));
        write_line("random_trials," + std::to_string(random_trials));
        write_line("seed," + std::to_string(seed));
        write_line("p,k,w_k,period,message,expected,got,ms,prune_skipped,prune_total,prune_rate");

        int total = 0;
        int passed = 0;
        for (const int k : k_values) {
            const int w_k = 1 << (p - 1 - k);
            const int period = MetaPBS2::BitExtractPeriodFromChapterBit(p, k);
            int k_total = 0;
            int k_passed = 0;
            double k_ms = 0.0;
            MetaPBS2::BlindRotatePruneStats k_stats{};
            for (const auto m : samples) {
                TFHEpp::TLWE<P> ct{};
                TFHEpp::tlweSymEncrypt<P>(
                    ct, EncodeMessage<P>(m, p), P::α, sk.key.get<P>());

                MetaPBS2::GapMSBOptions options;
                options.p = p;
                options.k = k;
                options.enable_periodic_pruning = true;

                MetaPBS2::BlindRotatePruneStats stats{};
                const auto start = std::chrono::steady_clock::now();
                const auto out = MetaPBS2::GapMSB<brP>(
                    ct, *bk, trkeys, cfg, options, &stats);
                const double ms = MsSince(start);

                const int got = DecodeSign<P>(out, sk);
                const int expected = m >= (std::uint64_t{1} << (p - 1)) ? 1 : 0;
                total++;
                k_total++;
                k_ms += ms;
                k_stats.skipped += stats.skipped;
                k_stats.total += stats.total;
                if (got == expected) {
                    passed++;
                    k_passed++;
                }

                std::ostringstream row;
                row << p << "," << k << "," << w_k << "," << period
                    << "," << m << "," << expected << "," << got << ","
                    << std::fixed << std::setprecision(2) << ms << ","
                    << stats.skipped << "," << stats.total << ","
                    << std::setprecision(6) << stats.prune_rate();
                write_line(row.str());
            }
            const double failure_rate =
                k_total == 0 ? 0.0 : 1.0 - static_cast<double>(k_passed) / k_total;
            const double avg_ms = k_total == 0 ? 0.0 : k_ms / k_total;
            std::ostringstream summary;
            summary << "k_summary," << p << "," << k
                    << ",passed=" << k_passed
                    << ",total=" << k_total
                    << ",failure_rate=" << std::setprecision(6) << failure_rate
                    << ",avg_ms=" << avg_ms
                    << ",first_round_prune_rate=" << k_stats.prune_rate();
            write_line(summary.str());
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
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_k_impact failed: " << e.what() << "\n";
        return 1;
    }
}
