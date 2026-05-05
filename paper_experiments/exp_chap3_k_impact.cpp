#include "algorithms.hpp"

#include <chrono>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
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

}  // namespace

int main(int argc, char** argv) {
    int p = 12;
    int radius = 4;
    std::vector<int> k_values{7, 8, 9, 10, 11};

    try {
        for (int i = 1; i < argc; i++) {
            const std::string arg = argv[i];
            if (arg == "--k" && i + 1 < argc) {
                k_values = {std::stoi(argv[++i])};
            } else if (arg == "--sweep") {
                k_values = {7, 8, 9, 10, 11};
            } else if (arg == "--radius" && i + 1 < argc) {
                radius = std::stoi(argv[++i]);
            } else {
                throw std::invalid_argument(
                    "usage: exp_chap3_k_impact [--sweep|--k K] [--radius R]");
            }
        }

        const auto cfg = PaperReview::Chapter3MetaPBSConfig();
        p = MetaPBS2::MessagePrecisionFromPowerOfTwoModulus(cfg.t);
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

        const auto samples = BoundarySamples(p, radius);
        std::cout << "=== Chapter 3 GapMSB zeroed-bit k impact ===\n";
        std::cout << "method=encrypt m as TFHE lvl2 TLWE, run GapMSB_k, decrypt MSB sign bit, compare with plaintext bit_0(m)\n";
        std::cout << "security=" << PaperReview::SecuritySummary() << "\n";
        std::cout << "p,k,w_k,period,message,expected,got,ms,prune_skipped,prune_total,prune_rate\n";

        int total = 0;
        int passed = 0;
        for (const int k : k_values) {
            const int w_k = 1 << (p - 1 - k);
            const int period = MetaPBS2::BitExtractPeriodFromChapterBit(p, k);
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
                if (got == expected) passed++;

                std::cout << p << "," << k << "," << w_k << "," << period
                          << "," << m << "," << expected << "," << got << ","
                          << std::fixed << std::setprecision(2) << ms << ","
                          << stats.skipped << "," << stats.total << ","
                          << std::setprecision(6) << stats.prune_rate() << "\n";
            }
        }
        std::cout << "summary,passed,total,accuracy\n";
        std::cout << "summary," << passed << "," << total << ","
                  << std::setprecision(6)
                  << (total == 0 ? 0.0 : static_cast<double>(passed) / total)
                  << "\n";
        return passed == total ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap3_k_impact failed: " << e.what() << "\n";
        return 1;
    }
}
