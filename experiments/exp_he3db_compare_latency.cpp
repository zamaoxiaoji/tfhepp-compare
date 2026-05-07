#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"

namespace {

struct Args {
    int argc;
    char** argv;

    std::string get(const std::string& key, const std::string& fallback) const {
        for (int i = 1; i + 1 < argc; i++) {
            if (argv[i] == key) return argv[i + 1];
        }
        return fallback;
    }

    bool has(const std::string& key) const {
        for (int i = 1; i < argc; i++) {
            if (argv[i] == key) return true;
        }
        return false;
    }

    int get_int(const std::string& key, int fallback) const {
        return std::stoi(get(key, std::to_string(fallback)));
    }

    std::uint64_t get_u64(const std::string& key, std::uint64_t fallback) const {
        return static_cast<std::uint64_t>(std::stoull(get(key, std::to_string(fallback))));
    }
};

std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> out;
    std::size_t start = 0;
    while (start <= s.size()) {
        const std::size_t comma = s.find(',', start);
        const auto token = s.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!token.empty()) out.push_back(std::stoi(token));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return out;
}

std::vector<std::string> parse_ops(const std::string& s) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start <= s.size()) {
        const std::size_t comma = s.find(',', start);
        const auto token = s.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!token.empty()) out.push_back(token);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return out;
}

struct Stats {
    double avg = 0;
    double median = 0;
    double min = 0;
    double max = 0;
    double stddev = 0;
};

Stats stats(std::vector<double> xs) {
    Stats s;
    if (xs.empty()) return s;
    std::sort(xs.begin(), xs.end());
    s.min = xs.front();
    s.max = xs.back();
    s.median = xs[xs.size() / 2];
    if (xs.size() % 2 == 0) s.median = (xs[xs.size() / 2 - 1] + xs[xs.size() / 2]) / 2.0;
    s.avg = std::accumulate(xs.begin(), xs.end(), 0.0) / static_cast<double>(xs.size());
    double sq = 0;
    for (const auto x : xs) sq += (x - s.avg) * (x - s.avg);
    s.stddev = std::sqrt(sq / static_cast<double>(xs.size()));
    return s;
}

bool expected(std::uint64_t a, std::uint64_t b, const std::string& op) {
    if (op == "lt") return a < b;
    if (op == "gt") return a > b;
    if (op == "le") return a <= b;
    if (op == "ge") return a >= b;
    if (op == "eq") return a == b;
    throw std::invalid_argument("unsupported op: " + op);
}

template <typename P>
void run_for_precision(int L, int runs, int warmup, std::uint64_t seed,
                       const std::vector<std::string>& ops,
                       HEDB::TFHESecretKey& sk, HEDB::TFHEEvalKey& ek) {
    using T = typename P::T;
    const int scale_bits = std::numeric_limits<T>::digits - L - 1;
    const double scale = std::ldexp(1.0, scale_bits);
    const std::uint64_t max_value = L == 64 ? ~std::uint64_t{0} : ((std::uint64_t{1} << L) - 1);

    for (const auto& op : ops) {
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint64_t> dist(0, max_value);

        for (int i = 0; i < warmup; i++) {
            const auto a64 = dist(rng);
            const auto b64 = dist(rng);
            auto ca = TFHEpp::tlweSymInt32Encrypt<P>(static_cast<T>(a64), P::α, scale, sk.key.get<P>());
            auto cb = TFHEpp::tlweSymInt32Encrypt<P>(static_cast<T>(b64), P::α, scale, sk.key.get<P>());
            HEDB::TLWELvl1 cres;
            constexpr bool result_type = LOGIC;
            if (op == "lt") HEDB::less_than<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "gt") HEDB::greater_than<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "le") HEDB::less_than_equal<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "ge") HEDB::greater_than_equal<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "eq") HEDB::equal<P>(ca, cb, cres, L, ek, result_type);
        }

        std::vector<double> latencies;
        int correct = 0;
        int fail = 0;
        for (int i = 0; i < runs; i++) {
            const auto a64 = dist(rng);
            const auto b64 = dist(rng);
            auto ca = TFHEpp::tlweSymInt32Encrypt<P>(static_cast<T>(a64), P::α, scale, sk.key.get<P>());
            auto cb = TFHEpp::tlweSymInt32Encrypt<P>(static_cast<T>(b64), P::α, scale, sk.key.get<P>());
            HEDB::TLWELvl1 cres;
            constexpr bool result_type = LOGIC;
            const auto start = std::chrono::steady_clock::now();
            if (op == "lt") HEDB::less_than<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "gt") HEDB::greater_than<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "le") HEDB::less_than_equal<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "ge") HEDB::greater_than_equal<P>(ca, cb, cres, L, ek, result_type);
            else if (op == "eq") HEDB::equal<P>(ca, cb, cres, L, ek, result_type);
            const auto end = std::chrono::steady_clock::now();
            const auto got = TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(cres, sk.key.lvl1);
            const bool exp = expected(a64, b64, op);
            if (static_cast<bool>(got) == exp) correct++;
            else {
                fail++;
                if (fail <= 20) {
                    std::cerr << "# he3db_failure L=" << L << " op=" << op
                              << " a=" << a64 << " b=" << b64
                              << " expected=" << exp << " got=" << got << "\n";
                }
            }
            latencies.push_back(std::chrono::duration<double, std::milli>(end - start).count());
        }

        const auto s = stats(latencies);
        std::cout << "HE3DB," << L << "," << op << "," << runs << "," << warmup
                  << "," << correct << "," << fail << "," << std::fixed
                  << std::setprecision(4) << s.avg << "," << s.median << ","
                  << s.min << "," << s.max << "," << s.stddev << ","
                  << seed << "\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Args args{argc, argv};
    const std::string bits_arg =
        args.has("--bits") ? args.get("--bits", "") : args.get("--L", "4,8,16,24,32");
    const std::string ops_arg =
        args.has("--ops") ? args.get("--ops", "") : args.get("--op", "lt,gt,le,ge,eq");
    const auto bits = parse_int_list(bits_arg);
    const auto ops = parse_ops(ops_arg);
    const int runs = args.get_int("--runs", 1);
    const int warmup = args.get_int("--warmup", 0);
    const auto seed = args.get_u64("--seed", 123456789);

    std::cout << "# experiment: he3db_compare_latency\n";
    std::cout << "# bits: " << bits_arg << "\n";
    std::cout << "# ops: " << ops_arg << "\n";
    std::cout << "# runs: " << runs << "\n";
    std::cout << "# warmup: " << warmup << "\n";
    std::cout << "# seed: " << seed << "\n";
    std::cout << "# timing_scope: comparison_call_only; keygen/encryption/decryption excluded\n";
    std::cout << "scheme,L,op,runs,warmup,correct,fail,latency_ms_avg,"
                 "latency_ms_median,latency_ms_min,latency_ms_max,"
                 "latency_ms_std,seed\n";

    HEDB::TFHESecretKey sk;
    HEDB::TFHEEvalKey ek;
    ek.emplacebkfft<HEDB::Lvl01>(sk);
    ek.emplacebkfft<HEDB::Lvl02>(sk);
    ek.emplaceiksk<HEDB::Lvl20>(sk);
    ek.emplaceiksk<HEDB::Lvl10>(sk);
    ek.emplaceiksk<HEDB::Lvl21>(sk);

    for (const int L : bits) {
        if (L <= 0 || L > 32) {
            std::cerr << "# invalid L=" << L << "\n";
            continue;
        }
        if (L <= 9) run_for_precision<HEDB::Lvl1>(L, runs, warmup, seed, ops, sk, ek);
        else run_for_precision<HEDB::Lvl2>(L, runs, warmup, seed, ops, sk, ek);
    }
}
