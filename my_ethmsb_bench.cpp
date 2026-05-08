#include "my_ethmsb_fixed.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using my_ethmsb::BOOL_ONE;
using my_ethmsb::P;
using my_ethmsb::Torus;
using my_ethmsb::my_lvl22pbsparam;

struct Options {
    std::vector<int> bits{8, 16, 24, 32};
    int kappa = 5;
    int trials = 10;
    uint64_t seed = 0x4554484d5342554cULL;
    bool print_keygen_timing = true;
};

struct PlainPair {
    Torus a;
    Torus b;
};

struct EncPair {
    PlainPair plain;
    TFHEpp::TLWE<P> a;
    TFHEpp::TLWE<P> b;
};

std::string capture_cmd(const char* cmd)
{
    char buffer[256];
    std::string out;
    FILE* pipe = popen(cmd, "r");
    if (pipe == nullptr) return "unavailable";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) out += buffer;
    pclose(pipe);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r'))
        out.pop_back();
    return out.empty() ? "unavailable" : out;
}

std::vector<int> parse_bits_csv(const std::string& value)
{
    std::vector<int> bits;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, ',')) {
        const int bit = std::stoi(item);
        if (bit != 8 && bit != 16 && bit != 24 && bit != 32)
            throw std::invalid_argument("--bits supports only 8,16,24,32");
        bits.push_back(bit);
    }
    if (bits.empty()) throw std::invalid_argument("--bits needs a value");
    return bits;
}

Options parse_args(const int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* what) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(what);
            return argv[++i];
        };
        if (arg == "--bits") {
            opt.bits = parse_bits_csv(need_value("--bits needs csv"));
        }
        else if (arg == "--all") {
            opt.bits = {8, 16, 24, 32};
        }
        else if (arg == "--kappa") {
            opt.kappa = std::stoi(need_value("--kappa needs 3/4/5"));
            if (opt.kappa != 3 && opt.kappa != 4 && opt.kappa != 5)
                throw std::invalid_argument("--kappa supports only 3/4/5");
        }
        else if (arg == "--trials") {
            opt.trials = std::stoi(need_value("--trials needs N"));
        }
        else if (arg == "--seed") {
            opt.seed = std::stoull(need_value("--seed needs integer"));
        }
        else if (arg == "--no-keygen-timing") {
            opt.print_keygen_timing = false;
        }
        else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    if (opt.trials < 1) throw std::invalid_argument("--trials must be >= 1");
    return opt;
}

std::string hex64(const Torus value)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << value
       << std::dec;
    return os.str();
}

Torus phase_of(const TFHEpp::TLWE<P>& ct, const TFHEpp::SecretKey& sk)
{
    return TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
}

std::vector<PlainPair> targeted_pairs_for_bits(const int bits)
{
    const Torus top = Torus{1} << bits;
    return {{0, 0},
            {0, 1},
            {1, 0},
            {0, top - 1},
            {top - 1, 0},
            {(Torus{1} << (bits - 1)) - 1, Torus{1} << (bits - 1)},
            {Torus{1} << (bits - 1), (Torus{1} << (bits - 1)) - 1},
            {top - 1, top - 1}};
}

int run_targeted_correctness(
    const int bits, const int kappa, const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft)
{
    int failures = 0;
    for (const PlainPair& pair : targeted_pairs_for_bits(bits)) {
        const auto ca =
            my_ethmsb::encrypt_unsigned_for_compare(pair.a, bits, sk);
        const auto cb =
            my_ethmsb::encrypt_unsigned_for_compare(pair.b, bits, sk);
        alignas(64) TFHEpp::TLWE<P> lt_ct;
        alignas(64) TFHEpp::TLWE<P> gt_ct;
        alignas(64) TFHEpp::TLWE<P> eq_ct;
        my_ethmsb::lt(lt_ct, ca, cb, bits, kappa, bkfft);
        my_ethmsb::gt(gt_ct, ca, cb, bits, kappa, bkfft);
        my_ethmsb::eq(eq_ct, ca, cb, bits, kappa, bkfft);
        const bool got_lt = my_ethmsb::decrypt_arith_bit(lt_ct, sk);
        const bool got_gt = my_ethmsb::decrypt_arith_bit(gt_ct, sk);
        const bool got_eq = my_ethmsb::decrypt_arith_bit(eq_ct, sk);
        if (got_lt != (pair.a < pair.b) || got_gt != (pair.a > pair.b) ||
            got_eq != (pair.a == pair.b)) {
            ++failures;
            std::cerr << "BENCH_TARGETED_FAIL bits=" << bits
                      << " k=" << bits + 1 << " kappa=" << kappa
                      << " a=" << pair.a << " b=" << pair.b
                      << " got_lt_gt_eq=" << got_lt << got_gt << got_eq
                      << " expected_lt_gt_eq=" << (pair.a < pair.b)
                      << (pair.a > pair.b) << (pair.a == pair.b)
                      << " phase_lt=" << hex64(phase_of(lt_ct, sk))
                      << " phase_gt=" << hex64(phase_of(gt_ct, sk))
                      << " phase_eq=" << hex64(phase_of(eq_ct, sk)) << '\n';
        }
    }
    return failures;
}

std::vector<EncPair> make_inputs(const int bits, const int trials,
                                 const uint64_t seed,
                                 const TFHEpp::SecretKey& sk)
{
    std::mt19937_64 rng(seed + static_cast<uint64_t>(bits) * 1315423911ULL);
    std::uniform_int_distribution<Torus> dist(0, (Torus{1} << bits) - 1);
    std::vector<EncPair> inputs;
    inputs.reserve(trials);
    for (int i = 0; i < trials; ++i) {
        PlainPair pair{dist(rng), dist(rng)};
        inputs.push_back(
            {pair, my_ethmsb::encrypt_unsigned_for_compare(pair.a, bits, sk),
             my_ethmsb::encrypt_unsigned_for_compare(pair.b, bits, sk)});
    }
    return inputs;
}

struct Summary {
    double mean = 0;
    double median = 0;
    double p95 = 0;
    double min = 0;
    double max = 0;
};

Summary summarize(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    Summary s;
    s.min = values.front();
    s.max = values.back();
    s.mean = std::accumulate(values.begin(), values.end(), 0.0) /
             static_cast<double>(values.size());
    s.median = values[values.size() / 2];
    s.p95 = values[std::min(values.size() - 1,
                            static_cast<size_t>(0.95 * values.size()))];
    return s;
}

template <class Fn>
void run_op(const int bits, const int kappa, const std::string& op,
            const int pbs_count, const std::vector<EncPair>& inputs,
            const TFHEpp::SecretKey& sk, Fn&& fn)
{
    std::vector<double> times;
    times.reserve(inputs.size());
    int failures = 0;
    for (const EncPair& input : inputs) {
        alignas(64) TFHEpp::TLWE<P> out;
        const auto start = std::chrono::steady_clock::now();
        fn(out, input.a, input.b);
        const auto end = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration<double, std::milli>(end - start)
                            .count());
        const bool actual = my_ethmsb::decrypt_arith_bit(out, sk);
        bool expected = false;
        if (op == "lt")
            expected = input.plain.a < input.plain.b;
        else if (op == "gt")
            expected = input.plain.a > input.plain.b;
        else if (op == "eq")
            expected = input.plain.a == input.plain.b;
        failures += actual == expected ? 0 : 1;
    }
    const Summary s = summarize(std::move(times));
    std::cout << bits << ',' << (bits + 1) << ',' << kappa << ','
              << pbs_count << ',' << op << ',' << inputs.size() << ','
              << s.mean << ',' << s.median << ',' << s.p95 << ',' << s.min
              << ',' << s.max << ',' << failures << '\n';
}

}  // namespace

int main(int argc, char** argv)
{
    try {
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        const Options opt = parse_args(argc, argv);
        std::cout << "# git_describe=" << capture_cmd("git describe --tags --always")
                  << '\n';
        std::cout << "# git_rev_parse=" << capture_cmd("git rev-parse HEAD")
                  << '\n';
        std::cout << "# USE_CONCRETE="
#ifdef USE_CONCRETE
                  << "ON\n";
#else
                  << "OFF\n";
#endif
        std::cout << "# USE_FFTW3="
#ifdef USE_FFTW3
                  << "ON\n";
#else
                  << "OFF\n";
#endif

        const auto keygen_start = std::chrono::steady_clock::now();
        TFHEpp::SecretKey sk;
        auto bkfft = std::make_unique_for_overwrite<
            TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>>();
        TFHEpp::bkfftgen<my_lvl22pbsparam>(*bkfft, sk.key.get<P>(),
                                           sk.key.get<P>());
        const auto keygen_end = std::chrono::steady_clock::now();
        if (opt.print_keygen_timing) {
            std::cout << "# keygen_ms="
                      << std::chrono::duration<double, std::milli>(
                             keygen_end - keygen_start)
                             .count()
                      << '\n';
        }

        std::cout
            << "bits,k,kappa,pbs_count,op,trials,mean_ms,median_ms,p95_ms,"
               "min_ms,max_ms,failures\n";
        for (const int bits : opt.bits) {
            const int targeted_failures =
                run_targeted_correctness(bits, opt.kappa, sk, *bkfft);
            if (targeted_failures != 0) {
                std::cerr << "targeted correctness failed for bits=" << bits
                          << ", benchmark latency not claimed\n";
                return 1;
            }
            const std::vector<EncPair> inputs =
                make_inputs(bits, opt.trials, opt.seed, sk);
            const int msb_pbs_count =
                my_ethmsb::ethmsb_pbs_count(bits + 1, opt.kappa);
            run_op(bits, opt.kappa, "lt", msb_pbs_count, inputs, sk,
                   [&](TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
                       const TFHEpp::TLWE<P>& b) {
                       my_ethmsb::lt(out, a, b, bits, opt.kappa, *bkfft);
                   });
            run_op(bits, opt.kappa, "gt", msb_pbs_count, inputs, sk,
                   [&](TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
                       const TFHEpp::TLWE<P>& b) {
                       my_ethmsb::gt(out, a, b, bits, opt.kappa, *bkfft);
                   });
            run_op(bits, opt.kappa, "eq", 2 * msb_pbs_count, inputs, sk,
                   [&](TFHEpp::TLWE<P>& out, const TFHEpp::TLWE<P>& a,
                       const TFHEpp::TLWE<P>& b) {
                       my_ethmsb::eq(out, a, b, bits, opt.kappa, *bkfft);
                   });
        }
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: benchmark\n";
        std::cerr << "what_failed: benchmark harness threw an exception\n";
        std::cerr << "commands_run:\n  ./build/my_ethmsb_bench\n";
        std::cerr << "relevant_files:\n  my_ethmsb_fixed.hpp\n  "
                     "my_ethmsb_bench.cpp\n";
        std::cerr << "exact_error_or_failure:\n" << ex.what() << '\n';
        std::cerr << "tfhepp_version:\n";
        std::cerr << "  git_describe: "
                  << capture_cmd("git describe --tags --always") << '\n';
        std::cerr << "  git_rev_parse: " << capture_cmd("git rev-parse HEAD")
                  << '\n';
        std::cerr << "params:\n";
        std::cerr << "  sizeof_T: " << sizeof(P::T) << '\n';
        std::cerr << "  P_n: " << P::n << '\n';
        std::cerr << "  P_k: " << P::k << '\n';
        std::cerr << "  P_l: " << P::l << '\n';
        std::cerr << "  P_l_a: " << P::lₐ << '\n';
        std::cerr << "  Bgbit: " << P::Bgbit << '\n';
        std::cerr << "  alpha: " << std::setprecision(18) << P::α << '\n';
        std::cerr << "  key_min_max: " << P::key_value_min << ','
                  << P::key_value_max << '\n';
        std::cerr << "test_case_if_any:\n  none\n";
        std::cerr << "api_signatures_needed:\n";
        std::cerr << "  GateBootstrappingTLWE2TLWEFFT / bkfftgen signatures "
                     "are in include/gatebootstrapping.hpp and "
                     "include/evalkeygens.hpp\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
