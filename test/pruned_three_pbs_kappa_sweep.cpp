#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../comparison/comparison.h"

using namespace tfhepp_compare;

namespace
{
    struct Options {
        int                   trials = 20;
        std::vector<uint32_t> bits = {16, 32};
        std::vector<uint32_t> kappas = {3, 4, 5, 6, 7, 8};
        bool                  fast = true;
    };

    std::vector<uint32_t> split_u32(const std::string &s)
    {
        std::vector<uint32_t> out;
        std::stringstream     ss(s);
        std::string           item;
        while (std::getline(ss, item, ','))
            if (!item.empty()) out.push_back(static_cast<uint32_t>(std::stoul(item)));
        return out;
    }

    Options parse_options(int argc, char **argv)
    {
        Options opts;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto need = [&](const char *name) -> const char * {
                if (i + 1 >= argc) {
                    std::cerr << "missing value for " << name << "\n";
                    std::exit(2);
                }
                return argv[++i];
            };
            if (arg == "--trials")
                opts.trials = std::stoi(need("--trials"));
            else if (arg == "--bits")
                opts.bits = split_u32(need("--bits"));
            else if (arg == "--kappa-list")
                opts.kappas = split_u32(need("--kappa-list"));
            else if (arg == "--standard")
                opts.fast = false;
            else if (arg == "--fast")
                opts.fast = true;
            else {
                std::cerr << "unknown option: " << arg << "\n";
                std::exit(2);
            }
        }
        if (opts.trials < 1) opts.trials = 1;
        return opts;
    }

    template <class P>
    TFHEpp::TLWE<P> encrypt_message(typename P::T value, uint32_t plain_bits,
                                    const TFHESecretKey &sk)
    {
        const uint32_t scale_bits =
            std::numeric_limits<typename P::T>::digits - plain_bits - 1;
        return tlweSymInt32Encrypt<P>(value, P::α, std::pow(2., scale_bits),
                                      sk.key.get<P>());
    }

    template <class P>
    TFHEpp::TLWE<P> subtract(const TFHEpp::TLWE<P> &a,
                             const TFHEpp::TLWE<P> &b)
    {
        TFHEpp::TLWE<P> out;
        for (size_t i = 0; i <= P::k * P::n; i++) out[i] = a[i] - b[i];
        return out;
    }

    std::string append_step(const std::string &prefix, uint32_t p)
    {
        return prefix.empty() ? std::to_string(p)
                              : prefix + "->" + std::to_string(p);
    }

    std::string lvl1_schedule(uint32_t p, uint32_t kappa,
                              uint32_t gap_parent_bits)
    {
        std::string out;
        for (;;) {
            out = append_step(out, p);
            if (p <= kappa || (gap_parent_bits != 0 && p <= kappa + 3))
                return out + "F";
            gap_parent_bits = p;
            p -= kappa;
        }
    }

    std::string lvl2_schedule(uint32_t p, uint32_t kappa)
    {
        std::string out;
        uint32_t    gap_parent_bits = 0;
        for (;;) {
            if (gap_parent_bits != 0 && p <= kappa + 9)
                return append_step(out, p) + "L1:" +
                       lvl1_schedule(p, kappa, gap_parent_bits);
            if (p <= 6) return append_step(out, p) + "F";
            if (p <= std::min<uint32_t>(10, kappa + 4))
                return append_step(out, p) + "L1:" +
                       lvl1_schedule(p, kappa, gap_parent_bits);
            out = append_step(out, p);
            gap_parent_bits = p;
            p -= kappa;
        }
    }

    struct TrialInputLvl2 {
        Lvl2::T          p0 = 0;
        Lvl2::T          p1 = 0;
        TFHEpp::TLWE<Lvl2> c0 = {};
        TFHEpp::TLWE<Lvl2> c1 = {};
    };

    void run_lvl2(uint32_t plain_bits, const Options &opts)
    {
        TFHESecretKey sk;
        TFHEEvalKey   ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplacebkfft<Lvl02>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        ek.emplaceiksk<Lvl20>(sk);
        ek.emplaceiksk<Lvl21>(sk);
        const auto micro_pack = three_pbs::GenerateFastB2AEvalKeyPack(sk, true);

        std::default_random_engine engine(0x6b617070u ^
                                          (plain_bits * 0x9e3779b9u));
        std::uniform_int_distribution<Lvl2::T> message(
            0, (Lvl2::T(1) << (plain_bits - 1)) - 1);

        std::vector<TrialInputLvl2> inputs;
        inputs.reserve(opts.trials);
        for (int t = 0; t < opts.trials; ++t) {
            TrialInputLvl2 in;
            in.p0 = message(engine);
            in.p1 = message(engine);
            in.c0 = encrypt_message<Lvl2>(in.p0, plain_bits, sk);
            in.c1 = encrypt_message<Lvl2>(in.p1, plain_bits, sk);
            inputs.push_back(in);
        }

        for (uint32_t kappa : opts.kappas) {
            uint32_t errors = 0;
            bool     first_fail = false;
            Lvl2::T  first_p0 = 0;
            Lvl2::T  first_p1 = 0;
            Lvl1::T  first_decoded = 0;
            double   total_ms = 0.0;

            for (const TrialInputLvl2 &in : inputs) {
                const auto sub = subtract<Lvl2>(in.c0, in.c1);
                TLWELvl1   out;
                const auto start = std::chrono::steady_clock::now();
                if (opts.fast)
                    three_pbs::HomMSBWithKappa(out, sub, plain_bits + 1, kappa,
                                               ek, micro_pack, LOGIC);
                else
                    three_pbs::HomMSBWithKappa(out, sub, plain_bits + 1, kappa,
                                               ek, LOGIC);
                const auto end = std::chrono::steady_clock::now();
                total_ms +=
                    std::chrono::duration<double, std::milli>(end - start)
                        .count();

                const Lvl1::T decoded =
                    TFHEpp::tlweSymDecrypt<Lvl1>(out, sk.key.lvl1);
                const Lvl1::T expected = static_cast<Lvl1::T>(in.p0 < in.p1);
                if (decoded != expected) {
                    errors++;
                    if (!first_fail) {
                        first_fail = true;
                        first_p0 = in.p0;
                        first_p1 = in.p1;
                        first_decoded = decoded;
                    }
                }
            }

            const double accuracy =
                100.0 * static_cast<double>(opts.trials - errors) /
                static_cast<double>(opts.trials);
            std::cout << plain_bits << "," << kappa << ","
                      << (opts.fast ? "fast" : "standard") << ","
                      << opts.trials << "," << std::fixed
                      << std::setprecision(3) << (total_ms / opts.trials)
                      << "," << errors << "," << std::setprecision(2)
                      << accuracy << "," << lvl2_schedule(plain_bits + 1, kappa);
            if (first_fail)
                std::cout << ",p0=" << static_cast<uint64_t>(first_p0)
                          << " p1=" << static_cast<uint64_t>(first_p1)
                          << " decoded=" << static_cast<uint64_t>(first_decoded);
            std::cout << "\n";
        }
    }
} // namespace

int main(int argc, char **argv)
{
    const Options opts = parse_options(argc, argv);
    std::cout << "bits,kappa,mode,trials,avg_ms,errors,accuracy_percent,"
                 "schedule,first_failure\n";
    for (uint32_t bits : opts.bits) {
        if (bits > 10 && bits <= 32)
            run_lvl2(bits, opts);
        else
            std::cerr << "skip bits=" << bits << " reason=only_lvl2_11_32\n";
    }
    return 0;
}
