#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "../comparison/bitextract_qhalf.h"

using namespace tfhepp_compare;
namespace be = tfhepp_compare::bitextract_qhalf;

namespace
{
    struct Options {
        uint32_t              plain_bits = 5;
        std::vector<uint32_t> lsb_bits = {0, 1, 2, 3};
        bool                  pruned = true;
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
            if (arg == "--plain-bits")
                opts.plain_bits = static_cast<uint32_t>(std::stoul(need("--plain-bits")));
            else if (arg == "--lsb-list")
                opts.lsb_bits = split_u32(need("--lsb-list"));
            else if (arg == "--no-prune")
                opts.pruned = false;
            else if (arg == "--pruned")
                opts.pruned = true;
            else {
                std::cerr << "unknown option: " << arg << "\n";
                std::exit(2);
            }
        }
        return opts;
    }

    TLWELvl1 encrypt_message(uint32_t m, uint32_t plain_bits,
                             const TFHESecretKey &sk)
    {
        const Lvl1::T delta =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - plain_bits);
        return TFHEpp::tlweSymEncrypt<Lvl1>(Lvl1::T(m) * delta, Lvl1::α,
                                            sk.key.lvl1);
    }

    void Boolean2Weight(TLWELvl1 &res, const TLWELvl1 &boolean_bit,
                        Lvl1::T weight, const TFHEEvalKey &ek)
    {
        constexpr Lvl1::T offset =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 6);
        TLWELvl1 tlweoffset = boolean_bit;
        tlweoffset[Lvl1::k * Lvl1::n] += offset;

        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
            res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<Lvl1>(weight));
        res[Lvl1::k * Lvl1::n] += weight;
    }

    uint32_t decrypt_message(const TLWELvl1 &ct, uint32_t plain_bits,
                             const TFHESecretKey &sk)
    {
        const double scale =
            std::ldexp(1.0, std::numeric_limits<Lvl1::T>::digits - plain_bits);
        return tlweSymInt32Decrypt<Lvl1>(ct, scale, sk.key.lvl1);
    }
} // namespace

int main(int argc, char **argv)
{
    const Options opts = parse_options(argc, argv);
    if (opts.plain_bits < 2 || opts.plain_bits > 10) {
        std::cerr << "plain_bits must be in [2,10]\n";
        return 2;
    }

    TFHESecretKey sk;
    TFHEEvalKey   ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);

    std::cout << "plain_bits,lsb_bit,avg_ms,errors\n";

    const uint32_t cases = uint32_t(1) << opts.plain_bits;
    for (uint32_t lsb_bit : opts.lsb_bits) {
        if (lsb_bit + 1 >= opts.plain_bits) {
            std::cerr << "skip lsb_bit=" << lsb_bit
                      << " reason=MSB_or_out_of_range\n";
            continue;
        }

        const uint32_t window_k = opts.plain_bits - 1 - lsb_bit;
        const Lvl1::T  delta =
            Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - opts.plain_bits);
        const Lvl1::T weight = delta << lsb_bit;
        const Lvl1::T b2a_weight = weight >> 1;

        uint32_t errors = 0;
        double   total_ms = 0.0;

        for (uint32_t m = 0; m < cases; ++m) {
            const TLWELvl1 in = encrypt_message(m, opts.plain_bits, sk);
            TLWELvl1       bit;
            be::Trace      trace;

            const auto start = std::chrono::steady_clock::now();
            if (opts.pruned)
                be::BitExtractQHalf_Pruned_Lvl1_CenteredCell(
                    bit, in, opts.plain_bits, window_k, ek, &trace);
            else
                be::BitExtractQHalf_NoPrune_Lvl1_CenteredCell(
                    bit, in, opts.plain_bits, window_k, ek, &trace);

            TLWELvl1 weight_bit;
            Boolean2Weight(weight_bit, bit, b2a_weight, ek);

            TLWELvl1 cleared;
            for (size_t i = 0; i <= Lvl1::k * Lvl1::n; ++i)
                cleared[i] = in[i] - weight_bit[i];
            const auto end = std::chrono::steady_clock::now();
            total_ms +=
                std::chrono::duration<double, std::milli>(end - start).count();

            const uint32_t decoded = decrypt_message(cleared, opts.plain_bits, sk);
            const uint32_t expected = m & ~(uint32_t(1) << lsb_bit);
            if (decoded != expected) errors++;
        }

        std::cout << opts.plain_bits << "," << lsb_bit << "," << std::fixed
                  << std::setprecision(3) << (total_ms / cases) << ","
                  << errors << "\n";
    }

    return 0;
}
