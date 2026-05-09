#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "my_ethmsb_pruned_fused_final.hpp"

namespace {

using namespace my_ethmsb_pruned_fused_final;

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

Context make_context(const std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    const auto keys = my_ethmsb_params::make_h3compat_keys(rng);
    Context ctx;
    ctx.key2 = keys.key2;
    ctx.key0 = keys.key0;
    ctx.iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    return ctx;
}

TFHEpp::TLWE<P2> encrypt_m(const Torus m, const int p, const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode_unsigned_phase(m, p), ctx.key2);
}

bool out_bit(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_value(phase_l2(ct, ctx.key2), BOOL_ONE);
}

std::vector<Torus> cases(const int p, const int k)
{
    std::set<Torus> vals = {0,
                            1,
                            (Torus{1} << (p - 1)) - 1,
                            Torus{1} << (p - 1),
                            (Torus{1} << (p - 1)) + 1,
                            (Torus{1} << p) - 1};
    const Torus W = bit_weight_msb_index(p, k);
    for (const Torus v : {W - 1, W, W + 1, 2 * W - 1, 2 * W, 2 * W + 1})
        if (v < (Torus{1} << p)) vals.insert(v);
    return {vals.begin(), vals.end()};
}

}  // namespace

int main(int argc, char** argv)
{
    int p = 9;
    int k = 1;
    std::uint64_t seed = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--p" && i + 1 < argc)
            p = std::stoi(argv[++i]);
        else if (a == "--k" && i + 1 < argc)
            k = std::stoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc)
            seed = std::stoull(argv[++i]);
        else
            throw std::runtime_error("unknown argument: " + a);
    }

    const Context ctx = make_context(seed);
    std::cout << "mode,p,k,m,expected,three_pbs_actual,fused_actual,"
                 "three_pbs_pass,fused_pass\n";
    int failures_3pbs = 0;
    int failures_fused = 0;
    for (const Torus m : cases(p, k)) {
        const bool expected = ((m >> (p - 1)) & 1) != 0;
        const TFHEpp::TLWE<P2> ct = encrypt_m(m, p, ctx);

        TFHEpp::TLWE<P2> out3;
        PrunedGapStats stats3;
        gapmsb_pruned_bitk_l2_to_l2(out3, ct, p, k, BOOL_ONE, *ctx.iksk,
                                    *ctx.bkfft, &stats3);
        const bool a3 = out_bit(out3, ctx);

        TFHEpp::TLWE<P2> outf;
        FusedGapStats statsf;
        gapmsb_pruned_bitk_fused_final_l2_to_l2(
            outf, ct, p, k, BOOL_ONE, FusedInputRoundingMode::DeltaHalf,
            *ctx.iksk, *ctx.bkfft, &statsf);
        const bool af = out_bit(outf, ctx);

        if (a3 != expected) ++failures_3pbs;
        if (af != expected) ++failures_fused;
        std::cout << "real_bitextract," << p << ',' << k << ',' << m << ','
                  << expected << ',' << a3 << ',' << af << ','
                  << (a3 == expected) << ',' << (af == expected) << "\n";
    }
    std::cerr << "FUSED_VS_CONVERSION_AB failures_3pbs=" << failures_3pbs
              << " failures_fused=" << failures_fused << "\n";
    return 0;
}
