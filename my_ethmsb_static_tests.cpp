#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"

namespace {

using namespace my_ethmsb_params;
using P2 = my_h3_lvl2param;
using P0 = my_h3_lvl0param;
using KS20 = my_h3_lvl20param;
using BR02 = my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Options {
    std::string impl = "h3compat_l20_l02";
    int kappa = 5;
    int exhaustive_k = 16;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

Options parse(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--impl" && i + 1 < argc)
            opt.impl = argv[++i];
        else if (a == "--kappa" && i + 1 < argc)
            opt.kappa = std::atoi(argv[++i]);
        else if ((a == "--exhaustive-k" || a == "--max-trivial-k") &&
                 i + 1 < argc)
            opt.exhaustive_k = std::atoi(argv[++i]);
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (opt.impl != "h3compat_l20_l02")
        throw std::runtime_error(
            "my_ethmsb_static_tests formal target supports h3compat_l20_l02 only");
    return opt;
}

Torus delta(const int k) { return my_ethmsb::delta(k); }

Torus encode(const Torus m, const int k)
{
    return static_cast<Torus>(Wide{m} * Wide{delta(k)});
}

bool closer_to_one(const Torus phase, const Torus out_value)
{
    return my_ethmsb::torus_abs_centered(phase - out_value) <
           my_ethmsb::torus_abs_centered(phase);
}

Context make_context()
{
    std::mt19937_64 rng(0);
    const auto keys = make_h3compat_keys(rng);
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

Torus phase(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return TFHEpp::tlweSymPhase<P2>(ct, ctx.key2);
}

void pbs(TFHEpp::TLWE<P2>& out, TFHEpp::TLWE<P2> in, const Torus offset,
         const Torus out_value, const Context& ctx)
{
    my_ethmsb::add_const_inplace<P2>(in, offset);
    TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, *ctx.iksk);
    const Torus half = out_value / 2;
    TFHEpp::Polynomial<P2> tv;
    tv.fill(Torus{0} - half);
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out, in0, *ctx.bkfft, tv);
    my_ethmsb::add_const_inplace<P2>(out, half);
}

void ethmsb(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& ct, const int k,
            const int kappa, const Torus out_value, const Context& ctx)
{
    if (k <= kappa) {
        pbs(out, ct, my_ethmsb::base_offset_for_current_layer(k), out_value,
            ctx);
        return;
    }
    TFHEpp::TLWE<P2> shifted;
    my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);
    const int suffix_bits = k - kappa;
    const Torus guard_value =
        my_ethmsb::guard_value_for_parent_scale(k, kappa);
    TFHEpp::TLWE<P2> guard;
    ethmsb(guard, shifted, suffix_bits, kappa, guard_value, ctx);
    TFHEpp::TLWE<P2> guarded;
    my_ethmsb::sub<P2>(guarded, ct, guard);
    pbs(out, guarded, my_ethmsb::gap_offset_for_current_layer(k, kappa),
        out_value, ctx);
}

bool semantic_ethmsb(const Torus m, const int k, const int kappa,
                     Torus& min_distance)
{
    if (k <= kappa) {
        const Torus after =
            encode(m, k) + my_ethmsb::base_offset_for_current_layer(k);
        min_distance = std::min(
            min_distance,
            my_ethmsb::torus_abs_centered(after - (Torus{1} << 63)));
        return (after & (Torus{1} << 63)) != 0;
    }
    const int suffix_bits = k - kappa;
    const Torus suffix_plain = m & ((Torus{1} << suffix_bits) - 1);
    const bool guard = semantic_ethmsb(suffix_plain, suffix_bits, kappa,
                                       min_distance);
    const Torus guard_weight = my_ethmsb::guard_weight(k, kappa);
    const Torus guarded_plain = m - (guard ? guard_weight : Torus{0});
    const Torus final_offset = my_ethmsb::gap_offset_for_current_layer(k, kappa);
    const Torus after = encode(guarded_plain, k) + final_offset;
    min_distance = std::min(
        min_distance,
        my_ethmsb::torus_abs_centered(after - (Torus{1} << 63)));
    return (after & (Torus{1} << 63)) != 0;
}

int exhaustive_semantic(const Options& opt)
{
    int failures = 0;
    Torus min_distance = std::numeric_limits<Torus>::max();
    const Torus switch_band = (Torus{1} << (64 - (P2::nbit + 1))) / 8;
    for (int k = 1; k <= opt.exhaustive_k; ++k) {
        for (Torus m = 0; m < (Torus{1} << k); ++m) {
            const bool actual = semantic_ethmsb(m, k, opt.kappa, min_distance);
            const bool expected = ((m >> (k - 1)) & 1) != 0;
            if (actual != expected) {
                std::cout << "H3COMPAT_EXHAUSTIVE_SEMANTIC_FAIL k=" << k
                          << " m=" << m << " expected=" << expected
                          << " actual=" << actual << "\n";
                return 1;
            }
        }
        std::cout << "H3COMPAT_EXHAUSTIVE_SEMANTIC_DONE k=" << k
                  << " failures_so_far=" << failures << "\n";
    }
    std::cout << "H3COMPAT_EXHAUSTIVE_SEMANTIC_RESULT failures=" << failures
              << " min_distance=" << min_distance
              << " switch_band=" << switch_band << "\n";
    return failures;
}

int representative_pbs(const Options& opt)
{
    const Context ctx = make_context();
    int failures = 0;
    for (int k = 1; k <= opt.exhaustive_k; ++k) {
        std::vector<Torus> cases = {0,
                                    1,
                                    (Torus{1} << (k - 1)) - 1,
                                    Torus{1} << (k - 1),
                                    (Torus{1} << k) - 1};
        std::sort(cases.begin(), cases.end());
        cases.erase(std::unique(cases.begin(), cases.end()), cases.end());
        for (const Torus m : cases) {
            TFHEpp::TLWE<P2> in = {};
            in[P2::k * P2::n] = encode(m, k);
            TFHEpp::TLWE<P2> out;
            ethmsb(out, in, k, opt.kappa, my_ethmsb::BOOL_ONE, ctx);
            const Torus ph = phase(out, ctx);
            const bool actual = closer_to_one(ph, my_ethmsb::BOOL_ONE);
            const bool expected = ((m >> (k - 1)) & 1) != 0;
            if (actual != expected) {
                std::cout << "H3COMPAT_REPRESENTATIVE_PBS_FAIL k=" << k
                          << " m=" << m << " expected=" << expected
                          << " actual=" << actual << " phase=" << ph << "\n";
                return 1;
            }
        }
    }
    std::cout << "H3COMPAT_REPRESENTATIVE_PBS_RESULT failures=" << failures
              << "\n";
    return failures;
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        std::cout << "MY_ETHMSB_STATIC_TESTS impl=h3compat_l20_l02 kappa="
                  << opt.kappa << " exhaustive_k=" << opt.exhaustive_k
                  << "\n";
        int failures = exhaustive_semantic(opt);
        if (failures == 0) failures += representative_pbs(opt);
        std::cout << "MY_ETHMSB_STATIC_TEST_RESULT failures=" << failures
                  << "\n";
        return failures == 0 ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: small_k_exhaustive\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
