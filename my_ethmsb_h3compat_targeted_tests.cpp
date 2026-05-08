#include <cstdint>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
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

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct Options {
    int kappa = 5;
    int random = 0;
    std::uint64_t seed = 0;
    bool adversarial = false;
    std::string csv_path;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

Torus delta(const int k) { return my_ethmsb::delta(k); }

Torus encode(const Torus m, const int k)
{
    return static_cast<Torus>(Wide{m} * Wide{delta(k)});
}

Torus phase(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return TFHEpp::tlweSymPhase<P2>(ct, ctx.key2);
}

bool closer_to_one(const Torus ph, const Torus out_value)
{
    return my_ethmsb::torus_abs_centered(ph - out_value) <
           my_ethmsb::torus_abs_centered(ph);
}

TFHEpp::TLWE<P2> encrypt_m(const Torus m, const int k, const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(static_cast<P2::T>(encode(m, k)),
                                      ctx.key2);
}

void pbs(TFHEpp::TLWE<P2>& out, TFHEpp::TLWE<P2> in, const int input_bits,
         const Torus offset, const Torus out_value, const Context& ctx)
{
    (void)input_bits;
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
        pbs(out, ct, k, delta(k) / 2, out_value, ctx);
        return;
    }
    TFHEpp::TLWE<P2> shifted;
    my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);
    const int suffix_bits = k - kappa;
    const int w = k - kappa - 1;
    const Torus guard_weight = Torus{1} << w;
    const Torus guard_value =
        static_cast<Torus>(Wide{delta(k)} * Wide{guard_weight});
    TFHEpp::TLWE<P2> guard;
    ethmsb(guard, shifted, suffix_bits, kappa, guard_value, ctx);
    TFHEpp::TLWE<P2> guarded;
    my_ethmsb::sub<P2>(guarded, ct, guard);
    const Torus final_offset =
        static_cast<Torus>((Wide{(Torus{1} << w) + Torus{1}} *
                            Wide{delta(k)}) /
                           Wide{2});
    pbs(out, guarded, k, final_offset, out_value, ctx);
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

int test_ethmsb_targeted(const Context& ctx, const int kappa)
{
    int failures = 0;
    for (const int k : {9, 17, 25, 33}) {
        const std::vector<Torus> cases = {
            0,
            1,
            (Torus{1} << (k - 1)) - 1,
            Torus{1} << (k - 1),
            (Torus{1} << (k - 1)) + 1,
            (Torus{1} << k) - 1,
        };
        for (const Torus m : cases) {
            const auto ct = encrypt_m(m, k, ctx);
            TFHEpp::TLWE<P2> out;
            ethmsb(out, ct, k, kappa, my_ethmsb::BOOL_ONE, ctx);
            const Torus ph = phase(out, ctx);
            const bool actual = closer_to_one(ph, my_ethmsb::BOOL_ONE);
            const bool expected = ((m >> (k - 1)) & 1) != 0;
            std::cout << "H3COMPAT_ETHMSB k=" << k << " kappa=" << kappa
                      << " m=" << m << " expected=" << expected
                      << " actual=" << actual << " phase=" << hex64(ph)
                      << " pass=" << (expected == actual) << "\n";
            if (expected != actual) {
                ++failures;
                return failures;
            }
        }
    }
    return failures;
}

TFHEpp::TLWE<P2> encrypt_compare_value(const Torus x, const int t,
                                       const Context& ctx)
{
    return encrypt_m(x, t + 1, ctx);
}

int test_compare_targeted(const Context& ctx, const int kappa)
{
    int failures = 0;
    for (const int t : {8, 16, 24, 32}) {
        const Torus maxv = (Torus{1} << t) - 1;
        const Torus mid = Torus{1} << (t - 1);
        std::set<std::pair<Torus, Torus>> pairset = {
            {0, 0}, {0, 1}, {1, 0}, {0, maxv}, {maxv, 0},
            {mid - 1, mid}, {mid, mid - 1}, {maxv, maxv},
        };
        const std::vector<Torus> xs = {0,
                                       1,
                                       2,
                                       3,
                                       mid - 2,
                                       mid - 1,
                                       mid,
                                       mid + 1,
                                       maxv - 3,
                                       maxv - 2,
                                       maxv - 1,
                                       maxv};
        for (const Torus x : xs) {
            if (x > maxv) continue;
            pairset.insert({x, x});
            if (x + 1 <= maxv) {
                pairset.insert({x, x + 1});
                pairset.insert({x + 1, x});
            }
            pairset.insert({0, x});
            pairset.insert({x, 0});
            pairset.insert({maxv, x});
            pairset.insert({x, maxv});
        }
        for (const auto& [a0, b0] : pairset) {
            const auto a = encrypt_compare_value(a0, t, ctx);
            const auto b = encrypt_compare_value(b0, t, ctx);
            TFHEpp::TLWE<P2> diff_lt;
            TFHEpp::TLWE<P2> diff_gt;
            my_ethmsb::sub<P2>(diff_lt, a, b);
            my_ethmsb::sub<P2>(diff_gt, b, a);
            TFHEpp::TLWE<P2> lt_ct;
            TFHEpp::TLWE<P2> gt_ct;
            ethmsb(lt_ct, diff_lt, t + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
            ethmsb(gt_ct, diff_gt, t + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
            const Torus lt_phase = phase(lt_ct, ctx);
            const Torus gt_phase = phase(gt_ct, ctx);
            const bool lt = closer_to_one(lt_phase, my_ethmsb::BOOL_ONE);
            const bool gt = closer_to_one(gt_phase, my_ethmsb::BOOL_ONE);
            const bool le = !gt;
            const bool ge = !lt;
            const bool eq = !lt && !gt;
            const bool neq = lt || gt;
            const bool exp_lt = a0 < b0;
            const bool exp_gt = a0 > b0;
            const bool exp_le = a0 <= b0;
            const bool exp_ge = a0 >= b0;
            const bool exp_eq = a0 == b0;
            const bool exp_neq = a0 != b0;
            std::cout << "H3COMPAT_CMP t=" << t << " a=" << a0 << " b=" << b0
                      << " lt=" << lt << " gt=" << gt << " le=" << le
                      << " ge=" << ge << " eq=" << eq << " neq=" << neq
                      << " expected_lt=" << exp_lt
                      << " expected_gt=" << exp_gt
                      << " expected_le=" << exp_le
                      << " expected_ge=" << exp_ge
                      << " expected_eq=" << exp_eq
                      << " expected_neq=" << exp_neq
                      << " pass="
                      << (lt == exp_lt && gt == exp_gt && le == exp_le &&
                          ge == exp_ge && eq == exp_eq && neq == exp_neq)
                      << "\n";
            if (lt != exp_lt || gt != exp_gt || le != exp_le ||
                ge != exp_ge || eq != exp_eq || neq != exp_neq) {
                std::cout << "H3COMPAT_CMP_PHASES t=" << t << " a=" << a0
                          << " b=" << b0 << " lt_phase=" << hex64(lt_phase)
                          << " gt_phase=" << hex64(gt_phase) << "\n";
                ++failures;
                return failures;
            }
        }
    }
    return failures;
}

int test_compare_random(const Context& ctx, const int kappa, const int trials,
                        const std::uint64_t seed)
{
    if (trials <= 0) return 0;
    std::mt19937_64 rng(seed);
    for (const int t : {8, 16, 24, 32}) {
        const Torus maxv = (Torus{1} << t) - 1;
        std::uniform_int_distribution<Torus> dist(0, maxv);
        for (int i = 0; i < trials; ++i) {
            const Torus a0 = dist(rng);
            const Torus b0 = dist(rng);
            const auto a = encrypt_compare_value(a0, t, ctx);
            const auto b = encrypt_compare_value(b0, t, ctx);
            TFHEpp::TLWE<P2> diff_lt;
            TFHEpp::TLWE<P2> diff_gt;
            my_ethmsb::sub<P2>(diff_lt, a, b);
            my_ethmsb::sub<P2>(diff_gt, b, a);
            TFHEpp::TLWE<P2> lt_ct;
            TFHEpp::TLWE<P2> gt_ct;
            ethmsb(lt_ct, diff_lt, t + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
            ethmsb(gt_ct, diff_gt, t + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
            const bool lt =
                closer_to_one(phase(lt_ct, ctx), my_ethmsb::BOOL_ONE);
            const bool gt =
                closer_to_one(phase(gt_ct, ctx), my_ethmsb::BOOL_ONE);
            const bool eq = !lt && !gt;
            const bool exp_lt = a0 < b0;
            const bool exp_gt = a0 > b0;
            const bool exp_eq = a0 == b0;
            if (lt != exp_lt || gt != exp_gt || eq != exp_eq) {
                std::cout << "H3COMPAT_RANDOM_FAIL t=" << t << " i=" << i
                          << " a=" << a0 << " b=" << b0 << " lt=" << lt
                          << " gt=" << gt << " eq=" << eq
                          << " expected_lt=" << exp_lt
                          << " expected_gt=" << exp_gt
                          << " expected_eq=" << exp_eq << "\n";
                return 1;
            }
        }
        std::cout << "H3COMPAT_RANDOM_DONE t=" << t
                  << " trials=" << trials << " failures=0\n";
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--kappa" && i + 1 < argc)
            opt.kappa = std::stoi(argv[++i]);
        else if (a == "--random" && i + 1 < argc)
            opt.random = std::stoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc)
            opt.seed = std::stoull(argv[++i]);
        else if (a == "--adversarial")
            opt.adversarial = true;
        else if (a == "--csv" && i + 1 < argc)
            opt.csv_path = argv[++i];
    }
    std::ofstream csv_file;
    std::streambuf* old_cout = nullptr;
    if (!opt.csv_path.empty()) {
        csv_file.open(opt.csv_path);
        old_cout = std::cout.rdbuf(csv_file.rdbuf());
    }
    const auto ctx = make_context();
    int failures = 0;
    failures += test_ethmsb_targeted(ctx, opt.kappa);
    if (failures == 0) failures += test_compare_targeted(ctx, opt.kappa);
    if (failures == 0)
        failures += test_compare_random(ctx, opt.kappa, opt.random, opt.seed);
    if (failures != 0) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage=comparison\n";
        std::cerr << "what_failed=h3compat full ETHMSB/comparison targeted mismatch\n";
        std::cerr << "commands_run:\n";
        std::cerr << "  ./build-128/my_ethmsb_h3compat_targeted_tests --kappa "
                  << opt.kappa << "\n";
        std::cerr << "params:\n";
        std::cerr << "  h3compat_lvl0_T_bits="
                  << std::numeric_limits<P0::T>::digits << "\n";
        std::cerr << "  lvl2_T_bits=" << std::numeric_limits<P2::T>::digits
                  << "\n";
        std::cerr << "END_NEED_INFO\n";
    }
    std::cout << "MY_ETHMSB_H3COMPAT_TARGETED_RESULT failures=" << failures
              << "\n";
    if (old_cout != nullptr) std::cout.rdbuf(old_cout);
    return failures == 0 ? 0 : 1;
}
