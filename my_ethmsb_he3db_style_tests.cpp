#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "my_ethmsb_he3db_style.hpp"

namespace {

using namespace my_ethmsb_he3db_style;

struct Options {
    int kappa = 5;
    int random_trials = 100;
    int max_trivial_k = 16;
};

struct Keys {
    TFHEpp::SecretKey sk;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk20;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft02;
};

Torus phase_l2(const TFHEpp::TLWE<P2>& ct, const TFHEpp::SecretKey& sk)
{
    return TFHEpp::tlweSymPhase<P2>(ct, sk.key.get<P2>());
}

bool closer_to_one(const Torus phase, const Torus out_value = BOOL_ONE_L2)
{
    return my_ethmsb::torus_abs_centered(phase - out_value) <
           my_ethmsb::torus_abs_centered(phase);
}

Options parse_args(const int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--kappa" && i + 1 < argc) opt.kappa = std::atoi(argv[++i]);
        else if (a == "--random" && i + 1 < argc)
            opt.random_trials = std::atoi(argv[++i]);
        else if (a == "--max-trivial-k" && i + 1 < argc)
            opt.max_trivial_k = std::atoi(argv[++i]);
    }
    return opt;
}

Keys make_keys()
{
    Keys keys;
    keys.iksk20 =
        std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    keys.bkfft02 =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*keys.iksk20, keys.sk);
    TFHEpp::bkfftgen<BR02>(*keys.bkfft02, keys.sk);
    return keys;
}

void print_need_info(const std::string& stage, const std::string& what,
                     const std::string& detail)
{
    std::cout << "BEGIN_NEED_INFO\n";
    std::cout << "stage: " << stage << "\n";
    std::cout << "what_failed: " << what << "\n";
    std::cout << "commands_run:\n";
    std::cout << "  cmake --build build-128 --target my_ethmsb_he3db_style_tests\n";
    std::cout << "  ./build-128/my_ethmsb_he3db_style_tests\n";
    std::cout << "exact_error_or_failure:\n" << detail << "\n";
    std::cout << "tfhepp_version:\n";
    std::cout << "  git_describe: v10\n";
    std::cout << "  git_rev_parse: d3653ef5608c11c60ddb2b754503f82d1d10f78f\n";
    std::cout << "he3db_paths_checked:\n";
    std::cout << "  ../HE3DB/src/HEDB/comparison/extract_msb.*\n";
    std::cout << "  ../HE3DB/src/HEDB/comparison/HomCompare.h\n";
    std::cout << "  ../HE3DB/src/HEDB/comparison/tfhepp_utils.*\n";
    std::cout << "  ../HE3DB/src/HEDB/comparison/operators.*\n";
    std::cout << "  ../HE3DB/src/HEDB/utils/types.h\n";
    std::cout << "params:\n";
    std::cout << "  P2_T_bits: " << std::numeric_limits<P2::T>::digits
              << "\n";
    std::cout << "  P2_n: " << P2::n << "\n";
    std::cout << "  P2_k: " << P2::k << "\n";
    std::cout << "  P2_l: " << P2::l << "\n";
    std::cout << "  P2_Bgbit: " << P2::Bgbit << "\n";
    std::cout << "  P2_alpha: " << P2::α << "\n";
    std::cout << "  P0_T_bits: " << std::numeric_limits<P0::T>::digits
              << "\n";
    std::cout << "  P0_n: " << P0::n << "\n";
    std::cout << "  P0_alpha: " << P0::α << "\n";
    std::cout << "api_signatures_needed:\n";
    std::cout << "  include/gatebootstrapping.hpp: GateBootstrappingTLWE2TLWEFFT / BRModSwitch / BlindRotate\n";
    std::cout << "  include/keyswitch.hpp: IdentityKeySwitch\n";
    std::cout << "  include/evalkeygens.hpp: bkfftgen / ikskgen\n";
    std::cout << "END_NEED_INFO\n";
}

int test_orientation(const Keys& keys, const int kappa)
{
    (void)kappa;
    for (const int input_bits : {3, 4, 5, 6}) {
        for (const Torus out_value : {BOOL_ONE_L2, Torus{1} << 58}) {
            const Torus lim = Torus{1} << input_bits;
            for (Torus m = 0; m < lim; ++m) {
                for (const bool encrypted : {false, true}) {
                    TFHEpp::TLWE<P2> in =
                        encrypted
                            ? TFHEpp::tlweSymEncrypt<P2>(
                                  encode_unsigned(m, input_bits),
                                  keys.sk.key.get<P2>())
                            : trivial_constant_l2(
                                  encode_unsigned(m, input_bits));
                    alignas(64) TFHEpp::TLWE<P2> out;
                    pbs_msb_value_l2_to_l2(
                        out, in, input_bits, delta(input_bits) / 2, out_value,
                        *keys.iksk20, *keys.bkfft02);
                    const Torus ph = phase_l2(out, keys.sk);
                    const bool actual = closer_to_one(ph, out_value);
                    const bool expected = ((m >> (input_bits - 1)) & 1) != 0;
                    if (actual != expected) {
                        std::cout << "HE3DB_STYLE_ORIENTATION_FAIL input_bits="
                                  << input_bits << " out_value=" << out_value
                                  << " encrypted=" << encrypted << " m=" << m
                                  << " expected=" << expected
                                  << " actual=" << actual
                                  << " raw_phase=" << ph << "\n";
                        return 1;
                    }
                }
            }
        }
    }
    std::cout << "HE3DB_STYLE_ORIENTATION_DONE failures=0\n";
    return 0;
}

int test_ethmsb(const Keys& keys, const int kappa, const int max_trivial_k)
{
    for (int k = 1; k <= max_trivial_k; ++k) {
        const Torus lim = Torus{1} << k;
        for (Torus m = 0; m < lim; ++m) {
            const auto in = trivial_constant_l2(encode_unsigned(m, k));
            alignas(64) TFHEpp::TLWE<P2> out;
            ethmsb_value_l2_to_l2(out, in, k, kappa, BOOL_ONE_L2,
                                  *keys.iksk20, *keys.bkfft02);
            const Torus ph = phase_l2(out, keys.sk);
            const bool actual = closer_to_one(ph);
            const bool expected = ((m >> (k - 1)) & 1) != 0;
            if (actual != expected) {
                std::cout << "HE3DB_STYLE_ETHMSB_TRIVIAL_FAIL k=" << k
                          << " kappa=" << kappa << " m=" << m
                          << " expected=" << expected
                          << " actual=" << actual << " raw_phase=" << ph
                          << "\n";
                return 1;
            }
        }
        std::cout << "HE3DB_STYLE_ETHMSB_TRIVIAL_DONE k=" << k
                  << " failures_so_far=0\n";
    }

    for (const int k : {9, 17, 25, 33}) {
        const std::vector<Torus> cases = {0,
                                          1,
                                          (Torus{1} << (k - 1)) - 1,
                                          Torus{1} << (k - 1),
                                          (Torus{1} << (k - 1)) + 1,
                                          (Torus{1} << k) - 1};
        for (const Torus m : cases) {
            const auto in = TFHEpp::tlweSymEncrypt<P2>(encode_unsigned(m, k),
                                                       keys.sk.key.get<P2>());
            alignas(64) TFHEpp::TLWE<P2> out;
            ethmsb_value_l2_to_l2(out, in, k, kappa, BOOL_ONE_L2,
                                  *keys.iksk20, *keys.bkfft02);
            const Torus ph = phase_l2(out, keys.sk);
            const bool actual = closer_to_one(ph);
            const bool expected = ((m >> (k - 1)) & 1) != 0;
            if (actual != expected) {
                std::cout << "HE3DB_STYLE_ETHMSB_ENC_FAIL k=" << k
                          << " kappa=" << kappa << " m=" << m
                          << " expected=" << expected
                          << " actual=" << actual << " raw_phase=" << ph
                          << "\n";
                return 1;
            }
        }
    }
    for (const auto& [k, m] :
         std::array<std::pair<int, Torus>, 2>{{{15, 16}, {15, 512}}}) {
        const auto in = TFHEpp::tlweSymEncrypt<P2>(encode_unsigned(m, k),
                                                   keys.sk.key.get<P2>());
        alignas(64) TFHEpp::TLWE<P2> out;
        ethmsb_value_l2_to_l2(out, in, k, kappa, BOOL_ONE_L2, *keys.iksk20,
                              *keys.bkfft02);
        const Torus ph = phase_l2(out, keys.sk);
        const bool actual = closer_to_one(ph);
        const bool expected = ((m >> (k - 1)) & 1) != 0;
        if (actual != expected) {
            std::cout << "HE3DB_STYLE_COUNTEREXAMPLE_FAIL k=" << k
                      << " kappa=" << kappa << " m=" << m
                      << " expected=" << expected << " actual=" << actual
                      << " raw_phase=" << ph << "\n";
            return 1;
        }
    }
    std::cout << "HE3DB_STYLE_ETHMSB_DONE kappa=" << kappa
              << " failures=0\n";
    return 0;
}

int test_comparison(const Keys& keys, const int kappa, const int random_trials)
{
    std::mt19937_64 rng(12345);
    for (const int t : {8, 16, 24, 32}) {
        const Torus maxv = (Torus{1} << t) - 1;
        const Torus mid = Torus{1} << (t - 1);
        std::vector<std::pair<Torus, Torus>> pairs = {
            {0, 0},       {0, 1},        {1, 0},       {0, maxv},
            {maxv, 0},    {mid - 1, mid}, {mid, mid - 1}, {maxv, maxv}};
        std::uniform_int_distribution<Torus> dist(0, maxv);
        for (int i = 0; i < random_trials; ++i)
            pairs.emplace_back(dist(rng), dist(rng));
        for (const auto& [a0, b0] : pairs) {
            const auto a = encrypt_unsigned_for_compare_l2(a0, t, keys.sk);
            const auto b = encrypt_unsigned_for_compare_l2(b0, t, keys.sk);
            alignas(64) TFHEpp::TLWE<P2> lt_ct;
            alignas(64) TFHEpp::TLWE<P2> gt_ct;
            alignas(64) TFHEpp::TLWE<P2> eq_ct;
            lt(lt_ct, a, b, t, kappa, *keys.iksk20, *keys.bkfft02);
            gt(gt_ct, a, b, t, kappa, *keys.iksk20, *keys.bkfft02);
            eq(eq_ct, a, b, t, kappa, *keys.iksk20, *keys.bkfft02);
            const bool actual_lt = decrypt_arith_bit_l2(lt_ct, keys.sk);
            const bool actual_gt = decrypt_arith_bit_l2(gt_ct, keys.sk);
            const bool actual_eq = decrypt_arith_bit_l2(eq_ct, keys.sk);
            if (actual_lt != (a0 < b0) || actual_gt != (a0 > b0) ||
                actual_eq != (a0 == b0)) {
                std::cout << "HE3DB_STYLE_COMPARISON_FAIL t=" << t
                          << " k=" << (t + 1) << " kappa=" << kappa
                          << " a=" << a0 << " b=" << b0
                          << " expected_lt=" << (a0 < b0)
                          << " actual_lt=" << actual_lt
                          << " expected_gt=" << (a0 > b0)
                          << " actual_gt=" << actual_gt
                          << " expected_eq=" << (a0 == b0)
                          << " actual_eq=" << actual_eq
                          << " raw_phase_lt=" << phase_l2(lt_ct, keys.sk)
                          << " raw_phase_gt=" << phase_l2(gt_ct, keys.sk)
                          << " raw_phase_eq=" << phase_l2(eq_ct, keys.sk)
                          << "\n";
                return 1;
            }
        }
        std::cout << "HE3DB_STYLE_COMPARISON_DONE t=" << t
                  << " random_trials=" << random_trials << " failures=0\n";
    }
    return 0;
}

int run_suite(const Keys& keys, const Options& opt, const int kappa,
              std::string& failed_stage)
{
    int failures = 0;
    failures += test_orientation(keys, kappa);
    if (failures != 0) failed_stage = "pbs_orientation";
    if (failures == 0)
        failures += test_ethmsb(keys, kappa, opt.max_trivial_k);
    if (failures != 0 && failed_stage.empty())
        failed_stage = "ethmsb_correctness";
    if (failures == 0)
        failures += test_comparison(keys, kappa, opt.random_trials);
    if (failures != 0 && failed_stage.empty()) failed_stage = "comparison";
    std::cout << "HE3DB_STYLE_SUITE_RESULT kappa=" << kappa
              << " failures=" << failures << "\n";
    return failures;
}

}  // namespace

int main(int argc, char** argv)
{
    const Options opt = parse_args(argc, argv);
    std::cout << "MY_ETHMSB_HE3DB_STYLE_TESTS kappa=" << opt.kappa
              << " random=" << opt.random_trials
              << " max_trivial_k=" << opt.max_trivial_k << "\n";
    auto keys = make_keys();
    std::string failed_stage;
    const int primary_failures = run_suite(keys, opt, opt.kappa, failed_stage);
    if (primary_failures != 0 && opt.kappa == 5) {
        std::cout << "HE3DB_STYLE_KAPPA5_FAILED_TRYING_KAPPA4\n";
        Options fallback = opt;
        fallback.kappa = 4;
        std::string fallback_stage;
        (void)run_suite(keys, fallback, 4, fallback_stage);
    }
    if (primary_failures != 0)
        print_need_info(failed_stage.empty() ? "comparison" : failed_stage,
                        "HE3DB-style correctness failed",
                        "See first HE3DB_STYLE_*_FAIL line above.");
    std::cout << "MY_ETHMSB_HE3DB_STYLE_TEST_RESULT failures="
              << primary_failures << "\n";
    return primary_failures == 0 ? 0 : 1;
}
