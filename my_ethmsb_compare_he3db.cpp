#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "HEDB/comparison/HomCompare.h"
#include "my_ethmsb_he3db_style.hpp"

namespace {

using P2 = my_ethmsb_he3db_style::P2;
using P1 = TFHEpp::lvl1param;
using Torus = my_ethmsb_he3db_style::Torus;

struct Options {
    int kappa = 5;
    int random_trials = 100;
};

Options parse_args(const int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--kappa" && i + 1 < argc) opt.kappa = std::atoi(argv[++i]);
        else if (a == "--random" && i + 1 < argc)
            opt.random_trials = std::atoi(argv[++i]);
    }
    return opt;
}

bool decrypt_l2_bit(const TFHEpp::TLWE<P2>& ct, const TFHEpp::SecretKey& sk)
{
    return my_ethmsb_he3db_style::decrypt_arith_bit_l2(ct, sk);
}

bool decrypt_l1_logic(const TFHEpp::TLWE<P1>& ct, const TFHEpp::SecretKey& sk)
{
    return TFHEpp::tlweSymDecrypt<P1>(ct, sk.key.get<P1>());
}

void make_he3db_eval_key(TFHEpp::EvalKey& ek, const TFHEpp::SecretKey& sk)
{
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    ek.emplaceiksk<TFHEpp::lvl21param>(sk);
}

int compare_one(const int t, const Torus a0, const Torus b0, const int kappa,
                TFHEpp::SecretKey& sk, TFHEpp::EvalKey& ek)
{
    auto a = my_ethmsb_he3db_style::encrypt_unsigned_for_compare_l2(a0, t, sk);
    auto b = my_ethmsb_he3db_style::encrypt_unsigned_for_compare_l2(b0, t, sk);

    HEDB::TLWELvl1 he_lt;
    HEDB::TLWELvl1 he_gt;
    HEDB::TLWELvl1 he_eq;
    HEDB::less_than<P2>(a, b, he_lt, t, ek, LOGIC);
    HEDB::greater_than<P2>(a, b, he_gt, t, ek, LOGIC);
    HEDB::equal<P2>(a, b, he_eq, t, ek, LOGIC);

    alignas(64) TFHEpp::TLWE<P2> my_lt;
    alignas(64) TFHEpp::TLWE<P2> my_gt;
    alignas(64) TFHEpp::TLWE<P2> my_eq;
    my_ethmsb_he3db_style::lt(my_lt, a, b, t, kappa, *ek.iksklvl20,
                              *ek.bkfftlvl02);
    my_ethmsb_he3db_style::gt(my_gt, a, b, t, kappa, *ek.iksklvl20,
                              *ek.bkfftlvl02);
    my_ethmsb_he3db_style::eq(my_eq, a, b, t, kappa, *ek.iksklvl20,
                              *ek.bkfftlvl02);

    const bool he_lt_b = decrypt_l1_logic(he_lt, sk);
    const bool he_gt_b = decrypt_l1_logic(he_gt, sk);
    const bool he_eq_b = decrypt_l1_logic(he_eq, sk);
    const bool my_lt_b = decrypt_l2_bit(my_lt, sk);
    const bool my_gt_b = decrypt_l2_bit(my_gt, sk);
    const bool my_eq_b = decrypt_l2_bit(my_eq, sk);

    std::cout << "COMPARE_HE3DB t=" << t << " a=" << a0 << " b=" << b0
              << " HE3DB_lt=" << he_lt_b << " ETHMSB_lt=" << my_lt_b
              << " HE3DB_gt=" << he_gt_b << " ETHMSB_gt=" << my_gt_b
              << " HE3DB_eq=" << he_eq_b << " ETHMSB_eq=" << my_eq_b
              << "\n";

    const bool expected_lt = a0 < b0;
    const bool expected_gt = a0 > b0;
    const bool expected_eq = a0 == b0;
    if (he_lt_b != my_lt_b || he_gt_b != my_gt_b || he_eq_b != my_eq_b ||
        my_lt_b != expected_lt || my_gt_b != expected_gt ||
        my_eq_b != expected_eq) {
        std::cout << "COMPARE_HE3DB_FAIL t=" << t << " a=" << a0
                  << " b=" << b0 << " expected_lt=" << expected_lt
                  << " expected_gt=" << expected_gt
                  << " expected_eq=" << expected_eq << "\n";
        return 1;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    const Options opt = parse_args(argc, argv);
    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    make_he3db_eval_key(ek, sk);

    std::mt19937_64 rng(67890);
    int failures = 0;
    for (const int t : {8, 16, 24, 32}) {
        const Torus maxv = (Torus{1} << t) - 1;
        const Torus mid = Torus{1} << (t - 1);
        std::vector<std::pair<Torus, Torus>> pairs = {
            {0, 0},       {0, 1},        {1, 0},       {0, maxv},
            {maxv, 0},    {mid - 1, mid}, {mid, mid - 1}, {maxv, maxv}};
        std::uniform_int_distribution<Torus> dist(0, maxv);
        for (int i = 0; i < opt.random_trials; ++i)
            pairs.emplace_back(dist(rng), dist(rng));
        for (const auto& [a, b] : pairs) {
            failures += compare_one(t, a, b, opt.kappa, sk, ek);
            if (failures != 0) break;
        }
        if (failures != 0) break;
    }

    if (failures != 0) {
        std::cout << "BEGIN_NEED_INFO\n";
        std::cout << "stage: he3db_link\n";
        std::cout << "what_failed: HE3DB baseline and ETHMSB logical outputs differ or correctness failed\n";
        std::cout << "commands_run:\n";
        std::cout << "  cmake --build build-128 --target my_ethmsb_compare_he3db\n";
        std::cout << "  ./build-128/my_ethmsb_compare_he3db\n";
        std::cout << "exact_error_or_failure:\n";
        std::cout << "  See COMPARE_HE3DB_FAIL line above.\n";
        std::cout << "tfhepp_version:\n";
        std::cout << "  git_describe: v10\n";
        std::cout << "  git_rev_parse: d3653ef5608c11c60ddb2b754503f82d1d10f78f\n";
        std::cout << "he3db_paths_checked:\n";
        std::cout << "  ../HE3DB/src/HEDB/comparison/HomCompare.h\n";
        std::cout << "  ../HE3DB/src/HEDB/comparison/extract_msb.cpp\n";
        std::cout << "  ../HE3DB/src/HEDB/comparison/tfhepp_utils.cpp\n";
        std::cout << "  ../HE3DB/src/HEDB/comparison/operators.cpp\n";
        std::cout << "api_signatures_needed:\n";
        std::cout << "  include/keyswitch.hpp: IdentityKeySwitch\n";
        std::cout << "  include/gatebootstrapping.hpp: GateBootstrappingTLWE2TLWEFFT\n";
        std::cout << "END_NEED_INFO\n";
    }
    std::cout << "MY_ETHMSB_COMPARE_HE3DB_RESULT failures=" << failures
              << "\n";
    return failures == 0 ? 0 : 1;
}
