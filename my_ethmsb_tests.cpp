#include "my_ethmsb_fixed.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

using my_ethmsb::BOOL_ONE;
using my_ethmsb::P;
using my_ethmsb::Torus;
using my_ethmsb::my_lvl22pbsparam;

struct Options {
    std::vector<int> kappas{3, 4, 5};
    int random_tests = 1000;
    int probe_reps = 100;
    int max_exhaustive_k = 12;
    bool exhaustive16 = false;
    bool skip_exhaustive = false;
    bool skip_comparison = false;
    bool skip_random = false;
    bool skip_probe = false;
};

std::string capture_cmd(const char* cmd)
{
    std::array<char, 256> buffer{};
    std::string out;
    FILE* pipe = popen(cmd, "r");
    if (pipe == nullptr) return "unavailable";
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) !=
           nullptr)
        out += buffer.data();
    pclose(pipe);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r'))
        out.pop_back();
    return out.empty() ? "unavailable" : out;
}

std::string hex64(const Torus value)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << value
       << std::dec;
    return os.str();
}

std::string first_cpu_flags()
{
    std::ifstream in("/proc/cpuinfo");
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("flags", 0) == 0) return line;
    }
    return "unavailable";
}

void print_params_block()
{
    std::cout << "sizeof_T=" << sizeof(P::T) << '\n';
    std::cout << "P_n=" << P::n << '\n';
    std::cout << "P_k=" << P::k << '\n';
    std::cout << "P_l=" << P::l << '\n';
    std::cout << "P_l_a=" << P::lₐ << '\n';
    std::cout << "Bgbit=" << P::Bgbit << '\n';
    std::cout << "Bg_a_bit=" << P::Bgₐbit << '\n';
    std::cout << "alpha=" << std::setprecision(18) << P::α << '\n';
    std::cout << "key_value_min=" << P::key_value_min << '\n';
    std::cout << "key_value_max=" << P::key_value_max << '\n';
    std::cout << "mu=" << P::μ << '\n';
}

void print_environment()
{
    std::cout << "=== ETHMSB fixed test environment ===\n";
    std::cout << "git_describe=" << capture_cmd("git describe --tags --always")
              << '\n';
    std::cout << "git_rev_parse=" << capture_cmd("git rev-parse HEAD")
              << '\n';
    print_params_block();
    std::cout << "macro_USE_CONCRETE="
#ifdef USE_CONCRETE
              << "ON\n";
#else
              << "OFF\n";
#endif
    std::cout << "macro_USE_FFTW3="
#ifdef USE_FFTW3
              << "ON\n";
#else
              << "OFF\n";
#endif
    std::cout << "macro_USE_MKL="
#ifdef USE_MKL
              << "ON\n";
#else
              << "OFF\n";
#endif
    std::cout << "macro_USE_KEY_BUNDLE="
#ifdef USE_KEY_BUNDLE
              << "ON\n";
#else
              << "OFF\n";
#endif
    std::cout << "macro_USE_TERNARY_CMUX="
#ifdef USE_TERNARY_CMUX
              << "ON\n";
#else
              << "OFF\n";
#endif
    const std::string flags = first_cpu_flags();
    std::cout << "cpu_flags_first=" << flags << '\n';
    std::cout << "cpu_has_avx2="
              << (flags.find(" avx2 ") != std::string::npos ? "yes" : "no")
              << '\n';
    std::cout << "cpu_has_avx512f="
              << (flags.find(" avx512f ") != std::string::npos ? "yes"
                                                                 : "no")
              << '\n';
}

void print_api_signatures()
{
    std::cout << "api_signatures_needed:\n";
    std::cout << "  template<class P> void bkfftgen("
                 "BootstrappingKeyFFT<P>&, const Key<typename P::domainP>&, "
                 "const Key<typename P::targetP>&)\n";
    std::cout << "  template<class P> void GateBootstrappingTLWE2TLWEFFT("
                 "TLWE<typename P::targetP>&, const TLWE<typename "
                 "P::domainP>&, const BootstrappingKeyFFT<P>&, const "
                 "Polynomial<typename P::targetP>&)\n";
    std::cout << "  template<class P> using TLWE = std::array<typename P::T, "
                 "P::k * P::n + 1>\n";
    std::cout << "  template<class P> using BootstrappingKeyFFT = "
                 "std::array<BootstrappingKeyElementFFT<P>, "
                 "P::domainP::k * P::domainP::n / P::Addends>\n";
}

void print_need_info(const std::string& stage, const std::string& what_failed,
                     const std::string& commands_run,
                     const std::string& exact_error,
                     const std::string& test_case = "none")
{
    std::cout << "BEGIN_NEED_INFO\n";
    std::cout << "stage: " << stage << '\n';
    std::cout << "what_failed: " << what_failed << '\n';
    std::cout << "commands_run:\n" << commands_run;
    std::cout << "relevant_files:\n";
    std::cout << "  my_ethmsb_fixed.hpp\n";
    std::cout << "  my_ethmsb_tests.cpp\n";
    std::cout << "  include/evalkeygens.hpp\n";
    std::cout << "  include/gatebootstrapping.hpp\n";
    std::cout << "exact_error_or_failure:\n" << exact_error << '\n';
    std::cout << "tfhepp_version:\n";
    std::cout << "  git_describe: "
              << capture_cmd("git describe --tags --always") << '\n';
    std::cout << "  git_rev_parse: " << capture_cmd("git rev-parse HEAD")
              << '\n';
    std::cout << "params:\n";
    std::cout << "  sizeof_T: " << sizeof(P::T) << '\n';
    std::cout << "  P_n: " << P::n << '\n';
    std::cout << "  P_k: " << P::k << '\n';
    std::cout << "  P_l: " << P::l << '\n';
    std::cout << "  P_l_a: " << P::lₐ << '\n';
    std::cout << "  Bgbit: " << P::Bgbit << '\n';
    std::cout << "  alpha: " << std::setprecision(18) << P::α << '\n';
    std::cout << "  key_min_max: " << P::key_value_min << ','
              << P::key_value_max << '\n';
    std::cout << "test_case_if_any:\n" << test_case << '\n';
    print_api_signatures();
    std::cout << "END_NEED_INFO\n";
}

Options parse_args(const int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(name);
            return argv[++i];
        };
        if (arg == "--kappa") {
            opt.kappas = {std::stoi(need_value("--kappa needs 3/4/5"))};
        }
        else if (arg == "--random") {
            opt.random_tests = std::stoi(need_value("--random needs N"));
        }
        else if (arg == "--probe-reps") {
            opt.probe_reps = std::stoi(need_value("--probe-reps needs N"));
        }
        else if (arg == "--max-exhaustive-k") {
            opt.max_exhaustive_k =
                std::stoi(need_value("--max-exhaustive-k needs N"));
        }
        else if (arg == "--exhaustive16") {
            opt.exhaustive16 = true;
            opt.max_exhaustive_k = std::max(opt.max_exhaustive_k, 16);
        }
        else if (arg == "--skip-exhaustive") {
            opt.skip_exhaustive = true;
        }
        else if (arg == "--skip-random") {
            opt.skip_random = true;
        }
        else if (arg == "--skip-comparison") {
            opt.skip_comparison = true;
        }
        else if (arg == "--skip-probe") {
            opt.skip_probe = true;
        }
        else if (arg == "--quick") {
            opt.random_tests = 2;
            opt.probe_reps = 2;
            opt.max_exhaustive_k = 3;
        }
        else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    return opt;
}

int run_encode_decode_tests()
{
    int failures = 0;
    const std::vector<int> ks{1, 2, 3, 4, 5, 8, 9, 16, 17, 24, 25, 32, 33};
    for (const int k : ks) {
        const Torus top = Torus{1} << k;
        const std::vector<Torus> values{0,
                                        1 % top,
                                        (Torus{1} << (k - 1)) - 1,
                                        Torus{1} << (k - 1),
                                        top - 1};
        for (const Torus m : values) {
            const Torus phase = my_ethmsb::encode_unsigned(m, k);
            const Torus decoded = my_ethmsb::decode_unsigned_phase(phase, k);
            if (decoded != m) {
                ++failures;
                std::cout << "ENCODE_DECODE_FAIL k=" << k << " m=" << m
                          << " phase=" << hex64(phase)
                          << " decoded=" << decoded << '\n';
            }
        }
    }
    std::cout << "ENCODE_DECODE failures=" << failures << '\n';
    return failures;
}

bool expected_msb(const Torus m, const int k)
{
    return ((m >> (k - 1)) & Torus{1}) != 0;
}

Torus phase_of(const TFHEpp::TLWE<P>& ct, const TFHEpp::SecretKey& sk)
{
    return TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
}

int run_orientation_tests(const TFHEpp::SecretKey& sk,
                          const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>&
                              bkfft)
{
    int failures = 0;
    bool need_info_printed = false;
    for (const int kappa : {3, 4, 5}) {
        for (const int input_bits : {kappa, std::min(kappa + 1, 6)}) {
            const Torus maxv = (Torus{1} << input_bits) - 1;
            const std::vector<Torus> ms{0,
                                        (Torus{1} << (input_bits - 1)) - 1,
                                        Torus{1} << (input_bits - 1),
                                        maxv};
            for (const Torus m : ms) {
                const Torus phase = my_ethmsb::encode_unsigned(m, input_bits);
                auto ct = my_ethmsb::trivial_constant<P>(phase);
                alignas(64) TFHEpp::TLWE<P> out;
                my_ethmsb::pbs_msb_value<my_lvl22pbsparam>(
                    out, ct, input_bits, my_ethmsb::delta(input_bits) / 2,
                    BOOL_ONE, bkfft);
                const Torus raw = phase_of(out, sk);
                const bool actual =
                    my_ethmsb::decrypt_arith_bit(out, sk, BOOL_ONE);
                const bool expected = expected_msb(m, input_bits);
                std::cout << "ORIENTATION kappa=" << kappa
                          << " input_bits=" << input_bits << " m=" << m
                          << " raw_phase=" << hex64(raw)
                          << " expected=" << expected
                          << " actual=" << actual << '\n';
                if (actual != expected) {
                    ++failures;
                    if (!need_info_printed) {
                        std::ostringstream tc;
                        tc << "  k: " << input_bits << "\n";
                        tc << "  kappa: " << kappa << "\n";
                        tc << "  m_or_a_b: " << m << "\n";
                        tc << "  expected: " << expected << "\n";
                        tc << "  actual: " << actual << "\n";
                        tc << "  raw_phase: " << hex64(raw) << "\n";
                        print_need_info(
                            "pbs_orientation",
                            "pbs_msb_value orientation is not the requested "
                            "0/out_value threshold",
                            "  ./build/my_ethmsb_tests --quick\n",
                            "orientation mismatch", tc.str());
                        need_info_printed = true;
                    }
                }
            }
        }
    }
    std::cout << "ORIENTATION failures=" << failures << '\n';
    return failures;
}

int check_ethmsb_case(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft, const int k,
    const int kappa, const Torus m, bool& need_info_printed)
{
    auto ct =
        TFHEpp::tlweSymEncrypt<P>(my_ethmsb::encode_unsigned(m, k),
                                  sk.key.get<P>());
    alignas(64) TFHEpp::TLWE<P> out;
    my_ethmsb::ethmsb_value<my_lvl22pbsparam>(out, ct, k, kappa, BOOL_ONE,
                                              bkfft);
    const bool actual = my_ethmsb::decrypt_arith_bit(out, sk, BOOL_ONE);
    const bool expected = expected_msb(m, k);
    if (actual == expected) return 0;
    const Torus raw = phase_of(out, sk);
    std::cout << "ETHMSB_FAIL k=" << k << " kappa=" << kappa << " m=" << m
              << " expected=" << expected << " actual=" << actual
              << " raw_phase=" << hex64(raw)
              << " pbs_count=" << my_ethmsb::ethmsb_pbs_count(k, kappa)
              << '\n';
    if (!need_info_printed) {
        std::ostringstream tc;
        tc << "  k: " << k << "\n";
        tc << "  kappa: " << kappa << "\n";
        tc << "  m_or_a_b: " << m << "\n";
        tc << "  expected: " << expected << "\n";
        tc << "  actual: " << actual << "\n";
        tc << "  raw_phase: " << hex64(raw) << "\n";
        print_need_info("ethmsb_correctness",
                        "ethmsb_value produced a wrong arithmetic MSB bit",
                        "  ./build/my_ethmsb_tests\n", "ETHMSB mismatch",
                        tc.str());
        need_info_printed = true;
    }
    return 1;
}

std::vector<Torus> targeted_m_values(const int k, const int kappa)
{
    const Torus maxv = (Torus{1} << k) - 1;
    std::vector<Torus> ms{0,
                          1,
                          (Torus{1} << (k - 1)) - 1,
                          Torus{1} << (k - 1),
                          (Torus{1} << (k - 1)) + 1,
                          maxv};
    if (kappa + 1 < k) {
        const Torus lower = (Torus{1} << kappa) - 1;
        ms.push_back(lower);
        ms.push_back((Torus{1} << kappa) | lower);
    }
    if (k == 15 && kappa == 5) {
        ms.push_back(16);
        ms.push_back(512);
    }
    std::sort(ms.begin(), ms.end());
    ms.erase(std::unique(ms.begin(), ms.end()), ms.end());
    return ms;
}

int run_ethmsb_tests(const TFHEpp::SecretKey& sk,
                     const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft,
                     const Options& opt, const int kappa)
{
    int failures = 0;
    bool need_info_printed = false;
    if (!opt.skip_exhaustive) {
        for (int k = 1; k <= opt.max_exhaustive_k; ++k) {
            if (k > 12 && !opt.exhaustive16) continue;
            const Torus limit = Torus{1} << k;
            for (Torus m = 0; m < limit; ++m)
                failures +=
                    check_ethmsb_case(sk, bkfft, k, kappa, m,
                                      need_info_printed);
            std::cout << "ETHMSB_EXHAUSTIVE_DONE k=" << k
                      << " kappa=" << kappa << " failures_so_far="
                      << failures << '\n';
        }
    }

    for (const int k : {8, 9, 15, 16, 17, 24, 25, 32, 33}) {
        for (const Torus m : targeted_m_values(k, kappa)) {
            failures += check_ethmsb_case(sk, bkfft, k, kappa, m,
                                          need_info_printed);
        }
    }
    std::cout << "ETHMSB kappa=" << kappa << " failures=" << failures << '\n';
    return failures;
}

struct PairCheck {
    int t;
    Torus a;
    Torus b;
};

std::vector<PairCheck> targeted_pairs_for_bits(const int t)
{
    const Torus top = Torus{1} << t;
    return {{t, 0, 0},
            {t, 0, 1},
            {t, 1, 0},
            {t, 0, top - 1},
            {t, top - 1, 0},
            {t, (Torus{1} << (t - 1)) - 1, Torus{1} << (t - 1)},
            {t, Torus{1} << (t - 1), (Torus{1} << (t - 1)) - 1},
            {t, top - 1, top - 1}};
}

int check_pair_full(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft,
    const PairCheck& tc, const int kappa, bool& need_info_printed)
{
    const auto ca = my_ethmsb::encrypt_unsigned_for_compare(tc.a, tc.t, sk);
    const auto cb = my_ethmsb::encrypt_unsigned_for_compare(tc.b, tc.t, sk);
    alignas(64) TFHEpp::TLWE<P> lt_ct;
    alignas(64) TFHEpp::TLWE<P> gt_ct;
    alignas(64) TFHEpp::TLWE<P> le_ct;
    alignas(64) TFHEpp::TLWE<P> ge_ct;
    alignas(64) TFHEpp::TLWE<P> eq_ct;
    alignas(64) TFHEpp::TLWE<P> neq_ct;
    my_ethmsb::lt(lt_ct, ca, cb, tc.t, kappa, bkfft);
    my_ethmsb::gt(gt_ct, ca, cb, tc.t, kappa, bkfft);
    my_ethmsb::le(le_ct, ca, cb, tc.t, kappa, bkfft);
    my_ethmsb::ge(ge_ct, ca, cb, tc.t, kappa, bkfft);
    my_ethmsb::eq(eq_ct, ca, cb, tc.t, kappa, bkfft);
    my_ethmsb::neq(neq_ct, ca, cb, tc.t, kappa, bkfft);

    const bool got_lt = my_ethmsb::decrypt_arith_bit(lt_ct, sk);
    const bool got_gt = my_ethmsb::decrypt_arith_bit(gt_ct, sk);
    const bool got_le = my_ethmsb::decrypt_arith_bit(le_ct, sk);
    const bool got_ge = my_ethmsb::decrypt_arith_bit(ge_ct, sk);
    const bool got_eq = my_ethmsb::decrypt_arith_bit(eq_ct, sk);
    const bool got_neq = my_ethmsb::decrypt_arith_bit(neq_ct, sk);

    const bool exp_lt = tc.a < tc.b;
    const bool exp_gt = tc.a > tc.b;
    const bool exp_le = tc.a <= tc.b;
    const bool exp_ge = tc.a >= tc.b;
    const bool exp_eq = tc.a == tc.b;
    const bool exp_neq = tc.a != tc.b;
    const bool ok = got_lt == exp_lt && got_gt == exp_gt &&
                    got_le == exp_le && got_ge == exp_ge &&
                    got_eq == exp_eq && got_neq == exp_neq;
    if (ok) return 0;

    const Torus phase_lt = phase_of(lt_ct, sk);
    const Torus phase_gt = phase_of(gt_ct, sk);
    const Torus phase_eq = phase_of(eq_ct, sk);
    std::cout << "COMPARISON_FAIL t=" << tc.t << " k=" << (tc.t + 1)
              << " kappa=" << kappa << " a=" << tc.a << " b=" << tc.b
              << " expected_lt_gt_le_ge_eq_neq=" << exp_lt << exp_gt
              << exp_le << exp_ge << exp_eq << exp_neq
              << " actual_lt_gt_le_ge_eq_neq=" << got_lt << got_gt << got_le
              << got_ge << got_eq << got_neq
              << " phase_lt=" << hex64(phase_lt)
              << " phase_gt=" << hex64(phase_gt)
              << " phase_eq=" << hex64(phase_eq)
              << " pbs_count="
              << my_ethmsb::ethmsb_pbs_count(tc.t + 1, kappa) << '\n';
    if (!need_info_printed) {
        std::ostringstream info;
        info << "  t: " << tc.t << "\n";
        info << "  k: " << (tc.t + 1) << "\n";
        info << "  kappa: " << kappa << "\n";
        info << "  m_or_a_b: " << tc.a << "," << tc.b << "\n";
        info << "  expected: " << exp_lt << exp_gt << exp_le << exp_ge
             << exp_eq << exp_neq << "\n";
        info << "  actual: " << got_lt << got_gt << got_le << got_ge
             << got_eq << got_neq << "\n";
        info << "  raw_phase: lt=" << hex64(phase_lt)
             << " gt=" << hex64(phase_gt) << " eq=" << hex64(phase_eq)
             << "\n";
        print_need_info(
            "comparison",
            "unsigned comparison result did not match the clear comparison",
            "  ./build/my_ethmsb_tests\n", "comparison mismatch",
            info.str());
        need_info_printed = true;
    }
    return 1;
}

int check_pair_lt_gt(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft,
    const PairCheck& tc, const int kappa, bool& need_info_printed)
{
    const auto ca = my_ethmsb::encrypt_unsigned_for_compare(tc.a, tc.t, sk);
    const auto cb = my_ethmsb::encrypt_unsigned_for_compare(tc.b, tc.t, sk);
    alignas(64) TFHEpp::TLWE<P> lt_ct;
    alignas(64) TFHEpp::TLWE<P> gt_ct;
    my_ethmsb::lt(lt_ct, ca, cb, tc.t, kappa, bkfft);
    my_ethmsb::gt(gt_ct, ca, cb, tc.t, kappa, bkfft);
    const bool got_lt = my_ethmsb::decrypt_arith_bit(lt_ct, sk);
    const bool got_gt = my_ethmsb::decrypt_arith_bit(gt_ct, sk);
    const bool got_le = !got_gt;
    const bool got_ge = !got_lt;
    const bool got_neq = got_lt || got_gt;
    const bool got_eq = !got_neq;

    const bool exp_lt = tc.a < tc.b;
    const bool exp_gt = tc.a > tc.b;
    const bool exp_le = tc.a <= tc.b;
    const bool exp_ge = tc.a >= tc.b;
    const bool exp_eq = tc.a == tc.b;
    const bool exp_neq = tc.a != tc.b;
    if (got_lt == exp_lt && got_gt == exp_gt && got_le == exp_le &&
        got_ge == exp_ge && got_eq == exp_eq && got_neq == exp_neq)
        return 0;

    const Torus phase_lt = phase_of(lt_ct, sk);
    const Torus phase_gt = phase_of(gt_ct, sk);
    std::cout << "COMPARISON_RANDOM_FAIL t=" << tc.t << " k=" << (tc.t + 1)
              << " kappa=" << kappa << " a=" << tc.a << " b=" << tc.b
              << " expected_lt_gt_le_ge_eq_neq=" << exp_lt << exp_gt
              << exp_le << exp_ge << exp_eq << exp_neq
              << " actual_lt_gt_le_ge_eq_neq=" << got_lt << got_gt << got_le
              << got_ge << got_eq << got_neq
              << " phase_lt=" << hex64(phase_lt)
              << " phase_gt=" << hex64(phase_gt)
              << " phase_eq=derived"
              << " pbs_count="
              << my_ethmsb::ethmsb_pbs_count(tc.t + 1, kappa) << '\n';
    if (!need_info_printed) {
        print_need_info("comparison",
                        "random unsigned comparison result mismatch",
                        "  ./build/my_ethmsb_tests\n",
                        "random comparison mismatch");
        need_info_printed = true;
    }
    return 1;
}

int run_comparison_tests(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft,
    const Options& opt, const int kappa)
{
    int failures = 0;
    bool need_info_printed = false;
    for (const int t : {8, 16, 24, 32}) {
        for (const PairCheck& tc : targeted_pairs_for_bits(t)) {
            failures +=
                check_pair_full(sk, bkfft, tc, kappa, need_info_printed);
        }
        std::cout << "COMPARISON_TARGETED_DONE t=" << t
                  << " kappa=" << kappa << " failures_so_far=" << failures
                  << '\n';
    }

    if (!opt.skip_random) {
        std::mt19937_64 rng(0x4554484d5342554cULL + kappa);
        for (const int t : {8, 16, 24, 32}) {
            std::uniform_int_distribution<Torus> dist(0, (Torus{1} << t) - 1);
            for (int i = 0; i < opt.random_tests; ++i) {
                const PairCheck tc{t, dist(rng), dist(rng)};
                failures += check_pair_lt_gt(sk, bkfft, tc, kappa,
                                             need_info_printed);
            }
            std::cout << "COMPARISON_RANDOM_DONE t=" << t
                      << " kappa=" << kappa
                      << " trials=" << opt.random_tests
                      << " failures_so_far=" << failures << '\n';
        }
    }
    std::cout << "COMPARISON kappa=" << kappa << " failures=" << failures
              << '\n';
    return failures;
}

void print_distance_stats(const std::vector<Torus>& values, const int k,
                          const int kappa, const Torus m)
{
    std::vector<Torus> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    const auto pick = [&](double q) -> Torus {
        if (sorted.empty()) return 0;
        const size_t idx = std::min(sorted.size() - 1,
                                    static_cast<size_t>(q * sorted.size()));
        return sorted[idx];
    };
    std::cout << "NOISE_PROBE k=" << k << " kappa=" << kappa << " m=" << m
              << " min=" << sorted.front()
              << " median=" << pick(0.50)
              << " p95=" << pick(0.95)
              << " max=" << sorted.back() << '\n';
}

int run_noise_probing(
    const TFHEpp::SecretKey& sk,
    const TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>& bkfft,
    const Options& opt, const int kappa)
{
    int failures = 0;
    const int k = 33;
    const std::vector<Torus> ms{0,
                                1,
                                (Torus{1} << 32) - 1,
                                Torus{1} << 32,
                                (Torus{1} << 33) - 1};
    for (const Torus m : ms) {
        std::vector<Torus> distances;
        distances.reserve(opt.probe_reps);
        for (int i = 0; i < opt.probe_reps; ++i) {
            auto ct =
                TFHEpp::tlweSymEncrypt<P>(my_ethmsb::encode_unsigned(m, k),
                                          sk.key.get<P>());
            alignas(64) TFHEpp::TLWE<P> out;
            my_ethmsb::ethmsb_value<my_lvl22pbsparam>(out, ct, k, kappa,
                                                      BOOL_ONE, bkfft);
            const Torus phase = phase_of(out, sk);
            const bool expected = expected_msb(m, k);
            const Torus target = expected ? BOOL_ONE : 0;
            const bool actual = my_ethmsb::decrypt_arith_bit(out, sk);
            failures += actual == expected ? 0 : 1;
            distances.push_back(my_ethmsb::torus_abs_centered(phase - target));
        }
        print_distance_stats(distances, k, kappa, m);
    }
    std::cout << "NOISE_PROBE kappa=" << kappa << " failures=" << failures
              << '\n';
    return failures;
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse_args(argc, argv);
        print_environment();

        int failures = run_encode_decode_tests();
        if (failures != 0) return 1;

        const auto keygen_start = std::chrono::steady_clock::now();
        TFHEpp::SecretKey sk;
        auto bkfft = std::make_unique_for_overwrite<
            TFHEpp::BootstrappingKeyFFT<my_lvl22pbsparam>>();
        TFHEpp::bkfftgen<my_lvl22pbsparam>(*bkfft, sk.key.get<P>(),
                                           sk.key.get<P>());
        const auto keygen_end = std::chrono::steady_clock::now();
        std::cout << "KEYGEN_MS="
                  << std::chrono::duration<double, std::milli>(keygen_end -
                                                               keygen_start)
                         .count()
                  << '\n';

        failures += run_orientation_tests(sk, *bkfft);
        if (failures != 0) return 1;

        int kappa5_failures = 0;
        for (const int kappa : opt.kappas) {
            const int before = failures;
            failures += run_ethmsb_tests(sk, *bkfft, opt, kappa);
            if (!opt.skip_comparison)
                failures += run_comparison_tests(sk, *bkfft, opt, kappa);
            if (!opt.skip_probe)
                failures += run_noise_probing(sk, *bkfft, opt, kappa);
            if (kappa == 5) kappa5_failures = failures - before;
        }

        if (kappa5_failures != 0 &&
            std::find(opt.kappas.begin(), opt.kappas.end(), 5) !=
                opt.kappas.end()) {
            std::cout << "KAPPA5_FAILED_AUTO_FALLBACK_BEGIN\n";
            for (const int fallback : {4, 3}) {
                Options quick = opt;
                quick.kappas = {fallback};
                quick.random_tests = 0;
                quick.skip_random = true;
                quick.skip_probe = true;
                quick.skip_exhaustive = true;
                int fallback_failures =
                    run_ethmsb_tests(sk, *bkfft, quick, fallback) +
                    run_comparison_tests(sk, *bkfft, quick, fallback);
                std::cout << "KAPPA_FALLBACK kappa=" << fallback
                          << " failures=" << fallback_failures << '\n';
            }
            std::cout << "KAPPA5_FAILED_AUTO_FALLBACK_END\n";
        }

        std::cout << "MY_ETHMSB_TEST_RESULT failures=" << failures << '\n';
        return failures == 0 ? 0 : 1;
    }
    catch (const std::exception& ex) {
        print_need_info("build", "test harness threw an exception",
                        "  ./build/my_ethmsb_tests\n", ex.what());
        return 1;
    }
}
