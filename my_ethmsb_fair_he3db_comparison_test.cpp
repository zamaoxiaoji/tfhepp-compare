#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "evalkeygens.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"

namespace {

using namespace my_ethmsb_params;
using Clock = std::chrono::system_clock;
using P0 = my_h3_lvl0param;
using P2 = my_h3_lvl2param;
using KS20 = my_h3_lvl20param;
using BR02 = my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

constexpr Torus BOOL_ONE_L2 = Torus{1} << 61;

struct Options {
    std::vector<int> bits = {16, 32};
    int trials = 100;
    int kappa = 5;
    std::uint64_t seed = 0;
    std::string range_mode = "he3db-code";
    std::string csv_path;
};

struct Context {
    TFHEpp::Key<P0> key0;
    TFHEpp::Key<P2> key2;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk20;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft02;
};

struct OpStats {
    std::uint64_t total_ms = 0;
    int errors = 0;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

std::vector<int> parse_bits(const std::string& s)
{
    if (s == "all") return {8, 16, 24, 32};
    if (s == "he3db-lvl2") return {16, 32};
    if (s == "he3db") return {4, 8, 16, 32};
    std::vector<int> out;
    std::size_t pos = 0;
    while (pos < s.size()) {
        const std::size_t comma = s.find(',', pos);
        out.push_back(std::atoi(s.substr(pos, comma - pos).c_str()));
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
}

Options parse_args(const int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--bits")
            opt.bits = parse_bits(need("--bits"));
        else if (a == "--trials")
            opt.trials = std::atoi(need("--trials"));
        else if (a == "--kappa")
            opt.kappa = std::atoi(need("--kappa"));
        else if (a == "--seed")
            opt.seed = std::strtoull(need("--seed"), nullptr, 10);
        else if (a == "--range")
            opt.range_mode = need("--range");
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (opt.trials <= 0) throw std::runtime_error("--trials must be positive");
    if (opt.kappa <= 0) throw std::runtime_error("--kappa must be positive");
    if (opt.range_mode != "he3db-code" && opt.range_mode != "comment" &&
        opt.range_mode != "full")
        throw std::runtime_error("--range must be he3db-code, comment, or full");
    for (const int b : opt.bits) {
        if (b < 1 || b > 32)
            throw std::runtime_error("--bits values must be in [1,32]");
    }
    return opt;
}

template <class P>
TFHEpp::TLWE<P> tlwe_int_encrypt(const typename P::T p,
                                 const std::uint32_t scale_bits,
                                 const TFHEpp::Key<P>& key)
{
    const double scale = std::pow(2.0, static_cast<double>(scale_bits));
    return TFHEpp::tlweSymEncrypt<P>(
        static_cast<typename P::T>(static_cast<double>(p) * scale), key);
}

bool decrypt_arith_l2(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    const P2::T phase = TFHEpp::tlweSymPhase<P2>(ct, ctx.key2);
    return my_ethmsb::torus_abs_centered(phase - BOOL_ONE_L2) <
           my_ethmsb::torus_abs_centered(phase);
}

Torus delta(const int k) { return my_ethmsb::delta(k); }

void strict_pbs(TFHEpp::TLWE<P2>& out, TFHEpp::TLWE<P2> in,
                const Torus offset, const Torus out_value, const Context& ctx)
{
    my_ethmsb::add_const_inplace<P2>(in, offset);
    TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, *ctx.iksk20);
    const Torus half = out_value / 2;
    TFHEpp::Polynomial<P2> tv;
    tv.fill(Torus{0} - half);
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out, in0, *ctx.bkfft02, tv);
    my_ethmsb::add_const_inplace<P2>(out, half);
}

void strict_ethmsb(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& in,
                   const int k, const int kappa, const Torus out_value,
                   const Context& ctx)
{
    if (k <= kappa) {
        strict_pbs(out, in, my_ethmsb::base_offset_for_current_layer(k),
                   out_value, ctx);
        return;
    }
    TFHEpp::TLWE<P2> shifted;
    my_ethmsb::scalar_mul_pow2<P2>(shifted, in, kappa);
    const int suffix_bits = k - kappa;
    const Torus guard_value =
        my_ethmsb::guard_value_for_parent_scale(k, kappa);
    TFHEpp::TLWE<P2> guard;
    strict_ethmsb(guard, shifted, suffix_bits, kappa, guard_value, ctx);
    TFHEpp::TLWE<P2> guarded;
    my_ethmsb::sub<P2>(guarded, in, guard);
    strict_pbs(out, guarded, my_ethmsb::gap_offset_for_current_layer(k, kappa),
               out_value, ctx);
}

void strict_lt(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int bits, const int kappa,
               const Context& ctx)
{
    TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, a, b);
    strict_ethmsb(out, diff, bits + 1, kappa, BOOL_ONE_L2, ctx);
}

void strict_gt(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int bits, const int kappa,
               const Context& ctx)
{
    TFHEpp::TLWE<P2> diff;
    my_ethmsb::sub<P2>(diff, b, a);
    strict_ethmsb(out, diff, bits + 1, kappa, BOOL_ONE_L2, ctx);
}

void strict_ge(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int bits, const int kappa,
               const Context& ctx)
{
    TFHEpp::TLWE<P2> is_lt;
    strict_lt(is_lt, a, b, bits, kappa, ctx);
    my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(BOOL_ONE_L2),
                       is_lt);
}

void strict_le(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int bits, const int kappa,
               const Context& ctx)
{
    TFHEpp::TLWE<P2> is_gt;
    strict_gt(is_gt, a, b, bits, kappa, ctx);
    my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(BOOL_ONE_L2),
                       is_gt);
}

void strict_eq(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
               const TFHEpp::TLWE<P2>& b, const int bits, const int kappa,
               const Context& ctx)
{
    TFHEpp::TLWE<P2> is_lt;
    TFHEpp::TLWE<P2> is_gt;
    TFHEpp::TLWE<P2> is_neq;
    strict_lt(is_lt, a, b, bits, kappa, ctx);
    strict_gt(is_gt, a, b, bits, kappa, ctx);
    my_ethmsb::add<P2>(is_neq, is_lt, is_gt);
    my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(BOOL_ONE_L2),
                       is_neq);
}

Context make_context(const std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    const auto keys = make_h3compat_keys(rng);
    Context ctx;
    ctx.key0 = keys.key0;
    ctx.key2 = keys.key2;
    ctx.iksk20 =
        std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft02 =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk20, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft02, ctx.key0, ctx.key2);
    return ctx;
}

Torus message_max(const int bits, const std::string& range_mode)
{
    if (range_mode == "full") return (Torus{1} << bits) - 1;
    if (range_mode == "comment")
        return bits == 1 ? 0 : ((Torus{1} << (bits - 1)) - 1);
    // Mirrors the effective C++ parse of HE3DB comparison_test.cpp:
    // (1 << (plain_bits - 1) - 1) == 1 << (plain_bits - 2).
    return bits <= 1 ? 0 : (Torus{1} << (bits - 2));
}

void csv_header(std::ostream& out)
{
    out << "impl,bits,k,op,trials,error_count,total_ms,mean_ms,"
           "timing_style,input_level,result_level,params_profile,range_mode,"
           "scale_bits,kappa\n";
}

void csv_row(std::ostream& out, const int bits, const std::string& op,
             const OpStats& s, const int trials,
             const std::string& range_mode, const int kappa)
{
    out << "ethmsb_h3compat_strict," << bits << ',' << (bits + 1) << ','
        << op << ',' << trials << ',' << s.errors << ',' << s.total_ms << ','
        << (static_cast<double>(s.total_ms) / trials)
        << ",he3db_comparison_test_system_clock_ms,lvl2,lvl2_arithmetic,"
           "my_h3compat_equal_to_HE3DB_128bit_params,"
        << range_mode << ','
        << (std::numeric_limits<P2::T>::digits - bits - 1) << ',' << kappa
        << '\n';
}

void run_bits(const int bits, const Options& opt, const Context& ctx,
              std::ostream* csv)
{
    static const std::vector<std::string> op_names = {"gt", "ge", "lt", "le",
                                                      "eq"};
    std::vector<OpStats> stats(5);
    const std::uint32_t scale_bits =
        std::numeric_limits<P2::T>::digits - bits - 1;
    const Torus max_msg = message_max(bits, opt.range_mode);
    std::mt19937_64 rng(opt.seed + static_cast<std::uint64_t>(bits) * 0x9e37U);
    std::uniform_int_distribution<Torus> dist(0, max_msg);

    std::cout << "------ ETHMSB strict with HE3DB comparison_test timing card ------\n";
    std::cout << "Test Time : " << opt.trials << "\n";
    std::cout << "Plain bits : " << bits << "\n";
    std::cout << "Input level : lvl2\n";
    std::cout << "Range mode : " << opt.range_mode << " max=" << max_msg
              << "\n";
    std::cout << "Timing card : system_clock now; duration_cast<milliseconds>; "
                 "timed region is exactly one comparison call\n";

    for (int test = 0; test < opt.trials; ++test) {
        const Torus p0 = dist(rng);
        const Torus p1 = dist(rng);
        const bool exp[5] = {p0 > p1, p0 >= p1, p0 < p1, p0 <= p1,
                             p0 == p1};
        const auto c0 =
            tlwe_int_encrypt<P2>(static_cast<P2::T>(p0), scale_bits, ctx.key2);
        const auto c1 =
            tlwe_int_encrypt<P2>(static_cast<P2::T>(p1), scale_bits, ctx.key2);

        TFHEpp::TLWE<P2> out;
        Clock::time_point start;
        Clock::time_point end;

        start = Clock::now();
        strict_gt(out, c0, c1, bits, opt.kappa, ctx);
        end = Clock::now();
        stats[0].total_ms +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        if (decrypt_arith_l2(out, ctx) != exp[0]) ++stats[0].errors;

        start = Clock::now();
        strict_ge(out, c0, c1, bits, opt.kappa, ctx);
        end = Clock::now();
        stats[1].total_ms +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        if (decrypt_arith_l2(out, ctx) != exp[1]) ++stats[1].errors;

        start = Clock::now();
        strict_lt(out, c0, c1, bits, opt.kappa, ctx);
        end = Clock::now();
        stats[2].total_ms +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        if (decrypt_arith_l2(out, ctx) != exp[2]) ++stats[2].errors;

        start = Clock::now();
        strict_le(out, c0, c1, bits, opt.kappa, ctx);
        end = Clock::now();
        stats[3].total_ms +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        if (decrypt_arith_l2(out, ctx) != exp[3]) ++stats[3].errors;

        start = Clock::now();
        strict_eq(out, c0, c1, bits, opt.kappa, ctx);
        end = Clock::now();
        stats[4].total_ms +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        if (decrypt_arith_l2(out, ctx) != exp[4]) ++stats[4].errors;
    }

    std::cout << "Test greater than, greater than or equal, less than, "
                 "less than or equal and equal to\n";
    for (std::size_t i = 0; i < op_names.size(); ++i) {
        std::cout << "impl=ethmsb_h3compat_strict op=" << op_names[i]
                  << " Error time" << i << " : " << stats[i].errors << "\n";
        std::cout << "impl=ethmsb_h3compat_strict op=" << op_names[i]
                  << " Time" << i << " : "
                  << static_cast<double>(stats[i].total_ms) / opt.trials
                  << "ms\n";
        if (csv)
            csv_row(*csv, bits, op_names[i], stats[i], opt.trials,
                    opt.range_mode, opt.kappa);
    }
}

void print_params(const Options& opt)
{
    std::cout << "PARAMS_PROFILE my_h3compat_equal_to_HE3DB_128bit_params\n";
    std::cout << "P0_T_bits=" << std::numeric_limits<P0::T>::digits
              << " P0_n=" << P0::n << " P0_alpha=" << P0::alpha << "\n";
    std::cout << "P2_T_bits=" << std::numeric_limits<P2::T>::digits
              << " P2_n=" << P2::n << " P2_alpha=" << P2::alpha << "\n";
    std::cout << "KS20_t=" << KS20::t << " KS20_basebit=" << KS20::basebit
              << " BR02=P0_to_P2 BOOL_ONE=" << hex64(BOOL_ONE_L2)
              << " kappa=" << opt.kappa << "\n";
    std::cout << "HE3DB_SOURCE_REFERENCE=../HE3DB/test/comparison_test.cpp\n";
    std::cout << "NOTE fair same-level comparison is HE3DB lvl2 path: bits 16 and 32.\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse_args(argc, argv);
        print_params(opt);
        auto ctx = make_context(opt.seed);
        std::ofstream csv_file;
        std::ostream* csv = nullptr;
        if (!opt.csv_path.empty()) {
            csv_file.open(opt.csv_path);
            if (!csv_file)
                throw std::runtime_error("failed to open csv: " +
                                         opt.csv_path);
            csv = &csv_file;
            csv_header(*csv);
        }
        for (const int bits : opt.bits) run_bits(bits, opt, ctx, csv);
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage=he3db_timing_card\n";
        std::cerr << "what_failed=" << e.what() << "\n";
        std::cerr << "commands_run:\n"
                  << "  ./build-ethmsb-release/"
                     "my_ethmsb_fair_he3db_comparison_test --bits he3db-lvl2 "
                     "--trials 100 --seed 0\n";
        std::cerr << "params:\n";
        std::cerr << "  P2_T_bits=" << std::numeric_limits<P2::T>::digits
                  << "\n";
        std::cerr << "  P2_n=" << P2::n << "\n";
        std::cerr << "  P2_alpha=" << P2::alpha << "\n";
        std::cerr << "  P0_T_bits=" << std::numeric_limits<P0::T>::digits
                  << "\n";
        std::cerr << "  P0_n=" << P0::n << "\n";
        std::cerr << "  P0_alpha=" << P0::alpha << "\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
