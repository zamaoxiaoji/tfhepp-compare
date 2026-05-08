#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
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
using Clock = std::chrono::steady_clock;
using P2 = my_h3_lvl2param;
using P0 = my_h3_lvl0param;
using KS20 = my_h3_lvl20param;
using BR02 = my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Options {
    std::string impl = "h3compat_l20_l02";
    std::vector<int> bits = {8, 16, 24, 32};
    std::vector<std::string> ops = {"lt", "gt", "eq"};
    int kappa = 5;
    int trials = 3;
    int warmup = 0;
    std::uint64_t seed = 0;
    std::string csv_path;
    std::string fixed_inputs = "random";
    bool include_keygen_time = false;
    std::string margin_cert_min_hex = "external_gate_not_supplied";
    std::string switch_band_width_hex;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct PairCT {
    Torus a0 = 0;
    Torus b0 = 0;
    TFHEpp::TLWE<P2> a;
    TFHEpp::TLWE<P2> b;
};

struct Stats {
    double mean = 0;
    double median = 0;
    double p95 = 0;
    double stddev = 0;
    double min = 0;
    double max = 0;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

Torus ring_step_l2()
{
    return Torus{1} << (64 - (P2::nbit + 1));
}

Torus default_switch_band_width()
{
    return ring_step_l2() / 8;
}

std::vector<std::string> parse_ops(const std::string& s)
{
    if (s == "all") return {"lt", "gt", "eq", "le", "ge", "neq"};
    std::vector<std::string> out;
    std::size_t pos = 0;
    while (pos < s.size()) {
        const std::size_t comma = s.find(',', pos);
        out.push_back(s.substr(pos, comma - pos));
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
}

Options parse_args(const int argc, char** argv)
{
    Options opt;
    opt.switch_band_width_hex = hex64(default_switch_band_width());
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--impl")
            opt.impl = need("--impl");
        else if (a == "--ops")
            opt.ops = parse_ops(need("--ops"));
        else if (a == "--bits") {
            const std::string b = need("--bits");
            if (b == "all")
                opt.bits = {8, 16, 24, 32};
            else
                opt.bits = {std::atoi(b.c_str())};
        }
        else if (a == "--all")
            opt.bits = {8, 16, 24, 32};
        else if (a == "--kappa")
            opt.kappa = std::atoi(need("--kappa"));
        else if (a == "--trials")
            opt.trials = std::atoi(need("--trials"));
        else if (a == "--warmup")
            opt.warmup = std::atoi(need("--warmup"));
        else if (a == "--seed")
            opt.seed = std::strtoull(need("--seed"), nullptr, 10);
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else if (a == "--fixed-inputs")
            opt.fixed_inputs = need("--fixed-inputs");
        else if (a == "--include-keygen-time") {
            const std::string v = need("--include-keygen-time");
            opt.include_keygen_time = (v == "yes");
        }
        else if (a == "--margin-cert-min-hex")
            opt.margin_cert_min_hex = need("--margin-cert-min-hex");
        else if (a == "--switch-band-width-hex")
            opt.switch_band_width_hex = need("--switch-band-width-hex");
        else if (a == "--compare-he3db") {
            // Accepted for compatibility with older scripts. This smoke target
            // intentionally benchmarks only the certified strict variant.
        }
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (opt.impl != "h3compat_l20_l02" &&
        opt.impl != "ethmsb_h3compat_strict")
        throw std::runtime_error(
            "benchmark smoke supports only h3compat_l20_l02/ethmsb_h3compat_strict");
    if (opt.trials <= 0) throw std::runtime_error("--trials must be positive");
    if (opt.warmup < 0) throw std::runtime_error("--warmup must be nonnegative");
    if (opt.fixed_inputs != "random" && opt.fixed_inputs != "targeted")
        throw std::runtime_error("--fixed-inputs must be random or targeted");
    return opt;
}

Stats summarize(std::vector<double> xs)
{
    std::sort(xs.begin(), xs.end());
    Stats s;
    s.min = xs.front();
    s.max = xs.back();
    s.median = xs[xs.size() / 2];
    s.p95 = xs[std::min(xs.size() - 1, (xs.size() * 95) / 100)];
    s.mean = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    double var = 0;
    for (const double x : xs) var += (x - s.mean) * (x - s.mean);
    s.stddev = std::sqrt(var / xs.size());
    return s;
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

Context make_context(const std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
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

TFHEpp::TLWE<P2> encrypt_compare_value(const Torus x, const int bits,
                                       const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode(x, bits + 1), ctx.key2);
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

void eval_op(TFHEpp::TLWE<P2>& out, const TFHEpp::TLWE<P2>& a,
             const TFHEpp::TLWE<P2>& b, const int bits, const int kappa,
             const std::string& op, const Context& ctx)
{
    if (op == "lt" || op == "gt" || op == "le" || op == "ge") {
        TFHEpp::TLWE<P2> diff;
        if (op == "lt" || op == "ge")
            my_ethmsb::sub<P2>(diff, a, b);
        else
            my_ethmsb::sub<P2>(diff, b, a);
        TFHEpp::TLWE<P2> strict_ct;
        ethmsb(strict_ct, diff, bits + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
        if (op == "le" || op == "ge")
            my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(
                                        my_ethmsb::BOOL_ONE),
                                strict_ct);
        else
            out = strict_ct;
        return;
    }
    if (op == "eq" || op == "neq") {
        TFHEpp::TLWE<P2> lt_diff;
        TFHEpp::TLWE<P2> gt_diff;
        TFHEpp::TLWE<P2> lt_ct;
        TFHEpp::TLWE<P2> gt_ct;
        TFHEpp::TLWE<P2> neq_ct;
        my_ethmsb::sub<P2>(lt_diff, a, b);
        my_ethmsb::sub<P2>(gt_diff, b, a);
        ethmsb(lt_ct, lt_diff, bits + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
        ethmsb(gt_ct, gt_diff, bits + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
        my_ethmsb::add<P2>(neq_ct, lt_ct, gt_ct);
        if (op == "neq")
            out = neq_ct;
        else
            my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(
                                        my_ethmsb::BOOL_ONE),
                                neq_ct);
        return;
    }
    throw std::runtime_error("unknown op: " + op);
}

bool decrypt_bool(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return closer_to_one(phase(ct, ctx), my_ethmsb::BOOL_ONE);
}

bool expected_op(const Torus a, const Torus b, const std::string& op)
{
    if (op == "lt") return a < b;
    if (op == "gt") return a > b;
    if (op == "eq") return a == b;
    if (op == "le") return a <= b;
    if (op == "ge") return a >= b;
    if (op == "neq") return a != b;
    throw std::runtime_error("unknown op: " + op);
}

std::vector<std::pair<Torus, Torus>> targeted_pairs(const int bits)
{
    const Torus maxv = (Torus{1} << bits) - 1;
    const Torus mid = Torus{1} << (bits - 1);
    return {{0, 0},
            {0, 1},
            {1, 0},
            {0, maxv},
            {maxv, 0},
            {mid - 1, mid},
            {mid, mid - 1},
            {maxv, maxv}};
}

bool targeted_correctness(const int bits, const int kappa, const Context& ctx)
{
    for (const auto& [a0, b0] : targeted_pairs(bits)) {
        const auto a = encrypt_compare_value(a0, bits, ctx);
        const auto b = encrypt_compare_value(b0, bits, ctx);
        for (const std::string op : {"lt", "gt", "eq", "le", "ge", "neq"}) {
            TFHEpp::TLWE<P2> out;
            eval_op(out, a, b, bits, kappa, op, ctx);
            if (decrypt_bool(out, ctx) != expected_op(a0, b0, op))
                return false;
        }
    }
    return true;
}

std::vector<PairCT> make_pairs(const int bits, const int count,
                               const std::string& fixed_inputs,
                               std::mt19937_64& rng, const Context& ctx)
{
    const Torus maxv = (Torus{1} << bits) - 1;
    std::uniform_int_distribution<Torus> dist(0, maxv);
    const auto targets = targeted_pairs(bits);
    std::vector<PairCT> pairs;
    pairs.reserve(count);
    for (int i = 0; i < count; ++i) {
        PairCT p;
        if (fixed_inputs == "targeted") {
            p.a0 = targets[i % targets.size()].first;
            p.b0 = targets[i % targets.size()].second;
        }
        else {
            p.a0 = dist(rng);
            p.b0 = dist(rng);
        }
        p.a = encrypt_compare_value(p.a0, bits, ctx);
        p.b = encrypt_compare_value(p.b0, bits, ctx);
        pairs.push_back(p);
    }
    return pairs;
}

std::string build_type()
{
#ifdef NDEBUG
    return "Release";
#else
    return "Debug";
#endif
}

std::string compiler_string()
{
#ifdef __VERSION__
    return __VERSION__;
#else
    return "unknown";
#endif
}

std::string cpu_model()
{
    std::ifstream in("/proc/cpuinfo");
    std::string line;
    while (std::getline(in, line)) {
        const std::string key = "model name";
        if (line.rfind(key, 0) == 0) {
            const std::size_t colon = line.find(':');
            if (colon != std::string::npos) {
                std::string v = line.substr(colon + 1);
                while (!v.empty() && v.front() == ' ') v.erase(v.begin());
                std::replace(v.begin(), v.end(), ',', ';');
                return v;
            }
        }
    }
    return "unknown";
}

std::string env_or_unknown(const char* name)
{
    const char* v = std::getenv(name);
    return v == nullptr ? "unknown" : v;
}

void emit_row(std::ostream& out, const int bits, const int kappa,
              const int pbs_count, const std::string& op, const int trials,
              const int warmup, const Stats& stats, const int failures,
              const double keygen_ms, const std::string& margin_cert_min_hex,
              const std::string& switch_band_width_hex)
{
    out << "ethmsb_h3compat_strict," << bits << ',' << (bits + 1) << ','
        << kappa << ',' << pbs_count << ',' << op << ',' << trials << ','
        << warmup << ',' << stats.mean << ',' << stats.median << ','
        << stats.p95 << ',' << stats.stddev << ',' << stats.min << ','
        << stats.max << ',' << failures << ',' << keygen_ms
        << ",yes,yes," << build_type() << ',' << compiler_string() << ','
        << cpu_model() << ',' << env_or_unknown("ETHMSB_TASKSET_CORE") << ','
#ifdef USE_CONCRETE
        << "ON"
#else
        << "OFF"
#endif
        << ','
#ifdef USE_FFTW3
        << "ON"
#else
        << "OFF"
#endif
        << ','
#ifdef USE_MKL
        << "ON"
#else
        << "OFF"
#endif
        << ",h3compat_l20_l02," << margin_cert_min_hex << ','
        << switch_band_width_hex << '\n';
}

void print_need_info(const int bits)
{
    std::cerr << "BEGIN_NEED_INFO\n";
    std::cerr << "stage: bench_gate\n";
    std::cerr << "what_failed: targeted correctness failed before benchmark smoke\n";
    std::cerr << "commands_run:\n  ./build-128/my_ethmsb_bench_he3db_style --trials 3 --bits all --impl h3compat_l20_l02\n";
    std::cerr << "first_failure:\n";
    std::cerr << "  impl=ethmsb_h3compat_strict\n";
    std::cerr << "  stage=comparison\n";
    std::cerr << "  bits=" << bits << "\n";
    std::cerr << "params:\n";
    std::cerr << "  P2_T_bits=64\n  P2_n=" << P2::n
              << "\n  P2_alpha=" << P2::α << "\n";
    std::cerr << "  P0_T_bits=32\n  P0_n=" << P0::n
              << "\n  P0_alpha=" << P0::α << "\n";
    std::cerr << "  lvl20_t=" << KS20::t
              << "\n  lvl20_basebit=" << KS20::basebit << "\n";
    std::cerr << "  lvl02_domain=my_h3_lvl0param\n";
    std::cerr << "  lvl02_target=my_h3_lvl2param\n";
    std::cerr << "vendor_diff:\n";
    std::cerr << "  check with: git diff -- include/params/128bit.hpp include/params/concrete.hpp\n";
    std::cerr << "END_NEED_INFO\n";
}

}  // namespace

int main(int argc, char** argv)
{
    const Options opt = parse_args(argc, argv);

    const auto keygen_begin = Clock::now();
    const Context ctx = make_context(opt.seed);
    const auto keygen_end = Clock::now();
    const double keygen_ms =
        std::chrono::duration<double, std::milli>(keygen_end - keygen_begin)
            .count();

    std::ofstream csv_file;
    std::ostream* out = &std::cout;
    if (!opt.csv_path.empty()) {
        csv_file.open(opt.csv_path);
        out = &csv_file;
    }

    *out << "# BENCHMARK_SMOKE_ONLY=0\n";
    *out << "# impl=ethmsb_h3compat_strict\n";
    *out << "# keygen_ms=" << keygen_ms << "\n";
    *out << "# margin_cert_min_hex=" << opt.margin_cert_min_hex << "\n";
    *out << "# switch_band_width_hex=" << opt.switch_band_width_hex << "\n";
    *out << "# P0_T_bits=32 P0_n=" << P0::n << " P0_alpha=" << P0::α
         << "\n";
    *out << "# P2_T_bits=64 P2_n=" << P2::n << " P2_alpha=" << P2::α
         << "\n";
    *out << "# lvl20_t=" << KS20::t
         << " lvl20_basebit=" << KS20::basebit << "\n";
    *out << "# lvl02_domain=my_h3_lvl0param lvl02_target=my_h3_lvl2param\n";
    *out << "# BOOL_ONE=" << hex64(my_ethmsb::BOOL_ONE)
         << " kappa=" << opt.kappa << " delta33="
         << hex64(my_ethmsb::delta(33)) << "\n";
#ifdef USE_CONCRETE
    *out << "# USE_CONCRETE=ON\n";
#else
    *out << "# USE_CONCRETE=OFF\n";
#endif
    *out << "impl,bits,k,kappa,pbs_count,op,trials,warmup,mean_ms,median_ms,p95_ms,std_ms,min_ms,max_ms,failures,keygen_ms,encrypt_excluded,decrypt_excluded,build_type,compiler,cpu_model,taskset_core,USE_CONCRETE,USE_FFTW3,USE_MKL,params_profile,min_safe_margin_hex,switch_band_width_hex\n";

    std::mt19937_64 rng(opt.seed ^ 0x9e3779b97f4a7c15ULL);
    for (const int bits : opt.bits) {
        if (!targeted_correctness(bits, opt.kappa, ctx)) {
            print_need_info(bits);
            return 2;
        }

        const auto pairs = make_pairs(bits, opt.trials + opt.warmup,
                                      opt.fixed_inputs, rng, ctx);
        const int primitive_count =
            my_ethmsb::ethmsb_pbs_count(bits + 1, opt.kappa);
        for (const std::string& op : opt.ops) {
            std::vector<double> times;
            times.reserve(opt.trials);
            int failures = 0;
            for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
                const auto& p = pairs[i];
                TFHEpp::TLWE<P2> out_ct;
                const auto begin = Clock::now();
                eval_op(out_ct, p.a, p.b, bits, opt.kappa, op, ctx);
                const auto end = Clock::now();
                if (i >= opt.warmup) {
                    times.push_back(
                        std::chrono::duration<double, std::milli>(end - begin)
                            .count());
                    failures += decrypt_bool(out_ct, ctx) !=
                                expected_op(p.a0, p.b0, op);
                }
            }
            const int op_pbs = (op == "eq" || op == "neq")
                                   ? 2 * primitive_count
                                   : primitive_count;
            emit_row(*out, bits, opt.kappa, op_pbs, op, opt.trials,
                     opt.warmup, summarize(times), failures, keygen_ms,
                     opt.margin_cert_min_hex, opt.switch_band_width_hex);
        }
    }
    return 0;
}
