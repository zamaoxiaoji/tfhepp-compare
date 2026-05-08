#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"

namespace {

using Clock = std::chrono::steady_clock;
using namespace my_ethmsb_params;
using P2 = my_h3_lvl2param;
using P0 = my_h3_lvl0param;
using KS20 = my_h3_lvl20param;
using BR02 = my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Options {
    std::string impl = "h3compat_l20_l02";
    std::vector<int> bits = {8, 16, 24, 32};
    int kappa = 5;
    std::vector<std::string> ops = {"lt", "gt"};
    int trials_per_bit = 20;
    std::uint64_t seed = 0;
    int max_seconds = 300;
    std::string checkpoint = "random_ethmsb.jsonl";
    bool resume = false;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
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
    return {std::stoi(s)};
}

std::vector<std::string> parse_ops(const std::string& s)
{
    if (s == "all") return {"lt", "gt", "eq"};
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

Options parse(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--impl")
            opt.impl = need("--impl");
        else if (a == "--bits")
            opt.bits = parse_bits(need("--bits"));
        else if (a == "--kappa")
            opt.kappa = std::stoi(need("--kappa"));
        else if (a == "--ops")
            opt.ops = parse_ops(need("--ops"));
        else if (a == "--trials-per-bit")
            opt.trials_per_bit = std::stoi(need("--trials-per-bit"));
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
        else if (a == "--max-seconds")
            opt.max_seconds = std::stoi(need("--max-seconds"));
        else if (a == "--checkpoint")
            opt.checkpoint = need("--checkpoint");
        else if (a == "--resume")
            opt.resume = true;
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    if (opt.impl != "h3compat_l20_l02")
        throw std::runtime_error("random_resumable currently supports h3compat_l20_l02 only");
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

Torus phase(const TFHEpp::TLWE<P2>& ct, const Context& ctx)
{
    return TFHEpp::tlweSymPhase<P2>(ct, ctx.key2);
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

TFHEpp::TLWE<P2> encrypt_value(const Torus x, const int t, const Context& ctx)
{
    return TFHEpp::tlweSymEncrypt<P2>(encode(x, t + 1), ctx.key2);
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

bool eval_op(const Context& ctx, const int bits, const int kappa,
             const std::string& op, const Torus a0, const Torus b0,
             Torus& out_phase)
{
    const auto a = encrypt_value(a0, bits, ctx);
    const auto b = encrypt_value(b0, bits, ctx);
    TFHEpp::TLWE<P2> diff;
    if (op == "lt" || op == "eq")
        my_ethmsb::sub<P2>(diff, a, b);
    else if (op == "gt")
        my_ethmsb::sub<P2>(diff, b, a);
    else
        throw std::runtime_error("unknown op");

    TFHEpp::TLWE<P2> out;
    if (op == "eq") {
        TFHEpp::TLWE<P2> lt_ct;
        TFHEpp::TLWE<P2> gt_ct;
        TFHEpp::TLWE<P2> diff_gt;
        ethmsb(lt_ct, diff, bits + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
        my_ethmsb::sub<P2>(diff_gt, b, a);
        ethmsb(gt_ct, diff_gt, bits + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
        TFHEpp::TLWE<P2> neq_ct;
        my_ethmsb::add<P2>(neq_ct, lt_ct, gt_ct);
        my_ethmsb::sub<P2>(out, my_ethmsb::trivial_constant<P2>(my_ethmsb::BOOL_ONE),
                           neq_ct);
    }
    else {
        ethmsb(out, diff, bits + 1, kappa, my_ethmsb::BOOL_ONE, ctx);
    }
    out_phase = phase(out, ctx);
    return closer_to_one(out_phase, my_ethmsb::BOOL_ONE);
}

std::set<std::string> load_done(const std::string& path)
{
    std::set<std::string> done;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        const auto find_num = [&](const std::string& key) -> std::string {
            const std::string pat = "\"" + key + "\":";
            const std::size_t p = line.find(pat);
            if (p == std::string::npos) return "";
            std::size_t b = p + pat.size();
            if (line[b] == '"') ++b;
            std::size_t e = b;
            while (e < line.size() && line[e] != ',' && line[e] != '"') ++e;
            return line.substr(b, e - b);
        };
        const std::string bits = find_num("bits");
        const std::string index = find_num("index");
        const std::string op = find_num("op");
        if (!bits.empty() && !index.empty() && !op.empty())
            done.insert(bits + ":" + op + ":" + index);
    }
    return done;
}

void print_failure(const int bits, const int index, const std::string& op,
                   const Torus a, const Torus b, const bool expected,
                   const bool actual, const Torus phase)
{
    std::cerr << "BEGIN_NEED_INFO\n";
    std::cerr << "stage: random_resumable\n";
    std::cerr << "what_failed: random comparison mismatch\n";
    std::cerr << "commands_run:\n  ./build-128/my_ethmsb_random_resumable\n";
    std::cerr << "first_failure:\n";
    std::cerr << "  impl=h3compat_l20_l02\n";
    std::cerr << "  stage=comparison\n";
    std::cerr << "  m_or_a_b=(" << a << "," << b << ")\n";
    std::cerr << "  root_k=" << (bits + 1) << "\n";
    std::cerr << "  kappa=5\n";
    std::cerr << "  expected=" << expected << "\n";
    std::cerr << "  actual=" << actual << "\n";
    std::cerr << "  output_phase=" << hex64(phase) << "\n";
    std::cerr << "  index=" << index << "\n";
    std::cerr << "  op=" << op << "\n";
    std::cerr << "END_NEED_INFO\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const auto begin_all = Clock::now();
        const auto ctx = make_context(opt.seed);
        const std::set<std::string> done =
            opt.resume ? load_done(opt.checkpoint) : std::set<std::string>{};
        std::ofstream out(opt.checkpoint, std::ios::app);
        int completed = 0;
        for (const int bits : opt.bits) {
            const Torus maxv = (Torus{1} << bits) - 1;
            for (const std::string& op : opt.ops) {
                for (int index = 0; index < opt.trials_per_bit; ++index) {
                    const std::string key =
                        std::to_string(bits) + ":" + op + ":" +
                        std::to_string(index);
                    if (done.contains(key)) continue;
                    const auto elapsed_s =
                        std::chrono::duration<double>(Clock::now() - begin_all)
                            .count();
                    if (elapsed_s > opt.max_seconds) {
                        std::cout << "RANDOM_RESUMABLE_TIMEOUT completed="
                                  << completed << " checkpoint="
                                  << opt.checkpoint << "\n";
                        return 0;
                    }
                    std::mt19937_64 rng(opt.seed ^
                                        (static_cast<std::uint64_t>(bits) << 32) ^
                                        (static_cast<std::uint64_t>(index) * 0x9e3779b97f4a7c15ULL) ^
                                        static_cast<std::uint64_t>(op[0]));
                    std::uniform_int_distribution<Torus> dist(0, maxv);
                    const Torus a = dist(rng);
                    const Torus b = dist(rng);
                    const bool expected =
                        op == "lt" ? a < b : op == "gt" ? a > b : a == b;
                    Torus out_phase = 0;
                    const auto t0 = Clock::now();
                    const bool actual =
                        eval_op(ctx, bits, opt.kappa, op, a, b, out_phase);
                    const auto t1 = Clock::now();
                    const double elapsed_ms =
                        std::chrono::duration<double, std::milli>(t1 - t0)
                            .count();
                    const bool pass = expected == actual;
                    out << "{\"seed\":" << opt.seed << ",\"index\":" << index
                        << ",\"bits\":" << bits << ",\"a\":" << a
                        << ",\"b\":" << b << ",\"op\":\"" << op
                        << "\",\"expected\":" << expected
                        << ",\"actual\":" << actual << ",\"pass\":" << pass
                        << ",\"elapsed_ms\":" << elapsed_ms << "}\n";
                    out.flush();
                    ++completed;
                    if (!pass) {
                        print_failure(bits, index, op, a, b, expected, actual,
                                      out_phase);
                        return 1;
                    }
                }
                std::cout << "RANDOM_RESUMABLE_DONE bits=" << bits
                          << " op=" << op
                          << " trials=" << opt.trials_per_bit << "\n";
            }
        }
        std::cout << "RANDOM_RESUMABLE_RESULT failures=0 completed="
                  << completed << " checkpoint=" << opt.checkpoint << "\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: random_resumable\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n  ./build-128/my_ethmsb_random_resumable\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
    return 0;
}
