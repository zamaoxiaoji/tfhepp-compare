#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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

using namespace my_ethmsb_params;
using P2 = my_h3_lvl2param;
using P0 = my_h3_lvl0param;
using KS20 = my_h3_lvl20param;
using BR02 = my_h3_lvl02param;
using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Options {
    int kappa = 5;
    int max_k = 33;
    std::string mode = "trivial";
    std::uint64_t seed = 0;
    int trials_per_case = 3;
    std::string csv_path;
    Torus switch_band_width = 0;
    Torus output_noise_margin = 0;
};

struct Context {
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;
};

struct CertRow {
    std::uint64_t cert_id = 0;
    int root_k = 0;
    int depth = 0;
    std::string stage;
    int current_k = 0;
    int kappa = 0;
    Torus out_value = 0;
    Torus offset = 0;
    Torus plaintext = 0;
    Torus phase_before = 0;
    Torus phase_after = 0;
    Torus distance_to_half = 0;
    Torus measured_switch_band = 0;
    Torus safe_margin = 0;
    bool expected = false;
    bool actual = false;
    Torus output_phase = 0;
    bool pass = false;
};

struct Summary {
    bool valid_failure = false;
    bool unsafe = false;
    Torus min_distance = std::numeric_limits<Torus>::max();
    Torus min_safe_margin = std::numeric_limits<Torus>::max();
    CertRow first_bad;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

Options parse(int argc, char** argv)
{
    Options opt;
    opt.switch_band_width = (Torus{1} << (64 - (P2::nbit + 1))) / 8;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--kappa")
            opt.kappa = std::stoi(need("--kappa"));
        else if (a == "--max-k")
            opt.max_k = std::stoi(need("--max-k"));
        else if (a == "--mode")
            opt.mode = need("--mode");
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
        else if (a == "--trials-per-case")
            opt.trials_per_case = std::stoi(need("--trials-per-case"));
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else if (a == "--switch-band-width-hex")
            opt.switch_band_width =
                static_cast<Torus>(std::strtoull(
                    need("--switch-band-width-hex"), nullptr, 0));
        else if (a == "--empirical-output-noise-margin-hex")
            opt.output_noise_margin =
                static_cast<Torus>(std::strtoull(
                    need("--empirical-output-noise-margin-hex"), nullptr, 0));
        else if (a == "--impl") {
            const std::string impl = need("--impl");
            if (impl != "h3compat_l20_l02")
                throw std::runtime_error("margin_cert currently supports h3compat_l20_l02 only");
        }
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    return opt;
}

Torus delta(const int k) { return my_ethmsb::delta(k); }

Torus mask_for(const int k) { return (Torus{1} << k) - 1; }

Torus encode(const Torus m, const int k)
{
    return static_cast<Torus>(Wide{m} * Wide{delta(k)});
}

Torus dist_to_half(const Torus phase_after_offset)
{
    return my_ethmsb::torus_abs_centered(phase_after_offset - (Torus{1} << 63));
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

TFHEpp::TLWE<P2> make_input(const Torus m, const int k, const Context& ctx,
                            const std::string& mode)
{
    if (mode == "trivial") {
        TFHEpp::TLWE<P2> ct = {};
        ct[P2::k * P2::n] = encode(m, k);
        return ct;
    }
    if (mode == "encrypted")
        return TFHEpp::tlweSymEncrypt<P2>(encode(m, k), ctx.key2);
    throw std::runtime_error("unknown --mode");
}

void emit_header(std::ostream& out)
{
    out << "cert_id,root_k,depth,stage,current_k,kappa,out_value_hex,offset_hex,plaintext_m_at_this_node,phase_before_offset_hex,phase_after_offset_hex,distance_to_half_hex,measured_switch_band_hex,safe_margin_hex,expected,actual,output_phase_hex,pass\n";
}

void emit_row(std::ostream& out, const CertRow& r)
{
    out << r.cert_id << ',' << r.root_k << ',' << r.depth << ',' << r.stage
        << ',' << r.current_k << ',' << r.kappa << ',' << hex64(r.out_value)
        << ',' << hex64(r.offset) << ',' << r.plaintext << ','
        << hex64(r.phase_before) << ',' << hex64(r.phase_after) << ','
        << hex64(r.distance_to_half) << ',' << hex64(r.measured_switch_band)
        << ',' << hex64(r.safe_margin) << ',' << r.expected << ','
        << r.actual << ',' << hex64(r.output_phase) << ',' << r.pass << '\n';
}

TFHEpp::TLWE<P2> pbs_trace(std::ostream& out, const Context& ctx,
                           TFHEpp::TLWE<P2> in, const int root_k,
                           const int depth, const std::string& stage,
                           const int current_k, const int kappa,
                           const Torus plaintext, const Torus offset,
                           const Torus out_value,
                           const Torus measured_switch_band,
                           const Torus safety_budget,
                           std::uint64_t& cert_id, Summary& summary)
{
    const Torus before = phase(in, ctx);
    my_ethmsb::add_const_inplace<P2>(in, offset);
    const Torus after = phase(in, ctx);
    TFHEpp::TLWE<P0> in0;
    TFHEpp::IdentityKeySwitch<KS20>(in0, in, *ctx.iksk);
    const Torus half = out_value / 2;
    TFHEpp::Polynomial<P2> tv;
    tv.fill(Torus{0} - half);
    TFHEpp::TLWE<P2> out_ct;
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out_ct, in0, *ctx.bkfft, tv);
    my_ethmsb::add_const_inplace<P2>(out_ct, half);
    const Torus out_phase = phase(out_ct, ctx);

    CertRow row;
    row.cert_id = cert_id++;
    row.root_k = root_k;
    row.depth = depth;
    row.stage = stage;
    row.current_k = current_k;
    row.kappa = kappa;
    row.out_value = out_value;
    row.offset = offset;
    row.plaintext = plaintext & mask_for(current_k);
    row.phase_before = before;
    row.phase_after = after;
    row.distance_to_half = dist_to_half(after);
    row.measured_switch_band = measured_switch_band;
    row.safe_margin = row.distance_to_half > safety_budget
                          ? row.distance_to_half - safety_budget
                          : 0;
    row.expected =
        ((encode(row.plaintext, current_k) + offset) & (Torus{1} << 63)) != 0;
    row.actual = closer_to_one(out_phase, out_value);
    row.output_phase = out_phase;
    row.pass = row.expected == row.actual;
    emit_row(out, row);

    summary.min_distance = std::min(summary.min_distance, row.distance_to_half);
    summary.min_safe_margin =
        std::min(summary.min_safe_margin, row.safe_margin);
    if ((!row.pass || row.safe_margin == 0) && !summary.valid_failure &&
        !summary.unsafe) {
        summary.first_bad = row;
    }
    if (!row.pass) summary.valid_failure = true;
    if (row.safe_margin == 0) summary.unsafe = true;
    return out_ct;
}

TFHEpp::TLWE<P2> ethmsb_trace(std::ostream& out, const Context& ctx,
                              const TFHEpp::TLWE<P2>& ct, const int root_k,
                              const int current_k, const int kappa,
                              const Torus plaintext, const Torus out_value,
                              const int depth, const std::string& stage,
                              const Torus measured_switch_band,
                              const Torus safety_budget,
                              std::uint64_t& cert_id, Summary& summary)
{
    if (current_k <= kappa) {
        return pbs_trace(out, ctx, ct, root_k, depth, stage, current_k, kappa,
                         plaintext, delta(current_k) / 2, out_value,
                         measured_switch_band, safety_budget, cert_id,
                         summary);
    }

    TFHEpp::TLWE<P2> shifted;
    my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);
    const int suffix_bits = current_k - kappa;
    const Torus suffix_plain = plaintext & mask_for(suffix_bits);
    const int w = current_k - kappa - 1;
    const Torus guard_weight = Torus{1} << w;
    const Torus guard_value =
        static_cast<Torus>(Wide{delta(current_k)} * Wide{guard_weight});
    TFHEpp::TLWE<P2> guard = ethmsb_trace(
        out, ctx, shifted, root_k, suffix_bits, kappa, suffix_plain,
        guard_value, depth + 1, "guard", measured_switch_band, safety_budget,
        cert_id, summary);

    TFHEpp::TLWE<P2> guarded;
    my_ethmsb::sub<P2>(guarded, ct, guard);
    const bool expected_guard =
        ((suffix_plain >> (suffix_bits - 1)) & Torus{1}) != 0;
    const Torus guarded_plain =
        plaintext - (expected_guard ? guard_weight : Torus{0});
    const Torus final_offset =
        static_cast<Torus>((Wide{(Torus{1} << w) + Torus{1}} *
                            Wide{delta(current_k)}) /
                           Wide{2});
    return pbs_trace(out, ctx, guarded, root_k, depth, "final", current_k,
                     kappa, guarded_plain, final_offset, out_value,
                     measured_switch_band, safety_budget, cert_id, summary);
}

std::vector<Torus> boundary_candidates(const int k)
{
    std::set<Torus> vals;
    vals.insert(0);
    vals.insert(1);
    vals.insert((Torus{1} << (k - 1)) - 1);
    vals.insert(Torus{1} << (k - 1));
    vals.insert((Torus{1} << (k - 1)) + 1);
    vals.insert((Torus{1} << k) - 1);
    return {vals.begin(), vals.end()};
}

std::vector<Torus> candidates_for_k(const int k, const int kappa)
{
    std::set<Torus> vals;
    for (const Torus v : boundary_candidates(k)) vals.insert(v);
    if (k > kappa) {
        const int w = k - kappa - 1;
        const Torus W = Torus{1} << w;
        const Torus half = Torus{1} << (k - 1);
        for (const Torus v : {half - 1 - W, half - W, half - 1, half,
                              half + W - 1, half + W, W - 1, W, W + 1}) {
            if (v < (Torus{1} << k)) vals.insert(v);
        }
    }
    for (int suffix_bits = k - kappa; suffix_bits > 0;
         suffix_bits -= kappa) {
        const int prefix_bits = k - suffix_bits;
        if (prefix_bits <= 0) continue;
        const int use_prefix_bits = std::min(kappa, prefix_bits);
        std::set<Torus> prefixes = {0, 1,
                                    (Torus{1} << use_prefix_bits) - 1};
        for (const Torus s : boundary_candidates(suffix_bits)) {
            for (const Torus prefix : prefixes) {
                const Torus m = (prefix << suffix_bits) | s;
                if (m < (Torus{1} << k)) vals.insert(m);
            }
        }
    }
    return {vals.begin(), vals.end()};
}

void print_need_info(const Summary& summary)
{
    const CertRow& r = summary.first_bad;
    std::cerr << "BEGIN_NEED_INFO\n";
    std::cerr << "stage: margin_cert\n";
    std::cerr << "what_failed: "
              << (summary.valid_failure ? "VALID reachable case failed"
                                        : "reachable case has non-positive safe margin")
              << "\n";
    std::cerr << "commands_run:\n  ./build-128/my_ethmsb_margin_cert\n";
    std::cerr << "first_failure:\n";
    std::cerr << "  impl=h3compat_l20_l02\n";
    std::cerr << "  root_k=" << r.root_k << "\n";
    std::cerr << "  current_k=" << r.current_k << "\n";
    std::cerr << "  depth=" << r.depth << "\n";
    std::cerr << "  stage=" << r.stage << "\n";
    std::cerr << "  kappa=" << r.kappa << "\n";
    std::cerr << "  m_or_a_b=" << r.plaintext << "\n";
    std::cerr << "  expected=" << r.expected << "\n";
    std::cerr << "  actual=" << r.actual << "\n";
    std::cerr << "  phase_before_offset=" << hex64(r.phase_before) << "\n";
    std::cerr << "  phase_after_offset=" << hex64(r.phase_after) << "\n";
    std::cerr << "  distance_to_half=" << hex64(r.distance_to_half) << "\n";
    std::cerr << "  measured_switch_band=" << hex64(r.measured_switch_band)
              << "\n";
    std::cerr << "  safe_margin=" << hex64(r.safe_margin) << "\n";
    std::cerr << "  output_phase=" << hex64(r.output_phase) << "\n";
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
    std::cerr << "  run: git diff -- include/params/128bit.hpp include/params/concrete.hpp\n";
    std::cerr << "END_NEED_INFO\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const Context ctx = make_context(opt.seed);
        std::ofstream file;
        std::ostream* out = &std::cout;
        if (!opt.csv_path.empty()) {
            file.open(opt.csv_path);
            out = &file;
        }
        emit_header(*out);

        const Torus measured_switch_band = opt.switch_band_width;
        const Torus safety_budget =
            opt.switch_band_width + opt.output_noise_margin;
        Summary summary;
        std::uint64_t cert_id = 0;
        for (int k = 1; k <= opt.max_k; ++k) {
            for (const Torus m : candidates_for_k(k, opt.kappa)) {
                const int trials =
                    opt.mode == "encrypted" ? opt.trials_per_case : 1;
                for (int trial = 0; trial < trials; ++trial) {
                    const auto ct = make_input(m, k, ctx, opt.mode);
                    ethmsb_trace(*out, ctx, ct, k, k, opt.kappa, m,
                                  my_ethmsb::BOOL_ONE, 0, "base",
                                  measured_switch_band, safety_budget, cert_id,
                                  summary);
                    if (summary.valid_failure || summary.unsafe) {
                        print_need_info(summary);
                        return 1;
                    }
                }
            }
            std::cerr << "MARGIN_CERT_DONE root_k=" << k
                      << " min_distance=" << hex64(summary.min_distance)
                      << " min_safe_margin=" << hex64(summary.min_safe_margin)
                      << "\n";
        }
        std::cerr << "MARGIN_CERT_RESULT failures=0 min_distance="
                  << hex64(summary.min_distance) << " min_safe_margin="
                  << hex64(summary.min_safe_margin) << "\n";
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: margin_cert\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n  ./build-128/my_ethmsb_margin_cert\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
    return 0;
}
