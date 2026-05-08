#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "my_ethmsb_fixed.hpp"
#include "my_he3db_compat_params.hpp"
#include "tlwe.hpp"

namespace {

using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Options {
    std::string impl = "all";
    std::string mode = "all";
    bool stop_after_fail = false;
    Torus noise_budget = 0;
};

struct PbsParamCase {
    int k;
    int kappa;
    Torus offset;
    Torus out_value;
};

struct FirstFailure {
    bool seen = false;
    std::string impl;
    std::string mode;
    int k = 0;
    int kappa = 0;
    std::string m_or_phase;
    bool expected = false;
    bool actual = false;
    Torus phase_in = 0;
    Torus phase_out = 0;
    Torus offset = 0;
    Torus out_value = 0;
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
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> char* {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing ") + name);
            return argv[++i];
        };
        if (a == "--impl")
            opt.impl = need("--impl");
        else if (a == "--mode")
            opt.mode = need("--mode");
        else if (a == "--valid-only")
            opt.mode = "valid-trivial";
        else if (a == "--stop-after-fail")
            opt.stop_after_fail = true;
        else if (a == "--noise-budget")
            opt.noise_budget = static_cast<Torus>(std::stoull(need("--noise-budget")));
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    return opt;
}

Torus delta(const int k) { return my_ethmsb::delta(k); }

Torus encode(const Torus m, const int k)
{
    return static_cast<Torus>(Wide{m} * Wide{delta(k)});
}

bool expected_from_phase(const Torus phase, const Torus offset)
{
    return ((phase + offset) & (Torus{1} << 63)) != 0;
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

std::vector<PbsParamCase> build_param_cases()
{
    std::vector<PbsParamCase> cases;
    for (const int k : {1, 2, 3, 4, 5}) {
        for (const Torus out_value : {my_ethmsb::BOOL_ONE, Torus{1} << 58}) {
            cases.push_back({k, 0, delta(k) / 2, out_value});
        }
    }
    for (const int k : {6, 7, 8, 9, 12, 17, 25, 33}) {
        for (const int kappa : {4, 5}) {
            if (k <= kappa) continue;
            const int w = k - kappa - 1;
            const Torus offset =
                static_cast<Torus>((Wide{(Torus{1} << w) + Torus{1}} *
                                    Wide{delta(k)}) /
                                   Wide{2});
            for (const Torus out_value :
                 {my_ethmsb::BOOL_ONE, Torus{1} << 58}) {
                cases.push_back({k, kappa, offset, out_value});
            }
        }
    }
    return cases;
}

Torus ring_step_l2(const int nbit) { return Torus{1} << (64 - (nbit + 1)); }

std::vector<Torus> threshold_phases(const PbsParamCase& pc,
                                    const Torus ring_step)
{
    const Torus threshold_input = (Torus{1} << 63) - pc.offset;
    const Torus quarter = std::max<Torus>(1, ring_step / 4);
    const Torus half = std::max<Torus>(1, ring_step / 2);
    std::set<Torus> values;
    for (const Torus d : {Torus{1}, Torus{2}, quarter, half, ring_step}) {
        values.insert(threshold_input - d);
        values.insert(threshold_input + d);
    }
    values.insert(threshold_input);
    return {values.begin(), values.end()};
}

std::vector<Torus> encrypted_messages(const PbsParamCase& pc)
{
    std::set<Torus> values;
    values.insert(0);
    values.insert(1);
    values.insert((Torus{1} << (pc.k - 1)) - 1);
    values.insert(Torus{1} << (pc.k - 1));
    values.insert((Torus{1} << pc.k) - 1);
    if (pc.kappa > 0) {
        const int w = pc.k - pc.kappa - 1;
        const Torus guard_weight = Torus{1} << w;
        const Torus half = Torus{1} << (pc.k - 1);
        if (half > guard_weight) values.insert(half - guard_weight - 1);
        values.insert(half);
    }
    return {values.begin(), values.end()};
}

std::vector<std::pair<Torus, std::string>> valid_phase_cases(
    const PbsParamCase& pc, const Torus noise_budget)
{
    std::vector<std::pair<Torus, std::string>> out;
    for (const Torus m : encrypted_messages(pc)) {
        const Torus base = encode(m, pc.k);
        out.emplace_back(base, std::to_string(m) + ":noise0");
        if (noise_budget != 0) {
            out.emplace_back(base + noise_budget,
                             std::to_string(m) + ":noise+" +
                                 std::to_string(noise_budget));
            out.emplace_back(base - noise_budget,
                             std::to_string(m) + ":noise-" +
                                 std::to_string(noise_budget));
        }
    }
    return out;
}

void emit_valid_row(const std::string& impl, const std::string& mode,
                    const PbsParamCase& pc, const std::string& label,
                    const bool expected, const bool actual,
                    const Torus phase_in, const Torus phase_out)
{
    std::cout << impl << ",VALID_ENCODING_CASE," << mode << "," << pc.k
              << "," << pc.kappa
              << "," << hex64(pc.offset) << "," << hex64(pc.out_value)
              << "," << label << "," << expected << "," << actual << ","
              << hex64(phase_in) << "," << hex64(phase_out) << ","
              << hex64(phase_in + pc.offset) << ","
              << hex64(dist_to_half(phase_in + pc.offset))
              << ",NA,NA,NA," << (expected == actual) << "\n";
}

void emit_calibration_row(const std::string& impl, const PbsParamCase& pc,
                          const Torus raw_phase, const bool expected,
                          const bool actual, const Torus phase_in,
                          const Torus phase_out, const Torus ring_step,
                          const Torus p0_phase, const Torus p0_lift)
{
    const Torus after = phase_in + pc.offset;
    const Torus distance = dist_to_half(after);
    const std::string classification =
        distance <= ring_step / 2 ? "inside_switch_band"
                                  : "outside_switch_band";
    std::cout << impl << ",THRESHOLD_CALIBRATION_CASE,trivial," << pc.k
              << "," << pc.kappa << "," << hex64(pc.offset) << ","
              << hex64(pc.out_value) << "," << hex64(raw_phase) << ","
              << expected << "," << actual << "," << hex64(phase_in) << ","
              << hex64(phase_out) << "," << hex64(after) << ","
              << hex64(distance) << "," << classification << ","
              << hex64(p0_phase) << "," << hex64(p0_lift) << ",1\n";
}

void remember_failure(FirstFailure& first, const std::string& impl,
                      const std::string& mode, const PbsParamCase& pc,
                      const std::string& label, const bool expected,
                      const bool actual, const Torus phase_in,
                      const Torus phase_out)
{
    if (first.seen || expected == actual) return;
    first = {true, impl, mode, pc.k, pc.kappa, label, expected, actual,
             phase_in, phase_out, pc.offset, pc.out_value};
}

template <class P2, class P0, class KS20, class BR02>
struct L20L02Context {
    std::string impl;
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;

    Torus phase(const TFHEpp::TLWE<P2>& ct) const
    {
        return static_cast<Torus>(TFHEpp::tlweSymPhase<P2>(ct, key2));
    }

    TFHEpp::TLWE<P2> trivial(const Torus body) const
    {
        TFHEpp::TLWE<P2> ct = {};
        ct[P2::k * P2::n] = static_cast<typename P2::T>(body);
        return ct;
    }

    TFHEpp::TLWE<P2> encrypt(const Torus body) const
    {
        return TFHEpp::tlweSymEncrypt<P2>(static_cast<typename P2::T>(body),
                                          key2);
    }

    TFHEpp::TLWE<P2> pbs(TFHEpp::TLWE<P2> in, const PbsParamCase& pc) const
    {
        my_ethmsb::add_const_inplace<P2>(in, pc.offset);
        TFHEpp::TLWE<P0> in0;
        TFHEpp::IdentityKeySwitch<KS20>(in0, in, *iksk);
        const Torus half = pc.out_value / 2;
        TFHEpp::Polynomial<P2> tv;
        tv.fill(Torus{0} - half);
        TFHEpp::TLWE<P2> out;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out, in0, *bkfft, tv);
        my_ethmsb::add_const_inplace<P2>(out, half);
        return out;
    }

    std::pair<Torus, Torus> p0_after_offset(TFHEpp::TLWE<P2> in,
                                            const PbsParamCase& pc) const
    {
        my_ethmsb::add_const_inplace<P2>(in, pc.offset);
        TFHEpp::TLWE<P0> in0;
        TFHEpp::IdentityKeySwitch<KS20>(in0, in, *iksk);
        const auto p0 = TFHEpp::tlweSymPhase<P0>(in0, key0);
        constexpr int p0_bits = std::numeric_limits<typename P0::T>::digits;
        const Torus lift = static_cast<Torus>(p0) << (64 - p0_bits);
        return {static_cast<Torus>(p0), lift};
    }

    Torus ring_step() const { return ring_step_l2(P2::nbit); }
};

template <class P2, class BR22>
struct DirectContext {
    std::string impl;
    TFHEpp::Key<P2> key2;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR22>> bkfft;

    Torus phase(const TFHEpp::TLWE<P2>& ct) const
    {
        return static_cast<Torus>(TFHEpp::tlweSymPhase<P2>(ct, key2));
    }

    TFHEpp::TLWE<P2> trivial(const Torus body) const
    {
        TFHEpp::TLWE<P2> ct = {};
        ct[P2::k * P2::n] = static_cast<typename P2::T>(body);
        return ct;
    }

    TFHEpp::TLWE<P2> encrypt(const Torus body) const
    {
        return TFHEpp::tlweSymEncrypt<P2>(static_cast<typename P2::T>(body),
                                          key2);
    }

    TFHEpp::TLWE<P2> pbs(TFHEpp::TLWE<P2> in, const PbsParamCase& pc) const
    {
        my_ethmsb::add_const_inplace<P2>(in, pc.offset);
        const Torus half = pc.out_value / 2;
        TFHEpp::Polynomial<P2> tv;
        tv.fill(Torus{0} - half);
        TFHEpp::TLWE<P2> out;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR22>(out, in, *bkfft, tv);
        my_ethmsb::add_const_inplace<P2>(out, half);
        return out;
    }

    std::pair<Torus, Torus> p0_after_offset(const TFHEpp::TLWE<P2>&,
                                            const PbsParamCase&) const
    {
        return {0, 0};
    }

    Torus ring_step() const { return ring_step_l2(P2::nbit); }
};

template <class Context>
int run_context(const Context& ctx, const Options& opt, FirstFailure& first)
{
    int failures = 0;
    for (const auto& pc : build_param_cases()) {
        if (opt.mode == "all" || opt.mode == "trivial") {
            for (const Torus raw_phase : threshold_phases(pc, ctx.ring_step())) {
                const auto in = ctx.trivial(raw_phase);
                const Torus phase_in = ctx.phase(in);
                const auto out = ctx.pbs(in, pc);
                const Torus phase_out = ctx.phase(out);
                const bool expected = expected_from_phase(raw_phase, pc.offset);
                const bool actual = closer_to_one(phase_out, pc.out_value);
                const auto [p0_phase, p0_lift] =
                    ctx.p0_after_offset(in, pc);
                emit_calibration_row(ctx.impl, pc, raw_phase, expected, actual,
                                     phase_in, phase_out, ctx.ring_step(),
                                     p0_phase, p0_lift);
            }
        }
        if (opt.mode == "all" || opt.mode == "encrypted") {
            for (const Torus m : encrypted_messages(pc)) {
                const Torus body = encode(m, pc.k);
                const auto in = ctx.encrypt(body);
                const Torus phase_in = ctx.phase(in);
                const auto out = ctx.pbs(in, pc);
                const Torus phase_out = ctx.phase(out);
                const bool expected = expected_from_phase(body, pc.offset);
                const bool actual = closer_to_one(phase_out, pc.out_value);
                const std::string label = std::to_string(m);
                emit_valid_row(ctx.impl, "encrypted", pc, label, expected,
                               actual, phase_in, phase_out);
                failures += expected != actual;
                remember_failure(first, ctx.impl, "encrypted", pc, label,
                                 expected, actual, phase_in, phase_out);
                if (opt.stop_after_fail && expected != actual) return failures;
            }
        }
        if (opt.mode == "all" || opt.mode == "valid-trivial") {
            for (const auto& [raw_phase, label] :
                 valid_phase_cases(pc, opt.noise_budget)) {
                const auto in = ctx.trivial(raw_phase);
                const Torus phase_in = ctx.phase(in);
                const auto out = ctx.pbs(in, pc);
                const Torus phase_out = ctx.phase(out);
                const bool expected = expected_from_phase(raw_phase, pc.offset);
                const bool actual = closer_to_one(phase_out, pc.out_value);
                emit_valid_row(ctx.impl, "valid-trivial", pc, label, expected,
                               actual, phase_in, phase_out);
                failures += expected != actual;
                remember_failure(first, ctx.impl, "valid-trivial", pc, label,
                                 expected, actual, phase_in, phase_out);
                if (opt.stop_after_fail && expected != actual) return failures;
            }
        }
    }
    return failures;
}

int run_direct(const Options& opt, FirstFailure& first)
{
    using namespace my_ethmsb_params;
    using P2 = my_v10_lvl2param;
    using BR22 = my_v10_lvl22param;
    std::mt19937_64 rng(0);
    DirectContext<P2, BR22> ctx;
    ctx.impl = "direct_v10_l22";
    ctx.key2 = make_v10_direct_keys(rng).key2;
    ctx.bkfft = std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR22>>();
    TFHEpp::bkfftgen<BR22>(*ctx.bkfft, ctx.key2, ctx.key2);
    return run_context(ctx, opt, first);
}

int run_v10_l20_l02(const Options& opt, FirstFailure& first)
{
    using P2 = TFHEpp::lvl2param;
    using P0 = TFHEpp::lvl0param;
    using KS20 = TFHEpp::lvl20param;
    using BR02 = TFHEpp::lvl02param;
    TFHEpp::SecretKey sk;
    L20L02Context<P2, P0, KS20, BR02> ctx;
    ctx.impl = "v10_l20_l02";
    ctx.key2 = sk.key.get<P2>();
    ctx.key0 = sk.key.get<P0>();
    ctx.iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft = std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    return run_context(ctx, opt, first);
}

int run_h3compat(const Options& opt, FirstFailure& first)
{
    using namespace my_ethmsb_params;
    using P2 = my_h3_lvl2param;
    using P0 = my_h3_lvl0param;
    using KS20 = my_h3_lvl20param;
    using BR02 = my_h3_lvl02param;
    std::mt19937_64 rng(0);
    auto keys = make_h3compat_keys(rng);
    L20L02Context<P2, P0, KS20, BR02> ctx;
    ctx.impl = "h3compat_l20_l02";
    ctx.key2 = keys.key2;
    ctx.key0 = keys.key0;
    ctx.iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft = std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    return run_context(ctx, opt, first);
}

void print_need_info(const FirstFailure& first)
{
    std::cerr << "BEGIN_NEED_INFO\n";
    std::cerr << "stage=pbs_matrix\n";
    std::cerr << "what_failed=PBS primitive matrix found a mismatch\n";
    std::cerr << "commands_run:\n";
    std::cerr << "  ./build-128/my_ethmsb_pbs_matrix_tests --stop-after-fail\n";
    std::cerr << "first_failure_summary:\n";
    std::cerr << "  impl=" << first.impl << "\n";
    std::cerr << "  mode=" << first.mode << "\n";
    std::cerr << "  k=" << first.k << "\n";
    std::cerr << "  kappa=" << first.kappa << "\n";
    std::cerr << "  m=" << first.m_or_phase << "\n";
    std::cerr << "  expected=" << first.expected << "\n";
    std::cerr << "  actual=" << first.actual << "\n";
    std::cerr << "  phase_in=" << hex64(first.phase_in) << "\n";
    std::cerr << "  phase_out=" << hex64(first.phase_out) << "\n";
    std::cerr << "params:\n";
    std::cerr << "  v10_lvl0_T_bits="
              << std::numeric_limits<typename TFHEpp::lvl0param::T>::digits
              << "\n";
    std::cerr << "  h3compat_lvl0_T_bits="
              << std::numeric_limits<
                     typename my_ethmsb_params::my_h3_lvl0param::T>::digits
              << "\n";
    std::cerr << "  lvl2_T_bits=64\n";
    std::cerr << "  lvl20_t_basebit=v10(" << TFHEpp::lvl20param::t << ","
              << TFHEpp::lvl20param::basebit << ") h3compat("
              << my_ethmsb_params::my_h3_lvl20param::t << ","
              << my_ethmsb_params::my_h3_lvl20param::basebit << ")\n";
    std::cerr << "  lvl02_target_l_Bgbit=v10(" << TFHEpp::lvl2param::l << ","
              << TFHEpp::lvl2param::Bgbit << ") h3compat("
              << my_ethmsb_params::my_h3_lvl2param::l << ","
              << my_ethmsb_params::my_h3_lvl2param::Bgbit << ")\n";
    std::cerr << "vendor_diff:\n";
    std::cerr << "  run: git diff -- include/params/128bit.hpp include/params/concrete.hpp\n";
    std::cerr << "api_signatures_needed:\n";
    std::cerr << "  rg -n \"IdentityKeySwitch|GateBootstrappingTLWE2TLWEFFT|bkfftgen|ikskgen\" include\n";
    std::cerr << "END_NEED_INFO\n";
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        std::cout << "impl,case_type,mode,k,kappa,offset_hex,out_value_hex,m_or_phase,expected,actual,phase_in_hex,phase_out_hex,phase_after_offset_hex,distance_to_half,classification,p0_phase_hex,p0_lift_hex,pass\n";
        FirstFailure first;
        int failures = 0;
        if (opt.impl == "all" || opt.impl == "direct_v10_l22") {
            failures += run_direct(opt, first);
            if (opt.stop_after_fail && first.seen) {
                print_need_info(first);
                return 1;
            }
        }
        if (opt.impl == "all" || opt.impl == "v10_l20_l02") {
            failures += run_v10_l20_l02(opt, first);
            if (opt.stop_after_fail && first.seen) {
                print_need_info(first);
                return 1;
            }
        }
        if (opt.impl == "all" || opt.impl == "h3compat_l20_l02") {
            failures += run_h3compat(opt, first);
        }
        if (first.seen) print_need_info(first);
        return failures == 0 ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage=pbs_matrix\n";
        std::cerr << "what_failed=" << ex.what() << "\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
}
