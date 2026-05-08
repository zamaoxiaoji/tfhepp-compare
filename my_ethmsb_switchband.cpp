#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
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

using Torus = std::uint64_t;
using Wide = unsigned __int128;

struct Options {
    std::string impl = "h3compat_l20_l02";
    std::vector<int> input_bits = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 17, 25, 33};
    std::string out_value = "bool";
    int repeat = 1;
    std::string csv_path;
};

struct PbsCase {
    int input_bits;
    Torus offset;
    Torus out_value;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

Torus delta(const int k) { return my_ethmsb::delta(k); }

Torus ring_step_from_nbit(const int nbit)
{
    return Torus{1} << (64 - (nbit + 1));
}

Torus choose_out_value(const int k, const std::string& mode)
{
    if (mode == "bool") return my_ethmsb::BOOL_ONE;
    if (mode != "guard") throw std::runtime_error("unknown --out-value");
    if (k > 5) {
        constexpr int kappa = 5;
        return my_ethmsb::guard_value_for_parent_scale(k, kappa);
    }
    return Torus{1} << 58;
}

PbsCase make_case(const int k, const std::string& out_value_mode)
{
    if (k <= 5)
        return {k, my_ethmsb::base_offset_for_current_layer(k),
                choose_out_value(k, out_value_mode)};
    constexpr int kappa = 5;
    return {k, my_ethmsb::gap_offset_for_current_layer(k, kappa),
            choose_out_value(k, out_value_mode)};
}

std::vector<int> parse_bits(const std::string& s)
{
    if (s == "all") return {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 17, 25, 33};
    std::vector<int> out;
    std::size_t pos = 0;
    while (pos < s.size()) {
        const std::size_t comma = s.find(',', pos);
        out.push_back(std::stoi(s.substr(pos, comma - pos)));
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
        else if (a == "--input-bits")
            opt.input_bits = parse_bits(need("--input-bits"));
        else if (a == "--all")
            opt.input_bits = parse_bits("all");
        else if (a == "--out-value")
            opt.out_value = need("--out-value");
        else if (a == "--repeat")
            opt.repeat = std::stoi(need("--repeat"));
        else if (a == "--csv")
            opt.csv_path = need("--csv");
        else
            throw std::runtime_error("unknown argument: " + a);
    }
    return opt;
}

bool closer_to_one(const Torus phase, const Torus out_value)
{
    return my_ethmsb::torus_abs_centered(phase - out_value) <
           my_ethmsb::torus_abs_centered(phase);
}

template <class P2, class P0, class KS20, class BR02>
struct L20L02Context {
    std::string impl;
    TFHEpp::Key<P2> key2;
    TFHEpp::Key<P0> key0;
    std::unique_ptr<TFHEpp::KeySwitchingKey<KS20>> iksk;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR02>> bkfft;

    TFHEpp::TLWE<P2> trivial(const Torus body) const
    {
        TFHEpp::TLWE<P2> ct = {};
        ct[P2::k * P2::n] = static_cast<typename P2::T>(body);
        return ct;
    }

    Torus phase(const TFHEpp::TLWE<P2>& ct) const
    {
        return TFHEpp::tlweSymPhase<P2>(ct, key2);
    }

    TFHEpp::TLWE<P2> pbs(TFHEpp::TLWE<P2> in, const PbsCase& pc) const
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

    std::pair<Torus, Torus> p0_lift_at_threshold(const PbsCase& pc) const
    {
        auto in = trivial((Torus{1} << 63) - pc.offset);
        my_ethmsb::add_const_inplace<P2>(in, pc.offset);
        TFHEpp::TLWE<P0> in0;
        TFHEpp::IdentityKeySwitch<KS20>(in0, in, *iksk);
        const auto p0 = TFHEpp::tlweSymPhase<P0>(in0, key0);
        constexpr int p0_bits = std::numeric_limits<typename P0::T>::digits;
        return {static_cast<Torus>(p0), static_cast<Torus>(p0) << (64 - p0_bits)};
    }

    Torus ring_step() const { return ring_step_from_nbit(P2::nbit); }
};

template <class P2, class BR22>
struct DirectContext {
    std::string impl;
    TFHEpp::Key<P2> key2;
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<BR22>> bkfft;

    TFHEpp::TLWE<P2> trivial(const Torus body) const
    {
        TFHEpp::TLWE<P2> ct = {};
        ct[P2::k * P2::n] = static_cast<typename P2::T>(body);
        return ct;
    }

    Torus phase(const TFHEpp::TLWE<P2>& ct) const
    {
        return TFHEpp::tlweSymPhase<P2>(ct, key2);
    }

    TFHEpp::TLWE<P2> pbs(TFHEpp::TLWE<P2> in, const PbsCase& pc) const
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

    std::pair<Torus, Torus> p0_lift_at_threshold(const PbsCase&) const
    {
        return {0, 0};
    }

    Torus ring_step() const { return ring_step_from_nbit(P2::nbit); }
};

template <class Context>
void measure_case(std::ostream& out, const Context& ctx, const PbsCase& pc,
                  const int repeat)
{
    const Torus ring_step = ctx.ring_step();
    const Torus half_step = ring_step / 2;
    const Torus scan_step = std::max<Torus>(1, ring_step / 8);
    const Torus threshold = (Torus{1} << 63) - pc.offset;
    Torus first_one = 0;
    Torus last_zero = 0;
    bool saw_one = false;
    bool saw_zero = false;

    for (int rep = 0; rep < repeat; ++rep) {
        for (int i = -32; i <= 32; ++i) {
            const Torus phase = threshold + static_cast<std::int64_t>(i) *
                                                static_cast<std::int64_t>(scan_step);
            const auto ct = ctx.trivial(phase);
            const auto out_ct = ctx.pbs(ct, pc);
            const bool actual = closer_to_one(ctx.phase(out_ct), pc.out_value);
            if (actual && !saw_one) {
                first_one = phase;
                saw_one = true;
            }
            if (!actual) {
                last_zero = phase;
                saw_zero = true;
            }
        }
    }

    const Torus width =
        saw_one && saw_zero && first_one > last_zero ? first_one - last_zero : 0;
    const auto [p0_phase, p0_lift] = ctx.p0_lift_at_threshold(pc);
    out << ctx.impl << ',' << pc.input_bits << ',' << hex64(pc.out_value) << ','
        << hex64(pc.offset) << ',' << hex64(ring_step) << ','
        << hex64(half_step) << ',' << (saw_one ? hex64(first_one) : "NA")
        << ',' << (saw_zero ? hex64(last_zero) : "NA") << ','
        << hex64(width) << '\n';
    out << "# threshold_p0 impl=" << ctx.impl
        << " input_bits=" << pc.input_bits << " p0_phase=" << hex64(p0_phase)
        << " p0_lift=" << hex64(p0_lift) << "\n";
    if (width > 4 * ring_step) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: switchband\n";
        std::cerr << "what_failed: switch band is wider than four ring steps\n";
        std::cerr << "commands_run:\n  ./build-128/my_ethmsb_switchband\n";
        std::cerr << "first_failure:\n  impl=" << ctx.impl
                  << "\n  current_k=" << pc.input_bits
                  << "\n  phase_before_offset=" << hex64(threshold)
                  << "\n  measured_switch_band=" << hex64(width) << "\n";
        std::cerr << "params:\n  P2_T_bits=64\n";
        std::cerr << "END_NEED_INFO\n";
    }
}

template <class Context>
void run_context(std::ostream& out, const Context& ctx, const Options& opt)
{
    for (const int k : opt.input_bits)
        measure_case(out, ctx, make_case(k, opt.out_value), opt.repeat);
}

void run_direct(std::ostream& out, const Options& opt)
{
    using namespace my_ethmsb_params;
    using P2 = my_v10_lvl2param;
    using BR22 = my_v10_lvl22param;
    std::mt19937_64 rng(0);
    DirectContext<P2, BR22> ctx;
    ctx.impl = "direct_v10_l22";
    ctx.key2 = make_v10_direct_keys(rng).key2;
    ctx.bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR22>>();
    TFHEpp::bkfftgen<BR22>(*ctx.bkfft, ctx.key2, ctx.key2);
    run_context(out, ctx, opt);
}

void run_v10_l20_l02(std::ostream& out, const Options& opt)
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
    ctx.bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    run_context(out, ctx, opt);
}

void run_h3compat(std::ostream& out, const Options& opt)
{
    using namespace my_ethmsb_params;
    using P2 = my_h3_lvl2param;
    using P0 = my_h3_lvl0param;
    using KS20 = my_h3_lvl20param;
    using BR02 = my_h3_lvl02param;
    std::mt19937_64 rng(0);
    const auto keys = make_h3compat_keys(rng);
    L20L02Context<P2, P0, KS20, BR02> ctx;
    ctx.impl = "h3compat_l20_l02";
    ctx.key2 = keys.key2;
    ctx.key0 = keys.key0;
    ctx.iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
    ctx.bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
    TFHEpp::ikskgen<KS20>(*ctx.iksk, ctx.key2, ctx.key0);
    TFHEpp::bkfftgen<BR02>(*ctx.bkfft, ctx.key0, ctx.key2);
    run_context(out, ctx, opt);
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        std::ofstream file;
        std::ostream* out = &std::cout;
        if (!opt.csv_path.empty()) {
            file.open(opt.csv_path);
            out = &file;
        }
        *out << "impl,input_bits,out_value_hex,offset_hex,ring_step_hex,ring_half_step_hex,first_one_phase_hex,last_zero_phase_hex,switch_band_width_hex\n";
        if (opt.impl == "direct_v10_l22")
            run_direct(*out, opt);
        else if (opt.impl == "v10_l20_l02")
            run_v10_l20_l02(*out, opt);
        else if (opt.impl == "h3compat_l20_l02")
            run_h3compat(*out, opt);
        else
            throw std::runtime_error("unknown --impl");
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage: switchband\n";
        std::cerr << "what_failed: " << ex.what() << "\n";
        std::cerr << "commands_run:\n  ./build-128/my_ethmsb_switchband\n";
        std::cerr << "vendor_diff:\n";
        std::cerr << "  run: git diff -- include/params/128bit.hpp include/params/concrete.hpp\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
    return 0;
}
