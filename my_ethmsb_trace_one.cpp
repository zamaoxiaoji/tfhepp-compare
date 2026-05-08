#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

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
    std::string mode = "encrypted";
    int k = 17;
    int kappa = 5;
    Torus m = 0;
    std::string out_value = "bool_one";
    std::uint64_t seed = 0;
};

std::string hex64(const Torus v)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
}

template <class T>
std::string hex_t(const T v)
{
    std::ostringstream os;
    os << "0x" << std::hex << +static_cast<std::uint64_t>(v);
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
        else if (a == "--k")
            opt.k = std::stoi(need("--k"));
        else if (a == "--kappa")
            opt.kappa = std::stoi(need("--kappa"));
        else if (a == "--m")
            opt.m = static_cast<Torus>(std::stoull(need("--m")));
        else if (a == "--out-value")
            opt.out_value = need("--out-value");
        else if (a == "--seed")
            opt.seed = std::stoull(need("--seed"));
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

Torus mask_for(const int k) { return k == 64 ? ~Torus{0} : ((Torus{1} << k) - 1); }

bool closer_to_one(const Torus phase, const Torus out_value)
{
    return my_ethmsb::torus_abs_centered(phase - out_value) <
           my_ethmsb::torus_abs_centered(phase);
}

template <class P2, class P0, class KS20, class BR02>
class L20L02Tracer {
public:
    using TLWE2 = TFHEpp::TLWE<P2>;
    using Key2 = TFHEpp::Key<P2>;
    using Key0 = TFHEpp::Key<P0>;

    L20L02Tracer(const Key2& key2, const Key0& key0,
                 const TFHEpp::KeySwitchingKey<KS20>& iksk20,
                 const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02)
        : key2_(key2), key0_(key0), iksk20_(iksk20), bkfft02_(bkfft02)
    {
    }

    TLWE2 trace_ethmsb(const TLWE2& ct, const int k, const int kappa,
                       const Torus m, const Torus out_value, const int depth)
    {
        const Torus input_phase = phase64(ct);
        const Torus expected_plain = m & mask_for(k);
        std::cout << "BEGIN_ETHMSB_NODE\n";
        std::cout << "depth=" << depth << "\n";
        std::cout << "k=" << k << "\n";
        std::cout << "kappa=" << kappa << "\n";
        std::cout << "out_value_hex=" << hex64(out_value) << "\n";
        std::cout << "input_phase_hex=" << hex64(input_phase) << "\n";
        std::cout << "input_centered_noise_est_hex="
                  << hex64(input_phase - encode(expected_plain, k)) << "\n";
        std::cout << "expected_msb=" << ((expected_plain >> (k - 1)) & 1) << "\n";

        if (k <= kappa) {
            std::cout << "case=base\n";
            TLWE2 out =
                trace_pbs(ct, k, my_ethmsb::base_offset_for_current_layer(k),
                          out_value, expected_plain, depth);
            std::cout << "END_ETHMSB_NODE\n";
            return out;
        }

        std::cout << "case=recursive\n";
        TLWE2 shifted;
        my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);
        const int suffix_bits = k - kappa;
        const Torus suffix_plain = expected_plain & mask_for(suffix_bits);
        const int w = k - kappa - 1;
        const Torus guard_weight = my_ethmsb::guard_weight(k, kappa);
        const Torus guard_value =
            my_ethmsb::guard_value_for_parent_scale(k, kappa);
        std::cout << "suffix_bits=" << suffix_bits << "\n";
        std::cout << "shifted_phase_hex=" << hex64(phase64(shifted)) << "\n";
        std::cout << "expected_suffix_plain=" << suffix_plain << "\n";
        std::cout << "guard_weight=" << guard_weight << "\n";
        std::cout << "guard_value_hex=" << hex64(guard_value) << "\n";

        TLWE2 guard =
            trace_ethmsb(shifted, suffix_bits, kappa, suffix_plain,
                         guard_value, depth + 1);
        const Torus guard_phase = phase64(guard);
        const bool guard_decoded = closer_to_one(guard_phase, guard_value);
        const bool expected_guard =
            ((suffix_plain >> (suffix_bits - 1)) & Torus{1}) != 0;
        std::cout << "guard_phase_hex=" << hex64(guard_phase) << "\n";
        std::cout << "guard_decoded=" << guard_decoded << "\n";
        std::cout << "expected_guard_bit=" << expected_guard << "\n";

        TLWE2 guarded;
        my_ethmsb::sub<P2>(guarded, ct, guard);
        const Torus guarded_phase = phase64(guarded);
        const Torus final_offset =
            my_ethmsb::gap_offset_for_current_layer(k, kappa);
        std::cout << "guarded_phase_hex=" << hex64(guarded_phase) << "\n";
        std::cout << "w=" << w << "\n";
        std::cout << "final_offset_hex=" << hex64(final_offset) << "\n";
        std::cout << "final_input_after_offset_phase_hex="
                  << hex64(guarded_phase + final_offset) << "\n";

        TLWE2 out = trace_pbs(guarded, k, final_offset, out_value,
                              expected_plain - (expected_guard ? guard_weight : 0),
                              depth);
        const Torus final_phase = phase64(out);
        std::cout << "final_output_phase_hex=" << hex64(final_phase) << "\n";
        std::cout << "final_decoded=" << closer_to_one(final_phase, out_value)
                  << "\n";
        std::cout << "END_ETHMSB_NODE\n";
        return out;
    }

private:
    Torus phase64(const TLWE2& ct) const
    {
        return static_cast<Torus>(TFHEpp::tlweSymPhase<P2>(ct, key2_));
    }

    TLWE2 trace_pbs(TLWE2 in, const int input_bits, const Torus offset,
                    const Torus out_value, const Torus expected_plain,
                    const int depth)
    {
        const Torus before = phase64(in);
        const Torus exact_before = encode(expected_plain & mask_for(input_bits),
                                          input_bits);
        const bool expected =
            ((exact_before + offset) & (Torus{1} << 63)) != 0;
        std::cout << "BEGIN_PBS_CALL\n";
        std::cout << "depth=" << depth << "\n";
        std::cout << "input_bits=" << input_bits << "\n";
        std::cout << "offset_hex=" << hex64(offset) << "\n";
        std::cout << "out_value_hex=" << hex64(out_value) << "\n";
        std::cout << "input_phase_before_offset_hex=" << hex64(before) << "\n";
        std::cout << "input_phase_after_offset_hex=" << hex64(before + offset)
                  << "\n";
        std::cout << "expected_after_offset_msb=" << expected << "\n";

        my_ethmsb::add_const_inplace<P2>(in, offset);
        TFHEpp::TLWE<P0> in0;
        TFHEpp::IdentityKeySwitch<KS20>(in0, in, iksk20_);
        const Torus half = out_value / 2;
        TFHEpp::Polynomial<P2> tv;
        tv.fill(Torus{0} - half);
        TLWE2 out;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR02>(out, in0, bkfft02_, tv);
        my_ethmsb::add_const_inplace<P2>(out, half);

        const Torus out_phase = phase64(out);
        std::cout << "output_phase_hex=" << hex64(out_phase) << "\n";
        std::cout << "decoded=" << closer_to_one(out_phase, out_value) << "\n";
        std::cout << "END_PBS_CALL\n";
        return out;
    }

    const Key2& key2_;
    const Key0& key0_;
    const TFHEpp::KeySwitchingKey<KS20>& iksk20_;
    const TFHEpp::BootstrappingKeyFFT<BR02>& bkfft02_;
};

template <class P2, class BR22>
class DirectTracer {
public:
    using TLWE2 = TFHEpp::TLWE<P2>;
    using Key2 = TFHEpp::Key<P2>;

    DirectTracer(const Key2& key2,
                 const TFHEpp::BootstrappingKeyFFT<BR22>& bkfft)
        : key2_(key2), bkfft_(bkfft)
    {
    }

    TLWE2 trace_ethmsb(const TLWE2& ct, const int k, const int kappa,
                       const Torus m, const Torus out_value, const int depth)
    {
        const Torus input_phase = phase64(ct);
        const Torus expected_plain = m & mask_for(k);
        std::cout << "BEGIN_ETHMSB_NODE\n";
        std::cout << "depth=" << depth << "\n";
        std::cout << "k=" << k << "\n";
        std::cout << "kappa=" << kappa << "\n";
        std::cout << "out_value_hex=" << hex64(out_value) << "\n";
        std::cout << "input_phase_hex=" << hex64(input_phase) << "\n";
        std::cout << "input_centered_noise_est_hex="
                  << hex64(input_phase - encode(expected_plain, k)) << "\n";
        std::cout << "expected_msb=" << ((expected_plain >> (k - 1)) & 1) << "\n";
        if (k <= kappa) {
            std::cout << "case=base\n";
            TLWE2 out =
                trace_pbs(ct, k, my_ethmsb::base_offset_for_current_layer(k),
                          out_value, expected_plain, depth);
            std::cout << "END_ETHMSB_NODE\n";
            return out;
        }

        std::cout << "case=recursive\n";
        TLWE2 shifted;
        my_ethmsb::scalar_mul_pow2<P2>(shifted, ct, kappa);
        const int suffix_bits = k - kappa;
        const Torus suffix_plain = expected_plain & mask_for(suffix_bits);
        const int w = k - kappa - 1;
        const Torus guard_weight = my_ethmsb::guard_weight(k, kappa);
        const Torus guard_value =
            my_ethmsb::guard_value_for_parent_scale(k, kappa);
        std::cout << "suffix_bits=" << suffix_bits << "\n";
        std::cout << "shifted_phase_hex=" << hex64(phase64(shifted)) << "\n";
        std::cout << "expected_suffix_plain=" << suffix_plain << "\n";
        std::cout << "guard_weight=" << guard_weight << "\n";
        std::cout << "guard_value_hex=" << hex64(guard_value) << "\n";
        TLWE2 guard =
            trace_ethmsb(shifted, suffix_bits, kappa, suffix_plain,
                         guard_value, depth + 1);
        const Torus guard_phase = phase64(guard);
        const bool guard_decoded = closer_to_one(guard_phase, guard_value);
        const bool expected_guard =
            ((suffix_plain >> (suffix_bits - 1)) & Torus{1}) != 0;
        std::cout << "guard_phase_hex=" << hex64(guard_phase) << "\n";
        std::cout << "guard_decoded=" << guard_decoded << "\n";
        std::cout << "expected_guard_bit=" << expected_guard << "\n";

        TLWE2 guarded;
        my_ethmsb::sub<P2>(guarded, ct, guard);
        const Torus guarded_phase = phase64(guarded);
        const Torus final_offset =
            my_ethmsb::gap_offset_for_current_layer(k, kappa);
        std::cout << "guarded_phase_hex=" << hex64(guarded_phase) << "\n";
        std::cout << "w=" << w << "\n";
        std::cout << "final_offset_hex=" << hex64(final_offset) << "\n";
        std::cout << "final_input_after_offset_phase_hex="
                  << hex64(guarded_phase + final_offset) << "\n";
        TLWE2 out = trace_pbs(guarded, k, final_offset, out_value,
                              expected_plain - (expected_guard ? guard_weight : 0),
                              depth);
        const Torus final_phase = phase64(out);
        std::cout << "final_output_phase_hex=" << hex64(final_phase) << "\n";
        std::cout << "final_decoded=" << closer_to_one(final_phase, out_value)
                  << "\n";
        std::cout << "END_ETHMSB_NODE\n";
        return out;
    }

private:
    Torus phase64(const TLWE2& ct) const
    {
        return static_cast<Torus>(TFHEpp::tlweSymPhase<P2>(ct, key2_));
    }

    TLWE2 trace_pbs(TLWE2 in, const int input_bits, const Torus offset,
                    const Torus out_value, const Torus expected_plain,
                    const int depth)
    {
        const Torus before = phase64(in);
        const Torus exact_before = encode(expected_plain & mask_for(input_bits),
                                          input_bits);
        const bool expected =
            ((exact_before + offset) & (Torus{1} << 63)) != 0;
        std::cout << "BEGIN_PBS_CALL\n";
        std::cout << "depth=" << depth << "\n";
        std::cout << "input_bits=" << input_bits << "\n";
        std::cout << "offset_hex=" << hex64(offset) << "\n";
        std::cout << "out_value_hex=" << hex64(out_value) << "\n";
        std::cout << "input_phase_before_offset_hex=" << hex64(before) << "\n";
        std::cout << "input_phase_after_offset_hex=" << hex64(before + offset)
                  << "\n";
        std::cout << "expected_after_offset_msb=" << expected << "\n";

        my_ethmsb::add_const_inplace<P2>(in, offset);
        const Torus half = out_value / 2;
        TFHEpp::Polynomial<P2> tv;
        tv.fill(Torus{0} - half);
        TLWE2 out;
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<BR22>(out, in, bkfft_, tv);
        my_ethmsb::add_const_inplace<P2>(out, half);
        const Torus out_phase = phase64(out);
        std::cout << "output_phase_hex=" << hex64(out_phase) << "\n";
        std::cout << "decoded=" << closer_to_one(out_phase, out_value) << "\n";
        std::cout << "END_PBS_CALL\n";
        return out;
    }

    const Key2& key2_;
    const TFHEpp::BootstrappingKeyFFT<BR22>& bkfft_;
};

template <class P>
TFHEpp::TLWE<P> make_input(const std::string& mode, const Torus m, const int k,
                           const TFHEpp::Key<P>& key)
{
    const typename P::T phase = static_cast<typename P::T>(encode(m, k));
    if (mode == "trivial") return my_ethmsb_params::trivial_tlwe<P>(phase);
    if (mode == "encrypted") return TFHEpp::tlweSymEncrypt<P>(phase, key);
    throw std::runtime_error("unknown mode: " + mode);
}

}  // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parse(argc, argv);
        const Torus out_value =
            opt.out_value == "guard" ? (delta(opt.k) << (opt.k - opt.kappa - 1))
                                     : my_ethmsb::BOOL_ONE;
        std::mt19937_64 rng(opt.seed);

        std::cout << "TRACE_CONFIG impl=" << opt.impl << " mode=" << opt.mode
                  << " k=" << opt.k << " kappa=" << opt.kappa
                  << " m=" << opt.m << " seed=" << opt.seed
                  << " out_value_hex=" << hex64(out_value) << "\n";

        if (opt.impl == "h3compat_l20_l02") {
            using namespace my_ethmsb_params;
            auto keys = make_h3compat_keys(rng);
            auto iksk = std::make_unique_for_overwrite<
                TFHEpp::KeySwitchingKey<my_h3_lvl20param>>();
            auto bkfft = std::make_unique_for_overwrite<
                TFHEpp::BootstrappingKeyFFT<my_h3_lvl02param>>();
            TFHEpp::ikskgen<my_h3_lvl20param>(*iksk, keys.key2, keys.key0);
            TFHEpp::bkfftgen<my_h3_lvl02param>(*bkfft, keys.key0, keys.key2);
            auto ct = make_input<my_h3_lvl2param>(opt.mode, opt.m, opt.k,
                                                  keys.key2);
            L20L02Tracer<my_h3_lvl2param, my_h3_lvl0param, my_h3_lvl20param,
                         my_h3_lvl02param>
                tracer(keys.key2, keys.key0, *iksk, *bkfft);
            auto out = tracer.trace_ethmsb(ct, opt.k, opt.kappa, opt.m,
                                           out_value, 0);
            const Torus phase =
                TFHEpp::tlweSymPhase<my_h3_lvl2param>(out, keys.key2);
            std::cout << "TRACE_RESULT decoded="
                      << closer_to_one(phase, out_value)
                      << " expected=" << ((opt.m >> (opt.k - 1)) & 1)
                      << " phase_hex=" << hex64(phase) << "\n";
        }
        else if (opt.impl == "v10_l20_l02") {
            using P2 = TFHEpp::lvl2param;
            using P0 = TFHEpp::lvl0param;
            using KS20 = TFHEpp::lvl20param;
            using BR02 = TFHEpp::lvl02param;
            TFHEpp::SecretKey sk;
            const auto key2 = sk.key.get<P2>();
            const auto key0 = sk.key.get<P0>();
            auto iksk = std::make_unique_for_overwrite<TFHEpp::KeySwitchingKey<KS20>>();
            auto bkfft = std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<BR02>>();
            TFHEpp::ikskgen<KS20>(*iksk, key2, key0);
            TFHEpp::bkfftgen<BR02>(*bkfft, key0, key2);
            auto ct = make_input<P2>(opt.mode, opt.m, opt.k, key2);
            L20L02Tracer<P2, P0, KS20, BR02> tracer(key2, key0, *iksk,
                                                    *bkfft);
            auto out = tracer.trace_ethmsb(ct, opt.k, opt.kappa, opt.m,
                                           out_value, 0);
            const Torus phase = TFHEpp::tlweSymPhase<P2>(out, key2);
            std::cout << "TRACE_RESULT decoded="
                      << closer_to_one(phase, out_value)
                      << " expected=" << ((opt.m >> (opt.k - 1)) & 1)
                      << " phase_hex=" << hex64(phase) << "\n";
        }
        else if (opt.impl == "direct_v10_l22") {
            using namespace my_ethmsb_params;
            auto keys = make_v10_direct_keys(rng);
            auto bkfft = std::make_unique_for_overwrite<
                TFHEpp::BootstrappingKeyFFT<my_v10_lvl22param>>();
            TFHEpp::bkfftgen<my_v10_lvl22param>(*bkfft, keys.key2, keys.key2);
            auto ct = make_input<my_v10_lvl2param>(opt.mode, opt.m, opt.k,
                                                   keys.key2);
            DirectTracer<my_v10_lvl2param, my_v10_lvl22param> tracer(
                keys.key2, *bkfft);
            auto out = tracer.trace_ethmsb(ct, opt.k, opt.kappa, opt.m,
                                           out_value, 0);
            const Torus phase =
                TFHEpp::tlweSymPhase<my_v10_lvl2param>(out, keys.key2);
            std::cout << "TRACE_RESULT decoded="
                      << closer_to_one(phase, out_value)
                      << " expected=" << ((opt.m >> (opt.k - 1)) & 1)
                      << " phase_hex=" << hex64(phase) << "\n";
        }
        else {
            throw std::runtime_error("unknown impl: " + opt.impl);
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "BEGIN_NEED_INFO\n";
        std::cerr << "stage=trace\n";
        std::cerr << "what_failed=" << ex.what() << "\n";
        std::cerr << "commands_run:\n  ./build-128/my_ethmsb_trace_one ...\n";
        std::cerr << "vendor_diff:\n";
        std::cerr << "  run: git diff -- include/params/128bit.hpp include/params/concrete.hpp\n";
        std::cerr << "api_signatures_needed:\n";
        std::cerr << "  rg -n \"IdentityKeySwitch|GateBootstrappingTLWE2TLWEFFT|bkfftgen|ikskgen\" include\n";
        std::cerr << "END_NEED_INFO\n";
        return 1;
    }
    return 0;
}
