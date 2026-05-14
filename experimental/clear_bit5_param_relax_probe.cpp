#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "cloudkey.hpp"
#include "evalkeygens.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "tlwe.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace experimental::clear_bit5
{
using Lvl0 = TFHEpp::lvl0param;
using Lvl1 = TFHEpp::lvl1param;
using Lvl10 = TFHEpp::lvl10param;
using TLWELvl0 = TFHEpp::TLWE<Lvl0>;
using TLWELvl1 = TFHEpp::TLWE<Lvl1>;
using TFHESecretKey = TFHEpp::SecretKey;

constexpr uint32_t kFixedP = 9;
constexpr uint32_t kFixedBitPosLsb1 = 5;
constexpr uint32_t kTorusBits = std::numeric_limits<Lvl1::T>::digits;
constexpr uint32_t kMessageCells = 1U << kFixedP;
constexpr uint32_t kRotSlots = 2U * Lvl1::n;

template <uint32_t Level, uint32_t Basebit>
struct RelaxedLvl1Param {
    static constexpr int32_t key_value_max = Lvl1::key_value_max;
    static constexpr int32_t key_value_min = Lvl1::key_value_min;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr uint32_t nbit = Lvl1::nbit;
    static constexpr uint32_t n = Lvl1::n;
    static constexpr uint32_t k = Lvl1::k;
    static constexpr uint32_t l = Level;
    static constexpr uint32_t lₐ = Level;
    static constexpr uint32_t Bgbit = Basebit;
    static constexpr uint32_t Bgₐbit = Basebit;
    static constexpr uint32_t Bg = 1U << Bgbit;
    static constexpr uint32_t Bgₐ = 1U << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist = Lvl1::errordist;
    static const inline double α = Lvl1::α;
    using T = typename Lvl1::T;
    static constexpr std::make_signed_t<T> μ = Lvl1::μ;
    static constexpr uint32_t plain_modulus = Lvl1::plain_modulus;
    static constexpr double Δ = Lvl1::Δ;
};

using Relax1 = RelaxedLvl1Param<1, 10>;
using Relax2 = RelaxedLvl1Param<1, 8>;
using Relax3 = RelaxedLvl1Param<1, 6>;

} // namespace experimental::clear_bit5

namespace TFHEpp
{
template <>
inline void TwistFFT<experimental::clear_bit5::Relax1>(
    Polynomial<experimental::clear_bit5::Relax1> &res,
    const PolynomialInFD<experimental::clear_bit5::Relax1> &a)
{
    fftplvl1.execute_direct_torus32(res.data(), a.data());
}

template <>
inline void TwistIFFT<experimental::clear_bit5::Relax1>(
    PolynomialInFD<experimental::clear_bit5::Relax1> &res,
    const Polynomial<experimental::clear_bit5::Relax1> &a)
{
    fftplvl1.execute_reverse_torus32(res.data(), a.data());
}

template <>
inline void TwistFFT<experimental::clear_bit5::Relax2>(
    Polynomial<experimental::clear_bit5::Relax2> &res,
    const PolynomialInFD<experimental::clear_bit5::Relax2> &a)
{
    fftplvl1.execute_direct_torus32(res.data(), a.data());
}

template <>
inline void TwistIFFT<experimental::clear_bit5::Relax2>(
    PolynomialInFD<experimental::clear_bit5::Relax2> &res,
    const Polynomial<experimental::clear_bit5::Relax2> &a)
{
    fftplvl1.execute_reverse_torus32(res.data(), a.data());
}

template <>
inline void TwistFFT<experimental::clear_bit5::Relax3>(
    Polynomial<experimental::clear_bit5::Relax3> &res,
    const PolynomialInFD<experimental::clear_bit5::Relax3> &a)
{
    fftplvl1.execute_direct_torus32(res.data(), a.data());
}

template <>
inline void TwistIFFT<experimental::clear_bit5::Relax3>(
    PolynomialInFD<experimental::clear_bit5::Relax3> &res,
    const Polynomial<experimental::clear_bit5::Relax3> &a)
{
    fftplvl1.execute_reverse_torus32(res.data(), a.data());
}
} // namespace TFHEpp

namespace experimental::clear_bit5
{

template <class TargetP>
struct Lvl01LikeParam {
    using domainP = Lvl0;
    using targetP = TargetP;
    static constexpr uint32_t Addends = 1;
};

template <class TargetP>
std::string TargetParamName()
{
    std::ostringstream out;
    if constexpr (std::is_same_v<TargetP, Lvl1>)
        out << "lvl1param";
    else if constexpr (std::is_same_v<TargetP, Relax1>)
        out << "experimental_lvl1_l1_b10";
    else if constexpr (std::is_same_v<TargetP, Relax2>)
        out << "experimental_lvl1_l1_b8";
    else if constexpr (std::is_same_v<TargetP, Relax3>)
        out << "experimental_lvl1_l1_b6";
    else
        out << "unknown_target";
    out << "(n=" << TargetP::n << ",nbit=" << TargetP::nbit
        << ",l=" << TargetP::l << ",Bgbit=" << TargetP::Bgbit
        << ",alpha=" << TargetP::α << ")";
    return out.str();
}

template <class TargetP>
TFHEpp::Key<TargetP> CopyLvl1Key(const TFHEpp::Key<Lvl1> &key)
{
    static_assert(std::tuple_size<TFHEpp::Key<TargetP>>::value ==
                  std::tuple_size<TFHEpp::Key<Lvl1>>::value);
    TFHEpp::Key<TargetP> out{};
    std::copy(key.begin(), key.end(), out.begin());
    return out;
}

template <class TargetP>
struct StageKeys {
    using PBS = Lvl01LikeParam<TargetP>;
    TFHEpp::Key<TargetP> target_key{};
    std::unique_ptr<TFHEpp::BootstrappingKeyFFT<PBS>> bkfft;

    void Generate(const TFHESecretKey &sk)
    {
        target_key = CopyLvl1Key<TargetP>(sk.key.lvl1);
        bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<PBS>>();
        TFHEpp::bkfftgen<PBS>(*bkfft, sk.key.lvl0, target_key);
    }
};

void SeedGenerator(uint64_t seed)
{
#if defined(USE_BLAKE3)
    TFHEpp::generator = BLAKE3PRNG::BLAKE3PRNG<uint64_t>(seed);
#else
    (void) seed;
#endif
}

enum class Profile {
    Baseline,
    FinalRelax1,
    FinalRelax2,
    FinalRelax3,
    ConversionRelax1,
    ConversionRelax2,
    BitExtractRelax1,
    AllRelax1,
    AllRelax2,
    Custom
};

struct KeyMaterial {
    std::unique_ptr<TFHESecretKey> sk;
    std::unique_ptr<TFHEpp::KeySwitchingKey<Lvl10>> iksk;
    StageKeys<Lvl1> baseline;
    StageKeys<Relax1> relax1;
    StageKeys<Relax2> relax2;
    StageKeys<Relax3> relax3;

    void Generate(uint64_t seed, const std::vector<Profile> &profiles)
    {
        SeedGenerator(seed);
        sk = std::make_unique<TFHESecretKey>();
        iksk = std::make_unique<TFHEpp::KeySwitchingKey<Lvl10>>();
        TFHEpp::ikskgen<Lvl10>(*iksk, sk->key.lvl1, sk->key.lvl0);
        bool need_baseline = false;
        bool need_relax1 = false;
        bool need_relax2 = false;
        bool need_relax3 = false;
        for (Profile p : profiles) {
            switch (p) {
                case Profile::Baseline:
                    need_baseline = true;
                    break;
                case Profile::FinalRelax1:
                case Profile::ConversionRelax1:
                case Profile::BitExtractRelax1:
                case Profile::AllRelax1:
                    need_baseline = true;
                    need_relax1 = true;
                    break;
                case Profile::FinalRelax2:
                case Profile::ConversionRelax2:
                    need_baseline = true;
                    need_relax2 = true;
                    break;
                case Profile::FinalRelax3:
                    need_baseline = true;
                    need_relax3 = true;
                    break;
                case Profile::AllRelax2:
                    need_relax1 = true;
                    need_relax2 = true;
                    break;
                case Profile::Custom:
                    need_baseline = true;
                    break;
            }
        }
        if (need_baseline) baseline.Generate(*sk);
        if (need_relax1) relax1.Generate(*sk);
        if (need_relax2) relax2.Generate(*sk);
        if (need_relax3) relax3.Generate(*sk);
    }
};

template <class TargetP>
const StageKeys<TargetP> &SelectStage(const KeyMaterial &keys)
{
    if constexpr (std::is_same_v<TargetP, Lvl1>)
        return keys.baseline;
    else if constexpr (std::is_same_v<TargetP, Relax1>)
        return keys.relax1;
    else if constexpr (std::is_same_v<TargetP, Relax2>)
        return keys.relax2;
    else if constexpr (std::is_same_v<TargetP, Relax3>)
        return keys.relax3;
    else
        static_assert(TFHEpp::false_v<typename TargetP::T>, "unknown target");
}

std::vector<std::string> Split(const std::string &s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

std::vector<int> SplitInts(const std::string &s)
{
    std::vector<int> out;
    for (const std::string &item : Split(s)) out.push_back(std::stoi(item));
    return out;
}

bool Contains(const std::vector<std::string> &xs, const std::string &x)
{
    return std::find(xs.begin(), xs.end(), x) != xs.end();
}

std::string JsonEscape(const std::string &s)
{
    std::ostringstream out;
    for (char c : s) {
        if (c == '"' || c == '\\') out << '\\';
        if (c == '\n')
            out << "\\n";
        else if (c != '"' && c != '\\')
            out << c;
    }
    return out.str();
}

std::string Hex32(Lvl1::T x)
{
    std::ostringstream out;
    out << "0x" << std::hex << std::setw(8) << std::setfill('0')
        << static_cast<uint32_t>(x);
    return out.str();
}

std::string HexSigned(int64_t x)
{
    return Hex32(static_cast<Lvl1::T>(x));
}

struct ClearBitSpec {
    uint32_t p = kFixedP;
    uint32_t bit_pos_lsb1 = kFixedBitPosLsb1;
    uint32_t j_lsb0 = 4;
    uint32_t k_msb0 = 4;
    uint32_t W = 16;
    uint32_t exact_period_cells = 32;
    Lvl1::T delta = Lvl1::T(1) << (kTorusBits - kFixedP);
    Lvl1::T delta_prime = Lvl1::T(32) * delta;
    Lvl1::T G = Lvl1::T(16) * delta;
    Lvl1::T margin_gap = Lvl1::T(17) * (delta >> 1);
};

ClearBitSpec MakeSpec(uint32_t p, uint32_t bit_pos_lsb1)
{
    if (p != kFixedP || bit_pos_lsb1 != kFixedBitPosLsb1)
        throw std::invalid_argument(
            "this experimental probe is fixed to --p 9 --bit-pos-lsb1 5");
    ClearBitSpec spec;
    spec.p = p;
    spec.bit_pos_lsb1 = bit_pos_lsb1;
    spec.j_lsb0 = bit_pos_lsb1 - 1;
    spec.k_msb0 = p - 1 - spec.j_lsb0;
    spec.W = 1U << spec.j_lsb0;
    spec.exact_period_cells = 2U * spec.W;
    spec.delta = Lvl1::T(1) << (kTorusBits - p);
    spec.delta_prime = Lvl1::T(spec.exact_period_cells) * spec.delta;
    spec.G = Lvl1::T(spec.W) * spec.delta;
    spec.margin_gap = Lvl1::T(spec.W + 1) * (spec.delta >> 1);
    return spec;
}

uint32_t CircularCellDistance(uint32_t a, uint32_t b, uint32_t modulus)
{
    const uint32_t hi = a > b ? a - b : b - a;
    return std::min(hi, modulus - hi);
}

uint32_t RotationToNearestCell(uint32_t rot, uint32_t p)
{
    const uint32_t cells = 1U << p;
    const uint32_t slots_per_cell = kRotSlots / cells;
    return ((rot + slots_per_cell / 2) / slots_per_cell) % cells;
}

uint32_t PhaseToCell(Lvl1::T phase, const ClearBitSpec &spec)
{
    return static_cast<uint32_t>((phase + (spec.delta >> 1)) / spec.delta) &
           ((1U << spec.p) - 1);
}

int64_t SignedTorusDiff(Lvl1::T actual, Lvl1::T expected)
{
    return static_cast<int32_t>(actual - expected);
}

uint64_t AbsTorusDiff(Lvl1::T actual, Lvl1::T expected)
{
    return static_cast<uint64_t>(
        std::llabs(static_cast<long long>(SignedTorusDiff(actual, expected))));
}

struct MarginAudit {
    ClearBitSpec spec;
    uint32_t reachable_count = 0;
    uint32_t unreachable_count = 0;
    uint32_t ambiguous_cell_count = 0;
    uint32_t min_opposite_label_distance_cells = 0;
    bool negacyclic_audit_pass = false;
    bool safe_coarse_margin = false;
    std::array<int, kMessageCells> final_labels{};
    std::map<uint32_t, std::set<int>> reachable_labels;
};

MarginAudit BuildMarginAudit(const ClearBitSpec &spec)
{
    MarginAudit audit;
    audit.spec = spec;

    for (uint32_t m = 0; m < (1U << spec.p); ++m) {
        const uint32_t bit = (m >> spec.j_lsb0) & 1U;
        const uint32_t clear = m - bit * spec.W;
        const int label = m >= (1U << (spec.p - 1)) ? 1 : 0;
        audit.reachable_labels[clear].insert(label);
    }

    audit.reachable_count = static_cast<uint32_t>(audit.reachable_labels.size());
    audit.unreachable_count = (1U << spec.p) - audit.reachable_count;

    uint32_t min_dist = 1U << 30;
    for (const auto &[a, la] : audit.reachable_labels) {
        for (const auto &[b, lb] : audit.reachable_labels) {
            if (a == b) continue;
            if (*la.begin() == *lb.begin()) continue;
            min_dist =
                std::min(min_dist, CircularCellDistance(a, b, 1U << spec.p));
        }
    }
    audit.min_opposite_label_distance_cells = min_dist;

    for (uint32_t c = 0; c < (1U << spec.p); ++c) {
        uint32_t best = 1U << 30;
        std::set<int> labels;
        for (const auto &[r, ls] : audit.reachable_labels) {
            const uint32_t dist = CircularCellDistance(c, r, 1U << spec.p);
            if (dist < best) {
                best = dist;
                labels.clear();
            }
            if (dist == best) labels.insert(ls.begin(), ls.end());
        }
        if (labels.size() != 1) {
            audit.ambiguous_cell_count++;
            audit.final_labels[c] = -1;
        }
        else {
            audit.final_labels[c] = *labels.begin();
        }
    }

    audit.negacyclic_audit_pass = true;
    const uint32_t half = 1U << (spec.p - 1);
    for (uint32_t c = 0; c < half; ++c) {
        const int a = audit.final_labels[c];
        const int b = audit.final_labels[(c + half) & ((1U << spec.p) - 1)];
        if (a < 0 || b < 0 || a == b) audit.negacyclic_audit_pass = false;
    }
    audit.safe_coarse_margin =
        audit.ambiguous_cell_count == 0 &&
        audit.min_opposite_label_distance_cells == spec.W + 1 &&
        audit.negacyclic_audit_pass;
    return audit;
}

template <class TargetP>
TFHEpp::Polynomial<TargetP> ConstantPoly(typename TargetP::T v)
{
    TFHEpp::Polynomial<TargetP> poly{};
    for (auto &x : poly) x = v;
    return poly;
}

template <class TargetP>
TFHEpp::Polynomial<TargetP> BitExtractQHalfPoly(const ClearBitSpec &spec)
{
    TFHEpp::Polynomial<TargetP> poly{};
    const uint32_t slots_per_cell = kRotSlots / (1U << spec.p);
    const typename TargetP::T qhalf =
        typename TargetP::T(1) << (std::numeric_limits<typename TargetP::T>::digits - 1);
    for (uint32_t i = 0; i < TargetP::n; ++i) {
        const uint32_t cell = (i / slots_per_cell) & ((1U << spec.p) - 1);
        const uint32_t bit = (cell >> spec.j_lsb0) & 1U;
        poly[i] = bit ? qhalf : typename TargetP::T(0);
    }
    return poly;
}

template <class TargetP>
TFHEpp::Polynomial<TargetP> FinalWeightedApproxPoly(const MarginAudit &audit)
{
    TFHEpp::Polynomial<TargetP> poly{};
    for (uint32_t i = 0; i < TargetP::n; ++i) {
        const uint32_t cell = RotationToNearestCell(i, audit.spec.p);
        const int label = audit.final_labels[cell];
        const auto pos = static_cast<typename TargetP::T>(TargetP::μ);
        const auto neg = static_cast<typename TargetP::T>(-TargetP::μ);
        poly[i] = label == 1 ? pos : neg;
    }
    return poly;
}

struct PbsTrace {
    uint32_t modswitch_index = 0;
    uint32_t selected_cell = 0;
    uint32_t skip_count = 0;
    uint32_t cmux_count = 0;
    uint32_t skip_period_rot = 0;
};

template <class P>
void GateBootstrappingTrace(
    TFHEpp::TLWE<typename P::targetP> &res,
    const TFHEpp::TLWE<typename P::domainP> &tlwe,
    const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
    const TFHEpp::Polynomial<typename P::targetP> &testvector,
    uint32_t p, uint32_t skip_period_rot, PbsTrace *trace)
{
    static_assert(P::Addends == 1, "experimental trace path expects Addends=1");
    using Target = typename P::targetP;
    alignas(64) TFHEpp::TRLWE<Target> acc{};
    TFHEpp::ModswitchTLWE<typename P::domainP> moded{};
    TFHEpp::BRModSwitch<P, 1>(moded, tlwe);
    const uint32_t body_index =
        static_cast<uint32_t>(moded[P::domainP::k * P::domainP::n]) %
        (2U * Target::n);
    TFHEpp::PolynomialMulByXai<Target>(acc[Target::k], testvector, body_index);

    uint32_t skip = 0;
    uint32_t cmux = 0;
    for (int i = 0; i < P::domainP::k * P::domainP::n; ++i) {
        const uint32_t a = static_cast<uint32_t>(moded[i]) % (2U * Target::n);
        if (a == 0) continue;
        if (skip_period_rot > 1 && (a % skip_period_rot) == 0) {
            skip++;
            continue;
        }
        cmux++;
        TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<P>(acc, bkfft[i], a);
    }
    TFHEpp::SampleExtractIndex<Target>(res, acc, 0);

    if (trace != nullptr) {
        const uint32_t mod_index = (2U * Target::n - body_index) % (2U * Target::n);
        trace->modswitch_index = mod_index;
        trace->selected_cell = RotationToNearestCell(mod_index, p);
        trace->skip_count = skip;
        trace->cmux_count = cmux;
        trace->skip_period_rot = skip_period_rot;
    }
}

template <class TargetP>
void BitExtractQHalfPBS(TLWELvl1 &out, const TLWELvl1 &in,
                        const ClearBitSpec &spec, const KeyMaterial &keys,
                        uint32_t period_cells, PbsTrace *trace)
{
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, in, *keys.iksk);
    const auto &stage = SelectStage<TargetP>(keys);
    const uint32_t slots_per_cell = kRotSlots / (1U << spec.p);
    const uint32_t skip_period_rot = period_cells >= (1U << spec.p)
                                         ? 0
                                         : period_cells * slots_per_cell;
    TFHEpp::TLWE<TargetP> tmp;
    GateBootstrappingTrace<Lvl01LikeParam<TargetP>>(
        tmp, lvl0, *stage.bkfft, BitExtractQHalfPoly<TargetP>(spec), spec.p,
        skip_period_rot, trace);
    out = tmp;
}

template <class TargetP>
void BoolQHalfToValueCenteredPBS(TLWELvl1 &out, const TLWELvl1 &bit_qhalf,
                                 const ClearBitSpec &spec,
                                 const KeyMaterial &keys, PbsTrace *trace)
{
    TLWELvl1 centered = bit_qhalf;
    centered[Lvl1::k * Lvl1::n] -=
        Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 2);

    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, centered, *keys.iksk);
    const auto &stage = SelectStage<TargetP>(keys);
    const typename TargetP::T A = static_cast<typename TargetP::T>(spec.G >> 1);

    TFHEpp::TLWE<TargetP> tmp;
    GateBootstrappingTrace<Lvl01LikeParam<TargetP>>(
        tmp, lvl0, *stage.bkfft, ConstantPoly<TargetP>(A), spec.p, 0, trace);
    tmp[TargetP::k * TargetP::n] += A;
    out = tmp;
}

template <class TargetP>
void WeightedApproxFinalPBS(TLWELvl1 &out, const TLWELvl1 &clear,
                            const MarginAudit &audit, const KeyMaterial &keys,
                            PbsTrace *trace)
{
    TLWELvl0 lvl0;
    TFHEpp::IdentityKeySwitch<Lvl10>(lvl0, clear, *keys.iksk);
    const auto &stage = SelectStage<TargetP>(keys);
    TFHEpp::TLWE<TargetP> tmp;
    GateBootstrappingTrace<Lvl01LikeParam<TargetP>>(
        tmp, lvl0, *stage.bkfft, FinalWeightedApproxPoly<TargetP>(audit),
        audit.spec.p, 0, trace);
    out = tmp;
}

void TLWESubAssign(TLWELvl1 &a, const TLWELvl1 &b)
{
    for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
}

void TLWEAddBodySigned(TLWELvl1 &a, int64_t offset)
{
    a[Lvl1::k * Lvl1::n] += static_cast<Lvl1::T>(offset);
}

TLWELvl1 TrivialTLWE(Lvl1::T phase)
{
    TLWELvl1 out{};
    out[Lvl1::k * Lvl1::n] = phase;
    return out;
}

TLWELvl1 MakeInputCipher(const std::string &input, Lvl1::T phase,
                         const TFHESecretKey &sk)
{
    if (input == "trivial") return TrivialTLWE(phase);
    if (input == "controlled-zero-noise")
        return TFHEpp::tlweSymEncrypt<Lvl1>(phase, 0.0, sk.key.lvl1);
    if (input == "encrypted-noisy")
        return TFHEpp::tlweSymEncrypt<Lvl1>(phase, Lvl1::α, sk.key.lvl1);
    throw std::invalid_argument("unknown --input " + input);
}

Lvl1::T Phase(const TLWELvl1 &ct, const TFHESecretKey &sk)
{
    return TFHEpp::tlweSymPhase<Lvl1>(ct, sk.key.lvl1);
}

uint32_t DecodeQHalfBit(const TLWELvl1 &ct, const TFHESecretKey &sk)
{
    const Lvl1::T phase = Phase(ct, sk);
    const Lvl1::T qhalf =
        Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 1);
    return AbsTorusDiff(phase, qhalf) < AbsTorusDiff(phase, 0) ? 1U : 0U;
}

uint32_t DecodeGuardBit(const TLWELvl1 &ct, const ClearBitSpec &spec,
                        const TFHESecretKey &sk)
{
    const Lvl1::T phase = Phase(ct, sk);
    return AbsTorusDiff(phase, spec.G) < AbsTorusDiff(phase, 0) ? 1U : 0U;
}

uint32_t DecodeLabel(const TLWELvl1 &ct, const TFHESecretKey &sk)
{
    const Lvl1::T phase = Phase(ct, sk);
    return static_cast<std::make_signed_t<Lvl1::T>>(phase) > 0 ? 1U : 0U;
}

Profile ParseProfile(const std::string &s)
{
    if (s == "baseline") return Profile::Baseline;
    if (s == "final_relax1" || s == "relax1") return Profile::FinalRelax1;
    if (s == "final_relax2" || s == "relax2") return Profile::FinalRelax2;
    if (s == "final_relax3" || s == "relax3") return Profile::FinalRelax3;
    if (s == "conversion_relax1") return Profile::ConversionRelax1;
    if (s == "conversion_relax2") return Profile::ConversionRelax2;
    if (s == "bitextract_relax1") return Profile::BitExtractRelax1;
    if (s == "all_relax1") return Profile::AllRelax1;
    if (s == "all_relax2") return Profile::AllRelax2;
    if (s == "custom") return Profile::Custom;
    throw std::invalid_argument("unknown --param-profile item " + s);
}

std::string ProfileName(Profile p)
{
    switch (p) {
        case Profile::Baseline: return "baseline";
        case Profile::FinalRelax1: return "final_relax1";
        case Profile::FinalRelax2: return "final_relax2";
        case Profile::FinalRelax3: return "final_relax3";
        case Profile::ConversionRelax1: return "conversion_relax1";
        case Profile::ConversionRelax2: return "conversion_relax2";
        case Profile::BitExtractRelax1: return "bitextract_relax1";
        case Profile::AllRelax1: return "all_relax1";
        case Profile::AllRelax2: return "all_relax2";
        case Profile::Custom: return "custom";
    }
    return "unknown";
}

struct ProfileInfo {
    std::string changed_stage;
    std::string changed_params;
    bool security_params_changed = false;
    bool correctness_params_changed = false;
    std::string expected_speed_effect;
};

ProfileInfo GetProfileInfo(Profile p)
{
    switch (p) {
        case Profile::Baseline:
            return {"none", "baseline lvl1param l=2 Bgbit=8", false, false,
                    "reference"};
        case Profile::FinalRelax1:
            return {"final", "final PBS target l=1 Bgbit=10", false, true,
                    "fewer final external products"};
        case Profile::FinalRelax2:
            return {"final", "final PBS target l=1 Bgbit=8", false, true,
                    "fewer final external products, more approximation error"};
        case Profile::FinalRelax3:
            return {"final", "final PBS target l=1 Bgbit=6", false, true,
                    "aggressive final external-product relaxation"};
        case Profile::ConversionRelax1:
            return {"conversion", "conversion PBS target l=1 Bgbit=10", false,
                    true, "fewer conversion external products"};
        case Profile::ConversionRelax2:
            return {"conversion", "conversion PBS target l=1 Bgbit=8", false,
                    true, "fewer conversion external products, lower precision"};
        case Profile::BitExtractRelax1:
            return {"bitextract", "BitExtract PBS target l=1 Bgbit=10", false,
                    true, "fewer BitExtract external products"};
        case Profile::AllRelax1:
            return {"final+conversion",
                    "final and conversion PBS target l=1 Bgbit=10", false,
                    true, "two post-clear stages use fewer external products"};
        case Profile::AllRelax2:
            return {"bitextract+conversion+final",
                    "all PBS targets relaxed; bitextract l=1 b10, conversion/final l=1 b8",
                    false, true, "ablation only; all stages use fewer levels"};
        case Profile::Custom:
            return {"none", "custom is parsed but mapped to baseline in this probe",
                    false, false, "none"};
    }
    return {};
}

template <class BitTarget, class ConvTarget, class FinalTarget>
std::string ParamValues()
{
    std::ostringstream out;
    out << "bitextract=" << TargetParamName<BitTarget>()
        << ";conversion=" << TargetParamName<ConvTarget>()
        << ";final=" << TargetParamName<FinalTarget>()
        << ";iks=lvl10param(t=" << Lvl10::t
        << ",basebit=" << Lvl10::basebit << ")";
    return out.str();
}

std::string ParamValuesForProfile(Profile p)
{
    switch (p) {
        case Profile::Baseline:
            return ParamValues<Lvl1, Lvl1, Lvl1>();
        case Profile::FinalRelax1:
            return ParamValues<Lvl1, Lvl1, Relax1>();
        case Profile::FinalRelax2:
            return ParamValues<Lvl1, Lvl1, Relax2>();
        case Profile::FinalRelax3:
            return ParamValues<Lvl1, Lvl1, Relax3>();
        case Profile::ConversionRelax1:
            return ParamValues<Lvl1, Relax1, Lvl1>();
        case Profile::ConversionRelax2:
            return ParamValues<Lvl1, Relax2, Lvl1>();
        case Profile::BitExtractRelax1:
            return ParamValues<Relax1, Lvl1, Lvl1>();
        case Profile::AllRelax1:
            return ParamValues<Lvl1, Relax1, Relax1>();
        case Profile::AllRelax2:
            return ParamValues<Relax1, Relax2, Relax2>();
        case Profile::Custom:
            return ParamValues<Lvl1, Lvl1, Lvl1>();
    }
    return "";
}

uint32_t PeriodCellsForMode(const std::string &mode, const ClearBitSpec &spec)
{
    if (mode == "exact") return spec.exact_period_cells;
    if (mode == "2x") return spec.exact_period_cells * 2;
    if (mode == "4x") return spec.exact_period_cells * 4;
    if (mode == "8x") return spec.exact_period_cells * 8;
    if (mode == "full") return 1U << spec.p;
    throw std::invalid_argument("unknown --period-mode " + mode);
}

struct Options {
    uint32_t p = kFixedP;
    uint32_t bit_pos_lsb1 = kFixedBitPosLsb1;
    std::vector<std::string> modes = {"margin-audit"};
    std::string pipeline = "conversion-pbs";
    std::string input = "controlled-zero-noise";
    std::string policy = "tfhepp_v10_poly";
    std::vector<std::string> period_modes = {"exact"};
    int preoffset_denominator = 8;
    std::vector<int> preoffset_r_list = {0};
    std::string param_scope = "all";
    std::vector<Profile> profiles = {Profile::Baseline};
    std::string m_generator = "all";
    int boundary_radius = 16;
    uint64_t seed_start = 0;
    uint32_t seeds = 1;
    bool summary_only = false;
    bool verbose = false;
    std::string jsonl_path;
    std::string summary_path;
};

Options ParseOptions(int argc, char **argv)
{
    Options opts;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need = [&](const char *name) -> const char * {
            if (i + 1 >= argc)
                throw std::invalid_argument(std::string("missing value for ") +
                                            name);
            return argv[++i];
        };
        if (arg == "--p")
            opts.p = static_cast<uint32_t>(std::stoul(need("--p")));
        else if (arg == "--bit-pos-lsb1")
            opts.bit_pos_lsb1 =
                static_cast<uint32_t>(std::stoul(need("--bit-pos-lsb1")));
        else if (arg == "--mode")
            opts.modes = Split(need("--mode"));
        else if (arg == "--pipeline")
            opts.pipeline = need("--pipeline");
        else if (arg == "--input")
            opts.input = need("--input");
        else if (arg == "--policy")
            opts.policy = need("--policy");
        else if (arg == "--period-mode")
            opts.period_modes = Split(need("--period-mode"));
        else if (arg == "--preoffset-denominator")
            opts.preoffset_denominator =
                std::stoi(need("--preoffset-denominator"));
        else if (arg == "--preoffset-r-list")
            opts.preoffset_r_list = SplitInts(need("--preoffset-r-list"));
        else if (arg == "--param-scope")
            opts.param_scope = need("--param-scope");
        else if (arg == "--param-profile") {
            opts.profiles.clear();
            for (const std::string &name : Split(need("--param-profile")))
                opts.profiles.push_back(ParseProfile(name));
        }
        else if (arg == "--m-generator")
            opts.m_generator = need("--m-generator");
        else if (arg == "--boundary-radius")
            opts.boundary_radius = std::stoi(need("--boundary-radius"));
        else if (arg == "--seed-start")
            opts.seed_start = std::stoull(need("--seed-start"));
        else if (arg == "--seeds")
            opts.seeds = static_cast<uint32_t>(std::stoul(need("--seeds")));
        else if (arg == "--summary-only")
            opts.summary_only = true;
        else if (arg == "--verbose")
            opts.verbose = true;
        else if (arg == "--jsonl")
            opts.jsonl_path = need("--jsonl");
        else if (arg == "--summary")
            opts.summary_path = need("--summary");
        else
            throw std::invalid_argument("unknown option: " + arg);
    }
    if (opts.preoffset_denominator <= 0)
        throw std::invalid_argument("--preoffset-denominator must be positive");
    if (opts.seeds == 0) opts.seeds = 1;
    if (opts.preoffset_r_list.empty()) opts.preoffset_r_list.push_back(0);
    if (opts.period_modes.empty()) opts.period_modes.push_back("exact");
    if (opts.profiles.empty()) opts.profiles.push_back(Profile::Baseline);
    return opts;
}

std::vector<uint32_t> GenerateMessages(const Options &opts,
                                       const ClearBitSpec &spec)
{
    std::set<uint32_t> ms;
    const uint32_t modulus = 1U << spec.p;
    auto add_window = [&](int center) {
        for (int d = -opts.boundary_radius; d <= opts.boundary_radius; ++d) {
            const int x = center + d;
            if (x >= 0 && x < static_cast<int>(modulus))
                ms.insert(static_cast<uint32_t>(x));
        }
    };

    if (opts.m_generator == "all") {
        for (uint32_t m = 0; m < modulus; ++m) ms.insert(m);
    }
    else if (opts.m_generator == "bit-boundaries") {
        for (uint32_t b = 0; b <= modulus; b += spec.W) add_window(b);
        add_window(1U << (spec.p - 1));
    }
    else if (opts.m_generator == "final-boundary") {
        add_window(0);
        add_window(1U << (spec.p - 1));
        add_window(modulus - 1);
    }
    else if (opts.m_generator == "random") {
        std::mt19937_64 rng(opts.seed_start);
        std::uniform_int_distribution<uint32_t> dist(0, modulus - 1);
        const uint32_t count = std::min<uint32_t>(modulus, 4U * opts.seeds + 64);
        for (uint32_t i = 0; i < count; ++i) ms.insert(dist(rng));
    }
    else {
        throw std::invalid_argument("unknown --m-generator " + opts.m_generator);
    }
    return {ms.begin(), ms.end()};
}

struct CaseResult {
    uint32_t p = kFixedP;
    uint32_t bit_pos_lsb1 = kFixedBitPosLsb1;
    uint32_t j_lsb0 = 4;
    uint32_t k_msb0 = 4;
    uint32_t W = 16;
    std::string Delta_hex;
    std::string Delta_prime_hex;
    std::string G_hex;
    std::string M_gap_hex;
    std::string pipeline;
    std::string input;
    std::string period_mode;
    uint32_t period_cells = 32;
    std::string param_scope;
    std::string param_profile;
    std::string param_values;
    int preoffset_denominator = 8;
    int preoffset_r = 0;
    std::string preoffset_hex;
    uint32_t m = 0;
    uint64_t seed = 0;
    uint32_t exact_bit = 0;
    uint32_t bit_actual = 0;
    int bit_selected_c = -1;
    int selected_minus_nominal = 0;
    bool conversion_pass = true;
    std::string guard_expected_hex;
    std::string guard_actual_hex;
    std::string guard_error_hex;
    uint32_t clear_cell_expected = 0;
    uint32_t clear_cell_actual = 0;
    int final_selected_c = -1;
    uint32_t expected_label = 0;
    uint32_t actual_label = 0;
    bool pass = false;
    std::string failure_class = "none";
    double elapsed_ms = 0.0;
    uint64_t skipped_cmux_count = 0;
    uint64_t cmux_total = 0;
};

std::string CaseJson(const CaseResult &r)
{
    std::ostringstream out;
    out << "{\"type\":\"case\""
        << ",\"p\":" << r.p
        << ",\"bit_pos_lsb1\":" << r.bit_pos_lsb1
        << ",\"j_lsb0\":" << r.j_lsb0
        << ",\"k_msb0\":" << r.k_msb0
        << ",\"W\":" << r.W
        << ",\"Delta_hex\":\"" << r.Delta_hex << "\""
        << ",\"Delta_prime_hex\":\"" << r.Delta_prime_hex << "\""
        << ",\"G_hex\":\"" << r.G_hex << "\""
        << ",\"M_gap_hex\":\"" << r.M_gap_hex << "\""
        << ",\"pipeline\":\"" << JsonEscape(r.pipeline) << "\""
        << ",\"input\":\"" << JsonEscape(r.input) << "\""
        << ",\"period_mode\":\"" << JsonEscape(r.period_mode) << "\""
        << ",\"period_cells\":" << r.period_cells
        << ",\"param_scope\":\"" << JsonEscape(r.param_scope) << "\""
        << ",\"param_profile\":\"" << JsonEscape(r.param_profile) << "\""
        << ",\"param_values\":\"" << JsonEscape(r.param_values) << "\""
        << ",\"preoffset_denominator\":" << r.preoffset_denominator
        << ",\"preoffset_r\":" << r.preoffset_r
        << ",\"preoffset_hex\":\"" << r.preoffset_hex << "\""
        << ",\"m\":" << r.m
        << ",\"seed\":" << r.seed
        << ",\"exact_bit\":" << r.exact_bit
        << ",\"bit_actual\":" << r.bit_actual
        << ",\"bit_selected_c\":" << r.bit_selected_c
        << ",\"selected_minus_nominal\":" << r.selected_minus_nominal
        << ",\"conversion_pass\":" << (r.conversion_pass ? "true" : "false")
        << ",\"guard_expected_hex\":\"" << r.guard_expected_hex << "\""
        << ",\"guard_actual_hex\":\"" << r.guard_actual_hex << "\""
        << ",\"guard_error_hex\":\"" << r.guard_error_hex << "\""
        << ",\"clear_cell_expected\":" << r.clear_cell_expected
        << ",\"clear_cell_actual\":" << r.clear_cell_actual
        << ",\"final_selected_c\":" << r.final_selected_c
        << ",\"expected_label\":" << r.expected_label
        << ",\"actual_label\":" << r.actual_label
        << ",\"pass\":" << (r.pass ? "true" : "false")
        << ",\"failure_class\":\"" << JsonEscape(r.failure_class) << "\""
        << ",\"elapsed_ms\":" << std::fixed << std::setprecision(6)
        << r.elapsed_ms
        << ",\"skipped_cmux_count\":" << r.skipped_cmux_count
        << ",\"cmux_total\":" << r.cmux_total << "}";
    return out.str();
}

int SignedCellDelta(uint32_t selected, uint32_t nominal, uint32_t modulus)
{
    int diff = static_cast<int>(selected) - static_cast<int>(nominal);
    const int half = static_cast<int>(modulus / 2);
    if (diff > half) diff -= static_cast<int>(modulus);
    if (diff < -half) diff += static_cast<int>(modulus);
    return diff;
}

template <class BitTarget, class ConvTarget, class FinalTarget>
CaseResult RunCaseWithTargets(Profile profile, const Options &opts,
                              const MarginAudit &audit,
                              const KeyMaterial &keys, uint64_t seed,
                              const std::string &period_mode,
                              uint32_t period_cells, int preoffset_r,
                              uint32_t m)
{
    const ClearBitSpec &spec = audit.spec;
    CaseResult r;
    r.p = spec.p;
    r.bit_pos_lsb1 = spec.bit_pos_lsb1;
    r.j_lsb0 = spec.j_lsb0;
    r.k_msb0 = spec.k_msb0;
    r.W = spec.W;
    r.Delta_hex = Hex32(spec.delta);
    r.Delta_prime_hex = Hex32(spec.delta_prime);
    r.G_hex = Hex32(spec.G);
    r.M_gap_hex = Hex32(spec.margin_gap);
    r.pipeline = opts.pipeline;
    r.input = opts.input;
    r.period_mode = period_mode;
    r.period_cells = period_cells;
    r.param_scope = opts.param_scope;
    r.param_profile = ProfileName(profile);
    r.param_values = ParamValues<BitTarget, ConvTarget, FinalTarget>();
    r.preoffset_denominator = opts.preoffset_denominator;
    r.preoffset_r = preoffset_r;
    const int64_t preoffset =
        (static_cast<int64_t>(preoffset_r) * static_cast<int64_t>(spec.delta)) /
        opts.preoffset_denominator;
    r.preoffset_hex = HexSigned(preoffset);
    r.m = m;
    r.seed = seed;

    const uint32_t bit = (m >> spec.j_lsb0) & 1U;
    const uint32_t m_clear = m - bit * spec.W;
    r.exact_bit = bit;
    r.expected_label = m >= (1U << (spec.p - 1)) ? 1U : 0U;
    r.clear_cell_expected = m_clear;
    r.guard_expected_hex = Hex32(bit ? spec.G : Lvl1::T(0));

    const Lvl1::T phase_m = static_cast<Lvl1::T>(m) * spec.delta;
    const Lvl1::T phase_clear = static_cast<Lvl1::T>(m_clear) * spec.delta;
    const Lvl1::T phase_bit =
        bit ? (Lvl1::T(1) << (std::numeric_limits<Lvl1::T>::digits - 1))
            : Lvl1::T(0);

    TLWELvl1 original =
        MakeInputCipher(opts.input, phase_m, *keys.sk);
    TLWELvl1 clear_ct{};
    TLWELvl1 guard{};
    PbsTrace bit_trace;
    PbsTrace conv_trace;
    PbsTrace final_trace;

    const auto t0 = std::chrono::steady_clock::now();
    if (opts.pipeline == "oracle-clear") {
        clear_ct = MakeInputCipher(opts.input, phase_clear, *keys.sk);
        r.bit_actual = bit;
        r.conversion_pass = true;
        r.guard_actual_hex = r.guard_expected_hex;
        r.guard_error_hex = Hex32(0);
    }
    else {
        TLWELvl1 bit_ct{};
        if (opts.pipeline == "oracle-bit-conversion") {
            bit_ct = MakeInputCipher(opts.input, phase_bit, *keys.sk);
            r.bit_actual = DecodeQHalfBit(bit_ct, *keys.sk);
        }
        else if (opts.pipeline == "conversion-pbs" ||
                 opts.pipeline == "full") {
            TLWELvl1 bit_input = original;
            TLWEAddBodySigned(bit_input, preoffset);
            BitExtractQHalfPBS<BitTarget>(bit_ct, bit_input, spec, keys,
                                          period_cells, &bit_trace);
            r.bit_actual = DecodeQHalfBit(bit_ct, *keys.sk);
            r.bit_selected_c = static_cast<int>(bit_trace.selected_cell);
            r.selected_minus_nominal =
                SignedCellDelta(bit_trace.selected_cell, m, 1U << spec.p);
            r.skipped_cmux_count += bit_trace.skip_count;
            r.cmux_total += bit_trace.skip_count + bit_trace.cmux_count;
        }
        else {
            throw std::invalid_argument("unknown --pipeline " + opts.pipeline);
        }

        BoolQHalfToValueCenteredPBS<ConvTarget>(guard, bit_ct, spec, keys,
                                                &conv_trace);
        const Lvl1::T guard_phase = Phase(guard, *keys.sk);
        const Lvl1::T guard_expected = bit ? spec.G : Lvl1::T(0);
        r.guard_actual_hex = Hex32(guard_phase);
        r.guard_error_hex =
            HexSigned(SignedTorusDiff(guard_phase, guard_expected));
        r.conversion_pass =
            DecodeGuardBit(guard, spec, *keys.sk) == bit;
        r.skipped_cmux_count += conv_trace.skip_count;
        r.cmux_total += conv_trace.skip_count + conv_trace.cmux_count;

        clear_ct = original;
        TLWESubAssign(clear_ct, guard);
    }

    r.clear_cell_actual = PhaseToCell(Phase(clear_ct, *keys.sk), spec);

    TLWELvl1 final_ct{};
    WeightedApproxFinalPBS<FinalTarget>(final_ct, clear_ct, audit, keys,
                                        &final_trace);
    r.actual_label = DecodeLabel(final_ct, *keys.sk);
    r.final_selected_c = static_cast<int>(final_trace.selected_cell);
    r.skipped_cmux_count += final_trace.skip_count;
    r.cmux_total += final_trace.skip_count + final_trace.cmux_count;
    const auto t1 = std::chrono::steady_clock::now();
    r.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    const bool bit_ok = r.bit_actual == r.exact_bit;
    const bool final_ok = r.actual_label == r.expected_label;
    r.pass = bit_ok && r.conversion_pass && final_ok;
    if (r.pass)
        r.failure_class = "none";
    else if (!bit_ok && final_ok)
        r.failure_class = "benign_bit_error";
    else if (!bit_ok)
        r.failure_class = "fatal_bit_error";
    else if (!r.conversion_pass && final_ok)
        r.failure_class = "conversion_error_absorbed_by_final";
    else if (!r.conversion_pass)
        r.failure_class = "final_failure_after_conversion_error";
    else
        r.failure_class = "final_failure_given_correct_guard";
    return r;
}

CaseResult RunCase(Profile profile, const Options &opts,
                   const MarginAudit &audit, const KeyMaterial &keys,
                   uint64_t seed, const std::string &period_mode,
                   uint32_t period_cells, int preoffset_r, uint32_t m)
{
    switch (profile) {
        case Profile::Baseline:
            return RunCaseWithTargets<Lvl1, Lvl1, Lvl1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::FinalRelax1:
            return RunCaseWithTargets<Lvl1, Lvl1, Relax1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::FinalRelax2:
            return RunCaseWithTargets<Lvl1, Lvl1, Relax2>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::FinalRelax3:
            return RunCaseWithTargets<Lvl1, Lvl1, Relax3>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::ConversionRelax1:
            return RunCaseWithTargets<Lvl1, Relax1, Lvl1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::ConversionRelax2:
            return RunCaseWithTargets<Lvl1, Relax2, Lvl1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::BitExtractRelax1:
            return RunCaseWithTargets<Relax1, Lvl1, Lvl1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::AllRelax1:
            return RunCaseWithTargets<Lvl1, Relax1, Relax1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::AllRelax2:
            return RunCaseWithTargets<Relax1, Relax2, Relax2>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
        case Profile::Custom:
            return RunCaseWithTargets<Lvl1, Lvl1, Lvl1>(
                profile, opts, audit, keys, seed, period_mode, period_cells,
                preoffset_r, m);
    }
    throw std::invalid_argument("unreachable profile");
}

struct AggKey {
    std::string pipeline;
    std::string input;
    std::string period_mode;
    std::string profile;
    int preoffset_r = 0;

    bool operator<(const AggKey &o) const
    {
        return std::tie(pipeline, input, period_mode, profile, preoffset_r) <
               std::tie(o.pipeline, o.input, o.period_mode, o.profile,
                        o.preoffset_r);
    }
};

struct Agg {
    uint64_t cases = 0;
    uint64_t failures = 0;
    uint64_t bit_failures = 0;
    uint64_t conversion_failures = 0;
    uint64_t final_failures_given_correct_guard = 0;
    uint64_t benign_bit_errors = 0;
    uint64_t fatal_bit_errors = 0;
    uint64_t skipped_cmux_count = 0;
    uint64_t cmux_total = 0;
    std::vector<double> timings;
    std::map<uint32_t, uint64_t> failures_by_m;
    std::map<uint64_t, uint64_t> failures_by_seed;
    int64_t max_abs_guard_error = 0;
    uint64_t fixed_vs_baseline = 0;
    uint64_t introduced_vs_baseline = 0;
};

std::pair<double, double> Wilson95(uint64_t failures, uint64_t cases)
{
    if (cases == 0) return {0.0, 0.0};
    const double n = static_cast<double>(cases);
    const double phat = static_cast<double>(failures) / n;
    constexpr double z = 1.959963984540054;
    const double denom = 1.0 + z * z / n;
    const double center = (phat + z * z / (2.0 * n)) / denom;
    const double half =
        z * std::sqrt((phat * (1.0 - phat) + z * z / (4.0 * n)) / n) /
        denom;
    return {std::max(0.0, center - half), std::min(1.0, center + half)};
}

double Percentile(std::vector<double> xs, double p)
{
    if (xs.empty()) return 0.0;
    std::sort(xs.begin(), xs.end());
    const double pos = p * static_cast<double>(xs.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    if (lo == hi) return xs[lo];
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - static_cast<double>(lo));
}

std::map<AggKey, Agg> AggregateResults(const std::vector<CaseResult> &results)
{
    std::map<AggKey, Agg> aggs;
    std::map<std::tuple<std::string, std::string, std::string, int, uint64_t,
                        uint32_t>,
             bool>
        baseline_pass;
    for (const CaseResult &r : results) {
        if (r.param_profile == "baseline") {
            baseline_pass[{r.pipeline, r.input, r.period_mode, r.preoffset_r,
                           r.seed, r.m}] = r.pass;
        }
    }

    for (const CaseResult &r : results) {
        AggKey key{r.pipeline, r.input, r.period_mode, r.param_profile,
                   r.preoffset_r};
        Agg &a = aggs[key];
        a.cases++;
        if (!r.pass) {
            a.failures++;
            a.failures_by_m[r.m]++;
            a.failures_by_seed[r.seed]++;
        }
        if (r.bit_actual != r.exact_bit) a.bit_failures++;
        if (!r.conversion_pass) a.conversion_failures++;
        if (r.failure_class == "final_failure_given_correct_guard")
            a.final_failures_given_correct_guard++;
        if (r.failure_class == "benign_bit_error") a.benign_bit_errors++;
        if (r.failure_class == "fatal_bit_error") a.fatal_bit_errors++;
        a.skipped_cmux_count += r.skipped_cmux_count;
        a.cmux_total += r.cmux_total;
        a.timings.push_back(r.elapsed_ms);
        a.max_abs_guard_error = std::max<int64_t>(
            a.max_abs_guard_error,
            std::llabs(static_cast<long long>(
                std::stoll(r.guard_error_hex.substr(2), nullptr, 16))));

        if (r.param_profile != "baseline") {
            auto it = baseline_pass.find({r.pipeline, r.input, r.period_mode,
                                          r.preoffset_r, r.seed, r.m});
            if (it != baseline_pass.end()) {
                if (!it->second && r.pass) a.fixed_vs_baseline++;
                if (it->second && !r.pass) a.introduced_vs_baseline++;
            }
        }
    }
    return aggs;
}

std::string ShortMap(const std::map<uint32_t, uint64_t> &m, size_t limit = 12)
{
    if (m.empty()) return "{}";
    std::ostringstream out;
    out << "{";
    size_t n = 0;
    for (const auto &[k, v] : m) {
        if (n++ > 0) out << " ";
        out << k << ":" << v;
        if (n >= limit && n < m.size()) {
            out << " ...";
            break;
        }
    }
    out << "}";
    return out.str();
}

std::string RenderSummary(const Options &opts, const MarginAudit &audit,
                          const std::vector<CaseResult> &results)
{
    (void) opts;
    std::ostringstream out;
    const ClearBitSpec &spec = audit.spec;
    out << "# Clear Bit5 Parameter Relaxation Probe\n\n";
    out << "## Experiment Purpose\n\n";
    out << "A 9-bit experiment that clears LSB-based bit 5 only and checks "
           "whether the post-clear effective Delta_prime margin can be spent "
           "on final PBS parameter relaxation. Results are probabilistic "
           "ablation data, not certified correctness.\n\n";

    out << "## Bit Definition\n\n";
    out << "- p: " << spec.p << "\n";
    out << "- bit_pos_lsb1: " << spec.bit_pos_lsb1 << "\n";
    out << "- j_lsb0: " << spec.j_lsb0 << "\n";
    out << "- k_msb0: " << spec.k_msb0 << "\n";
    out << "- W: " << spec.W << "\n";
    out << "- exact_period_cells: " << spec.exact_period_cells << "\n";
    out << "- Delta_hex: " << Hex32(spec.delta) << "\n";
    out << "- Delta_prime_hex: " << Hex32(spec.delta_prime) << "\n";
    out << "- G_hex: " << Hex32(spec.G) << "\n";
    out << "- M_gap_hex: " << Hex32(spec.margin_gap) << "\n\n";

    out << "## Pipeline\n\n";
    out << "BitExtractQHalf -> BoolQHalfToValueCenteredPBS -> clear -> "
           "WeightedApproxFinalPBS. BitExtract keeps 0/Q2 output. Conversion "
           "uses qhalf - Q/4, centered PBS to +/-G/2, then post-adds +G/2.\n\n";

    out << "## Parameter Profiles\n\n";
    out << "| profile | changed stage | changed params | security params changed | "
           "correctness params changed | expected speed effect |\n";
    out << "|---|---|---|---|---|---|\n";
    for (Profile p : {Profile::Baseline, Profile::FinalRelax1,
                      Profile::FinalRelax2, Profile::FinalRelax3,
                      Profile::ConversionRelax1, Profile::ConversionRelax2,
                      Profile::BitExtractRelax1, Profile::AllRelax1,
                      Profile::AllRelax2}) {
        const ProfileInfo info = GetProfileInfo(p);
        out << "| " << ProfileName(p) << " | " << info.changed_stage << " | "
            << info.changed_params << " | "
            << (info.security_params_changed ? "yes" : "no") << " | "
            << (info.correctness_params_changed ? "yes" : "no") << " | "
            << info.expected_speed_effect << " |\n";
    }
    out << "\n";

    out << "## Margin Audit\n\n";
    out << "| reachable_count | unreachable_count | ambiguous_cell_count | "
           "min_opposite_label_distance_cells | expected_min_gap_cells | "
           "expected_margin_cells | negacyclic_audit_pass | pass |\n";
    out << "|---:|---:|---:|---:|---:|---:|---|---|\n";
    out << "| " << audit.reachable_count << " | " << audit.unreachable_count
        << " | " << audit.ambiguous_cell_count << " | "
        << audit.min_opposite_label_distance_cells << " | " << (spec.W + 1)
        << " | 8.5 | " << (audit.negacyclic_audit_pass ? "yes" : "no")
        << " | " << (audit.safe_coarse_margin ? "pass" : "fail") << " |\n\n";
    out << "Boundary labels: c239=" << audit.final_labels[239]
        << ", c240=" << audit.final_labels[240]
        << ", c255=" << audit.final_labels[255]
        << ", c256=" << audit.final_labels[256]
        << ", c257=" << audit.final_labels[257] << ".\n\n";

    if (results.empty()) {
        out << "## Result Tables\n\n";
        out << "No PBS pipeline cases were run in this invocation.\n\n";
    }
    else {
        const auto aggs = AggregateResults(results);
        out << "## Result Tables\n\n";
        out << "| pipeline | input | period | preoffset_r | profile | cases | "
               "failures | failure_rate | Wilson_95_CI | bit_failures | "
               "conversion_failures | final_failures_given_correct_guard | "
               "benign_bit_errors | fatal_bit_errors | fixed_vs_baseline | "
               "introduced_vs_baseline | net_improvement | avg_ms | p50_ms | "
               "p95_ms | skipped_cmux | cmux_total | pruning_ratio | "
               "failures_by_m |\n";
        out << "|---|---|---|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n";
        for (const auto &[key, agg] : aggs) {
            const double failure_rate =
                agg.cases == 0 ? 0.0
                               : static_cast<double>(agg.failures) /
                                     static_cast<double>(agg.cases);
            const auto ci = Wilson95(agg.failures, agg.cases);
            const double avg =
                agg.timings.empty()
                    ? 0.0
                    : std::accumulate(agg.timings.begin(), agg.timings.end(),
                                      0.0) /
                          static_cast<double>(agg.timings.size());
            const double p50 = Percentile(agg.timings, 0.50);
            const double p95 = Percentile(agg.timings, 0.95);
            const double pruning_ratio =
                agg.cmux_total == 0
                    ? 0.0
                    : static_cast<double>(agg.skipped_cmux_count) /
                          static_cast<double>(agg.cmux_total);
            out << "| " << key.pipeline << " | " << key.input << " | "
                << key.period_mode << " | " << key.preoffset_r << " | "
                << key.profile << " | " << agg.cases << " | "
                << agg.failures << " | " << std::fixed << std::setprecision(6)
                << failure_rate << " | [" << ci.first << "," << ci.second
                << "] | " << agg.bit_failures << " | "
                << agg.conversion_failures << " | "
                << agg.final_failures_given_correct_guard << " | "
                << agg.benign_bit_errors << " | " << agg.fatal_bit_errors
                << " | " << agg.fixed_vs_baseline << " | "
                << agg.introduced_vs_baseline << " | "
                << static_cast<int64_t>(agg.fixed_vs_baseline) -
                       static_cast<int64_t>(agg.introduced_vs_baseline)
                << " | " << avg << " | " << p50 << " | " << p95 << " | "
                << agg.skipped_cmux_count << " | " << agg.cmux_total << " | "
                << pruning_ratio << " | " << ShortMap(agg.failures_by_m)
                << " |\n";
        }
        out << "\n";
    }

    out << "## Conclusion Template\n\n";
    out << "- Final PBS relaxation is a probabilistic candidate only if "
           "oracle-clear holdout and all-m smoke both beat baseline.\n";
    out << "- Conversion relaxation must be interpreted separately from final "
           "absorption of guard error.\n";
    out << "- BitExtract relaxation is expected to be boundary sensitive and is "
           "not justified by the post-clear Delta_prime margin.\n";
    out << "- Encrypted-noisy results, if absent, are a blocker for extending "
           "controlled-zero-noise conclusions to normal ciphertexts.\n";
    out << "- No result in this report is certified correctness.\n";
    return out.str();
}

void EnsureParentDir(const std::string &path)
{
    if (path.empty()) return;
    const std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
}

void WriteTextFile(const std::string &path, const std::string &text)
{
    if (path.empty()) return;
    EnsureParentDir(path);
    std::ofstream out(path);
    out << text;
}

void AppendJsonl(const std::string &path, const std::vector<CaseResult> &results)
{
    if (path.empty()) return;
    EnsureParentDir(path);
    std::ofstream out(path, std::ios::app);
    for (const CaseResult &r : results) out << CaseJson(r) << "\n";
}

std::string AuditJson(const MarginAudit &audit)
{
    const ClearBitSpec &s = audit.spec;
    std::ostringstream out;
    out << "{\"type\":\"margin_audit\""
        << ",\"p\":" << s.p
        << ",\"bit_pos_lsb1\":" << s.bit_pos_lsb1
        << ",\"j_lsb0\":" << s.j_lsb0
        << ",\"k_msb0\":" << s.k_msb0
        << ",\"W\":" << s.W
        << ",\"exact_period_cells\":" << s.exact_period_cells
        << ",\"Delta_hex\":\"" << Hex32(s.delta) << "\""
        << ",\"Delta_prime_hex\":\"" << Hex32(s.delta_prime) << "\""
        << ",\"G_hex\":\"" << Hex32(s.G) << "\""
        << ",\"gap_margin_hex\":\"" << Hex32(s.margin_gap) << "\""
        << ",\"reachable_count\":" << audit.reachable_count
        << ",\"unreachable_count\":" << audit.unreachable_count
        << ",\"ambiguous_cell_count\":" << audit.ambiguous_cell_count
        << ",\"min_opposite_label_distance_cells\":"
        << audit.min_opposite_label_distance_cells
        << ",\"expected_min_gap_cells\":" << (s.W + 1)
        << ",\"expected_margin_cells\":8.5"
        << ",\"final_boundary_cells\":256"
        << ",\"label_c239\":" << audit.final_labels[239]
        << ",\"label_c240\":" << audit.final_labels[240]
        << ",\"label_c255\":" << audit.final_labels[255]
        << ",\"label_c256\":" << audit.final_labels[256]
        << ",\"label_c257\":" << audit.final_labels[257]
        << ",\"negacyclic_audit_pass\":"
        << (audit.negacyclic_audit_pass ? "true" : "false")
        << ",\"safe_coarse_margin\":"
        << (audit.safe_coarse_margin ? "true" : "false") << "}";
    return out.str();
}

void WriteAuditJsonl(const std::string &path, const MarginAudit &audit)
{
    if (path.empty()) return;
    EnsureParentDir(path);
    std::ofstream out(path, std::ios::app);
    out << AuditJson(audit) << "\n";
}

std::vector<CaseResult> RunExperiment(const Options &opts,
                                      const MarginAudit &audit)
{
    std::vector<CaseResult> results;
    const std::vector<uint32_t> messages = GenerateMessages(opts, audit.spec);

    for (uint32_t si = 0; si < opts.seeds; ++si) {
        const uint64_t seed = opts.seed_start + si;
        KeyMaterial keys;
        keys.Generate(seed, opts.profiles);
        for (Profile profile : opts.profiles) {
            for (const std::string &period_mode : opts.period_modes) {
                const uint32_t period_cells =
                    PeriodCellsForMode(period_mode, audit.spec);
                for (int preoffset_r : opts.preoffset_r_list) {
                    for (uint32_t m : messages) {
                        CaseResult r = RunCase(profile, opts, audit, keys, seed,
                                               period_mode, period_cells,
                                               preoffset_r, m);
                        if (!opts.summary_only || opts.verbose)
                            std::cout << CaseJson(r) << "\n";
                        results.push_back(std::move(r));
                    }
                }
            }
        }
    }
    return results;
}

Options NormalizeModeOptions(Options opts)
{
    if (Contains(opts.modes, "oracle-final")) opts.pipeline = "oracle-clear";
    if (Contains(opts.modes, "conversion-sanity"))
        opts.pipeline = "oracle-bit-conversion";
    if (Contains(opts.modes, "all")) {
        opts.modes = {"margin-audit", "param-sweep"};
    }
    return opts;
}

int Main(int argc, char **argv)
{
    Options opts = NormalizeModeOptions(ParseOptions(argc, argv));
    const ClearBitSpec spec = MakeSpec(opts.p, opts.bit_pos_lsb1);
    const MarginAudit audit = BuildMarginAudit(spec);

    if (Contains(opts.modes, "margin-audit")) {
        std::cout << AuditJson(audit) << "\n";
        WriteAuditJsonl(opts.jsonl_path, audit);
    }

    const bool wants_pbs =
        Contains(opts.modes, "param-sweep") ||
        Contains(opts.modes, "full-pipeline") ||
        Contains(opts.modes, "oracle-final") ||
        Contains(opts.modes, "conversion-sanity");

    if (wants_pbs && !audit.safe_coarse_margin) {
        std::cerr << "margin-audit failed; refusing to run PBS pipeline\n";
        return 1;
    }

    std::vector<CaseResult> results;
    if (wants_pbs) {
        results = RunExperiment(opts, audit);
        AppendJsonl(opts.jsonl_path, results);
    }

    const std::string report = RenderSummary(opts, audit, results);
    if (!opts.summary_path.empty()) WriteTextFile(opts.summary_path, report);
    WriteTextFile("experimental/results/CLEAR_BIT5_PARAM_RELAX_REPORT.md",
                  report);

    if (opts.summary_only || !wants_pbs) std::cout << report;
    return 0;
}

} // namespace experimental::clear_bit5

int main(int argc, char **argv)
{
    try {
        return experimental::clear_bit5::Main(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << "\n";
        return 2;
    }
}
