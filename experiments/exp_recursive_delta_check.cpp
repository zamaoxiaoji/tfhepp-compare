#include "experiment_utils.hpp"

#include <array>

using namespace OursExperiments;

namespace {

template <class P>
typename P::T Phase(const TFHEpp::TLWE<P>& ct, const TFHEpp::SecretKey& sk) {
    return TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
}

template <typename T>
std::uint64_t DecodeArithmeticPhase(T phase, int p) {
    const int digits = std::numeric_limits<T>::digits;
    const T half_delta = T(1) << (digits - p - 1);
    const T rounded = phase + half_delta;
    const std::uint64_t mask = (std::uint64_t{1} << p) - 1;
    return (static_cast<std::uint64_t>(rounded >> (digits - p))) & mask;
}

template <typename T>
std::uint64_t TorusDistance(T a, T b) {
    const T ab = a - b;
    const T ba = b - a;
    return static_cast<std::uint64_t>(std::min(ab, ba));
}

struct Failure {
    std::string where;
    std::uint64_t a = 0;
    std::uint64_t b = 0;
    int expected = 0;
    int got = 0;
    std::string detail;
};

void RecordFailure(std::vector<Failure>& failures, const Failure& f) {
    if (failures.size() < 50) failures.push_back(f);
}

TFHEpp::TLWE<PIn> FinalGapMSBResidualAware(
    const TFHEpp::TLWE<PIn>& ct_work,
    int p_work,
    int k,
    MetaPBS2::ResidualInterval rho_work,
    OursRuntime& rt,
    MetaPBS2::BlindRotatePruneStats* stats = nullptr) {
    const int s = MetaPBS2::LSBIndexFromChapterBit(p_work, k) + 1;
    const auto weight_ct = MetaPBS2::ExtractWeightedLowWindowPBS_Lvl02<
        iksP_t, brP_logari>(
            ct_work, p_work, s,
            MetaPBS2::ArithmeticWeightForChapterBit<PIn::T>(p_work, k),
            *rt.iksk, *rt.bk_logari, stats);
    TFHEpp::TLWE<PIn> gap_ct{};
    for (std::size_t i = 0; i < gap_ct.size(); i++)
        gap_ct[i] = ct_work[i] - weight_ct[i];
    MetaPBS2::RecordLvl2ToLvl0GatePBSStats<iksP_t, brP_logari>(stats);
    if (stats) stats->pbs_count_final_msb++;
    return MetaPBS2::NaiveSignPBS_Lvl02<iksP_t, brP_logari>(
        gap_ct,
        MetaPBS2::FinalGapOffsetForResidual<PIn::T>(p_work, k, rho_work),
        *rt.iksk, *rt.bk_logari);
}

bool SupportedCurrentEncryptedPath(int p, const MetaPBS2::Algorithm1Config& cfg) {
    (void)cfg;
    return p <= 33;
}

void PrintUnsupported(int L, int p, const MetaPBS2::Algorithm1Config& cfg) {
    const int cfg_p = MetaPBS2::MessagePrecisionFromPowerOfTwoModulus(cfg.t);
    std::cout << "section,status,L,p_original,cfg_t,cfg_p,reason\n"
              << "support,unsupported," << L << ',' << p << ','
              << cfg.t << ',' << cfg_p
              << ",p_original_above_supported_33\n";
}

struct WorkState {
    TFHEpp::TLWE<PIn> ct_work{};
    int p_reduced = 0;
    int p_work = 0;
    std::uint64_t m_reduced = 0;
    std::uint64_t m_work = 0;
    MetaPBS2::ResidualInterval rho_work{};
};

WorkState ReduceOrZeroExtendToWork(
    const TFHEpp::TLWE<PIn>& ct,
    std::uint64_t m,
    int p,
    OursRuntime& rt,
    MetaPBS2::BlindRotatePruneStats* stats = nullptr,
    std::vector<MetaPBS2::RecursiveGapRoundOutput<PIn>>* rounds = nullptr) {
    const int p_work = NativePrecision(rt.cfg);
    if (p <= p_work) {
        return WorkState{
            .ct_work = ct,
            .p_reduced = p,
            .p_work = p_work,
            .m_reduced = m,
            .m_work = m << (p_work - p),
            .rho_work = MetaPBS2::ResidualInterval{.lo = 0.0, .hi = 0.0},
        };
    }
    const auto reduced = MetaPBS2::RecursiveGapReduceToNative<iksP_t, brP_logari>(
        ct, p, p_work, *rt.iksk, *rt.bk_logari, stats);
    if (rounds) *rounds = reduced.round_outputs;
    const auto sim = MetaPBS2::SimulateRecursiveGapReduce(
        m, p, p_work, DefaultK(p_work));
    const auto m_reduced = sim.m_work;
    return WorkState{
        .ct_work = reduced.ct_out,
        .p_reduced = reduced.p_out,
        .p_work = p_work,
        .m_reduced = m_reduced,
        .m_work = m_reduced,
        .rho_work = reduced.rho_out,
    };
}

bool CheckDeltaStateRows(int L, int trials, std::uint64_t seed,
                         OursRuntime& rt,
                         std::vector<Failure>& failures) {
    const int p = L + 1;
    const int p_work = WorkPrecisionForComparison(L, rt.cfg);
    std::mt19937_64 rng(seed ^ 0x44504c544143484bULL);
    std::uniform_int_distribution<std::uint64_t> dist(
        0, (std::uint64_t{1} << p) - 1);

    std::cout << "section,L,trial,round,p_in,p_out,delta_in,delta_out,"
                 "m_original,decoded_current,expected_current,ok\n";
    for (int i = 0; i < trials; i++) {
        const auto m = dist(rng);
        TFHEpp::TLWE<PIn> ct{};
        TFHEpp::tlweSymEncrypt<PIn>(
            ct, EncodeMessage<PIn>(m, p), PIn::α, rt.sk.key.get<PIn>());
        std::vector<MetaPBS2::RecursiveGapRoundOutput<PIn>> rounds;
        MetaPBS2::BlindRotatePruneStats stats{};
        const auto work =
            ReduceOrZeroExtendToWork(ct, m, p, rt, &stats, &rounds);
        const auto sim = MetaPBS2::SimulateRecursiveGapReduce(
            m, p, p_work, DefaultK(p_work));

        bool printed_round = false;
        for (const auto& round : rounds) {
            const auto delta_in =
                MetaPBS2::ReducerTorusScaleForPrecision<PIn::T>(round.p_in);
            const auto delta_out =
                MetaPBS2::ReducerTorusScaleForPrecision<PIn::T>(round.p_out);
            const auto decoded =
                DecodeArithmeticPhase(Phase<PIn>(round.ct_out, rt.sk), round.p_out);
            const auto& sim_round = sim.rounds.at(static_cast<std::size_t>(round.round));
            const auto expected = static_cast<std::uint64_t>(
                std::floor(static_cast<long double>(sim_round.q) +
                           static_cast<long double>(sim_round.rho_out) + 0.5L)) &
                ((std::uint64_t{1} << round.p_out) - 1);
            const bool ok = decoded == expected;
            printed_round = true;
            std::cout << "delta_state," << L << ',' << i
                      << ',' << round.round << ',' << round.p_in << ','
                      << round.p_out << ',' << delta_in << ',' << delta_out
                      << ',' << m << ',' << decoded << ',' << expected << ','
                      << (ok ? 1 : 0) << '\n';
            if (!ok) {
                RecordFailure(failures, Failure{
                    .where = "delta_state",
                    .a = m,
                    .expected = static_cast<int>(expected),
                    .got = static_cast<int>(decoded),
                    .detail = "precision reducer round did not decode at its output Delta"});
            }
        }

        const auto decoded_work =
            DecodeArithmeticPhase(Phase<PIn>(work.ct_work, rt.sk), p_work);
        const auto expected_work = static_cast<std::uint64_t>(
            std::floor(static_cast<long double>(sim.m_work) +
                       static_cast<long double>(sim.rho_work) + 0.5L)) &
            ((std::uint64_t{1} << p_work) - 1);
        const bool work_ok = decoded_work == expected_work;
        const int round_index =
            printed_round ? static_cast<int>(rounds.size()) : 0;
        std::cout << "delta_state," << L << ',' << i
                  << ',' << round_index << ',' << work.p_reduced << ','
                  << p_work << ','
                  << MetaPBS2::ReducerTorusScaleForPrecision<PIn::T>(work.p_reduced)
                  << ','
                  << MetaPBS2::ReducerTorusScaleForPrecision<PIn::T>(p_work)
                  << ',' << m << ',' << decoded_work << ',' << expected_work
                  << ',' << (work_ok ? 1 : 0) << '\n';
        if (!work_ok) {
            RecordFailure(failures, Failure{
                .where = "delta_state",
                .a = m,
                .expected = static_cast<int>(work.m_work),
                .got = static_cast<int>(decoded_work),
                .detail = "work ciphertext did not decode at native/current Delta"});
        }
    }
    return failures.empty();
}

bool CheckGapIntermediates(int L, int trials, std::uint64_t seed,
                           OursRuntime& rt,
                           std::vector<Failure>& failures) {
    const int p = L + 1;
    const int p_work = WorkPrecisionForComparison(L, rt.cfg);
    const int k = DefaultK(p_work);
    const auto params = MetaPBS2::MakeGapMSBRoundParams<PIn::T>(p, p_work, 0, k);
    const int lsb_k = MetaPBS2::LSBIndexFromChapterBit(p_work, k);
    const auto delta = params.current.delta;
    const auto half_delta = delta >> 1;
    std::mt19937_64 rng(seed ^ 0x474150494e544dULL);
    std::uniform_int_distribution<std::uint64_t> dist(
        0, (std::uint64_t{1} << p) - 1);

    std::cout << "section,L,trial,p_cur,k,w_k,m_original,m_cur,bit_k,"
                 "decoded_bit,expected_weight,decoded_weight,expected_gap,"
                 "decoded_gap,shifted_phase_ok,final_expected,final_got,ok,"
                 "pbs_count,cmux_total,cmux_skipped,first_round_period\n";
    for (int i = 0; i < trials; i++) {
        const auto m = dist(rng);
        TFHEpp::TLWE<PIn> ct{};
        TFHEpp::tlweSymEncrypt<PIn>(
            ct, EncodeMessage<PIn>(m, p), PIn::α, rt.sk.key.get<PIn>());

        MetaPBS2::BlindRotatePruneStats stats{};
        const auto work = ReduceOrZeroExtendToWork(ct, m, p, rt, &stats);
        auto bit_ct = MetaPBS2::BitExtractLowWindowBoolPruned<
            iksP_t, brP_logari>(
                work.ct_work, p_work, lsb_k + 1, *rt.iksk,
                *rt.bk_logari, &stats);
        const int decoded_bit = MetaPBS2::DecodeBinaryCout<PIn>(
            bit_ct, rt.sk.key.get<PIn>());
        const int bit = static_cast<int>((work.m_work >> lsb_k) & 1);

        auto weight_ct = MetaPBS2::BoolToWeightPBS_Lvl02<
            iksP_t, brP_logari>(
                bit_ct, params.A_k, *rt.iksk, *rt.bk_logari, &stats, false);
        TFHEpp::TLWE<PIn> gap_ct{};
        MetaPBS2::ClearChapterBitAssign<brP_meta>(gap_ct, work.ct_work, weight_ct);
        TFHEpp::TLWE<PIn> shifted_ct = gap_ct;
        const auto final_offset = MetaPBS2::FinalGapOffsetForResidual<PIn::T>(
            p_work, k, work.rho_work);
        shifted_ct[PIn::k * PIn::n] += final_offset;

        auto final_ct = FinalGapMSBResidualAware(
            work.ct_work, p_work, k, work.rho_work, rt, &stats);
        const int final_got = MetaPBS2::DecodeSignCout<PIn>(
            final_ct, rt.sk.key.get<PIn>());
        const int final_expected = static_cast<int>((m >> (p - 1)) & 1);

        const std::uint64_t expected_weight = bit ? params.w_k : 0;
        const std::uint64_t decoded_weight =
            DecodeArithmeticPhase(Phase<PIn>(weight_ct, rt.sk), p_work);
        const std::uint64_t expected_gap = work.m_work - expected_weight;
        const std::uint64_t decoded_gap =
            DecodeArithmeticPhase(Phase<PIn>(gap_ct, rt.sk), p_work);
        const auto expected_shifted_phase =
            static_cast<PIn::T>(expected_gap) * delta + final_offset;
        const bool shifted_ok =
            TorusDistance(Phase<PIn>(shifted_ct, rt.sk), expected_shifted_phase) <= half_delta;
        const bool ok = decoded_bit == bit &&
                        decoded_weight == expected_weight &&
                        decoded_gap == expected_gap &&
                        shifted_ok &&
                        final_got == final_expected;

        const auto counters = ToCounters(stats);
        std::cout << "gap_intermediate," << L << ',' << i << ','
                  << p_work << ',' << k << ',' << params.w_k << ','
                  << m << ',' << work.m_work << ',' << bit << ','
                  << decoded_bit << ',' << expected_weight << ','
                  << decoded_weight << ',' << expected_gap << ','
                  << decoded_gap << ',' << (shifted_ok ? 1 : 0) << ','
                  << final_expected << ',' << final_got << ','
                  << (ok ? 1 : 0) << ','
                  << counters.pbs_count << ',' << counters.total_cmux << ','
                  << counters.skipped_cmux << ','
                  << counters.first_round_slot_period << '\n';
        if (!ok) {
            RecordFailure(failures, Failure{
                .where = "gap_intermediate",
                .a = m,
                .expected = final_expected,
                .got = final_got,
                .detail = "bit/weight/gap/shifted/final invariant mismatch"});
        }
    }
    return failures.empty();
}

bool CheckBoundaryMSB(int L, int trials_per_case, std::uint64_t seed,
                      OursRuntime& rt,
                      std::vector<Failure>& failures) {
    const int p = L + 1;
    const int p_work = WorkPrecisionForComparison(L, rt.cfg);
    const int k = DefaultK(p_work);
    std::vector<std::pair<std::uint64_t, std::string>> cases;
    const std::int64_t threshold = std::int64_t{1} << (p - 1);
    for (std::int64_t d = -64; d <= 64; d++) {
        const auto v = threshold + d;
        if (v >= 0 && v < (std::int64_t{1} << p))
            cases.push_back({static_cast<std::uint64_t>(v),
                             d < 0 ? "lower_boundary" : "upper_boundary"});
    }
    for (auto v : {std::uint64_t{0}, std::uint64_t{1}, std::uint64_t{2},
                   (std::uint64_t{1} << p) - 3,
                   (std::uint64_t{1} << p) - 2,
                   (std::uint64_t{1} << p) - 1}) {
        cases.push_back({v, v < 3 ? "zero" : "max"});
    }

    std::cout << "section,L,m,category,trials,fail,failure_rate_percent\n";
    bool ok_all = true;
    for (const auto& [m, category] : cases) {
        int fail = 0;
        for (int t = 0; t < trials_per_case; t++) {
            TFHEpp::TLWE<PIn> ct{};
            TFHEpp::tlweSymEncrypt<PIn>(
                ct, EncodeMessage<PIn>(m, p), PIn::α, rt.sk.key.get<PIn>());
            const auto work = ReduceOrZeroExtendToWork(ct, m, p, rt);
            auto out = FinalGapMSBResidualAware(
                work.ct_work, p_work, k, work.rho_work, rt, nullptr);
            const int got = MetaPBS2::DecodeSignCout<PIn>(
                out, rt.sk.key.get<PIn>());
            const int expected = static_cast<int>((m >> (p - 1)) & 1);
            if (got != expected) {
                fail++;
                RecordFailure(failures, Failure{
                    .where = "boundary_msb",
                    .a = m,
                    .expected = expected,
                    .got = got,
                    .detail = category});
            }
        }
        const double rate = Percent(fail, trials_per_case);
        std::cout << "boundary_msb," << L << ',' << m << ',' << category << ','
                  << trials_per_case << ',' << fail << ','
                  << std::fixed << std::setprecision(4) << rate << '\n';
        if (rate > 0.1) ok_all = false;
    }
    (void)seed;
    return ok_all && failures.empty();
}

bool CheckUnsignedRegression(int L, int trials_per_case, std::uint64_t,
                             OursRuntime& rt,
                             std::vector<Failure>& failures) {
    const int p = L + 1;
    const int k = DefaultK(WorkPrecisionForComparison(L, rt.cfg));
    const std::uint64_t max_l = MaxValueForL(L);
    const std::uint64_t half = std::uint64_t{1} << (L - 1);
    const std::vector<std::pair<std::uint64_t, std::uint64_t>> cases = {
        {0, max_l},
        {max_l, 0},
        {0, 0},
        {max_l, max_l},
        {half - 1, half},
        {half, half - 1},
    };
    const std::array<const char*, 6> ops = {"lt", "le", "gt", "ge", "eq", "ne"};

    std::cout << "section,L,a,b,op,expected,trials,fail,failure_rate_percent\n";
    bool ok_all = true;
    for (const auto& [a, b] : cases) {
        for (const char* op_c : ops) {
            const std::string op(op_c);
            int fail = 0;
            const bool expected = ExpectedCompare(a, b, op);
            for (int t = 0; t < trials_per_case; t++) {
                auto result = RunOursCompare(rt, a, b, L, k, op);
                if (result.got != expected) {
                    fail++;
                    RecordFailure(failures, Failure{
                        .where = "unsigned_regression",
                        .a = a,
                        .b = b,
                        .expected = expected ? 1 : 0,
                        .got = result.got ? 1 : 0,
                        .detail = op});
                }
            }
            const double rate = Percent(fail, trials_per_case);
            std::cout << "unsigned_regression," << L << ',' << a << ',' << b
                      << ',' << op << ',' << (expected ? 1 : 0) << ','
                      << trials_per_case << ',' << fail << ','
                      << std::fixed << std::setprecision(4) << rate << '\n';
            if (fail > 0) ok_all = false;
        }
    }
    (void)p;
    return ok_all && failures.empty();
}

bool CheckNoiseSweep(int L, int trials_per_point, std::uint64_t seed,
                     OursRuntime& rt,
                     std::vector<Failure>& failures) {
    const int p = L + 1;
    const int p_work = WorkPrecisionForComparison(L, rt.cfg);
    const int k = DefaultK(p_work);
    const std::array<double, 5> multipliers = {0.5, 1.0, 1.5, 2.0, 3.0};
    std::mt19937_64 rng(seed ^ 0x4e4f495345ULL);
    std::uniform_int_distribution<std::uint64_t> dist(
        0, (std::uint64_t{1} << p) - 1);
    std::cout << "section,L,noise_multiplier,trials,fail,failure_rate_percent\n";
    bool ok_all = true;
    for (double mult : multipliers) {
        int fail = 0;
        for (int i = 0; i < trials_per_point; i++) {
            const auto m = dist(rng);
            TFHEpp::TLWE<PIn> ct{};
            TFHEpp::tlweSymEncrypt<PIn>(
                ct, EncodeMessage<PIn>(m, p), PIn::α * mult, rt.sk.key.get<PIn>());
            const auto work = ReduceOrZeroExtendToWork(ct, m, p, rt);
            auto out = FinalGapMSBResidualAware(
                work.ct_work, p_work, k, work.rho_work, rt, nullptr);
            const int got = MetaPBS2::DecodeSignCout<PIn>(
                out, rt.sk.key.get<PIn>());
            const int expected = static_cast<int>((m >> (p - 1)) & 1);
            if (got != expected) {
                fail++;
                RecordFailure(failures, Failure{
                    .where = "noise_sweep",
                    .a = work.m_work,
                    .expected = expected,
                    .got = got,
                    .detail = std::to_string(mult)});
            }
        }
        const double rate = Percent(fail, trials_per_point);
        std::cout << "noise_sweep," << L << ',' << mult << ','
                  << trials_per_point << ',' << fail << ','
                  << std::fixed << std::setprecision(4) << rate << '\n';
        if (rate > 0.1) ok_all = false;
    }
    return ok_all && failures.empty();
}

void PrintScheduleAblationRows(int L, int p, const MetaPBS2::Algorithm1Config& cfg) {
    const int p_work = WorkPrecisionForComparison(L, cfg);
    const int k = DefaultK(p_work);
    const auto params = MetaPBS2::MakeGapMSBRoundParams<PIn::T>(p, p_work, 0, k);
    std::cout << "section,L,schedule,p_final,k_final,Delta_final,trials,fail,"
                 "failure_rate_percent,latency_ms_avg,latency_ms_median,status\n";
    std::cout << "schedule_ablation," << L << ','
              << PrecisionScheduleString(L, cfg) << ','
              << p_work << ',' << k << ',' << params.current.delta
              << ",0,0,0.0000,0,0,implemented_default\n";
}

void PrintPerfCounterRow(int L, int p, std::uint64_t seed, OursRuntime& rt) {
    const int p_work = WorkPrecisionForComparison(L, rt.cfg);
    const int k = DefaultK(p_work);
    const auto params = MetaPBS2::MakeGapMSBRoundParams<PIn::T>(p, p_work, 0, k);
    std::mt19937_64 rng(seed ^ 0x50455246ULL);
    std::uniform_int_distribution<std::uint64_t> dist(
        0, (std::uint64_t{1} << p) - 1);
    TFHEpp::TLWE<PIn> ct{};
    MetaPBS2::BlindRotatePruneStats stats{};
    const auto m = dist(rng);
    TFHEpp::tlweSymEncrypt<PIn>(
        ct, EncodeMessage<PIn>(m, p), PIn::α, rt.sk.key.get<PIn>());
    const auto start = std::chrono::steady_clock::now();
    const auto work = ReduceOrZeroExtendToWork(ct, m, p, rt, &stats);
    auto out = FinalGapMSBResidualAware(
        work.ct_work, p_work, k, work.rho_work, rt, &stats);
    const double ms = MsSince(start);
    (void)out;
    const auto counters = ToCounters(stats);
    std::cout << "section,L,p_original,p_final,k_final,schedule,pbs_count_total,"
                 "pbs_count_recursive,pbs_count_bit_extract,pbs_count_bool_to_weight,"
                 "pbs_count_final_msb,key_switch_count_total,cmux_total,cmux_skipped,"
                 "pruning_ratio_percent,latency_ms,delta_current,offset\n";
    std::cout << "perf_counters," << L << ',' << p << ',' << work.p_reduced << ','
              << k << ',' << PrecisionScheduleString(L, rt.cfg) << ','
              << counters.pbs_count << ',' << counters.pbs_count_reducer << ','
              << counters.pbs_count_bit_extract << ','
              << counters.pbs_count_bool_to_weight << ','
              << counters.pbs_count_final_msb << ','
              << counters.key_switch_count << ','
              << counters.total_cmux << ',' << counters.skipped_cmux << ','
              << std::fixed << std::setprecision(4)
              << Percent(counters.skipped_cmux, counters.total_cmux) << ','
              << ms << ',' << params.current.delta << ',' << params.offset << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 11);
    const int trials = args.GetInt("--trials", 200000);
    const int boundary_trials = args.GetInt("--boundary-trials", 20000);
    const int noise_trials = args.GetInt("--noise-trials", 50000);
    const auto seed = args.GetU64("--seed", 123456789ULL);
    const bool verbose = args.Has("--verbose");

    const int p = L + 1;
    const auto cfg_probe = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = WorkPrecisionForComparison(L, cfg_probe);
    const int k = DefaultK(p_work);
    PrintExperimentHeader("recursive_delta_check", L, p, k, trials, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg_probe) << "\n"
              << "# boundary_trials: " << boundary_trials << "\n"
              << "# noise_trials: " << noise_trials << "\n"
              << "# note: Delta_prime must be produced by PBS re-encoding, "
                 "not decoder reinterpretation\n";

    if (!SupportedCurrentEncryptedPath(p, cfg_probe)) {
        PrintUnsupported(L, p, cfg_probe);
        std::cout << std::flush;
        std::cerr << "FAILED_UNSUPPORTED_DELTA_PRIME: current cfg.t="
                  << cfg_probe.t << " supports p="
                  << MetaPBS2::MessagePrecisionFromPowerOfTwoModulus(cfg_probe.t)
                  << ", but requested p=" << p
                  << ". No true recursive PBS re-encoding schedule is currently "
                     "implemented, so this test refuses to reinterpret Delta'.\n";
        return 2;
    }

    OursRuntime rt;
    std::vector<Failure> failures;
    bool ok = true;
    ok = CheckDeltaStateRows(L, trials, seed, rt, failures) && ok;
    ok = CheckGapIntermediates(L, trials, seed, rt, failures) && ok;
    ok = CheckBoundaryMSB(L, boundary_trials, seed, rt, failures) && ok;
    ok = CheckUnsignedRegression(L, boundary_trials, seed, rt, failures) && ok;
    ok = CheckNoiseSweep(L, noise_trials, seed, rt, failures) && ok;
    PrintScheduleAblationRows(L, p, rt.cfg);
    PrintPerfCounterRow(L, p, seed, rt);

    if (!failures.empty()) {
        std::cerr << "first_failures,where,a,b,expected,got,detail\n";
        for (const auto& f : failures) {
            std::cerr << "failure," << f.where << ',' << f.a << ',' << f.b
                      << ',' << f.expected << ',' << f.got << ','
                      << f.detail << '\n';
        }
    }
    if (verbose)
        std::cerr << "recursive_delta_check_status=" << (ok ? "ok" : "failed") << '\n';
    return ok ? 0 : 1;
}
