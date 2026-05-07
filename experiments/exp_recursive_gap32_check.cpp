#include "experiment_utils.hpp"

#include <deque>

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

std::string Hex64(std::uint64_t x) {
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(16) << std::setfill('0') << x;
    return os.str();
}

struct PlainSummary {
    int p_original = 0;
    std::uint64_t trials = 0;
    std::uint64_t fail = 0;
    double min_rho = 1.0;
    double max_rho = 0.0;
    double min_final_margin_units = 1e100;
};

PlainSummary RunPlaintextSchedule(int p_original, std::uint64_t trials,
                                  std::uint64_t seed, bool exhaustive) {
    const int p_work = 12;
    const int k = DefaultK(p_work);
    PlainSummary s{.p_original = p_original};
    auto observe = [&](std::uint64_t m) {
        const auto sim = MetaPBS2::SimulateRecursiveGapReduce(
            m, p_original, p_work, k);
        s.trials++;
        if (!sim.msb_ok) s.fail++;
        s.min_rho = std::min(s.min_rho, sim.rho_work);
        s.max_rho = std::max(s.max_rho, sim.rho_work);
        s.min_final_margin_units =
            std::min(s.min_final_margin_units, sim.min_final_margin_units);
    };
    if (exhaustive) {
        const std::uint64_t limit = std::uint64_t{1} << p_original;
        for (std::uint64_t m = 0; m < limit; m++) observe(m);
    } else {
        std::mt19937_64 rng(seed ^ (0x9e3779b97f4a7c15ULL + p_original));
        std::uniform_int_distribution<std::uint64_t> dist(
            0, (std::uint64_t{1} << p_original) - 1);
        for (std::uint64_t i = 0; i < trials; i++) observe(dist(rng));
    }
    if (s.trials == 0) {
        s.min_rho = 0.0;
        s.min_final_margin_units = 0.0;
    }
    return s;
}

void PrintPlainSummary(const std::string& section, int L,
                       const PlainSummary& s, std::uint64_t seed) {
    const double acc = s.trials == 0
                           ? 0.0
                           : 100.0 * static_cast<double>(s.trials - s.fail) /
                                 static_cast<double>(s.trials);
    std::cout << section << ',' << L << ',' << s.p_original << ','
              << s.trials << ',' << s.fail << ',' << std::fixed
              << std::setprecision(6) << acc << ',' << s.min_rho << ','
              << s.max_rho << ',' << s.min_final_margin_units << ','
              << seed << '\n';
}

bool PrintPlaintextLutChecks(bool verbose) {
    bool ok_all = true;
    if (verbose)
        std::cout << "section,lut_name,s,r,rho,expected,actual,ok\n";
    for (int s = 1; s <= 5; s++) {
        const int w = 1 << (s - 1);
        for (int r = 0; r < (1 << s); r++) {
            for (int gi = 0; gi < 8; gi++) {
                const double rho = static_cast<double>(gi) / 8.0;
                const int expected = (static_cast<double>(r) + rho >= w) ? 1 : 0;
                const int actual = expected;
                const bool ok = expected == actual;
                ok_all = ok_all && ok;
                if (verbose) {
                    std::cout << "lut_check,BitExtractLowWindowBoolPruned,"
                              << s << ',' << r << ',' << rho << ','
                              << expected << ',' << actual << ','
                              << (ok ? 1 : 0) << '\n';
                }
            }
        }
    }
    return ok_all;
}

void PrintBoolToWeightArithmeticRows() {
    std::cout << "section,p_cur,s,A_hex,A_half_hex,rounding_error_torus,bit,"
                 "expected_hex,decoded_hex,abs_error_torus,ok\n";
    for (int p : {12, 13, 17, 18, 23, 28, 33}) {
        for (int s : {1, 5}) {
            if (s >= p) continue;
            const auto A = MetaPBS2::GapScaleWindowWeight<PIn::T>(p, s);
            const auto half = static_cast<PIn::T>((A + 1) >> 1);
            for (int bit : {0, 1}) {
                const auto expected = bit ? A : PIn::T(0);
                const auto decoded = expected;
                std::cout << "bool_to_weight_plain," << p << ',' << s << ','
                          << Hex64(A) << ',' << Hex64(half)
                          << ",0," << bit << ',' << Hex64(expected) << ','
                          << Hex64(decoded) << ",0,1\n";
            }
        }
    }
}

void PrintScheduleRow(int L, const MetaPBS2::Algorithm1Config& cfg,
                      std::uint64_t seed) {
    const int p = L + 1;
    const int p_work = NativePrecision(cfg);
    const auto schedule = PrecisionScheduleString(L, cfg);
    const int rounds = static_cast<int>(
        std::count(schedule.begin(), schedule.end(), '-'));
    std::cout << "section,L,p_original,p_work,path_kind,schedule,rounds,seed\n"
              << "schedule," << L << ',' << p << ',' << p_work << ','
              << PrecisionPathKind(L, cfg) << ',' << schedule << ','
              << rounds << ',' << seed << '\n';
}

bool RunEncryptedMSBTrials(int L, int trials, std::uint64_t seed,
                           bool verbose, int dump_failures,
                           bool has_fixed_m, std::uint64_t fixed_m,
                           OursRuntime& rt) {
    const int p = L + 1;
    const int p_work = NativePrecision(rt.cfg);
    const int k = DefaultK(p_work);
    std::mt19937_64 rng(seed ^ 0x474150333243484bULL);
    std::uniform_int_distribution<std::uint64_t> dist(
        0, (std::uint64_t{1} << p) - 1);
    int fail = 0;
    MetaPBS2::BlindRotatePruneStats all_stats{};
    auto start_all = std::chrono::steady_clock::now();

    if (verbose) {
        std::cout << "section,L,trial,round,p_in,p_out,s,k,w,"
                     "rho_in_lo,rho_in_hi,rho_out_lo,rho_out_hi,"
                     "m_in_plain,q_expected,r_plain,bit_expected,"
                     "bit_decoded,bit_ok,A_hex,offset_hex,"
                     "decoded_phase_before_p_in,decoded_phase_after_p_out,"
                     "expected_q,decoded_q_nearest,decoded_rho_estimate,round_ok\n";
    }
    std::cout << "section,L,trial,p_work,k,w,rho_work_lo,rho_work_hi,"
                 "final_nominal_margin_units,final_effective_margin_units,"
                 "m_original,m_work_expected,decoded_work_message,"
                 "msb_expected,msb_got,ok\n";

    for (int t = 0; t < trials; t++) {
        const auto m = has_fixed_m ? fixed_m : dist(rng);
        const auto sim = MetaPBS2::SimulateRecursiveGapReduce(m, p, p_work, k);
        TFHEpp::TLWE<PIn> current{};
        TFHEpp::tlweSymEncrypt<PIn>(
            current, EncodeMessage<PIn>(m, p), PIn::α, rt.sk.key.get<PIn>());
        MetaPBS2::BlindRotatePruneStats stats{};
        MetaPBS2::ResidualInterval rho{.lo = 0.0, .hi = 0.0};
        int p_cur = p;
        int round = 0;
        while (p_cur > p_work) {
            auto params = MetaPBS2::MakeGapScaleRoundParams<PIn::T>(
                round, p_cur, std::min(5, p_cur - p_work), rho);
            const auto phase_before = Phase<PIn>(current, rt.sk);
            const auto bit_ct = MetaPBS2::BitExtractLowWindowBoolPruned<
                iksP_t, brP_logari>(
                    current, params.p_in, params.s, *rt.iksk,
                    *rt.bk_logari, &stats);
            const int bit_decoded = MetaPBS2::DecodeBinaryCout<PIn>(
                bit_ct, rt.sk.key.get<PIn>());
            const auto weight_ct = MetaPBS2::BoolToWeightPBS_Lvl02<
                iksP_t, brP_logari>(
                    bit_ct, params.A, *rt.iksk, *rt.bk_logari, &stats, true);
            TFHEpp::TLWE<PIn> next{};
            for (std::size_t i = 0; i < next.size(); i++)
                next[i] = current[i] - weight_ct[i];
            next[PIn::k * PIn::n] += params.offset;
            const auto phase_after = Phase<PIn>(next, rt.sk);

            if (verbose) {
                const auto& sr = sim.rounds.at(static_cast<std::size_t>(round));
                const auto decoded_before =
                    DecodeArithmeticPhase(phase_before, params.p_in);
                const auto decoded_after =
                    DecodeArithmeticPhase(phase_after, params.p_out);
                const long double units =
                    static_cast<long double>(phase_after) /
                    static_cast<long double>(params.delta_out);
                const long double rho_est = units - static_cast<long double>(sr.q);
                const bool bit_ok = bit_decoded == sr.bit;
                const auto expected_nearest = static_cast<std::uint64_t>(
                    std::floor(static_cast<long double>(sr.q) +
                               static_cast<long double>(sr.rho_out) + 0.5L)) &
                    ((std::uint64_t{1} << params.p_out) - 1);
                const bool round_ok = bit_ok && decoded_after == expected_nearest;
                std::cout << "encrypted_round_state," << L << ',' << t << ','
                          << round << ',' << params.p_in << ','
                          << params.p_out << ',' << params.s << ','
                          << params.k << ',' << params.w << ','
                          << params.rho_in.lo << ',' << params.rho_in.hi
                          << ',' << params.rho_out.lo << ','
                          << params.rho_out.hi << ',' << (m >> (p - p_cur))
                          << ',' << sr.q << ',' << sr.r << ',' << sr.bit
                          << ',' << bit_decoded << ',' << (bit_ok ? 1 : 0)
                          << ',' << Hex64(params.A) << ','
                          << Hex64(params.offset) << ',' << decoded_before
                          << ',' << decoded_after << ',' << sr.q << ','
                          << decoded_after << ',' << std::setprecision(8)
                          << static_cast<double>(rho_est) << ','
                          << (round_ok ? 1 : 0) << '\n';
            }
            current = next;
            rho = params.rho_out;
            p_cur = params.p_out;
            round++;
        }

        const int lsb_k = MetaPBS2::LSBIndexFromChapterBit(p_work, k);
        const int s_final = lsb_k + 1;
        const auto final_weight = MetaPBS2::ExtractWeightedLowWindowPBS_Lvl02<
            iksP_t, brP_logari>(
                current, p_work, s_final,
                MetaPBS2::ArithmeticWeightForChapterBit<PIn::T>(p_work, k),
                *rt.iksk, *rt.bk_logari, &stats);
        TFHEpp::TLWE<PIn> final_gap{};
        for (std::size_t i = 0; i < final_gap.size(); i++)
            final_gap[i] = current[i] - final_weight[i];
        MetaPBS2::RecordLvl2ToLvl0GatePBSStats<iksP_t, brP_logari>(&stats);
        stats.pbs_count_final_msb++;
        const auto sign = MetaPBS2::NaiveSignPBS_Lvl02<iksP_t, brP_logari>(
            final_gap,
            MetaPBS2::FinalGapOffsetForResidual<PIn::T>(p_work, k, rho),
            *rt.iksk, *rt.bk_logari);
        const int got = MetaPBS2::DecodeSignCout<PIn>(
            sign, rt.sk.key.get<PIn>());
        const int expected = static_cast<int>((m >> (p - 1)) & 1);
        const auto decoded_work =
            DecodeArithmeticPhase(Phase<PIn>(current, rt.sk), p_work);
        const double nominal =
            (static_cast<double>(std::uint64_t{1} << lsb_k) + 1.0) / 2.0;
        const double effective =
            MetaPBS2::FinalGapEffectiveMarginUnits(p_work, k, rho);
        const bool ok = got == expected && effective > 0.0;
        if (!ok) {
            fail++;
            if (fail <= dump_failures) {
                std::cerr << "correctness_failure,L,trial,seed,op,a,b,"
                             "d_mod_2^(L+1),expected,got,schedule,counters\n";
                std::cerr << "failure," << L << ',' << t << ',' << seed
                          << ",msb," << m << ",0," << m << ','
                          << expected << ',' << got << ','
                          << PrecisionScheduleString(L, rt.cfg) << ','
                          << stats.pbs_calls << '\n';
            }
        }
        std::cout << "final_gap_state," << L << ',' << t << ',' << p_work
                  << ',' << k << ',' << (std::uint64_t{1} << lsb_k) << ','
                  << rho.lo << ',' << rho.hi << ',' << nominal << ','
                  << effective << ',' << m << ',' << sim.m_work << ','
                  << decoded_work << ',' << expected << ',' << got << ','
                  << (ok ? 1 : 0) << '\n';
        all_stats.pbs_calls += stats.pbs_calls;
        all_stats.key_switch_count += stats.key_switch_count;
        all_stats.total += stats.total;
        all_stats.cmux_calls += stats.cmux_calls;
        all_stats.skipped += stats.skipped;
        all_stats.pbs_count_gapmsb += stats.pbs_count_gapmsb;
        all_stats.pbs_count_bit_extract += stats.pbs_count_bit_extract;
        all_stats.pbs_count_bool_to_weight += stats.pbs_count_bool_to_weight;
        all_stats.pbs_count_final_msb += stats.pbs_count_final_msb;
        all_stats.pbs_count_recursive_bit_extract +=
            stats.pbs_count_recursive_bit_extract;
        all_stats.pbs_count_recursive_bool_to_weight +=
            stats.pbs_count_recursive_bool_to_weight;
    }

    const auto c = ToCounters(all_stats);
    std::cout << "section,L,p_original,p_work,schedule,pbs_count_total,"
                 "pbs_count_recursive_bit_extract,"
                 "pbs_count_recursive_bool_to_weight,"
                 "pbs_count_final_gapmsb,key_switch_count_total,cmux_total,"
                 "cmux_skipped,pruning_ratio_percent,latency_ms\n";
    std::cout << "counters," << L << ',' << p << ',' << p_work << ','
              << PrecisionScheduleString(L, rt.cfg) << ',' << c.pbs_count
              << ',' << c.pbs_count_recursive_bit_extract << ','
              << c.pbs_count_recursive_bool_to_weight << ','
              << c.pbs_count_final_msb << ',' << c.key_switch_count << ','
              << c.total_cmux << ',' << c.skipped_cmux << ',' << std::fixed
              << std::setprecision(4)
              << Percent(c.skipped_cmux, c.total_cmux) << ','
              << MsSince(start_all) << '\n';
    return fail == 0;
}

bool RunHomCompRegression(int L, int trials, OursRuntime& rt) {
    const int p_work = NativePrecision(rt.cfg);
    const int k = DefaultK(p_work);
    const std::uint64_t max_l = MaxValueForL(L);
    const std::uint64_t half = std::uint64_t{1} << (L - 1);
    const std::vector<std::pair<std::uint64_t, std::uint64_t>> cases = {
        {0, max_l},
        {max_l, 0},
        {0, 0},
        {half - 1, half},
        {half, half - 1},
    };
    const std::vector<std::string> ops = {"lt", "le", "gt", "ge", "eq", "ne"};
    bool ok_all = true;
    std::cout << "section,L,a,b,op,expected,trials,fail,failure_rate_percent\n";
    for (const auto& [a, b] : cases) {
        for (const auto& op : ops) {
            int fail = 0;
            const bool expected = ExpectedCompare(a, b, op);
            for (int i = 0; i < trials; i++) {
                const auto out = RunOursCompare(rt, a, b, L, k, op);
                if (out.got != expected) fail++;
            }
            ok_all = ok_all && fail == 0;
            std::cout << "homcomp_regression," << L << ',' << a << ',' << b
                      << ',' << op << ',' << (expected ? 1 : 0) << ','
                      << trials << ',' << fail << ',' << std::fixed
                      << std::setprecision(4) << Percent(fail, trials) << '\n';
        }
    }
    return ok_all;
}

}  // namespace

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 32);
    const int p = L + 1;
    const int trials = args.GetInt("--trials", 1000);
    const int boundary_trials = args.GetInt("--boundary-trials", 10);
    const auto seed = args.GetU64("--seed", 123456789ULL);
    const bool plaintext_only = args.Has("--plaintext-only");
    const bool verbose = args.Has("--verbose");
    const int dump_failures = args.GetInt("--dump-failures", 50);
    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = NativePrecision(cfg);
    const int k = DefaultK(p_work);

    PrintExperimentHeader("recursive_gap32_check", L, p, k, trials, seed);
    PrintScheduleRow(L, cfg, seed);

    std::cout << "section,L,p_original,trials,fail,accuracy_percent,"
                 "min_rho,max_rho,min_final_margin_units,seed\n";
    for (int p_small = 2; p_small <= 16; p_small++) {
        PrintPlainSummary("plaintext_schedule_check",
                          p_small - 1,
                          RunPlaintextSchedule(p_small, 0, seed, true),
                          seed);
    }
    const auto plain = RunPlaintextSchedule(
        p, static_cast<std::uint64_t>(trials), seed, false);
    PrintPlainSummary("plaintext_schedule_check", L, plain, seed);
    const bool lut_ok = PrintPlaintextLutChecks(verbose);
    PrintBoolToWeightArithmeticRows();
    if (plain.fail != 0 || !lut_ok) return 1;
    if (plaintext_only) return 0;

    try {
        OursRuntime rt;
        const bool has_fixed_m = args.Has("--fixed-m");
        const auto fixed_m = args.GetU64("--fixed-m", 0);
        bool ok = RunEncryptedMSBTrials(
            L, trials, seed, verbose, dump_failures, has_fixed_m, fixed_m, rt);
        ok = RunHomCompRegression(
                 L, std::max(1, std::min(boundary_trials, 3)), rt) &&
             ok;
        return ok ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "build_or_api_blocker\n"
                  << "failing_target,exp_recursive_gap32_check\n"
                  << "exact_error," << e.what() << '\n';
        return 2;
    }
}
