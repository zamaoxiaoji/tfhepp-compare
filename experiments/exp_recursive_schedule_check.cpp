#include "experiment_utils.hpp"

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

bool OriginalMSB(std::uint64_t m, int p) {
    return ((m >> (p - 1)) & 1) != 0;
}

std::uint64_t WorkMessage(std::uint64_t reduced, int p_reduced, int p_work) {
    return p_reduced < p_work ? (reduced << (p_work - p_reduced)) : reduced;
}

}  // namespace

int main(int argc, char** argv) {
    Args args{argc, argv};
    const bool plaintext_only = args.Has("--plaintext-only");
    const int L = args.GetInt("--L", 16);
    const int p_original = args.GetInt("--p", L + 1);
    const int trials = args.GetInt("--trials", plaintext_only ? 1000000 : 1000);
    const auto seed = args.GetU64("--seed", 123456789ULL);
    const bool verbose = args.Has("--verbose");

    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_native = NativePrecision(cfg);
    const auto schedule =
        MetaPBS2::MakeHE3DBStylePrecisionSchedule(p_original, p_native);
    const int p_work = p_native;
    const std::string schedule_text =
        MetaPBS2::ScheduleString(schedule, p_work);

    std::mt19937_64 rng(seed ^ 0x5343484544554c45ULL);
    const std::uint64_t max_m =
        p_original == 64 ? ~std::uint64_t{0}
                         : ((std::uint64_t{1} << p_original) - 1);
    std::uniform_int_distribution<std::uint64_t> dist(0, max_m);

    int fail = 0;
    std::vector<std::string> failures;

    if (plaintext_only) {
        for (int i = 0; i < trials; i++) {
            const auto m = dist(rng);
            int p_reduced = p_original;
            const auto reduced = MetaPBS2::SimulateHE3DBStylePrecisionReduce(
                m, p_original, p_reduced, p_native);
            const auto work = WorkMessage(reduced, p_reduced, p_work);
            const bool ok = OriginalMSB(m, p_original) == OriginalMSB(work, p_work);
            if (!ok) {
                fail++;
                if (failures.size() < 50) {
                    std::ostringstream ss;
                    ss << "m=" << m << " reduced=" << reduced
                       << " work=" << work;
                    failures.push_back(ss.str());
                }
            }
        }
        std::cout << "scheme,p_original,p_out,p_work,trials,fail,"
                     "accuracy_percent,schedule,seed\n";
        std::cout << "ours," << p_original << ',' << schedule.p_final << ','
                  << p_work << ',' << trials << ',' << fail << ','
                  << std::fixed << std::setprecision(4)
                  << Percent(trials - fail, trials) << ','
                  << schedule_text << ',' << seed << '\n';
        for (const auto& f : failures) std::cerr << "# failure_case " << f << "\n";
        return fail == 0 ? 0 : 1;
    }

    OursRuntime rt;
    if (verbose) {
        std::cout << "round,p_in,p_out,lut_kind,delta_in,delta_out\n";
        for (const auto& step : schedule.steps) {
            std::cout << "round," << step.p_in << ',' << step.p_out << ','
                      << step.kind << ','
                      << MetaPBS2::ReducerTorusScaleForPrecision<PIn::T>(step.p_in)
                      << ','
                      << MetaPBS2::ReducerTorusScaleForPrecision<PIn::T>(step.p_out)
                      << '\n';
        }
    }

    std::cout << "scheme,L,p_original,p_final,schedule,path_kind,trials,correct,"
                 "fail,accuracy_percent,seed\n";

    for (int i = 0; i < trials; i++) {
        const auto m = dist(rng);
        TFHEpp::TLWE<PIn> ct{};
        TFHEpp::tlweSymEncrypt<PIn>(
            ct, EncodeMessage<PIn>(m, p_original), PIn::α, rt.sk.key.get<PIn>());

        MetaPBS2::BlindRotatePruneStats stats{};
        std::vector<MetaPBS2::PrecisionReduceRoundOutput<PIn>> rounds;
        const auto reduced = MetaPBS2::HE3DBStylePrecisionReduceToNative<iksP_t, brP_logari>(
            ct, p_original, p_native, *rt.iksk, *rt.bk_logari, &stats, &rounds);

        bool ok = true;
        int removed = 0;
        if (verbose)
            std::cout << "section,L,trial,round,p_in,p_out,m_original,"
                         "expected_round_message,decoded_round_message,round_ok\n";
        for (const auto& round : rounds) {
            removed += round.p_in - round.p_out;
            const auto expected = m >> removed;
            const auto decoded = DecodeArithmeticPhase(Phase<PIn>(round.ct_out, rt.sk),
                                                       round.p_out);
            const bool round_ok = decoded == expected;
            ok = ok && round_ok;
            if (verbose) {
                std::cout << "round_state," << L << ',' << i << ','
                          << round.round << ',' << round.p_in << ','
                          << round.p_out << ',' << m << ',' << expected
                          << ',' << decoded << ',' << (round_ok ? 1 : 0)
                          << '\n';
            }
        }
        const auto reduced_msg = m >> (p_original - reduced.p_out);
        const auto expected_work = WorkMessage(reduced_msg, reduced.p_out, p_work);
        const auto decoded_work = DecodeArithmeticPhase(
            Phase<PIn>(reduced.ct_out, rt.sk), p_work);
        ok = ok && (decoded_work == expected_work) &&
             (OriginalMSB(m, p_original) == OriginalMSB(expected_work, p_work));
        if (!ok) {
            fail++;
            if (failures.size() < 50) {
                std::ostringstream ss;
                ss << "m=" << m << " expected_work=" << expected_work
                   << " decoded_work=" << decoded_work;
                failures.push_back(ss.str());
            }
        }
    }

    std::cout << "ours," << L << ',' << p_original << ','
              << schedule.p_final << ',' << schedule_text << ','
              << PrecisionPathKind(L, cfg) << ',' << trials << ','
              << (trials - fail) << ',' << fail << ','
              << std::fixed << std::setprecision(4)
              << Percent(trials - fail, trials) << ',' << seed << '\n';
    for (const auto& f : failures) std::cerr << "# failure_case " << f << "\n";
    return fail == 0 ? 0 : 1;
}
