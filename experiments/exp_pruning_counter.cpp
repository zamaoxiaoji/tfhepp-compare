#include "experiment_utils.hpp"

using namespace OursExperiments;

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 16);
    const int p = L + 1;
    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = WorkPrecisionForComparison(L, cfg);
    const int trials = args.GetInt("--trials", 1000);
    const std::uint64_t seed = args.GetU64("--seed", 123456789);
    const int k = ParseK(args.Get("--k", "default"), p_work);

    PrintExperimentHeader("pruning_counter", L, p, k, trials, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg) << "\n";
    std::cout << "scheme,L,p,k,trials,first_round_slot_period,"
                 "total_cmux_first_round,skipped_cmux_first_round,"
                 "pruning_ratio_first_round_percent,total_cmux_all,"
                 "skipped_cmux_all,pruning_ratio_all_percent,seed\n";

    try {
        if (p > 33)
            throw std::invalid_argument("p_original above supported L+1=33");
        if (k <= 0 || k >= p_work)
            throw std::invalid_argument("k must satisfy 1 <= k < p");

        OursRuntime rt;
        const int period = ExactBitExtractPeriod(p_work, k, rt.cfg);
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint64_t> dist(0, MaxValueForL(L));
        MetaPBS2::BlindRotatePruneStats all_stats{};
        bool any_fail = false;

        for (int i = 0; i < trials; i++) {
            const auto a = dist(rng);
            const auto b = dist(rng);
            const bool expected = a < b;
            const auto out = RunOursCompare(rt, a, b, L, k, "lt");
            const bool got = out.got;
            const auto& stats = out.stats;
            if (got != expected) {
                any_fail = true;
                std::cerr << "# failure_case a=" << a << " b=" << b << " expected="
                          << expected << " got=" << got << "\n";
            }
            all_stats.pbs_calls += stats.pbs_calls;
            all_stats.total += stats.total;
            all_stats.cmux_calls += stats.cmux_calls;
            all_stats.skipped += stats.skipped;
            for (std::size_t j = 0; j < stats.periods.size(); j++) {
                if (stats.periods[j] != period) continue;
                if (j < stats.total_by_pbs.size())
                    all_stats.total_by_pbs.push_back(stats.total_by_pbs[j]);
                if (j < stats.skipped_by_pbs.size())
                    all_stats.skipped_by_pbs.push_back(stats.skipped_by_pbs[j]);
                break;
            }
        }

        std::uint64_t first_total = 0;
        std::uint64_t first_skipped = 0;
        for (auto v : all_stats.total_by_pbs) first_total += v;
        for (auto v : all_stats.skipped_by_pbs) first_skipped += v;
        std::cout << "ours," << L << "," << p << "," << k << "," << trials
                  << "," << period << "," << first_total << ","
                  << first_skipped << "," << std::fixed << std::setprecision(4)
                  << Percent(first_skipped, first_total) << ","
                  << all_stats.total << "," << all_stats.skipped << ","
                  << Percent(all_stats.skipped, all_stats.total) << ","
                  << seed << "\n";
        return any_fail ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_pruning_counter invalid: " << e.what() << "\n";
        return 2;
    }
}
