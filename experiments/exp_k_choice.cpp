#include "experiment_utils.hpp"

using namespace OursExperiments;

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 16);
    const int p = L + 1;
    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = WorkPrecisionForComparison(L, cfg);
    const int trials = args.GetInt("--trials", 100000);
    const std::uint64_t seed = args.GetU64("--seed", 123456789);
    const auto k_values = ParseKList(args.Get("--k-list", "p-2,p-3,p-4,p-5,1"), p_work);

    PrintExperimentHeader("k_choice", L, p, DefaultK(p_work), trials, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg) << "\n";
    std::cout << "scheme,L,p,k,k_label,w_k,gap_width_in_delta,trials,fail,"
                 "failure_rate_percent,actual_first_round_slot_period,"
                 "total_cmux,skipped_cmux,pruning_ratio_percent,"
                 "latency_ms_avg,latency_ms_median,seed,status,reason\n";

    try {
        const bool kernel_supported = p <= 33;
        std::unique_ptr<OursRuntime> rt;
        if (kernel_supported) rt = std::make_unique<OursRuntime>();
        for (const auto& [k, label] : k_values) {
            std::string invalid_reason;
            if (k <= 0 || k >= p_work) invalid_reason = "k_out_of_work_precision_range";
            if (!kernel_supported) invalid_reason = "p_original_above_supported_33";

            const std::uint64_t w_k =
                (k > 0 && k < p_work && p_work - 1 - k < 63)
                    ? (std::uint64_t{1} << (p_work - 1 - k))
                    : 0;
            const std::uint64_t gap_width = w_k == 0 ? 0 : 2 * w_k + 1;
            if (!invalid_reason.empty()) {
                std::cout << "ours," << L << "," << p << "," << k << ","
                          << label << "," << w_k << "," << gap_width
                          << ",0,0,0.0000,0,0,0,0.0000,0.0000,0.0000,"
                          << seed << ",invalid," << invalid_reason << "\n";
                continue;
            }

            std::mt19937_64 rng(seed + static_cast<std::uint64_t>(k));
            std::uniform_int_distribution<std::uint64_t> dist(0, MaxValueForL(L));
            int fail = 0;
            std::vector<double> latencies;
            MetaPBS2::BlindRotatePruneStats all_stats{};
            int period = 0;
            try {
                period = ExactBitExtractPeriod(p_work, k, rt->cfg);
            } catch (const std::exception& e) {
                invalid_reason = e.what();
            }
            if (!invalid_reason.empty()) {
                std::cout << "ours," << L << "," << p << "," << k << ","
                          << label << "," << w_k << "," << gap_width
                          << ",0,0,0.0000,0,0,0,0.0000,0.0000,0.0000,"
                          << seed << ",invalid," << invalid_reason << "\n";
                continue;
            }

            for (int i = 0; i < trials; i++) {
                const auto a = dist(rng), b = dist(rng);
                const bool expected = a < b;
                const auto out = RunOursCompare(*rt, a, b, L, k, "lt");
                latencies.push_back(out.ms);
                all_stats.pbs_calls += out.stats.pbs_calls;
                all_stats.total += out.stats.total;
                all_stats.cmux_calls += out.stats.cmux_calls;
                all_stats.skipped += out.stats.skipped;
                if (out.got != expected) {
                    fail++;
                    if (fail <= 20)
                        std::cerr << "# failure_case k=" << k << " a=" << a
                                  << " b=" << b << " expected=" << expected
                                  << " got=" << out.got << "\n";
                }
            }
            const auto s = ComputeLatencyStats(latencies);
            std::cout << "ours," << L << "," << p << "," << k << ","
                      << label << "," << w_k << "," << gap_width << ","
                      << trials << "," << fail << "," << std::fixed
                      << std::setprecision(4) << Percent(fail, trials) << ","
                      << period << "," << all_stats.total << ","
                      << all_stats.skipped << ","
                      << Percent(all_stats.skipped, all_stats.total) << ","
                      << s.avg << "," << s.median << "," << seed
                      << ",ok,\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_k_choice invalid: " << e.what() << "\n";
        return 2;
    }
}
