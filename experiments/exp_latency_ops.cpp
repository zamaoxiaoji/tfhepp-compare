#include "experiment_utils.hpp"

using namespace OursExperiments;

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 11);
    const int p = L + 1;
    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = WorkPrecisionForComparison(L, cfg);
    const int runs = args.GetInt("--runs", 30);
    const int warmup = args.GetInt("--warmup", 5);
    const std::uint64_t seed = args.GetU64("--seed", 123456789);
    const auto ops = ParseOps(args.Get("--ops", "lt,eq"));
    const int k = ParseK(args.Get("--k", "default"), p_work);

    PrintExperimentHeader("latency_ops", L, p, k, runs, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg) << "\n"
              << "# warmup: " << warmup << "\n";
    std::cout << "scheme,L,p,k,op,runs,warmup,latency_ms_avg,"
                 "latency_ms_median,latency_ms_min,latency_ms_max,"
                 "latency_ms_std,pbs_count,key_switch_count,total_cmux,"
                 "skipped_cmux,pruning_ratio_percent,seed\n";

    try {
        OursRuntime rt;

        bool any_fail = false;
        for (const auto& op : ops) {
            std::mt19937_64 rng(seed + static_cast<std::uint64_t>(op[0]));
            std::uniform_int_distribution<std::uint64_t> dist(0, MaxValueForL(L));
            for (int i = 0; i < warmup; i++) {
                const auto a = dist(rng), b = dist(rng);
                const auto out = RunOursCompare(rt, a, b, L, k, op);
                if (out.got != ExpectedCompare(a, b, op))
                    throw std::runtime_error("warmup produced incorrect comparison output");
            }

            std::vector<double> latencies;
            MetaPBS2::BlindRotatePruneStats all_stats{};
            for (int i = 0; i < runs; i++) {
                const auto a = dist(rng), b = dist(rng);
                const bool expected = ExpectedCompare(a, b, op);
                const auto out = RunOursCompare(rt, a, b, L, k, op);
                latencies.push_back(out.ms);
                all_stats.pbs_calls += out.stats.pbs_calls;
                all_stats.total += out.stats.total;
                all_stats.cmux_calls += out.stats.cmux_calls;
                all_stats.skipped += out.stats.skipped;
                if (out.got != expected) {
                    any_fail = true;
                    std::cerr << "# failure_case op=" << op << " a=" << a
                              << " b=" << b << " expected=" << expected
                              << " got=" << out.got << "\n";
                }
            }
            const auto s = ComputeLatencyStats(latencies);
            const auto c = ToCounters(all_stats);
            std::cout << "ours," << L << "," << p << "," << k << "," << op
                      << "," << runs << "," << warmup << "," << std::fixed
                      << std::setprecision(4) << s.avg << "," << s.median
                      << "," << s.min << "," << s.max << "," << s.stddev
                      << "," << c.pbs_count << "," << c.key_switch_count
                      << "," << c.total_cmux << "," << c.skipped_cmux << ","
                      << Percent(c.skipped_cmux, c.total_cmux) << ","
                      << seed << "\n";
        }
        return any_fail ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_latency_ops invalid: " << e.what() << "\n";
        return 2;
    }
}
