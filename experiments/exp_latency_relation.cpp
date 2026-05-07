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
    const std::string op = args.Get("--op", "lt");
    const int k = ParseK(args.Get("--k", "default"), p_work);

    PrintExperimentHeader("latency_relation", L, p, k, runs, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg) << "\n"
              << "# warmup: " << warmup << "\n";
    std::cout << "scheme,L,p,k,op,runs,warmup,latency_ms_avg,"
                 "latency_ms_median,latency_ms_min,latency_ms_max,"
                 "latency_ms_std,seed\n";

    try {
        OursRuntime rt;

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint64_t> dist(0, MaxValueForL(L));
        for (int i = 0; i < warmup; i++) {
            const auto a = dist(rng), b = dist(rng);
            const auto out = RunOursCompare(rt, a, b, L, k, op);
            if (out.got != ExpectedCompare(a, b, op))
                throw std::runtime_error("warmup produced incorrect comparison output");
        }

        std::vector<double> latencies;
        bool any_fail = false;
        for (int i = 0; i < runs; i++) {
            const auto a = dist(rng), b = dist(rng);
            const bool expected = ExpectedCompare(a, b, op);
            const auto out = RunOursCompare(rt, a, b, L, k, op);
            latencies.push_back(out.ms);
            if (out.got != expected) {
                any_fail = true;
                std::cerr << "# failure_case a=" << a << " b=" << b
                          << " expected=" << expected << " got=" << out.got
                          << "\n";
            }
        }
        const auto s = ComputeLatencyStats(latencies);
        std::cout << "ours," << L << "," << p << "," << k << "," << op << ","
                  << runs << "," << warmup << "," << std::fixed
                  << std::setprecision(4) << s.avg << "," << s.median << ","
                  << s.min << "," << s.max << "," << s.stddev << ","
                  << seed << "\n";
        return any_fail ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_latency_relation invalid: " << e.what() << "\n";
        return 2;
    }
}
