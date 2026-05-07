#include "experiment_utils.hpp"

#include <deque>

using namespace OursExperiments;

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 11);
    const int p = L + 1;
    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = WorkPrecisionForComparison(L, cfg);
    const int trials = args.GetInt("--trials", 100000);
    const std::uint64_t seed = args.GetU64("--seed", 123456789);
    const auto ops = ParseOps(args.Get("--op", "all"));
    const int k = ParseK(args.Get("--k", "default"), p_work);

    PrintExperimentHeader("correctness", L, p, k, trials, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg) << "\n";
    std::cout << "scheme,L,p,k,op,trials,correct,fail,accuracy_percent,"
                 "failure_rate_percent,seed\n";

    try {
        OursRuntime rt;

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint64_t> dist(0, MaxValueForL(L));
        bool any_fail = false;

        for (const auto& op : ops) {
            int correct = 0;
            std::deque<std::string> failures;
            for (int i = 0; i < trials; i++) {
                const auto a = dist(rng);
                const auto b = dist(rng);
                const bool expected = ExpectedCompare(a, b, op);
                const auto out = RunOursCompare(rt, a, b, L, k, op);
                if (out.got == expected) {
                    correct++;
                } else {
                    any_fail = true;
                    if (failures.size() < 20) {
                        std::ostringstream ss;
                        ss << "# failure_case op=" << op << " a=" << a << " b=" << b
                           << " expected=" << expected << " got=" << out.got;
                        failures.push_back(ss.str());
                    }
                }
            }
            const int fail = trials - correct;
            for (const auto& f : failures) std::cerr << f << "\n";
            std::cout << "ours," << L << "," << p << "," << k << "," << op << ","
                      << trials << "," << correct << "," << fail << ","
                      << std::fixed << std::setprecision(4)
                      << Percent(correct, trials) << "," << Percent(fail, trials)
                      << "," << seed << "\n";
        }
        return any_fail ? 1 : 0;
    } catch (const std::exception& e) {
        for (const auto& op : ops)
            std::cout << "ours," << L << "," << p << "," << k << "," << op
                      << ",0,0,0,0.0000,0.0000," << seed << "\n";
        std::cerr << "exp_correctness invalid: " << e.what() << "\n";
        return 2;
    }
}
