#include "experiment_utils.hpp"

#include <deque>

using namespace OursExperiments;

int main(int argc, char** argv) {
    Args args{argc, argv};
    const int L = args.GetInt("--L", 16);
    const int p = L + 1;
    const auto cfg = PaperReview::Chapter3MetaPBSConfig();
    const int p_work = WorkPrecisionForComparison(L, cfg);
    const int radius = args.GetInt("--radius", 4);
    const int trials = args.GetInt("--trials-per-offset", 10000);
    const std::uint64_t seed = args.GetU64("--seed", 123456789);
    const std::string op = args.Get("--op", "lt");
    const int k = ParseK(args.Get("--k", "default"), p_work);
    const std::uint64_t T = std::uint64_t{1} << (L - 1);

    PrintExperimentHeader("boundary", L, p, k, trials, seed);
    std::cout << "# p_work: " << p_work << "\n"
              << "# path_kind: " << PrecisionPathKind(L, cfg) << "\n"
              << "# radius: " << radius << "\n";
    std::cout << "scheme,L,p,k,T,offset,a,b,op,trials,fail,"
                 "failure_rate_percent,seed\n";

    try {
        OursRuntime rt;

        bool any_fail = false;
        for (int offset = -radius; offset <= radius; offset++) {
            const auto b_signed = static_cast<long long>(T) + offset;
            if (b_signed < 0 || static_cast<std::uint64_t>(b_signed) > MaxValueForL(L))
                continue;
            const auto a = T;
            const auto b = static_cast<std::uint64_t>(b_signed);
            int fail = 0;
            std::deque<std::string> failures;
            for (int i = 0; i < trials; i++) {
                const bool expected = ExpectedCompare(a, b, op);
                const auto out = RunOursCompare(rt, a, b, L, k, op);
                if (out.got != expected) {
                    fail++;
                    any_fail = true;
                    if (failures.size() < 20) {
                        std::ostringstream ss;
                        ss << "# failure_case offset=" << offset << " a=" << a
                           << " b=" << b << " expected=" << expected
                           << " got=" << out.got;
                        failures.push_back(ss.str());
                    }
                }
            }
            for (const auto& f : failures) std::cerr << f << "\n";
            std::cout << "ours," << L << "," << p << "," << k << "," << T << ","
                      << offset << "," << a << "," << b << "," << op << ","
                      << trials << "," << fail << "," << std::fixed
                      << std::setprecision(4) << Percent(fail, trials) << ","
                      << seed << "\n";
        }
        return any_fail ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_boundary invalid: " << e.what() << "\n";
        return 2;
    }
}
