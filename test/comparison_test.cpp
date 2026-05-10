/**
 * @file comparison_test.cpp
 * @brief Time + accuracy benchmark for the two homomorphic comparison
 *        algorithms exposed by comparison/ :
 *
 *           1) ethmsb     — Algorithm 1 (ETHMSB + gap offset, samplepaper)
 *           2) three_pbs  — Algorithm 2 (Pruned BitExtract + B2A + ETHMSB,
 *                                          Chapter 3 of the thesis)
 *
 * Layout follows HE3DB's test/comparison_test.cpp: one Lvl1 entry point for
 * 1–10-bit inputs, one Lvl2 entry point for 11–32-bit inputs, with five
 * predicates (>, ≥, <, ≤, =) timed per algorithm.
 */
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../comparison/comparison.h"

using namespace tfhepp_compare;

namespace
{
    constexpr const char *kPredicateNames[5] = {">", ">=", "<", "<=", "=="};

    struct AlgoMetrics {
        std::vector<uint32_t> error_time = std::vector<uint32_t>(5, 0);
        std::vector<double>   total_time_ms = std::vector<double>(5, 0.0);
    };

    template <typename P, typename FnGT, typename FnGE, typename FnLT,
              typename FnLE, typename FnEQ>
    void run_one_trial(AlgoMetrics &m, FnGT &gt, FnGE &ge, FnLT &lt, FnLE &le,
                       FnEQ &eq, TFHEpp::TLWE<P> &c0, TFHEpp::TLWE<P> &c1,
                       typename P::T p0, typename P::T p1, uint32_t plain_bits,
                       TFHEEvalKey &ek, const TFHESecretKey &sk)
    {
        (void) plain_bits;
        const auto truth = std::array<typename P::T, 5>{
            (typename P::T)(p0 > p1),  (typename P::T)(p0 >= p1),
            (typename P::T)(p0 < p1),  (typename P::T)(p0 <= p1),
            (typename P::T)(p0 == p1),
        };
        // Always decode in LOGIC mode for simplicity; that is what every op is
        // configured to emit below.
        constexpr bool result_type = LOGIC;
        TLWELvl1 cres;
        std::chrono::system_clock::time_point start, end;

        // Each lambda already binds the right HomMSB; we just call them in turn.
        const auto exec = [&](auto fn, int idx) {
            start = std::chrono::system_clock::now();
            fn(c0, c1, cres, plain_bits, ek, result_type);
            end   = std::chrono::system_clock::now();
            m.total_time_ms[idx] +=
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                    .count();
            const typename Lvl1::T decoded =
                TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
            if (decoded != truth[idx]) m.error_time[idx] += 1;
        };

        exec(gt, 0);
        exec(ge, 1);
        exec(lt, 2);
        exec(le, 3);
        exec(eq, 4);
    }

    void print_results(const std::string &label, const AlgoMetrics &m,
                       int num_test)
    {
        std::cout << "  [" << label << "]\n";
        std::cout << "    op    avg_time(ms)   errors\n";
        for (size_t i = 0; i < 5; i++) {
            std::cout << "    " << std::setw(4) << std::left
                      << kPredicateNames[i] << "  " << std::setw(11) << std::right
                      << std::fixed << std::setprecision(2)
                      << (m.total_time_ms[i] / num_test) << "   "
                      << std::setw(6) << m.error_time[i] << "\n";
        }
    }

    void tlwelvl1_comparison_test(uint32_t plain_bits, int num_test)
    {
        std::cout << "\n=== Lvl1 comparison: plain_bits=" << plain_bits
                  << ", trials=" << num_test << " ===\n";
        using P = Lvl1;

        TFHESecretKey sk;
        TFHEEvalKey   ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplaceiksk<Lvl10>(sk);

        const uint32_t scale_bits =
            std::numeric_limits<P::T>::digits - plain_bits - 1;
        std::random_device                          seed_gen;
        std::default_random_engine                  engine(seed_gen());
        std::uniform_int_distribution<typename P::T> message(
            0, (typename P::T(1) << (plain_bits - 1)) - 1);

        AlgoMetrics ethmsb_m, three_pbs_m;
        for (int t = 0; t < num_test; t++) {
            const typename P::T p0 = message(engine);
            const typename P::T p1 = message(engine);
            TFHEpp::TLWE<P>     c0 = tlweSymInt32Encrypt<P>(
                p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
            TFHEpp::TLWE<P>     c1 = tlweSymInt32Encrypt<P>(
                p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

            // Algorithm 1: ETHMSB + offset
            {
                auto gt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { ethmsb::greater_than<P>(a, b, r, pb, e, rt); };
                auto ge = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    ethmsb::greater_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto lt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { ethmsb::less_than<P>(a, b, r, pb, e, rt); };
                auto le = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    ethmsb::less_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto eq = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { ethmsb::equal<P>(a, b, r, pb, e, rt); };
                run_one_trial<P>(ethmsb_m, gt, ge, lt, le, eq, c0, c1, p0, p1,
                                 plain_bits, ek, sk);
            }
            // Algorithm 2: Pruned 3-PBS
            {
                auto gt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::greater_than<P>(a, b, r, pb, e, rt);
                };
                auto ge = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::greater_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto lt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::less_than<P>(a, b, r, pb, e, rt);
                };
                auto le = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::less_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto eq = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { three_pbs::equal<P>(a, b, r, pb, e, rt); };
                run_one_trial<P>(three_pbs_m, gt, ge, lt, le, eq, c0, c1, p0,
                                 p1, plain_bits, ek, sk);
            }
        }

        print_results("ETHMSB+offset    (samplepaper)", ethmsb_m, num_test);
        print_results("Pruned 3-PBS     (Chapter 3)  ", three_pbs_m, num_test);
    }

    void tlwelvl2_comparison_test(uint32_t plain_bits, int num_test)
    {
        std::cout << "\n=== Lvl2 comparison: plain_bits=" << plain_bits
                  << ", trials=" << num_test << " ===\n";
        using P = Lvl2;

        TFHESecretKey sk;
        TFHEEvalKey   ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplacebkfft<Lvl02>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        ek.emplaceiksk<Lvl20>(sk);
        ek.emplaceiksk<Lvl21>(sk);

        const uint32_t scale_bits =
            std::numeric_limits<P::T>::digits - plain_bits - 1;
        std::random_device                          seed_gen;
        std::default_random_engine                  engine(seed_gen());
        std::uniform_int_distribution<typename P::T> message(
            0, (typename P::T(1) << (plain_bits - 1)) - 1);

        AlgoMetrics ethmsb_m, three_pbs_m;
        for (int t = 0; t < num_test; t++) {
            const typename P::T p0 = message(engine);
            const typename P::T p1 = message(engine);
            TFHEpp::TLWE<P>     c0 = tlweSymInt32Encrypt<P>(
                p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
            TFHEpp::TLWE<P>     c1 = tlweSymInt32Encrypt<P>(
                p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

            {
                auto gt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { ethmsb::greater_than<P>(a, b, r, pb, e, rt); };
                auto ge = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    ethmsb::greater_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto lt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { ethmsb::less_than<P>(a, b, r, pb, e, rt); };
                auto le = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    ethmsb::less_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto eq = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { ethmsb::equal<P>(a, b, r, pb, e, rt); };
                run_one_trial<P>(ethmsb_m, gt, ge, lt, le, eq, c0, c1, p0, p1,
                                 plain_bits, ek, sk);
            }
            {
                auto gt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::greater_than<P>(a, b, r, pb, e, rt);
                };
                auto ge = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::greater_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto lt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::less_than<P>(a, b, r, pb, e, rt);
                };
                auto le = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) {
                    three_pbs::less_than_equal<P>(a, b, r, pb, e, rt);
                };
                auto eq = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,
                             bool rt) { three_pbs::equal<P>(a, b, r, pb, e, rt); };
                run_one_trial<P>(three_pbs_m, gt, ge, lt, le, eq, c0, c1, p0,
                                 p1, plain_bits, ek, sk);
            }
        }

        print_results("ETHMSB+offset    (samplepaper)", ethmsb_m, num_test);
        print_results("Pruned 3-PBS     (Chapter 3)  ", three_pbs_m, num_test);
    }

} // namespace

int main(int argc, char **argv)
{
    int num_test = 50;
    if (argc >= 2) num_test = std::stoi(argv[1]);

    std::cout << "TFHEpp comparison: ETHMSB+offset vs Pruned 3-PBS\n";

    tlwelvl1_comparison_test(4, num_test);
    tlwelvl1_comparison_test(8, num_test);
    tlwelvl2_comparison_test(16, num_test);
    tlwelvl2_comparison_test(32, num_test);
    return 0;
}
