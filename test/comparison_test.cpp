/**
 * @file comparison_test.cpp
 * @brief Time + accuracy benchmark for the homomorphic comparison pipelines.
 *
 *   1) ethmsb     — ETHMSB + gap offset           (samplepaper.tex)
 *   2) three_pbs  — Pruned BitExtract + Micro-PBS B2A + ETHMSB (Chapter 3)
 *   3) HEDB       — HE3DB original HomMSB         (only when HE3DB sources are
 *                                                  available at ../HE3DB/src)
 *
 * Layout follows HE3DB's test/comparison_test.cpp: a Lvl1 entry point for
 * 1–10-bit inputs, a Lvl2 entry point for 11–32-bit inputs, with five
 * predicates (>, ≥, <, ≤, =) timed per algorithm.
 */
#include <chrono>
#include <array>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../comparison/comparison.h"

#if HAVE_HE3DB_ORIGINAL
// HE3DB's utils/types.h defines LOGIC/ARITHMETIC as preprocessor macros.
// Drop them after the include so they don't textually replace identifiers
// elsewhere (our tfhepp_compare::LOGIC is a constexpr bool with the same name).
#include "HEDB/comparison/comparison.h"
#undef LOGIC
#undef ARITHMETIC
#endif

using namespace tfhepp_compare;

namespace
{
    constexpr const char *kPredicateNames[5] = {">", ">=", "<", "<=", "=="};

    struct AlgoMetrics {
        std::vector<uint32_t> error_time   = std::vector<uint32_t>(5, 0);
        std::vector<double>   total_time_ms = std::vector<double>(5, 0.0);
        bool                  has_first_failure = false;
        uint32_t              first_plain_bits = 0;
        size_t                first_op = 0;
        uint64_t              first_p0 = 0;
        uint64_t              first_p1 = 0;
        uint64_t              first_expected = 0;
        uint64_t              first_decoded = 0;
    };

    template <typename T>
    uint64_t printable_uint(T v)
    {
        return static_cast<uint64_t>(v);
    }

    template <typename P>
    void record_failure(AlgoMetrics &m, size_t idx, uint32_t plain_bits,
                        typename P::T p0, typename P::T p1,
                        typename P::T expected, typename Lvl1::T decoded)
    {
        if (m.has_first_failure) return;
        m.has_first_failure = true;
        m.first_plain_bits = plain_bits;
        m.first_op = idx;
        m.first_p0 = printable_uint(p0);
        m.first_p1 = printable_uint(p1);
        m.first_expected = printable_uint(expected);
        m.first_decoded = printable_uint(decoded);
    }

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
        constexpr bool result_type = LOGIC;
        TLWELvl1                              cres;
        std::chrono::system_clock::time_point start, end;

        const auto exec = [&](auto fn, int idx) {
            start = std::chrono::system_clock::now();
            fn(c0, c1, cres, plain_bits, ek, result_type);
            end = std::chrono::system_clock::now();
            m.total_time_ms[idx] +=
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                    .count();
            const typename Lvl1::T decoded =
                TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
            if (decoded != truth[idx]) {
                m.error_time[idx] += 1;
                record_failure<P>(m, idx, plain_bits, p0, p1, truth[idx],
                                  decoded);
            }
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
                      << kPredicateNames[i] << "  " << std::setw(11)
                      << std::right << std::fixed << std::setprecision(2)
                      << (m.total_time_ms[i] / num_test) << "   "
                      << std::setw(6) << m.error_time[i] << "\n";
        }
        if (m.has_first_failure) {
            std::cout << "    first_failure plain_bits=" << m.first_plain_bits
                      << " op=" << kPredicateNames[m.first_op]
                      << " p0=" << m.first_p0 << " p1=" << m.first_p1
                      << " expected=" << m.first_expected
                      << " decoded=" << m.first_decoded << "\n";
        }
    }

    template <typename P>
    void run_three_pbs_micro_trial(
        AlgoMetrics &m, TFHEpp::TLWE<P> &c0, TFHEpp::TLWE<P> &c1,
        typename P::T p0, typename P::T p1, uint32_t plain_bits,
        TFHEEvalKey &ek, const three_pbs::FastB2AEvalKeyPack &micro_pack,
        const TFHESecretKey &sk)
    {
        const auto truth = std::array<typename P::T, 5>{
            (typename P::T)(p0 > p1),  (typename P::T)(p0 >= p1),
            (typename P::T)(p0 < p1),  (typename P::T)(p0 <= p1),
            (typename P::T)(p0 == p1),
        };

        constexpr bool result_type = LOGIC;
        const auto exec = [&](auto fn, int idx) {
            TLWELvl1 cres;
            const auto start = std::chrono::system_clock::now();
            fn(cres);
            const auto end = std::chrono::system_clock::now();
            m.total_time_ms[idx] +=
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                    .count();
            const typename Lvl1::T decoded =
                TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
            if (decoded != truth[idx]) {
                m.error_time[idx] += 1;
                record_failure<P>(m, idx, plain_bits, p0, p1, truth[idx],
                                  decoded);
            }
        };

        exec(
            [&](TLWELvl1 &res) {
                TFHEpp::TLWE<P> sub_tlwe;
                for (size_t i = 0; i <= P::k * P::n; i++)
                    sub_tlwe[i] = c1[i] - c0[i];
                three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek,
                                  micro_pack, result_type);
            },
            0);

        exec(
            [&](TLWELvl1 &res) {
                TFHEpp::TLWE<P> sub_tlwe;
                for (size_t i = 0; i <= P::k * P::n; i++)
                    sub_tlwe[i] = c0[i] - c1[i];
                three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek,
                                  micro_pack, LOGIC);
                HomNOT<Lvl1>(res, res);
                if (IS_ARITHMETIC(result_type))
                    tfhepp_compare::LOG_to_ARI(res, res, ek);
            },
            1);

        exec(
            [&](TLWELvl1 &res) {
                TFHEpp::TLWE<P> sub_tlwe;
                for (size_t i = 0; i <= P::k * P::n; i++)
                    sub_tlwe[i] = c0[i] - c1[i];
                three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek,
                                  micro_pack, result_type);
            },
            2);

        exec(
            [&](TLWELvl1 &res) {
                TFHEpp::TLWE<P> sub_tlwe;
                for (size_t i = 0; i <= P::k * P::n; i++)
                    sub_tlwe[i] = c1[i] - c0[i];
                three_pbs::HomMSB(res, sub_tlwe, plain_bits + 1, ek,
                                  micro_pack, LOGIC);
                HomNOT<Lvl1>(res, res);
                if (IS_ARITHMETIC(result_type))
                    tfhepp_compare::LOG_to_ARI(res, res, ek);
            },
            3);

        exec(
            [&](TLWELvl1 &res) {
                TLWELvl1 ge_tlwe, le_tlwe;
                TFHEpp::TLWE<P> sub_ge, sub_le;
                for (size_t i = 0; i <= P::k * P::n; i++) {
                    sub_ge[i] = c0[i] - c1[i];
                    sub_le[i] = c1[i] - c0[i];
                }
                three_pbs::HomMSB(ge_tlwe, sub_ge, plain_bits + 1, ek,
                                  micro_pack, LOGIC);
                HomNOT<Lvl1>(ge_tlwe, ge_tlwe);
                three_pbs::HomMSB(le_tlwe, sub_le, plain_bits + 1, ek,
                                  micro_pack, LOGIC);
                HomNOT<Lvl1>(le_tlwe, le_tlwe);
                HomAND(res, ge_tlwe, le_tlwe, ek, result_type);
            },
            4);
    }

    // Each per-algorithm block builds five lambdas that bind to that
    // algorithm's namespace, then dispatches one trial.
#define RUN_ALGO_TRIAL(NS, METRICS)                                          \
    do {                                                                      \
        auto gt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,        \
                     bool rt) { NS::greater_than<P>(a, b, r, pb, e, rt); };  \
        auto ge = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,        \
                     bool rt) {                                               \
            NS::greater_than_equal<P>(a, b, r, pb, e, rt);                   \
        };                                                                    \
        auto lt = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,        \
                     bool rt) { NS::less_than<P>(a, b, r, pb, e, rt); };     \
        auto le = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,        \
                     bool rt) {                                               \
            NS::less_than_equal<P>(a, b, r, pb, e, rt);                      \
        };                                                                    \
        auto eq = [](auto &a, auto &b, auto &r, uint32_t pb, auto &e,        \
                     bool rt) { NS::equal<P>(a, b, r, pb, e, rt); };         \
        run_one_trial<P>(METRICS, gt, ge, lt, le, eq, c0, c1, p0, p1,        \
                         plain_bits, ek, sk);                                \
    } while (0)

    void tlwelvl1_comparison_test(uint32_t plain_bits, int num_test)
    {
        std::cout << "\n=== Lvl1 comparison: plain_bits=" << plain_bits
                  << ", trials=" << num_test << " ===\n";
        using P = Lvl1;

        TFHESecretKey sk;
        TFHEEvalKey   ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        const auto micro_pack = three_pbs::GenerateFastB2AEvalKeyPack(sk, false);

        const uint32_t scale_bits =
            std::numeric_limits<P::T>::digits - plain_bits - 1;
        std::default_random_engine                  engine(0x5eed1000u ^
                                         (plain_bits * 0x9e3779b9u));
        std::uniform_int_distribution<typename P::T> message(
            0, (typename P::T(1) << (plain_bits - 1)) - 1);

        AlgoMetrics ethmsb_m, three_pbs_m;
#if HAVE_HE3DB_ORIGINAL
        AlgoMetrics he3db_m;
#endif
        for (int t = 0; t < num_test; t++) {
            const typename P::T p0 = message(engine);
            const typename P::T p1 = message(engine);
            TFHEpp::TLWE<P>     c0 = tlweSymInt32Encrypt<P>(
                p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
            TFHEpp::TLWE<P> c1 = tlweSymInt32Encrypt<P>(
                p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

            RUN_ALGO_TRIAL(ethmsb, ethmsb_m);
            run_three_pbs_micro_trial<P>(three_pbs_m, c0, c1, p0, p1,
                                         plain_bits, ek, micro_pack, sk);
#if HAVE_HE3DB_ORIGINAL
            RUN_ALGO_TRIAL(HEDB,      he3db_m);
#endif
        }

        print_results("ETHMSB+offset    (samplepaper)", ethmsb_m, num_test);
        print_results("Pruned 3-PBS optimized        ", three_pbs_m, num_test);
#if HAVE_HE3DB_ORIGINAL
        print_results("HE3DB HomMSB     (original)   ", he3db_m, num_test);
#endif
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
        const auto micro_pack = three_pbs::GenerateFastB2AEvalKeyPack(sk, true);

        const uint32_t scale_bits =
            std::numeric_limits<P::T>::digits - plain_bits - 1;
        std::default_random_engine                  engine(0x5eed2000u ^
                                         (plain_bits * 0x9e3779b9u));
        std::uniform_int_distribution<typename P::T> message(
            0, (typename P::T(1) << (plain_bits - 1)) - 1);

        AlgoMetrics ethmsb_m, three_pbs_m;
#if HAVE_HE3DB_ORIGINAL
        AlgoMetrics he3db_m;
#endif
        for (int t = 0; t < num_test; t++) {
            const typename P::T p0 = message(engine);
            const typename P::T p1 = message(engine);
            TFHEpp::TLWE<P>     c0 = tlweSymInt32Encrypt<P>(
                p0, P::α, std::pow(2., scale_bits), sk.key.get<P>());
            TFHEpp::TLWE<P> c1 = tlweSymInt32Encrypt<P>(
                p1, P::α, std::pow(2., scale_bits), sk.key.get<P>());

            RUN_ALGO_TRIAL(ethmsb, ethmsb_m);
            run_three_pbs_micro_trial<P>(three_pbs_m, c0, c1, p0, p1,
                                         plain_bits, ek, micro_pack, sk);
#if HAVE_HE3DB_ORIGINAL
            RUN_ALGO_TRIAL(HEDB,      he3db_m);
#endif
        }

        print_results("ETHMSB+offset    (samplepaper)", ethmsb_m, num_test);
        print_results("Pruned 3-PBS optimized        ", three_pbs_m, num_test);
#if HAVE_HE3DB_ORIGINAL
        print_results("HE3DB HomMSB     (original)   ", he3db_m, num_test);
#endif
    }

#undef RUN_ALGO_TRIAL

} // namespace

int main(int argc, char **argv)
{
    int num_test = 50;
    if (argc >= 2) num_test = std::stoi(argv[1]);

    std::cout << "TFHEpp comparison benchmark\n";
    std::cout << "  Algorithms: ETHMSB+offset, Pruned 3-PBS optimized";
#if HAVE_HE3DB_ORIGINAL
    std::cout << ", HE3DB HomMSB";
#else
    std::cout << "  (HE3DB original not available — set "
                 "../HE3DB/src to enable)";
#endif
    std::cout << "\n";

    tlwelvl1_comparison_test(4, num_test);
    tlwelvl1_comparison_test(5, num_test);
    tlwelvl1_comparison_test(8, num_test);
    tlwelvl2_comparison_test(16, num_test);
    tlwelvl2_comparison_test(32, num_test);
    return 0;
}
