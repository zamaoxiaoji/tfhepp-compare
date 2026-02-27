#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <string>

#include <metapbs.hpp>
#include <tfhe++.hpp>

namespace {

#define CHECK(cond)                                                       \
    do {                                                                  \
        if (!(cond)) {                                                    \
            std::cerr << "CHECK failed: " << #cond << " @ " << __FILE__   \
                      << ":" << __LINE__ << std::endl;                    \
            return 1;                                                     \
        }                                                                 \
    } while (0)

template <class DomainP>
std::int64_t phase_star(const TFHEpp::metapbs::LWE<DomainP, std::int64_t> &c,
                        const TFHEpp::Key<DomainP> &key)
{
    std::int64_t acc = c.c[DomainP::k * DomainP::n];
    for (int i = 0; i < DomainP::k * DomainP::n; i++)
        acc -= c.c[i] * static_cast<std::int64_t>(key[i]);
    return acc;
}

}  // namespace

int main(int argc, char **argv)
{
    std::uint64_t failrate_trials = 0;
    bool failrate_full_domain = false;
    for (int i = 1; i < argc; i++) {
        const std::string arg(argv[i]);
        if (arg == "--failrate" && i + 1 < argc) {
            failrate_trials = std::stoull(argv[++i]);
        }
        else if (arg == "--full-domain") {
            failrate_full_domain = true;
        }
    }

    using domainP = TFHEpp::lvl0param;
    using targetP = TFHEpp::lvl1param;
    using brP = TFHEpp::lvl01param;

    constexpr std::int64_t q_quo = 2 * targetP::n;  // 2N

    TFHEpp::SecretKey sk;

    std::random_device seed_gen;
    std::mt19937 rng(seed_gen());
    std::uniform_int_distribution<std::uint32_t> msg_dist(0, q_quo - 1);

    for (int test = 0; test < 100; test++) {
        const std::uint32_t m = msg_dist(rng);

        TFHEpp::TLWE<domainP> c{};
        TFHEpp::tlweSymIntEncrypt<domainP, static_cast<std::uint32_t>(q_quo)>(
            c, m, sk);

        const auto lifted = TFHEpp::metapbs::LiftTLWEToInt<domainP>(c);
        const auto [cquo, crem] = TFHEpp::metapbs::HomDivRem<domainP>(lifted, q_quo);

        const std::int64_t psi = phase_star<domainP>(lifted, sk.key.get<domainP>());
        const std::int64_t psi_quo = phase_star<domainP>(cquo, sk.key.get<domainP>());
        const std::int64_t psi_rem = phase_star<domainP>(crem, sk.key.get<domainP>());

        const std::int64_t q = lifted.modulus;
        const std::int64_t q_rem = q / q_quo;
        CHECK(psi == q_rem * psi_quo + psi_rem);

        for (const auto &x : crem.c) {
            CHECK(x >= -q_rem / 2);
            CHECK(x <= (q_rem - 1) / 2);
        }
    }

    std::cout << "HomDivRem: Passed" << std::endl;

    // BlindRotate sanity: if the domain key is all-zero, the bootstrapping key
    // encrypts 0 and the message should only be rotated by b (CMUXs still add
    // noise but should not change the plaintext).
    {
        TFHEpp::Key<brP::domainP> domainkey{};
        domainkey.fill(0);
        const TFHEpp::Key<brP::targetP> targetkey = sk.key.get<brP::targetP>();

        auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
        TFHEpp::bkfftgen<brP>(*bkfft, domainkey, targetkey);

        TFHEpp::Polynomial<brP::targetP> tv{};
        for (int i = 0; i < brP::targetP::n; i++) {
            tv[i] = (i & 1) ? static_cast<brP::targetP::T>(-brP::targetP::μ)
                            : static_cast<brP::targetP::T>(brP::targetP::μ);
        }

        TFHEpp::metapbs::LWE<brP::domainP> tlwe_quo{};
        tlwe_quo.modulus = q_quo;
        // Arbitrary coefficients; secret is zero so they should have no effect.
        for (int i = 0; i < brP::domainP::k * brP::domainP::n; i++)
            tlwe_quo.c[i] = static_cast<std::int64_t>((i * 17) % (2 * brP::targetP::n));
        tlwe_quo.c[brP::domainP::k * brP::domainP::n] = 123;

        const int bbar = TFHEpp::metapbs::mod2N<brP::targetP, std::int64_t>(
            -tlwe_quo.c[brP::domainP::k * brP::domainP::n]);

        TFHEpp::Polynomial<brP::targetP> expected_tv{};
        TFHEpp::PolynomialMulByXai<brP::targetP>(
            expected_tv, tv, static_cast<brP::targetP::T>(bbar));

        alignas(64) TFHEpp::TRLWE<brP::targetP> out_tv;
        TFHEpp::metapbs::BlindRotate<brP>(out_tv, tlwe_quo, *bkfft, tv);
        const auto dec_tv = TFHEpp::trlweSymDecrypt<brP::targetP>(out_tv, targetkey);
        for (int i = 0; i < brP::targetP::n; i++) {
            const bool exp =
                static_cast<std::make_signed_t<brP::targetP::T>>(expected_tv[i]) >
                0;
            CHECK(dec_tv[i] == exp);
        }

        alignas(64) TFHEpp::TRLWE<brP::targetP> tv_ct;
        TFHEpp::trlweSymEncrypt<brP::targetP>(tv_ct, tv, targetkey);

        alignas(64) TFHEpp::TRLWE<brP::targetP> out_ct;
        TFHEpp::metapbs::BlindRotate<brP>(out_ct, tlwe_quo, *bkfft, tv_ct);
        const auto dec_ct = TFHEpp::trlweSymDecrypt<brP::targetP>(out_ct, targetkey);
        for (int i = 0; i < brP::targetP::n; i++) {
            const bool exp =
                static_cast<std::make_signed_t<brP::targetP::T>>(expected_tv[i]) >
                0;
            CHECK(dec_ct[i] == exp);
        }
    }

    std::cout << "MetaBlindRotate: Passed" << std::endl;

    // TruncRepeat correctness on an integer torus message (stable even when the
    // expected coefficient is 0).
    {
        using P = targetP;
        const TFHEpp::Key<P> key = sk.key.get<P>();

        auto ahk = std::make_unique<TFHEpp::AnnihilateKey<P>>();
        TFHEpp::annihilatekeygen<P>(*ahk, sk);

        constexpr int T = 5;
        constexpr int B = 4;
        constexpr std::uint32_t plain_modulus = 1u << 8;  // must divide 2^32
        constexpr int shift = 32 - 8;                     // 2^32 / plain_modulus

        std::array<std::uint32_t, P::n> m_int{};
        for (int i = 0; i < P::n; i++) m_int[i] = (i * 7 + 3) & (plain_modulus - 1);

        TFHEpp::Polynomial<P> m{};
        for (int i = 0; i < P::n; i++)
            m[i] = static_cast<typename P::T>(m_int[i] << shift);

        TFHEpp::TRLWE<P> ct;
        TFHEpp::trlweSymEncrypt<P>(ct, m, key);

        TFHEpp::TRLWE<P> ct_tr;
        TFHEpp::metapbs::TruncRepeat<P>(ct_tr, ct, T, B, *ahk);

        const auto dec =
            TFHEpp::trlweSymIntDecrypt<P, plain_modulus>(ct_tr, key);

        auto neg_mod = [](std::uint32_t x) -> std::uint32_t {
            return x ? (plain_modulus - x) : 0u;
        };
        auto mul_by_Xai = [&](const std::array<std::uint32_t, P::n> &poly,
                              const int a) {
            std::array<std::uint32_t, P::n> res{};
            if (a == 0) return poly;
            if (a < static_cast<int>(P::n)) {
                for (int i = 0; i < a; i++) res[i] = neg_mod(poly[i - a + P::n]);
                for (int i = a; i < static_cast<int>(P::n); i++)
                    res[i] = poly[i - a];
            }
            else {
                const int aa = a - static_cast<int>(P::n);
                for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n];
                for (int i = aa; i < static_cast<int>(P::n); i++)
                    res[i] = neg_mod(poly[i - aa]);
            }
            return res;
        };

        // Plain reference: truncPad then multiply by sum_{k in [B]_sym} X^k.
        std::array<std::uint32_t, P::n> pad{};
        for (int d = -T; d <= T; d++) {
            const std::uint32_t cd =
                (d >= 0) ? m_int[d] : neg_mod(m_int[static_cast<int>(P::n) + d]);
            const int D = d * B;
            if (D >= 0)
                pad[D] = cd;
            else
                pad[static_cast<int>(P::n) + D] = neg_mod(cd);
        }

        const int kmin = TFHEpp::metapbs::sym_min(B);
        const int kmax = TFHEpp::metapbs::sym_max(B);
        std::array<std::uint32_t, P::n> exp{};
        for (int k = kmin; k <= kmax; k++) {
            const int a = TFHEpp::metapbs::mod2N<P, int>(k);
            const auto rot = mul_by_Xai(pad, a);
            for (int i = 0; i < P::n; i++)
                exp[i] = (exp[i] + rot[i]) & (plain_modulus - 1);
        }

        for (int i = 0; i < P::n; i++) CHECK(dec[i] == exp[i]);
    }

    std::cout << "TruncRepeat: Passed" << std::endl;

    // End-to-end Meta-PBS (Algorithm 1) for a small plaintext space (<= 5 bits).
    // Use a folded-identity negacyclic LUT:
    //   f(x) = x for x in [0,t/2), and f(x+t/2) = -f(x) (mod t).
    {
        using bkP = TFHEpp::lvlh1param;  // 32-bit domain for higher-precision Meta-PBS
        constexpr std::int64_t t = 32;  // 5-bit plaintext space
        static_assert((2 * bkP::targetP::n) % t == 0, "t must divide 2N");

        auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<bkP>>();
        TFHEpp::bkfftgen<bkP>(*bkfft, sk);

        auto ahk = std::make_unique<TFHEpp::AnnihilateKey<bkP::targetP>>();
        TFHEpp::annihilatekeygen<bkP::targetP>(*ahk, sk);

        std::vector<std::int64_t> f_half(t / 2);
        for (int i = 0; i < static_cast<int>(t / 2); i++) f_half[i] = i;
        const auto tv =
            TFHEpp::metapbs::GenerateTestVectorNegacyclic<bkP::targetP>(f_half, t);

        auto run_one = [&](const std::vector<int> &betas,
                           const std::vector<int> &Ts,
                           const char *tag) -> int {
            for (std::uint32_t m = 0; m < static_cast<std::uint32_t>(t); m++) {
                TFHEpp::TLWE<bkP::domainP> cin{};
                TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                         static_cast<std::uint32_t>(t)>(cin, m,
                                                                       sk);

                TFHEpp::TLWE<bkP::targetP> cout{};
                TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv, betas,
                                              Ts, t);

                const auto dec_u = TFHEpp::tlweSymIntDecrypt<
                    bkP::targetP, static_cast<std::uint32_t>(t)>(cout, sk);
                const std::int32_t dec =
                    static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                        dec_u);

                const std::int32_t exp =
                    (m < static_cast<std::uint32_t>(t / 2))
                        ? static_cast<std::int32_t>(m)
                        : -static_cast<std::int32_t>(m - (t / 2));

                if (dec != exp) {
                    constexpr std::int64_t Delta =
                        (static_cast<std::int64_t>(1) << 32) /
                        t;  // q_target / t
                    const auto phase_u = TFHEpp::tlweSymPhase<bkP::targetP>(
                        cout, sk.key.get<bkP::targetP>());
                    const std::int64_t phase =
                        static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                            phase_u);
                    const std::int64_t exp_torus =
                        static_cast<std::int64_t>(exp) * Delta;
                    const std::int64_t err = phase - exp_torus;
                    std::cerr << "MetaPBS mismatch(" << tag << "): m=" << m
                              << " dec=" << dec << " exp=" << exp
                              << " phase_err=" << err << " (~"
                              << (static_cast<double>(err) / Delta) << " ulp)"
                              << std::endl;
                    return 1;
                }
            }
            return 0;
        };

        // For t=32, the initial redundancy r0 = 2N/t = 64 is already large, so
        // K=0 (no TruncRepeat iterations) already works well.
        if (run_one({}, {}, "K=0") != 0) return 1;

        // Also test the iterative path (K=1) to exercise TruncRepeat+BlindRotate.
        if (run_one({4}, {63}, "K=1") != 0) return 1;

        // Paper claim sanity: start from a non-redundant LUT (t = 2N, so r0 = 1)
        // and use the toy-pattern coefficients 0, q/2, 0, q/2, ... .
        //
        // This corresponds to a negacyclic LUT over Z_{2N} with f(i)=0 for even i
        // and f(i)=N (i.e., t/2) for odd i, so q/t*f(i) is 0 or q/2.
        {
            constexpr std::int64_t t_nr =
                2 * static_cast<std::int64_t>(bkP::targetP::n);  // 2N
            static_assert(t_nr == 2 * bkP::targetP::n);
            static_assert((2 * bkP::targetP::n) / t_nr == 1,
                          "t=2N should give r0=1 (non-redundant LUT)");

            TFHEpp::Polynomial<bkP::targetP> tv_nr{};
            tv_nr.fill(0);
            constexpr auto half_q =
                static_cast<typename bkP::targetP::T>(1)
                << (std::numeric_limits<typename bkP::targetP::T>::digits - 1);
            for (int i = 0; i < static_cast<int>(bkP::targetP::n); i++)
                tv_nr[i] = (i & 1) ? half_q : 0;

            // Use a small number of iterations (K=2) to increase redundancy from
            // r0=1. (Toy example in the paper shows K=1 can already work for
            // certain noise bounds.)
            const std::vector<int> betas = {8, 8};
            const std::vector<int> Ts = {63, 63};  // (2T+1)*beta <= N

            // Only test messages within 5 bits (0..31) as requested.
            for (std::uint32_t m = 0; m < 32; m++) {
                TFHEpp::TLWE<bkP::domainP> cin{};
                TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                         static_cast<std::uint32_t>(t_nr)>(
                    cin, m, sk);

                TFHEpp::TLWE<bkP::targetP> cout{};
                TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv_nr,
                                              betas, Ts, t_nr);

                // Decode as a 1-bit output (closer to 0 or q/2).
                // Expected: 0 for even m, q/2 for odd m.
                const auto phase_u = TFHEpp::tlweSymPhase<bkP::targetP>(
                    cout, sk.key.get<bkP::targetP>());
                const bool exp_one = (m & 1) != 0;

                // Classify by circular distance on the torus:
                // closer to 0 iff phase_u in [0,q/4) U [3q/4,q),
                // closer to q/2 iff phase_u in [q/4, 3q/4).
                constexpr std::uint32_t q_over_4 = 1U << 30;       // 2^(32-2)
                constexpr std::uint32_t three_q_over_4 = 3U << 30; // 3*2^(32-2)
                const bool got_one =
                    (phase_u >= q_over_4) && (phase_u < three_q_over_4);
                if (got_one != exp_one) {
                    std::cerr << "Non-redundant TV mismatch: m=" << m
                              << " exp=" << (exp_one ? "q/2" : "0")
                              << " phase_u=" << phase_u << std::endl;
                    return 1;
                }
            }

            // Extra boundary-case stress: fix m = 2N-1 (2047 for N=1024) and run
            // multiple encryptions to sample randomness/noise.
            constexpr std::uint64_t fixed_reps = 200;
            constexpr std::uint32_t fixed_m =
                static_cast<std::uint32_t>(t_nr - 1);  // 2N-1
            for (std::uint64_t rep = 0; rep < fixed_reps; rep++) {
                TFHEpp::TLWE<bkP::domainP> cin{};
                TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                         static_cast<std::uint32_t>(t_nr)>(
                    cin, fixed_m, sk);

                TFHEpp::TLWE<bkP::targetP> cout{};
                TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv_nr,
                                              betas, Ts, t_nr);

                const auto phase_u = TFHEpp::tlweSymPhase<bkP::targetP>(
                    cout, sk.key.get<bkP::targetP>());
                const bool exp_one = (fixed_m & 1U) != 0;

                constexpr std::uint32_t q_over_4 = 1U << 30;
                constexpr std::uint32_t three_q_over_4 = 3U << 30;
                const bool got_one =
                    (phase_u >= q_over_4) && (phase_u < three_q_over_4);
                if (got_one != exp_one) {
                    std::cerr << "Non-redundant TV mismatch (fixed m=2N-1): rep="
                              << rep << " m=" << fixed_m
                              << " exp=" << (exp_one ? "q/2" : "0")
                              << " phase_u=" << phase_u << std::endl;
                    return 1;
                }
            }

            std::cout << "MetaPBS (t=2N, non-redundant TV 0/q2): Passed"
                      << std::endl;
        }

        if (failrate_trials > 0) {
            auto measure_folded_identity = [&](const std::vector<int> &betas,
                                               const std::vector<int> &Ts,
                                               const char *tag) {
                std::uniform_int_distribution<std::uint32_t> mdist(
                    0, static_cast<std::uint32_t>(t - 1));
                std::uint64_t fails = 0;
                for (std::uint64_t rep = 0; rep < failrate_trials; rep++) {
                    const std::uint32_t m = mdist(rng);
                    TFHEpp::TLWE<bkP::domainP> cin{};
                    TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                             static_cast<std::uint32_t>(t)>(
                        cin, m, sk);

                    TFHEpp::TLWE<bkP::targetP> cout{};
                    TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv,
                                                  betas, Ts, t);

                    const auto dec_u = TFHEpp::tlweSymIntDecrypt<
                        bkP::targetP, static_cast<std::uint32_t>(t)>(cout, sk);
                    const std::int32_t dec =
                        static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                            dec_u);

                    const std::int32_t exp =
                        (m < static_cast<std::uint32_t>(t / 2))
                            ? static_cast<std::int32_t>(m)
                            : -static_cast<std::int32_t>(m - (t / 2));

                    if (dec != exp) fails++;
                }
                const double rate = static_cast<double>(fails) /
                                    static_cast<double>(failrate_trials);
                std::cout << "FailRate(" << tag << ", t=32): " << fails << "/"
                          << failrate_trials << " = " << rate;
                if (fails == 0) {
                    std::cout << " (95% upper bound ~" << (3.0 / failrate_trials)
                              << ")";
                }
                std::cout << std::endl;
            };

            auto measure_nonredundant_alt_tv = [&]() {
                constexpr std::int64_t t_nr =
                    2 * static_cast<std::int64_t>(bkP::targetP::n);  // 2N

                // 0, q/2, 0, q/2, ... (non-redundant TV for t=2N)
                TFHEpp::Polynomial<bkP::targetP> tv_nr{};
                tv_nr.fill(0);
                constexpr auto half_q =
                    static_cast<typename bkP::targetP::T>(1)
                    << (std::numeric_limits<typename bkP::targetP::T>::digits - 1);
                for (int i = 0; i < static_cast<int>(bkP::targetP::n); i++)
                    tv_nr[i] = (i & 1) ? half_q : 0;

                const std::vector<int> betas = {8, 8};
                const std::vector<int> Ts = {63, 63};

                std::uniform_int_distribution<std::uint32_t> mdist(
                    0,
                    failrate_full_domain ? static_cast<std::uint32_t>(t_nr - 1)
                                         : 31U);
                std::uint64_t fails = 0;

                for (std::uint64_t rep = 0; rep < failrate_trials; rep++) {
                    const std::uint32_t m = mdist(rng);
                    TFHEpp::TLWE<bkP::domainP> cin{};
                    TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                             static_cast<std::uint32_t>(t_nr)>(
                        cin, m, sk);

                    TFHEpp::TLWE<bkP::targetP> cout{};
                    TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv_nr,
                                                  betas, Ts, t_nr);

                    const auto phase_u = TFHEpp::tlweSymPhase<bkP::targetP>(
                        cout, sk.key.get<bkP::targetP>());
                    const bool exp_one = (m & 1U) != 0;

                    constexpr std::uint32_t q_over_4 = 1U << 30;
                    constexpr std::uint32_t three_q_over_4 = 3U << 30;
                    const bool got_one =
                        (phase_u >= q_over_4) && (phase_u < three_q_over_4);
                    if (got_one != exp_one) fails++;
                }

                const double rate = static_cast<double>(fails) /
                                    static_cast<double>(failrate_trials);
                std::cout
                    << "FailRate(nonredundant 0/q2, t=2N, "
                    << (failrate_full_domain ? "m in [0,2N)" : "m in [0,32)")
                    << "): " << fails << "/" << failrate_trials << " = " << rate;
                if (fails == 0) {
                    std::cout << " (95% upper bound ~" << (3.0 / failrate_trials)
                              << ")";
                }
                std::cout << std::endl;
            };

            // Measure on the same key material as the correctness tests.
            measure_folded_identity({}, {}, "K=0");
            measure_folded_identity({4}, {63}, "K=1");
            measure_nonredundant_alt_tv();
        }
    }

    std::cout << "MetaPBS: Passed" << std::endl;
    return 0;
}
