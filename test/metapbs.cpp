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

        // Random LUTs (t=32): periodic and non-periodic.
        // These are still *negacyclic* LUTs (f(x+t/2) = -f(x)), but the first-half
        // values are chosen pseudo-randomly to exercise the generic path.
        auto check_random_lut_t32 =
            [&](const std::vector<std::int64_t> &f_half_lut,
                const char *tag) -> int {
            const auto tv_lut =
                TFHEpp::metapbs::GenerateTestVectorNegacyclic<bkP::targetP>(
                    f_half_lut, t);
            for (std::uint32_t m = 0; m < static_cast<std::uint32_t>(t); m++) {
                TFHEpp::TLWE<bkP::domainP> cin{};
                TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                         static_cast<std::uint32_t>(t)>(cin, m,
                                                                       sk);

                TFHEpp::TLWE<bkP::targetP> cout{};
                TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv_lut,
                                              {}, {}, t);

                const auto dec_u = TFHEpp::tlweSymIntDecrypt<
                    bkP::targetP, static_cast<std::uint32_t>(t)>(cout, sk);
                const std::int32_t dec =
                    static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                        dec_u);

                std::int64_t exp64 = 0;
                if (m < static_cast<std::uint32_t>(t / 2))
                    exp64 = f_half_lut[m];
                else
                    exp64 = -f_half_lut[m - static_cast<std::uint32_t>(t / 2)];
                const std::int32_t exp = static_cast<std::int32_t>(exp64);

                if (dec != exp) {
                    std::cerr << "Random LUT mismatch(" << tag << "): m=" << m
                              << " dec=" << dec << " exp=" << exp << std::endl;
                    return 1;
                }
            }
            std::cout << "MetaPBS (t=32, random LUT " << tag << "): Passed"
                      << std::endl;
            return 0;
        };

        {
            // Periodic LUT: pick a short random pattern and repeat it.
            std::mt19937 lut_rng(0x4D505253U);  // "MPRS"
            std::uniform_int_distribution<std::int32_t> vdist(
                -static_cast<std::int32_t>(t / 2 - 1),
                static_cast<std::int32_t>(t / 2 - 1));

            constexpr std::uint32_t period = 4;
            std::vector<std::int64_t> pat(period);
            for (std::uint32_t i = 0; i < period; i++) pat[i] = vdist(lut_rng);

            std::vector<std::int64_t> f_half_lut(t / 2);
            for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(t / 2); i++)
                f_half_lut[i] = pat[i % period];

            if (check_random_lut_t32(f_half_lut, "periodic") != 0) return 1;
        }
        {
            // Non-periodic LUT: each entry is sampled independently.
            std::mt19937 lut_rng(0x4D505254U);  // "MPRT"
            std::uniform_int_distribution<std::int32_t> vdist(
                -static_cast<std::int32_t>(t / 2 - 1),
                static_cast<std::int32_t>(t / 2 - 1));

            std::vector<std::int64_t> f_half_lut(t / 2);
            for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(t / 2); i++)
                f_half_lut[i] = vdist(lut_rng);

            if (check_random_lut_t32(f_half_lut, "non-periodic") != 0) return 1;
        }

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

            // Idea sanity: periodic-invariance pruning inside BlindRotate.
            // For the 2-periodic LUT (LSB), Rot_{a}(TV)=TV whenever a is even.
            {
                TFHEpp::Polynomial<bkP::targetP> tv_rot{};
                TFHEpp::PolynomialMulByXai<bkP::targetP>(
                    tv_rot, tv_nr,
                    static_cast<typename bkP::targetP::T>(2));  // X^2
                for (int i = 0; i < static_cast<int>(bkP::targetP::n); i++)
                    CHECK(tv_rot[i] == tv_nr[i]);

                TFHEpp::metapbs::BlindRotatePruneStats stats{};
                constexpr std::uint32_t prune_trials = 4;
                std::uniform_int_distribution<std::uint32_t> mdist_full(
                    0, static_cast<std::uint32_t>(t_nr - 1));
                for (std::uint32_t rep = 0; rep < prune_trials; rep++) {
                    const std::uint32_t m = mdist_full(rng);
                    TFHEpp::TLWE<bkP::domainP> cin{};
                    TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                             static_cast<std::uint32_t>(t_nr)>(
                        cin, m, sk);

                    const auto lifted =
                        TFHEpp::metapbs::LiftTLWEToInt<bkP::domainP>(cin);
                    auto [cquo, crem] = TFHEpp::metapbs::HomDivRem<bkP::domainP>(
                        lifted, static_cast<std::int64_t>(t_nr));
                    (void)crem;

                    alignas(64) TFHEpp::TRLWE<bkP::targetP> out_ref{};
                    TFHEpp::metapbs::BlindRotate<bkP>(out_ref, cquo, *bkfft,
                                                      tv_nr);

                    alignas(64) TFHEpp::TRLWE<bkP::targetP> out_pruned{};
                    TFHEpp::metapbs::BlindRotatePeriodic<bkP>(
                        out_pruned, cquo, *bkfft, tv_nr, /*period=*/2, &stats);

                    const auto phase_ref = TFHEpp::trlwePhase<bkP::targetP>(
                        out_ref, sk.key.get<bkP::targetP>());
                    const auto phase_pruned = TFHEpp::trlwePhase<bkP::targetP>(
                        out_pruned, sk.key.get<bkP::targetP>());

                    constexpr std::uint32_t q_over_4 = 1U << 30;
                    constexpr std::uint32_t three_q_over_4 = 3U << 30;
                    for (int i = 0; i < static_cast<int>(bkP::targetP::n); i++) {
                        const bool b0 =
                            (phase_ref[i] >= q_over_4) &&
                            (phase_ref[i] < three_q_over_4);
                        const bool b1 =
                            (phase_pruned[i] >= q_over_4) &&
                            (phase_pruned[i] < three_q_over_4);
                        CHECK(b0 == b1);
                    }
                }
                std::cout << "MetaBlindRotate (periodic-prune M=2): skipped "
                          << stats.skipped << "/" << stats.total << std::endl;
            }

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

            // Interface test: MetaPBSExtractBit2N(bit) for a few bit positions.
            {
                TFHEpp::metapbs::BlindRotatePruneStats stats{};
                std::vector<std::uint32_t> msgs;
                msgs.reserve(64);
                msgs.push_back(0);
                msgs.push_back(1);
                msgs.push_back(2);
                msgs.push_back(3);
                msgs.push_back(31);
                msgs.push_back(32);
                msgs.push_back(static_cast<std::uint32_t>(t_nr / 2 - 1));
                msgs.push_back(static_cast<std::uint32_t>(t_nr / 2));
                msgs.push_back(static_cast<std::uint32_t>(t_nr - 2));
                msgs.push_back(static_cast<std::uint32_t>(t_nr - 1));
                std::mt19937 msg_rng(0x42495455U);  // "BITU"
                std::uniform_int_distribution<std::uint32_t> mdist_full(
                    0, static_cast<std::uint32_t>(t_nr - 1));
                for (int i = 0; i < 16; i++) msgs.push_back(mdist_full(msg_rng));

                const std::array<int, 3> bits_to_test = {0, 3, 9};
                for (const int bit : bits_to_test) {
                    for (const std::uint32_t m : msgs) {
                        TFHEpp::TLWE<bkP::domainP> cin{};
                        TFHEpp::tlweSymIntEncrypt<
                            bkP::domainP, static_cast<std::uint32_t>(t_nr)>(
                            cin, m, sk);

                        TFHEpp::TLWE<bkP::targetP> cout{};
                        TFHEpp::metapbs::MetaPBSExtractBit2N<bkP>(
                            cout, cin, *bkfft, *ahk, bit, &stats);

                        const auto phase_u = TFHEpp::tlweSymPhase<bkP::targetP>(
                            cout, sk.key.get<bkP::targetP>());
                        constexpr std::uint32_t q_over_4 = 1U << 30;
                        constexpr std::uint32_t three_q_over_4 = 3U << 30;
                        const bool got_one =
                            (phase_u >= q_over_4) && (phase_u < three_q_over_4);
                        const bool exp_one = ((m >> bit) & 1U) != 0;
                        CHECK(got_one == exp_one);
                    }
                }

                std::cout << "MetaPBSExtractBit2N: Passed (pruned skipped "
                          << stats.skipped << "/" << stats.total << ")"
                          << std::endl;
            }

            // High-level interface: ExtractBitInPlace(bit indexed from MSB=0).
            //
            // For a 10-bit message m in [0,1023], extracting bit_msb returns an
            // LWE encrypting either 0 or (1<<bit_lsb) in Z_{2N} (2N=2048),
            // i.e., the bit is "put back" into its original position.
            {
                // High-precision weighting step:
                //   lvl1 (bit half-turn) -> lvlhalf (KS) -> lvl2 (PBS) -> lvl1 (KS)
                // to make decoding under 2N=2048 reliable.
                using ksToDomP = TFHEpp::lvl1hparam;   // lvl1 -> lvlhalf
                using weightBkP = TFHEpp::lvlh2param;  // lvlhalf -> lvl2
                using ksDownP = TFHEpp::lvl21param;    // lvl2 -> lvl1

                auto ksk1h = std::make_unique<TFHEpp::KeySwitchingKey<ksToDomP>>();
                TFHEpp::ikskgen<ksToDomP>(*ksk1h, sk);

                auto bkfft_h2 =
                    std::make_unique<TFHEpp::BootstrappingKeyFFT<weightBkP>>();
                TFHEpp::bkfftgen<weightBkP>(*bkfft_h2, sk);

                auto ksk21 = std::make_unique<TFHEpp::KeySwitchingKey<ksDownP>>();
                TFHEpp::ikskgen<ksDownP>(*ksk21, sk);

                std::vector<std::uint32_t> msgs;
                msgs.reserve(64);
                // Fixed edge cases within 10 bits.
                msgs.push_back(0);
                msgs.push_back(1);
                msgs.push_back(2);
                msgs.push_back(3);
                msgs.push_back(15);
                msgs.push_back(16);
                msgs.push_back(31);
                msgs.push_back(32);
                msgs.push_back(511);
                msgs.push_back(512);
                msgs.push_back(1023);

                std::mt19937 msg_rng(0x45584254U);  // "EXBT"
                std::uniform_int_distribution<std::uint32_t> mdist_10bit(0, 1023);
                for (int i = 0; i < 32; i++) msgs.push_back(mdist_10bit(msg_rng));

                for (const std::uint32_t m : msgs) {
                    TFHEpp::TLWE<bkP::domainP> cin{};
                    TFHEpp::tlweSymIntEncrypt<
                        bkP::domainP, static_cast<std::uint32_t>(t_nr)>(cin, m,
                                                                      sk);

                    for (int bit_msb = 0;
                         bit_msb < static_cast<int>(bkP::targetP::nbit);
                         bit_msb++) {
                        const int bit_lsb =
                            (static_cast<int>(bkP::targetP::nbit) - 1) - bit_msb;
                        const std::uint32_t weight = 1u << bit_lsb;

                        TFHEpp::TLWE<bkP::targetP> cout{};
                        TFHEpp::metapbs::ExtractBitInPlaceViaLvl2<
                            bkP, weightBkP, ksToDomP, ksDownP>(
                            cout, cin, *bkfft, *ahk, *ksk1h, *bkfft_h2, *ksk21,
                            bit_msb);

                        const auto dec_u = TFHEpp::tlweSymIntDecrypt<
                            bkP::targetP,
                            static_cast<std::uint32_t>(t_nr)>(cout, sk);
                        const std::int32_t dec =
                            static_cast<std::make_signed_t<
                                typename bkP::targetP::T>>(dec_u);
                        const std::int32_t exp =
                            ((m >> bit_lsb) & 1U)
                                ? static_cast<std::int32_t>(weight)
                                : 0;

                        if (dec != exp) {
                            std::cerr << "ExtractBitInPlace mismatch: m=" << m
                                      << " bit_msb=" << bit_msb
                                      << " bit_lsb=" << bit_lsb
                                      << " dec=" << dec << " exp=" << exp
                                      << std::endl;
                            // Debug the 3-stage pipeline:
                            //  (1) MetaPBSExtractBit2N -> {0,q/2}
                            //  (2) rescale by 2^{-nbit} and negate -> {0,Δ}
                            //  (3) multiply by weight -> {0,weight}
                            {
                                using TargetP = typename bkP::targetP;
                                constexpr int digits =
                                    std::numeric_limits<typename TargetP::T>::digits;
                                constexpr int delta_shift =
                                    digits - (static_cast<int>(TargetP::nbit) + 1);
                                const auto w_half =
                                    static_cast<typename TargetP::T>(
                                        static_cast<typename TargetP::T>(weight)
                                        << (delta_shift - 1));

                                TFHEpp::TLWE<TargetP> bit_ct_dbg{};
                                TFHEpp::metapbs::MetaPBSExtractBit2N<bkP>(
                                    bit_ct_dbg, cin, *bkfft, *ahk, bit_lsb);
                                const auto bit_phase_u =
                                    TFHEpp::tlweSymPhase<TargetP>(
                                        bit_ct_dbg, sk.key.get<TargetP>());
                                const std::int64_t bit_phase =
                                    static_cast<std::make_signed_t<typename TargetP::T>>(
                                        bit_phase_u);
                                const auto bit_dec_u = TFHEpp::tlweSymIntDecrypt<
                                    TargetP,
                                    static_cast<std::uint32_t>(t_nr)>(bit_ct_dbg,
                                                                     sk);
                                const std::int32_t bit_dec =
                                    static_cast<std::make_signed_t<typename TargetP::T>>(
                                        bit_dec_u);

                                std::cerr << "  dbg: weight=" << weight
                                          << " w_half(int)="
                                          << static_cast<std::int64_t>(w_half >>
                                                                       delta_shift)
                                          << "\n"
                                          << "  dbg: bit_ct phase=" << bit_phase
                                          << " dec_mod2N=" << bit_dec << "\n";
                            }
                            return 1;
                        }
                    }
                }

                // 16-bit input adaptation idea:
                // If a 10-bit value x is embedded into a 16-bit plaintext by
                // shifting left by (16-(nbit+1)) bits, then its torus encoding
                // matches the mod-2N (2N=2048) encoding:
                //   (x << embed_shift) * 2^(32-16) == x * 2^(32-(nbit+1)).
                // This lets MetaPBS "see" x (no rounding) while the ciphertext
                // is still a valid 16-bit TLWE integer ciphertext.
                {
                    constexpr std::uint32_t t16 = 1u << 16;
                    constexpr int embed_shift =
                        16 - (static_cast<int>(bkP::targetP::nbit) + 1);
                    static_assert(embed_shift >= 0,
                                  "16-bit embedding requires 16 >= nbit+1");

                    std::mt19937 rng16(0x31364254U);  // "16BT"
                    std::uniform_int_distribution<std::uint32_t> xdist(
                        0, (1u << bkP::targetP::nbit) - 1);

                    for (int bit_msb = 0;
                         bit_msb < static_cast<int>(bkP::targetP::nbit);
                         bit_msb++) {
                        const int bit_lsb =
                            (static_cast<int>(bkP::targetP::nbit) - 1) - bit_msb;
                        for (int rep = 0; rep < 100; rep++) {
                            const std::uint32_t x = xdist(rng16);
                            const std::uint32_t m16 = x << embed_shift;

                            TFHEpp::TLWE<bkP::domainP> cin16{};
                            TFHEpp::tlweSymIntEncrypt<bkP::domainP, t16>(cin16, m16,
                                                                        sk);

                            TFHEpp::TLWE<bkP::targetP> cout16{};
                            TFHEpp::metapbs::ExtractBitInPlaceViaLvl2<
                                bkP, weightBkP, ksToDomP, ksDownP>(
                                cout16, cin16, *bkfft, *ahk, *ksk1h, *bkfft_h2,
                                *ksk21, bit_msb);

                            const auto dec_u =
                                TFHEpp::tlweSymIntDecrypt<bkP::targetP,
                                                         static_cast<std::uint32_t>(
                                                             t_nr)>(cout16, sk);
                            const std::uint32_t dec = dec_u;
                            const std::uint32_t exp =
                                ((x >> bit_lsb) & 1U)
                                    ? (1u << bit_lsb)
                                    : 0u;

                            if (dec != exp) {
                                std::cerr
                                    << "ExtractBitInPlace(16-bit embed) mismatch:"
                                    << " x=" << x << " m16=" << m16
                                    << " bit_msb=" << bit_msb
                                    << " bit_lsb=" << bit_lsb << " dec=" << dec
                                    << " exp=" << exp << std::endl;
                                return 1;
                            }
                        }
                    }

                    std::cout << "ExtractBitInPlace (16-bit embedded): Passed"
                              << std::endl;
                }

                // Experiment / fail-rate measurement:
                // Start from a random 16-bit ciphertext at lvl1 (mod 2^16),
                // rescale it to the mod-2N encoding by a torus left shift,
                // then key-switch to lvlhalf to feed MetaPBS extraction,
                // and finally subtract the extracted (bit*weight) ciphertext to
                // "clear that bit" in Z_{2N}.
                //
                // This is a noisy pipeline because the lvl1->lvlhalf key switch
                // happens on an 11-bit plaintext space (2N=2048). Use --failrate
                // to quantify the error probability.
                if (failrate_trials > 0) {
                    constexpr std::uint32_t t16 = 1u << 16;
                    constexpr int embed_shift =
                        16 - (static_cast<int>(bkP::targetP::nbit) + 1);
                    static_assert(embed_shift >= 0,
                                  "16-bit embedding requires 16 >= nbit+1");

                    std::mt19937 rng_clr(0x434C5232U);  // "CLR2"
                    std::uniform_int_distribution<std::uint32_t> mdist16(0,
                                                                         t16 - 1);
                    std::uniform_int_distribution<int> bdist(
                        0, static_cast<int>(bkP::targetP::nbit) - 1);

                    std::uint64_t fails = 0;
                    constexpr std::uint64_t max_print = 5;

                    for (std::uint64_t rep = 0; rep < failrate_trials; rep++) {
                        const std::uint32_t m16 = mdist16(rng_clr);
                        const int bit_msb = bdist(rng_clr);
                        const int bit_lsb =
                            (static_cast<int>(bkP::targetP::nbit) - 1) - bit_msb;
                        const std::uint32_t weight = 1u << bit_lsb;

                        // ctI: random 16-bit plaintext, encrypted at lvl1.
                        TFHEpp::TLWE<bkP::targetP> ctI{};
                        TFHEpp::tlweSymIntEncrypt<bkP::targetP, t16>(ctI, m16, sk);

                        // ct1_lvl1: rescale 16-bit encoding -> mod 2N encoding
                        // (same integer value modulo 2N).
                        TFHEpp::TLWE<bkP::targetP> ct1_lvl1 = ctI;
                        for (auto &x : ct1_lvl1)
                            x = static_cast<typename bkP::targetP::T>(
                                x << embed_shift);

                        // Key switch to the MetaPBS extract domain (lvlhalf).
                        TFHEpp::TLWE<bkP::domainP> ct1_dom{};
                        TFHEpp::IdentityKeySwitch<ksToDomP>(ct1_dom, ct1_lvl1,
                                                           *ksk1h);

                        // ct2_lvl1: extracted bit put back in-place (returned at lvl1),
                        // encoded for modulus 2N.
                        TFHEpp::TLWE<bkP::targetP> ct2_lvl1{};
                        TFHEpp::metapbs::ExtractBitInPlaceViaLvl2<
                            bkP, weightBkP, ksToDomP, ksDownP>(
                            ct2_lvl1, ct1_dom, *bkfft, *ahk, *ksk1h, *bkfft_h2,
                            *ksk21, bit_msb);

                        // Subtract at lvl1 in the mod-2N encoding.
                        TFHEpp::TLWE<bkP::targetP> ct_clear = ct1_lvl1;
                        for (size_t i = 0; i < ct_clear.size(); i++)
                            ct_clear[i] -= ct2_lvl1[i];

                        const auto dec_u = TFHEpp::tlweSymIntDecrypt<
                            bkP::targetP, static_cast<std::uint32_t>(t_nr)>(
                            ct_clear, sk);
                        const std::int32_t dec =
                            static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                                dec_u);

                        const std::uint32_t m_mod2N =
                            m16 & (static_cast<std::uint32_t>(t_nr) - 1);
                        const std::uint32_t exp_u = m_mod2N & (~weight);
                        const std::int32_t exp =
                            (exp_u >= static_cast<std::uint32_t>(t_nr / 2))
                                ? static_cast<std::int32_t>(exp_u) -
                                      static_cast<std::int32_t>(t_nr)
                                : static_cast<std::int32_t>(exp_u);

                        if (dec != exp) {
                            fails++;
                            if (fails <= max_print) {
                                const auto ct1_lvl1_u = TFHEpp::tlweSymIntDecrypt<
                                    bkP::targetP,
                                    static_cast<std::uint32_t>(t_nr)>(ct1_lvl1,
                                                                     sk);
                                const std::int32_t ct1_lvl1_dec =
                                    static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                                        ct1_lvl1_u);

                                const auto ct1_dom_u = TFHEpp::tlweSymIntDecrypt<
                                    bkP::domainP,
                                    static_cast<std::uint32_t>(t_nr)>(ct1_dom,
                                                                     sk);
                                const std::int32_t ct1_dom_dec =
                                    static_cast<std::make_signed_t<typename bkP::domainP::T>>(
                                        ct1_dom_u);

                                const auto ct2_u = TFHEpp::tlweSymIntDecrypt<
                                    bkP::targetP,
                                    static_cast<std::uint32_t>(t_nr)>(ct2_lvl1,
                                                                     sk);
                                const std::int32_t ct2_dec =
                                    static_cast<std::make_signed_t<typename bkP::targetP::T>>(
                                        ct2_u);

                                std::cerr
                                    << "ClearBit mismatch: rep=" << rep
                                    << " m16=" << m16 << " bit_msb=" << bit_msb
                                    << " bit_lsb=" << bit_lsb << " dec=" << dec
                                    << " exp=" << exp << "\n"
                                    << "  dbg: ct1_lvl1(dec mod2N)="
                                    << ct1_lvl1_dec
                                    << " ct1_dom(dec mod2N)=" << ct1_dom_dec
                                    << " ct2(dec mod2N)=" << ct2_dec << std::endl;
                            }
                        }
                    }

                    const double rate = static_cast<double>(fails) /
                                        static_cast<double>(failrate_trials);
                    std::cout << "FailRate(ClearBit 16-bit->shift->KS->ExtractBitInPlace->sub): "
                              << fails << "/" << failrate_trials << " = " << rate;
                    if (fails == 0) {
                        std::cout << " (95% upper bound ~"
                                  << (3.0 / failrate_trials) << ")";
                    }
                    std::cout << std::endl;
                }

                std::cout << "ExtractBitInPlace: Passed" << std::endl;
            }

            // Random 1-bit LUTs (t=2N): periodic and non-periodic.
            // Coefficients are 0 or q/2, so decoding uses a wide threshold
            // (closer to 0 vs closer to q/2). This makes the test robust even
            // when decoding a full 11-bit output (mod 2N) would be too tight
            // under the current TFHEpp noise parameters.
            auto check_random_bit_lut_t2N =
                [&](const TFHEpp::Polynomial<bkP::targetP> &tv_bits,
                    const std::vector<std::uint8_t> &bits,
                    const char *tag) -> int {
                const std::uint32_t half =
                    static_cast<std::uint32_t>(t_nr / 2);  // N

                std::vector<std::uint32_t> msgs;
                msgs.reserve(64);
                // A few fixed edge cases.
                msgs.push_back(0);
                msgs.push_back(1);
                msgs.push_back(2);
                msgs.push_back(3);
                msgs.push_back(31);
                msgs.push_back(32);
                msgs.push_back(half - 1);
                msgs.push_back(half);
                msgs.push_back(half + 1);
                msgs.push_back(static_cast<std::uint32_t>(t_nr - 2));
                msgs.push_back(static_cast<std::uint32_t>(t_nr - 1));

                // Plus a small deterministic pseudo-random sample across Z_{2N}.
                std::mt19937 msg_rng(0x4D455441U);  // "META"
                std::uniform_int_distribution<std::uint32_t> mdist(
                    0, static_cast<std::uint32_t>(t_nr - 1));
                constexpr std::uint32_t random_msgs = 32;
                for (std::uint32_t i = 0; i < random_msgs; i++)
                    msgs.push_back(mdist(msg_rng));

                for (const std::uint32_t m : msgs) {
                    TFHEpp::TLWE<bkP::domainP> cin{};
                    TFHEpp::tlweSymIntEncrypt<bkP::domainP,
                                             static_cast<std::uint32_t>(t_nr)>(
                        cin, m, sk);

                    TFHEpp::TLWE<bkP::targetP> cout{};
                    TFHEpp::metapbs::MetaPBS<bkP>(cout, cin, *bkfft, *ahk, tv_bits,
                                                  betas, Ts, t_nr);

                    const auto phase_u = TFHEpp::tlweSymPhase<bkP::targetP>(
                        cout, sk.key.get<bkP::targetP>());

                    constexpr std::uint32_t q_over_4 = 1U << 30;
                    constexpr std::uint32_t three_q_over_4 = 3U << 30;
                    const bool got_one =
                        (phase_u >= q_over_4) && (phase_u < three_q_over_4);

                    // For outputs in {0, q/2}, the negacyclic extension satisfies
                    // f(x+N) = -f(x) = f(x) (since -q/2 == q/2 on the torus).
                    const bool exp_one = bits[m % half] != 0;

                    if (got_one != exp_one) {
                        std::cerr << "Random 1-bit LUT mismatch(" << tag << "): m="
                                  << m << " exp=" << (exp_one ? "q/2" : "0")
                                  << " phase_u=" << phase_u << std::endl;
                        return 1;
                    }
                }

                std::cout << "MetaPBS (t=2N, random 1-bit LUT " << tag
                          << "): Passed" << std::endl;
                return 0;
            };

            auto build_tv_from_bits =
                [&](const std::vector<std::uint8_t> &bits)
                    -> TFHEpp::Polynomial<bkP::targetP> {
                TFHEpp::Polynomial<bkP::targetP> tv_bits{};
                tv_bits.fill(0);
                for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(t_nr / 2);
                     i++)
                    tv_bits[i] = bits[i] ? half_q : 0;
                return tv_bits;
            };

            {
                // Periodic LUT: pick a short random bit-pattern and repeat it.
                const std::uint32_t half =
                    static_cast<std::uint32_t>(t_nr / 2);  // N
                std::mt19937 lut_rng(0x4D544231U);         // "MTB1"
                std::bernoulli_distribution bdist(0.5);

                constexpr std::uint32_t period = 16;
                std::vector<std::uint8_t> pat(period);
                for (std::uint32_t i = 0; i < period; i++)
                    pat[i] = bdist(lut_rng) ? 1 : 0;

                std::vector<std::uint8_t> bits(half);
                for (std::uint32_t i = 0; i < half; i++) bits[i] = pat[i % period];

                const auto tv_bits = build_tv_from_bits(bits);
                if (check_random_bit_lut_t2N(tv_bits, bits, "periodic") != 0)
                    return 1;
            }
            {
                // Non-periodic LUT: each entry is sampled independently.
                const std::uint32_t half =
                    static_cast<std::uint32_t>(t_nr / 2);  // N
                std::mt19937 lut_rng(0x4D544232U);         // "MTB2"
                std::bernoulli_distribution bdist(0.5);

                std::vector<std::uint8_t> bits(half);
                for (std::uint32_t i = 0; i < half; i++)
                    bits[i] = bdist(lut_rng) ? 1 : 0;

                const auto tv_bits = build_tv_from_bits(bits);
                if (check_random_bit_lut_t2N(tv_bits, bits, "non-periodic") != 0)
                    return 1;
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
