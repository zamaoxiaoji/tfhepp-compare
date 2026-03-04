// verify_offset_pbs.cpp
//
// Verify: encrypt a 9-bit number m (bit k = 0 in MSB-first order),
// then do a Sign PBS with offset = (2^k + 1)/2 * Δ,
// where Δ = 2^(q-p), and check whether decryption succeeds.
//
// The Sign PBS outputs true (1) iff the shifted phase (m*Δ + offset) encodes
// a negative value when viewed as a signed torus element.
//
// Build: cmake .. -DENABLE_TEST=ON && make -j verify_offset_pbs
// Run:   ./test/verify_offset_pbs

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include <tfhe++.hpp>

namespace {

// ─── Sign PBS with custom offset ────────────────────────────────────────────
//
// Replicates HE3DB's MSBGateBootstrapping with a caller-specified offset.
//
// Flow:
//   1. Add `offset` to ciphertext body
//   2. Key-switch lvl1 → lvl0
//   3. Blind-rotate lvl0 → lvl1 with test vector TV[i] = -μ (sign LUT)
//
// Output semantics (after blind rotate + sample extract):
//   The negacyclic blind rotate with TV = [ -μ, -μ, ..., -μ ] means:
//     If the modswitched phase b̄ ∈ [0, N):  output ≈ -μ → decrypt = false (0)
//     If b̄ ∈ [N, 2N):                       output ≈ +μ → decrypt = true  (1)
//
//   b̄ = 2N - round( phase * 2N / 2^q )   (mod 2N)
//
//   So b̄ ∈ [N, 2N) ⟺ round(phase * 2N / 2^q) ∈ [1, N]
//   ⟺ phase ∈ (0, 2^{q-1}]   approximately
//   i.e., the signed torus value (int32_t) phase > 0.
//
template <class iksP, class brP>
void SignPBS_CustomOffset(
    TFHEpp::TLWE<typename brP::targetP> &res,
    const TFHEpp::TLWE<typename iksP::domainP> &tlwe,
    typename iksP::domainP::T offset,
    const TFHEpp::KeySwitchingKey<iksP> &ksk,
    const TFHEpp::BootstrappingKeyFFT<brP> &bk)
{
    using DomainP = typename iksP::domainP;
    using TargetP = typename brP::targetP;
    using DomBRP  = typename brP::domainP;

    TFHEpp::TLWE<DomainP> tlwe_off = tlwe;
    tlwe_off[DomainP::k * DomainP::n] += offset;

    TFHEpp::TLWE<DomBRP> dom{};
    TFHEpp::IdentityKeySwitch<iksP>(dom, tlwe_off, ksk);

    typename TargetP::T mu = TargetP::μ;
    TFHEpp::Polynomial<TargetP> tv{};
    for (auto &p : tv) p = static_cast<typename TargetP::T>(-mu);

    TFHEpp::GateBootstrappingTLWE2TLWE<brP>(res, dom, bk, tv);
}

// Compute the expected Sign PBS output based on the phase analysis.
//
// The relevant question: is (int32_t)(m * Δ + offset) > 0 ?
//
// From the blind rotate analysis:
//   phase_2N = round( (m*Δ + offset) * 2N / 2^q ) mod 2N
//
// b̄ = 2N − phase_2N.
//
// TV = [−μ, …, −μ]. After rotation by X^{b̄}:
//   X^a for a < N → coeff[0] = −μ          → decrypt = false
//   X^a for a ≥ N → X^N = −1 → coeff[0] = +μ → decrypt = true
//
// b̄ ≥ N ⟺ 2N − phase_2N ≥ N ⟺ phase_2N ≤ N
// Since phase_2N = round(phase * 2N / 2^q) and 2N / 2^q = 1/2^{q-nbit-1},
// phase_2N ≤ N ⟺ phase ≤ 2^{q-1} (approx)
// i.e., (uint32_t)phase ≤ 2^31
// i.e., (int32_t)phase ≥ 0.
//
// BUT WAIT: b̄ = 0 is a special case (maps to identity rotation → −μ).
// So b̄ ≥ N means phase_2N ∈ [0, N], specifically phase_2N ∈ {1, ..., N}.
// b̄ < N (including b̄ = 0) means phase_2N ∈ {0, N+1, ..., 2N-1}.
//
// Actually: b̄ = 2N − phase_2N (mod 2N).  If phase_2N = 0 → b̄ = 0 (or 2N ≡ 0 mod 2N).
// b̄ = 0 → X^0 = 1 → TV unchanged → coeff[0] = −μ → false.
//
// So output = true  iff  phase_2N ∈ [1, N]  (approximately: 0 < phase < 2^{q-1})
//    output = false iff  phase_2N ∈ {0} ∪ [N+1, 2N-1]
//
// Easier approximation: output = true iff (int32_t)(m*Δ + offset) > 0.
//
bool compute_expected(uint32_t m, uint32_t offset)
{
    // m*Δ + offset as uint32_t, then interpret as signed.
    //
    // The Sign PBS with TV = all -μ outputs:
    //   true  (1) iff  signed torus phase < 0  (upper half: [2^31, 2^32))
    //   false (0) iff  signed torus phase ≥ 0  (lower half: [0, 2^31))
    //
    // This is because blind rotate with b̄ = 2N - phase_2N:
    //   b̄ ∈ [1, N] → coeff[0] = +μ → true   (phase_2N ∈ [N, 2N-1])
    //   b̄ = 0 or b̄ ∈ (N, 2N)  → coeff[0] = -μ → false (phase_2N ∈ {0, 1..N-1})
    //
    constexpr uint32_t Delta = 1u << 23;  // 2^(q-p) for q=32, p=9
    uint32_t phase = m * Delta + offset;
    return static_cast<int32_t>(phase) < 0;
}

}  // namespace

int main()
{
    // ─── Type aliases ─────────────────────────────────────────────────────
    using iksP   = TFHEpp::lvl10param;   // KS: lvl1 → lvl0
    using brP    = TFHEpp::lvl01param;   // BR: lvl0 → lvl1
    using PlainP = TFHEpp::lvl1param;    // Plaintext lives in lvl1 torus

    constexpr int q = std::numeric_limits<PlainP::T>::digits;  // 32
    constexpr int p = 9;   // 9-bit plaintext
    constexpr PlainP::T Delta = static_cast<PlainP::T>(1) << (q - p);  // 2^23
    constexpr int N = PlainP::n;  // 1024

    constexpr int num_trials = 100;

    // ─── Key generation ──────────────────────────────────────────────────
    std::cout << "[verify_offset_pbs] generating keys..." << std::flush;

    TFHEpp::SecretKey sk;

    auto ksk = std::make_unique<TFHEpp::KeySwitchingKey<iksP>>();
    TFHEpp::ikskgen<iksP>(*ksk, sk);

    auto bk = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    TFHEpp::bkfftgen<brP>(*bk, sk);

    std::cout << " done.\n\n";

    // ─── Parameters ──────────────────────────────────────────────────────
    std::cout << "Parameters:\n"
              << "  q (torus bits)     = " << q << "\n"
              << "  p (plain bits)     = " << p << "\n"
              << "  N (ring dim)       = " << N << "\n"
              << "  Delta = 2^" << (q - p) << " = " << Delta << "\n"
              << "  mu                 = 0x" << std::hex << PlainP::μ
              << std::dec << " (" << PlainP::μ << ")\n"
              << "  trials per k       = " << num_trials << "\n\n";

    // ─── Header ──────────────────────────────────────────────────────────
    std::cout << "Sign PBS output semantics:\n"
              << "  TV = all -mu.  After blind rotate:\n"
              << "    phase ∈ (0, 2^31] → +mu → decrypt = true  (1)\n"
              << "    phase ∈ (2^31, 0] → -mu → decrypt = false (0)\n"
              << "  With offset, effective phase = m*Delta + offset.\n\n";

    std::cout << std::setw(4) << "k"
              << std::setw(14) << "offset"
              << std::setw(10) << "correct"
              << std::setw(10) << "errors"
              << std::setw(10) << "trials"
              << std::setw(10) << "status"
              << "\n";
    std::cout << std::string(58, '-') << "\n";

    std::mt19937 rng(42);

    int total_errors = 0;

    for (int k = 0; k < p; k++) {
        // offset = (2^k + 1) / 2 * Delta = (2^k + 1) * 2^(q-p-1)
        const PlainP::T offset =
            static_cast<PlainP::T>((1u << k) + 1u) *
            (static_cast<PlainP::T>(1) << (q - p - 1));

        // Generate random p-bit numbers with bit k = 0 (MSB-first).
        // MSB-first bit k means: the (p-1-k)-th bit from LSB is 0.
        const int bit_pos = p - 1 - k;  // LSB-indexed position
        const uint32_t mask = 1u << bit_pos;
        const uint32_t max_val = (1u << p) - 1;  // 0..511

        std::vector<uint32_t> valid_vals;
        for (uint32_t v = 0; v <= max_val; v++) {
            if ((v & mask) == 0) valid_vals.push_back(v);
        }
        std::uniform_int_distribution<size_t> idx_dist(0, valid_vals.size() - 1);

        int errors = 0;
        int first_mismatches = 0;
        for (int t = 0; t < num_trials; t++) {
            uint32_t m = valid_vals[idx_dist(rng)];

            // Encrypt: plaintext = m * Delta
            TFHEpp::TLWE<PlainP> ct{};
            PlainP::T encoded = static_cast<PlainP::T>(m) * Delta;
            TFHEpp::tlweSymEncrypt<PlainP>(ct, encoded, sk);

            // Sign PBS with custom offset
            TFHEpp::TLWE<PlainP> res{};
            SignPBS_CustomOffset<iksP, brP>(res, ct, offset, *ksk, *bk);

            // Decrypt
            bool dec = TFHEpp::tlweSymDecrypt<PlainP>(res, sk);

            // Expected: phase = m * Delta + offset
            bool expected = compute_expected(m, offset);

            if (dec != expected) {
                errors++;
                if (first_mismatches < 3) {
                    uint32_t eff = static_cast<PlainP::T>(m) * Delta + offset;
                    std::cerr << "  [MISMATCH] k=" << k
                              << " m=" << m
                              << " (0b" << std::bitset<9>(m) << ")"
                              << " eff_phase=0x" << std::hex << eff << std::dec
                              << " (signed=" << static_cast<int32_t>(eff) << ")"
                              << " expected=" << expected
                              << " got=" << dec << "\n";
                    first_mismatches++;
                }
            }
        }

        total_errors += errors;

        std::cout << std::setw(4) << k
                  << "  0x" << std::hex << std::setw(8) << std::setfill('0')
                  << offset << std::dec << std::setfill(' ')
                  << std::setw(10) << (num_trials - errors)
                  << std::setw(10) << errors
                  << std::setw(6) << "/" << num_trials
                  << "    " << (errors == 0 ? "PASS" : "FAIL")
                  << "\n";
    }

    std::cout << "\nTotal errors: " << total_errors << " / " << (p * num_trials) << "\n";
    std::cout << "[verify_offset_pbs] done.\n";

    return total_errors > 0 ? 1 : 0;
}
