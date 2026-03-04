//
// 验证: 密文左移 vs 右移对明文的影响
// 左移 (<<k) = 乘以 2^k，在模运算下良定义 → 正确
// 右移 (>>k) = 整数除法，在模运算下不良定义 → 错误
//
#include <cassert>
#include <cstdio>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace TFHEpp;

// 对 TLWE 密文做左移 k 位 (所有分量 ×2^k)
template <class P>
void tlweLeftShift(TLWE<P> &res, const TLWE<P> &c, int k)
{
    for (size_t i = 0; i <= P::k * P::n; i++)
        res[i] = c[i] << k;
}

// 对 TLWE 密文做右移 k 位 (所有分量整数除以 2^k)
template <class P>
void tlweRightShift(TLWE<P> &res, const TLWE<P> &c, int k)
{
    for (size_t i = 0; i <= P::k * P::n; i++)
        res[i] = c[i] >> k;
}

int main()
{
    using P = lvl1param;
    constexpr uint32_t plain_modulus = 256;  // 8-bit 明文空间
    constexpr int num_test = 100;

    SecretKey sk;

    std::random_device rd;
    std::default_random_engine engine(rd());
    // 明文范围 [0, 127] (正数，避免符号问题)
    std::uniform_int_distribution<int32_t> msg_dist(0, 127);

    printf("=== TLWE Ciphertext Shift Test ===\n\n");
    printf("plain_modulus=%u, Δ = 2^32 / %u = 2^%d\n\n",
           plain_modulus, plain_modulus, 32 - 8);

    // ========== Test 1: 左移 ==========
    printf("--- Test 1: Left Shift (<<k) ---\n");
    for (int k = 1; k <= 4; k++) {
        int errors = 0;
        for (int t = 0; t < num_test; t++) {
            int32_t m = msg_dist(engine);
            int32_t expected = (m << k) % (int32_t)plain_modulus;  // 模明文空间

            TLWE<P> c, c_shifted;
            tlweSymIntEncrypt<P, plain_modulus>(c, (P::T)m, sk);
            tlweLeftShift<P>(c_shifted, c, k);

            // 解密: 左移 k 后, Δ 也变为 Δ·2^k, 对应 plain_modulus/2^k
            // 但实际上我们可以直接用 phase 来验证
            P::T phase_orig = tlweSymPhase<P>(c, sk.key.get<P>());
            P::T phase_shifted = tlweSymPhase<P>(c_shifted, sk.key.get<P>());

            // phase_shifted 应该 ≈ phase_orig << k
            P::T expected_phase = phase_orig << k;
            int32_t diff = (int32_t)(phase_shifted - expected_phase);
            if (std::abs(diff) > 1) errors++;
        }
        printf("  k=%d: errors = %d / %d  %s\n", k, errors, num_test,
               errors == 0 ? "✓ PASS" : "✗ FAIL");
    }

    // ========== Test 2: 右移 ==========
    printf("\n--- Test 2: Right Shift (>>k) ---\n");
    for (int k = 1; k <= 4; k++) {
        int errors = 0;
        int64_t total_phase_error = 0;
        for (int t = 0; t < num_test; t++) {
            int32_t m = msg_dist(engine);

            TLWE<P> c, c_shifted;
            tlweSymIntEncrypt<P, plain_modulus>(c, (P::T)m, sk);
            tlweRightShift<P>(c_shifted, c, k);

            P::T phase_orig = tlweSymPhase<P>(c, sk.key.get<P>());
            P::T phase_shifted = tlweSymPhase<P>(c_shifted, sk.key.get<P>());

            // 如果右移正确, phase_shifted 应该 ≈ phase_orig >> k
            P::T expected_phase = phase_orig >> k;
            int32_t diff = (int32_t)(phase_shifted - expected_phase);
            total_phase_error += std::abs(diff);
            // 允许很大的容差(因为截断误差可能累积)
            // phase 量级 ~ 2^32, 如果误差超过 2^24 就算错
            if (std::abs(diff) > (1 << 24)) errors++;
        }
        printf("  k=%d: large errors = %d / %d, avg phase error = %ld  %s\n",
               k, errors, num_test, (long)(total_phase_error / num_test),
               errors == 0 ? "? check avg" : "✗ FAIL");
    }

    // ========== Test 3: 详细展示一个具体例子 ==========
    printf("\n--- Test 3: Detailed Example ---\n");
    {
        int32_t m = 42;
        int k = 2;

        TLWE<P> c, c_left, c_right;
        tlweSymIntEncrypt<P, plain_modulus>(c, (P::T)m, sk);

        tlweLeftShift<P>(c_left, c, k);
        tlweRightShift<P>(c_right, c, k);

        P::T phase_orig  = tlweSymPhase<P>(c, sk.key.get<P>());
        P::T phase_left  = tlweSymPhase<P>(c_left, sk.key.get<P>());
        P::T phase_right = tlweSymPhase<P>(c_right, sk.key.get<P>());

        // 理论期望
        P::T expected_left  = phase_orig << k;  // = phase * 4
        P::T expected_right = phase_orig >> k;  // = phase / 4 (整数)

        int32_t left_err  = (int32_t)(phase_left - expected_left);
        int32_t right_err = (int32_t)(phase_right - expected_right);

        printf("  m = %d, k = %d\n", m, k);
        printf("  Original phase   = 0x%08X  (≈ Δ·%d = 0x%08X)\n",
               phase_orig, m, (uint32_t)((uint64_t)m * ((1ULL << 32) / plain_modulus)));
        printf("  Left shift phase = 0x%08X  (expected 0x%08X, err = %d)  %s\n",
               phase_left, expected_left, left_err,
               std::abs(left_err) <= 1 ? "✓" : "✗");
        printf("  Right shift phase= 0x%08X  (expected 0x%08X, err = %d)  %s\n",
               phase_right, expected_right, right_err,
               std::abs(right_err) <= (1 << 20) ? "~ok (small noise)" : "✗ WRONG");

        // 但是真正的问题: 右移后的密文能正确解密吗?
        // 用 plain_modulus 解密原始密文
        auto dec_orig = tlweSymIntDecrypt<P, plain_modulus>(c, sk);
        printf("\n  Original decrypt = %u (expected %d) %s\n",
               dec_orig, m, (int32_t)dec_orig == m ? "✓" : "✗");

        // 左移后用 plain_modulus/4 解密? 不，问题更本质:
        // 密文的 phase = Δ·m + e
        // 左移: phase' = 4Δm + 4e, 这是用 new_Δ = 4Δ 编码的 m，或用 Δ 编码的 4m
        // 右移: phase' ≠ Δm/4 + e/4 (因为右移不是模运算下的除法!)
        printf("\n  --- Interpretation ---\n");
        printf("  Left shift:  phase*4 is well-defined mod 2^32 (multiplication)\n");
        printf("  Right shift: phase/4 is NOT well-defined mod 2^32 (depends on representative)\n");
        printf("  Example: if phase wrapped to value V, then V>>2 gives top bits of V,\n");
        printf("  but the 'true' phase could be V + 2^32, and (V+2^32)>>2 = V>>2 + 2^30 ≠ V>>2\n");
    }

    // ========== Test 4: 直接展示发散性 ==========
    printf("\n--- Test 4: Right Shift Decryption Failure ---\n");
    {
        int errors = 0;
        int k = 1;
        // 加密 m, 右移, 然后尝试把它当作 Enc(m/2) 来解密
        for (int t = 0; t < num_test; t++) {
            int32_t m = msg_dist(engine);
            int32_t expected = m >> k;  // 明文右移

            TLWE<P> c, c_shifted;
            tlweSymIntEncrypt<P, plain_modulus>(c, (P::T)m, sk);
            tlweRightShift<P>(c_shifted, c, k);

            // 尝试用 plain_modulus*2 来解密 (因为 Δ 变为 Δ/2)
            constexpr uint32_t new_mod = plain_modulus * 2;
            auto dec = tlweSymIntDecrypt<P, new_mod>(c_shifted, sk);
            if ((int32_t)dec != expected) errors++;
        }
        printf("  Right shift >>%d then decrypt: errors = %d / %d  %s\n",
               k, errors, num_test,
               errors > num_test/4 ? "✗ CONFIRMED BROKEN" : "? unexpectedly few errors");
    }

    printf("\n=== Conclusion ===\n");
    printf("Left shift (<<k):  CORRECT - multiplication by 2^k is well-defined mod 2^32\n");
    printf("Right shift (>>k): BROKEN  - integer division is NOT well-defined mod 2^32\n");
    printf("  The phase (b - a·s) wraps around mod 2^32,\n");
    printf("  so (b mod 2^32) >> k ≠ (b_true >> k) mod 2^32\n");

    return 0;
}
