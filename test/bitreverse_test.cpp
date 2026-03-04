//
// 单次传统 Bootstrap 实现 5-bit 无符号整数的二进制反转
//
// 使用 GateBootstrappingTLWE2TLWE (标准 PBS)，不使用 MetaPBS。
//
// 核心思路:
//   bitrev(m) 不满足 negacyclic 约束 f(x+t/2) = -f(x)，
//   因此将明文空间翻倍: t' = 2t = 64，构造
//     f_half[m] = bitrev(m, 5) 对 m ∈ [0, 32)
//   输入始终在 [0, 32) 范围内，negacyclic 第二半不会被访问。
//
//   测试向量手动构造:
//     block = 2N / t' = 2048 / 64 = 32
//     对每个 i ∈ [0, t'/2) = [0, 32):
//       TV[i*block + j - block/2] = Δ * bitrev(i, 5)
//     其中 Δ = 2^32 / t' 是 torus 编码系数。
//
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#include <tfhe++.hpp>

namespace {

// 反转 val 的低 bits 位
uint32_t bitrev(uint32_t val, int bits)
{
    uint32_t result = 0;
    for (int i = 0; i < bits; i++) {
        result |= ((val >> i) & 1u) << (bits - 1 - i);
    }
    return result;
}

}  // namespace

int main()
{
    constexpr int B = 5;                       // 5-bit 输入
    constexpr uint32_t t = 1u << B;            // t = 32
    constexpr uint32_t t_prime = 2 * t;        // t' = 64 (翻倍绕过 negacyclic)

    // --- 标准 bootstrap 参数: lvl0 → lvl1 ---
    using bkP     = TFHEpp::lvl01param;
    using domainP = typename bkP::domainP;     // lvl0param
    using targetP = typename bkP::targetP;     // lvl1param

    constexpr int N    = static_cast<int>(targetP::n);    // 1024
    constexpr int twoN = 2 * N;                           // 2048

    static_assert(twoN % t_prime == 0, "t' 必须整除 2N");

    constexpr int block = twoN / static_cast<int>(t_prime);  // 32
    constexpr int half_block = block / 2;                     // 16

    // --- 生成密钥 ---
    std::cout << "=== 单次传统 Bootstrap 二进制反转测试 ===" << std::endl;
    std::cout << "参数: B=" << B << ", t=" << t << ", t'=" << t_prime
              << ", N=" << N << ", block=" << block << std::endl;
    std::cout << "使用: GateBootstrappingTLWE2TLWE (标准 PBS)" << std::endl;

    std::cout << "生成密钥..." << std::endl;
    TFHEpp::SecretKey sk;

    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<bkP>>();
    TFHEpp::bkfftgen<bkP>(*bkfft, sk);

    std::cout << "密钥生成完毕。" << std::endl;

    // --- 手动构造 negacyclic 测试向量 ---
    //
    // 编码: Δ = 2^32 / t'
    // 对每个消息值 i ∈ [0, t'/2 = 32):
    //   在多项式的第 i 个 block 内填入 Δ * bitrev(i, B)
    //   block 以中心对齐 (偏移 -half_block) 以容忍舍入误差
    //
    // negacyclic 自动保证 i ∈ [t'/2, t') 部分取反: f(i+32) = -f(i)

    TFHEpp::Polynomial<targetP> tv{};
    for (auto &x : tv) x = 0;

    for (int i = 0; i < static_cast<int>(t_prime / 2); i++) {
        // Torus 编码: bitrev(i, B) * 2^32 / t'
        const auto y = static_cast<typename targetP::T>(
            (static_cast<__int128_t>(bitrev(static_cast<uint32_t>(i), B))
             << std::numeric_limits<typename targetP::T>::digits) /
            t_prime);

        for (int j = 0; j < block; j++) {
            int e = i * block + j - half_block;
            e %= twoN;
            if (e < 0) e += twoN;
            if (e >= N)
                tv[e - N] -= y;
            else
                tv[e] += y;
        }
    }

    // --- 打印 bit reversal 表 ---
    std::cout << "\n--- 二进制反转表 (5-bit) ---" << std::endl;
    for (uint32_t m = 0; m < t; m++) {
        std::cout << "  ";
        for (int b = B - 1; b >= 0; b--)
            std::cout << ((m >> b) & 1);
        std::cout << " (" << m << ") -> ";
        uint32_t rev = bitrev(m, B);
        for (int b = B - 1; b >= 0; b--)
            std::cout << ((rev >> b) & 1);
        std::cout << " (" << rev << ")" << std::endl;
    }

    // --- 全域扫描测试 ---
    std::cout << "\n--- 全域扫描: m ∈ [0, " << t << ") ---" << std::endl;

    int errors = 0;
    for (uint32_t m = 0; m < t; m++) {
        // 在 lvl0 加密 m，明文模数 t'
        TFHEpp::TLWE<domainP> cin{};
        TFHEpp::tlweSymIntEncrypt<domainP, t_prime>(cin, m, sk);

        // 一次标准 PBS
        TFHEpp::TLWE<targetP> cout{};
        TFHEpp::GateBootstrappingTLWE2TLWE<bkP>(cout, cin, *bkfft, tv);

        // 在 lvl1 解密，明文模数 t'
        const auto dec = TFHEpp::tlweSymIntDecrypt<targetP, t_prime>(cout, sk);
        const uint32_t expected = bitrev(m, B);

        if (dec != expected) {
            std::cerr << "  ✗ m=" << m << " dec=" << dec
                      << " expected=" << expected << std::endl;
            errors++;
        } else {
            std::cout << "  ✓ m=" << m << " -> bitrev=" << dec << std::endl;
        }
    }

    std::cout << "\n=== 结果 ===" << std::endl;
    if (errors == 0) {
        std::cout << "全部 " << t << " 个输入通过! 传统 Bootstrap 二进制反转正确。"
                  << std::endl;
    } else {
        std::cerr << errors << "/" << t << " 个输入失败!" << std::endl;
        return 1;
    }

    // --- 随机测试 ---
    std::cout << "\n--- 随机测试 (50 次随机 5-bit 输入) ---" << std::endl;
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<uint32_t> dist(0, t - 1);

    int random_errors = 0;
    constexpr int random_trials = 50;
    for (int trial = 0; trial < random_trials; trial++) {
        const uint32_t m = dist(rng);

        TFHEpp::TLWE<domainP> cin{};
        TFHEpp::tlweSymIntEncrypt<domainP, t_prime>(cin, m, sk);

        TFHEpp::TLWE<targetP> cout{};
        TFHEpp::GateBootstrappingTLWE2TLWE<bkP>(cout, cin, *bkfft, tv);

        const auto dec = TFHEpp::tlweSymIntDecrypt<targetP, t_prime>(cout, sk);
        const uint32_t expected = bitrev(m, B);

        if (dec != expected) {
            std::cerr << "  ✗ trial=" << trial << " m=" << m
                      << " dec=" << dec << " expected=" << expected << std::endl;
            random_errors++;
        }
    }

    std::cout << "随机测试: " << (random_trials - random_errors) << "/"
              << random_trials << " 通过" << std::endl;

    if (random_errors > 0) {
        std::cerr << "随机测试有 " << random_errors << " 个失败!" << std::endl;
        return 1;
    }

    std::cout << "\n=== 全部测试通过 ===" << std::endl;
    return 0;
}
