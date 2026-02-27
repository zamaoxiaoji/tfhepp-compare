// TFHEpp/test/metapbs_packing_fullN.cpp
// Minimal Unit F3a (int-plaintext): verify TLWE2TRLWEPacking works for FULL N
// coefficients by encrypting integer coefficients mod M, packing, extracting,
// and integer-decrypting exactly.

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <tfhe++.hpp>
#include <vector>
static inline uint32_t norm_mod_M(uint32_t x_u32, uint32_t M)
{
    int32_t x = (int32_t)x_u32;  // interpret as signed
    int32_t r = x % (int32_t)M;
    if (r < 0) r += (int32_t)M;
    return (uint32_t)r;
}

int main()
{
    using P = TFHEpp::lvl1param;
    constexpr int N = P::n;

    // Choose small plaintext modulus for exact integer decrypt
    constexpr uint32_t M = 256;

    std::cout << "[metapbs_packing_fullN:int] N=" << N << " M=" << M << "\n";
    std::cout << "[metapbs_packing_fullN:int] building keys...\n";

    TFHEpp::SecretKey sk;

    // packing key
    TFHEpp::AnnihilateKey<P> ahk;
    TFHEpp::annihilatekeygen<P>(ahk, sk);

    std::mt19937_64 rng(0);
    std::uniform_int_distribution<uint32_t> dist(0, M - 1);

    // Encrypt each integer coefficient as TLWE<P, M>
    std::vector<uint32_t> plain(N);
    std::vector<TFHEpp::TLWE<P>> coeffs(N);

    for (int i = 0; i < N; i++) {
        plain[i] = dist(rng);
        TFHEpp::tlweSymIntEncrypt<P, M>(coeffs[i], (typename P::T)plain[i], sk);
    }

    // Pack to TRLWE
    TFHEpp::TRLWE<P> packed{};
    TFHEpp::TLWE2TRLWEPacking<P>(packed, coeffs, ahk);

    // Verify by extracting each coeff and integer-decrypting
    for (int i = 0; i < N; i++) {
        TFHEpp::TLWE<P> ci{};
        TFHEpp::SampleExtractIndex<P>(ci, packed, i);
        uint32_t got_raw = (uint32_t)TFHEpp::tlweSymIntDecrypt<P, M>(ci, sk);
        uint32_t got = norm_mod_M(got_raw, M);
        if (got != plain[i]) {
            std::cerr << "FAIL coeff " << i << " got=" << got
                      << " expect=" << plain[i] << "\n";
            return 1;
        }
    }

    std::cout
        << "PASS: TLWE2TRLWEPacking packs FULL N integer coeffs exactly (M="
        << M << ").\n";
    return 0;
}
