#include <algorithm>
#include <iostream>

#include "my_ethmsb_pruned_fused_final.hpp"

namespace {

using namespace my_ethmsb_pruned_fused_final;

void print_case(const int p, const int k, const Torus c)
{
    const Torus b = (c >> (p - 1 - k)) & 1;
    std::cout << p << ',' << k << ',' << c << ',' << b << ','
              << fused_final_truth_bit(p, k, c) << "\n";
}

}  // namespace

int main()
{
    const int p = 9;
    const int k = 1;
    bool ok = fused_truth_antisymmetry_ok(p, k);

    std::cout << "p,k,c,bit_k,truth_bit\n";
    for (const Torus c : {Torus{127}, Torus{128}, Torus{383}, Torus{384}})
        print_case(p, k, c);

    for (int pp : {4, 5, 6, 8, 9, 16, 17, 25, 33}) {
        for (int kk = 1; kk <= std::min(pp - 1, 11); ++kk) {
            if (!fused_truth_antisymmetry_ok(pp, kk)) {
                ok = false;
                std::cerr << "antisymmetry_failed p=" << pp << " k=" << kk
                          << "\n";
            }
        }
    }
    std::cerr << "FUSED_LUT_TRUTH_TABLE_AUDIT antisymmetry_ok=" << ok
              << "\n";
    return ok ? 0 : 1;
}
