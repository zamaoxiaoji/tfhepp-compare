#include <iostream>
#include <set>
#include <vector>

#include "my_ethmsb_pruned_fused_final.hpp"

namespace {

using namespace my_ethmsb_pruned_fused_final;

std::vector<int> k_choices(const int p)
{
    std::set<int> ks;
    if (p > 1) ks.insert(1);
    if (p > 1) ks.insert(p - 1);
    if (p > 2) ks.insert(p - 2);
    if (p > 3) ks.insert(p - 3);
    if (p > 4) ks.insert(p - 4);
    if (p > 5) ks.insert(p - 5);
    if (p > 11) ks.insert(11);
    return {ks.begin(), ks.end()};
}

}  // namespace

int main()
{
    std::cout << "p,k,W,guard_value_hex,period_idx,period_divides_N,"
                 "direct_mask_periodic_possible,reason\n";
    for (const int p : {4, 5, 6, 8, 9, 16, 17, 25, 33}) {
        for (const int k : k_choices(p)) {
            const auto info = bit_index_info(p, k);
            const Torus guard = guard_value_for_bit_index(p, k);
            const bool period_divides =
                info.rotation_period_index != 0 &&
                (BR_CYCLE % info.rotation_period_index) == 0;
            const bool possible = direct_guard_mask_periodic_possible(p, k);
            std::cout << p << ',' << k << ','
                      << bit_weight_msb_index(p, k) << ',' << hex64(guard)
                      << ',' << info.rotation_period_index << ','
                      << period_divides << ',' << possible << ','
                      << (possible
                              ? "guard_value_is_qhalf"
                              : "negacyclic_periodic_output_requires_self_"
                                "negating_0_or_qhalf")
                      << "\n";
        }
    }
    return 0;
}
