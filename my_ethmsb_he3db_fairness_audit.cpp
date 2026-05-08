#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

#include "my_he3db_compat_params.hpp"
#include "params.hpp"

namespace {

template <class P>
void print_lwe_param(const std::string& label)
{
    std::cout << label << "_T_bits="
              << std::numeric_limits<typename P::T>::digits << "\n";
    std::cout << label << "_n=" << P::n << "\n";
    std::cout << label << "_alpha=" << P::α << "\n";
}

template <class A, class B>
bool same_lwe_core()
{
    return std::numeric_limits<typename A::T>::digits ==
               std::numeric_limits<typename B::T>::digits &&
           A::n == B::n && A::α == B::α;
}

}  // namespace

int main()
{
    using namespace my_ethmsb_params;
    std::cout << "HE3DB_FAIRNESS_AUDIT\n";
    std::cout << "current_tfhepp_link=this_project_tfhepp\n";
    std::cout << "he3db_sources=../HE3DB/src/HEDB/comparison when target is linked by CMake\n";
    std::cout << "vendored_tfhepp=not_linked_by_my_ethmsb_targets\n";

    print_lwe_param<TFHEpp::lvl0param>("v10_lvl0");
    print_lwe_param<TFHEpp::lvl1param>("v10_lvl1");
    print_lwe_param<TFHEpp::lvl2param>("v10_lvl2");
    print_lwe_param<my_h3_lvl0param>("h3_lvl0");
    print_lwe_param<my_h3_lvl1param>("h3_lvl1");
    print_lwe_param<my_h3_lvl2param>("h3_lvl2");

    std::cout << "v10_lvl20_t=" << TFHEpp::lvl20param::t << "\n";
    std::cout << "v10_lvl20_basebit=" << TFHEpp::lvl20param::basebit << "\n";
    std::cout << "h3_lvl20_t=" << my_h3_lvl20param::t << "\n";
    std::cout << "h3_lvl20_basebit=" << my_h3_lvl20param::basebit << "\n";
    std::cout << "h3_lvl02_domain=my_h3_lvl0param\n";
    std::cout << "h3_lvl02_target=my_h3_lvl2param\n";

    const bool lvl0_same =
        same_lwe_core<TFHEpp::lvl0param, my_h3_lvl0param>();
    const bool lvl2_same =
        same_lwe_core<TFHEpp::lvl2param, my_h3_lvl2param>();
    const bool lvl20_same =
        TFHEpp::lvl20param::t == my_h3_lvl20param::t &&
        TFHEpp::lvl20param::basebit == my_h3_lvl20param::basebit;

    std::cout << "output_encoding_ethmsb_strict=arithmetic_0_or_BOOL_ONE_lvl2\n";
    std::cout << "output_encoding_he3db_baseline=logical_lvl1_in_HE3DB_wrappers\n";
    std::cout << "comparison_scale_ethmsb=k_plain_bits_plus_1\n";
    std::cout << "comparison_scale_he3db=plain_bits_plus_1_expected_by_HE3DB_tests\n";
    std::cout << "pbs_count_ethmsb_k9=2\n";
    std::cout << "pbs_count_ethmsb_k17=4\n";
    std::cout << "pbs_count_ethmsb_k25=5\n";
    std::cout << "pbs_count_ethmsb_k33=7\n";
    std::cout << "baseline_fairness="
              << (lvl0_same && lvl2_same && lvl20_same ? "passed"
                                                        : "failed")
              << "\n";
    if (!(lvl0_same && lvl2_same && lvl20_same)) {
        std::cout << "baseline_status=non_comparable_or_logical_only\n";
        std::cout << "reason=current_HE3DB_baseline_target_links_project_TFHEpp_params_not_h3compat_custom_params\n";
    }
    return 0;
}
