/**
 * @file exp_chap5_query5.cpp
 * @brief Chapter 5 实验: TPC-H Query 5 — 过滤 + 六表多路连接 + 分组聚合
 *
 * Query 5 语义:
 *   SELECT n_name, SUM(l_extendedprice * (1 - l_discount)) AS revenue
 *   FROM customer, orders, lineitem, supplier, nation, region
 *   WHERE c_custkey = o_custkey
 *     AND l_orderkey = o_orderkey
 *     AND l_suppkey = s_suppkey
 *     AND c_nationkey = s_nationkey
 *     AND s_nationkey = n_nationkey
 *     AND n_regionkey = r_regionkey
 *     AND r_name = ':1'
 *     AND o_orderdate >= date ':2'
 *     AND o_orderdate < date ':2' + interval '1' year
 *   GROUP BY n_name;
 *
 * 密态执行流程:
 *   1. 生成合成的 nation/region/customer/supplier/orders/lineitem 六表数据
 *   2. TFHE 加密过滤列 (region 等值, orderdate 范围)
 *   3. GapMSB 比较过滤
 *   4. TFHE→CKKS repack 过滤掩码
 *   5. CKKS 多路 LookupJoin:
 *      - orders JOIN customer → customer_nation
 *      - lineitem JOIN orders → order_date 掩码 + customer_nation
 *      - lineitem JOIN supplier → supplier_nation
 *      - supplier JOIN nation → asia 掩码
 *      - EqualityMask(customer_nation, supplier_nation) → same_nation 掩码
 *   6. CKKS: 组合掩码 × revenue → MatrixGroupBySum (nation)
 *
 * 参数:
 *   --row-exp E   数据行数 = 2^E, 对应论文 X 轴 (10 ~ 14)
 *
 * 输出:
 *   总查询延迟 (秒), 分组聚合结果, 对应论文 Y 轴
 */

#include "tpch_experiments.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

static int ParseRowExp(int argc, char** argv) {
    int row_exp = 10;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--row-exp" && i + 1 < argc) {
            row_exp = std::stoi(argv[++i]);
        } else {
            throw std::invalid_argument(
                "usage: exp_chap5_query5 --row-exp 10|11|12|13|14");
        }
    }
    if (row_exp < 1 || row_exp > 20)
        throw std::invalid_argument("row-exp must be between 1 and 20");
    return row_exp;
}

int main(int argc, char** argv) {
    try {
        const int row_exp = ParseRowExp(argc, argv);

        // 构建实验选项 — 仅变化行数, 其余参数使用与论文一致的默认值
        PaperReview::ExperimentOptions options;
        options.rows = std::size_t{1} << row_exp;

        std::cout << "=== Chapter 5: TPC-H Query 5 ===" << std::endl;
        std::cout << "rows = 2^" << row_exp
                  << " = " << options.rows << std::endl;

        // 执行真实同态加密实验:
        //   TFHE GapMSB 过滤 → repack → CKKS 六表多路 LookupJoin → GROUP BY (nation)
        auto result = PaperReview::RunTpchQ5Experiment(options);

        // 输出实验结果
        std::cout << "\n--- Timing ---" << std::endl;
        std::cout << "filter_time_s     = "
                  << result.filter_time_ms / 1000.0 << std::endl;
        std::cout << "aggregation_time_s= "
                  << result.aggregation_time_ms / 1000.0 << std::endl;
        std::cout << "total_query_time_s= "
                  << result.total_query_time_ms / 1000.0 << std::endl;

        std::cout << "\n--- Correctness ---" << std::endl;
        std::cout << "total_revenue_plain    = " << result.plain_scalar
                  << std::endl;
        std::cout << "total_revenue_encrypted= " << result.encrypted_scalar
                  << std::endl;
        std::cout << "predicate_errors       = " << result.predicate_errors
                  << std::endl;

        if (!result.plain_groups.empty()) {
            std::cout << "\n--- Group-by results (per nation) ---"
                      << std::endl;
            for (std::size_t g = 0; g < result.plain_groups.size(); g++) {
                std::cout << "  nation[" << g << "]  plain="
                          << result.plain_groups[g]
                          << "  encrypted=" << result.encrypted_groups[g]
                          << "  error="
                          << std::abs(result.plain_groups[g] -
                                      result.encrypted_groups[g])
                          << std::endl;
            }
            std::cout << "max_abs_error= " << result.abs_error << std::endl;
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query5 failed: " << e.what() << "\n";
        return 1;
    }
}
