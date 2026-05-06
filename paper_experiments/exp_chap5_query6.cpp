/**
 * @file exp_chap5_query6.cpp
 * @brief Chapter 5 实验: TPC-H Query 6 — 多条件区间过滤 + 算术聚合 (SUM)
 *
 * Query 6 语义:
 *   SELECT SUM(l_extendedprice * l_discount) AS revenue
 *   FROM lineitem
 *   WHERE l_shipdate >= :1
 *     AND l_shipdate < :1 + interval '1' year
 *     AND l_discount BETWEEN :2 - 0.01 AND :2 + 0.01
 *     AND l_quantity < :3;
 *
 * 密态执行流程 (仿照 HE3DB tpch_q6.cpp):
 *   1. 生成合成 lineitem 数据 (shipdate, discount, quantity, price)
 *   2. TFHE 加密各属性列 + 谓词阈值
 *   3. 对每条记录执行 GapMSB 比较 (5 个谓词), 通过 HomAND 组合多谓词
 *   4. 将 TFHE 过滤掩码 repack 到 CKKS 域
 *   5. CKKS: mask × (price * discount) + rotate-and-sum 聚合
 *
 * 参数:
 *   --row-exp E   数据行数 = 2^E, 对应论文 X 轴 (10 ~ 14)
 *
 * 输出:
 *   总查询延迟 (秒), 对应论文 Y 轴
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
                "usage: exp_chap5_query6 --row-exp 10|11|12|13|14");
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

        std::cout << "=== Chapter 5: TPC-H Query 6 ===" << std::endl;
        std::cout << "rows = 2^" << row_exp
                  << " = " << options.rows << std::endl;

        // 执行真实同态加密实验:
        //   TFHE GapMSB 多谓词过滤 → repack → CKKS SUM 聚合
        auto result = PaperReview::RunTpchQ6Experiment(options);

        // 输出实验结果
        std::cout << "\n--- Timing ---" << std::endl;
        std::cout << "filter_time_s     = "
                  << result.filter_time_ms / 1000.0 << std::endl;
        std::cout << "aggregation_time_s= "
                  << result.aggregation_time_ms / 1000.0 << std::endl;
        std::cout << "total_query_time_s= "
                  << result.total_query_time_ms / 1000.0 << std::endl;

        std::cout << "\n--- Correctness ---" << std::endl;
        std::cout << "plain_result    = " << result.plain_scalar << std::endl;
        std::cout << "encrypted_result= " << result.encrypted_scalar << std::endl;
        std::cout << "abs_error       = " << result.abs_error << std::endl;
        std::cout << "predicate_errors= " << result.predicate_errors << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exp_chap5_query6 failed: " << e.what() << "\n";
        return 1;
    }
}
