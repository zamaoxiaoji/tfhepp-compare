# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|12|1|exact|190|23|[0.0820,0.1751]|59|0|0|36|23|2|43|120840|89.380|1025|0|0:2;1:1;2044:1;2045:4;2046:2;2047:3;2048:2;2049:1;4093:1;4094:4;4095:2|probabilistic-candidate|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|12|1|exact|190|24|[0.0864,0.1811]|61|0|0|37|24|1|43|120840|84.787|1025|0|0:2;1:1;2044:1;2045:4;2046:2;2047:3;2048:2;2049:1;4093:1;4094:4;4095:3|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|12|1|exact|190|25|[0.0907,0.1870]|62|0|0|37|25|0|43|120840|89.238|1025|0|0:2;1:1;2044:1;2045:4;2046:2;2047:3;2048:2;2049:1;4093:2;4094:4;4095:3|ablation|
