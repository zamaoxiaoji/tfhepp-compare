# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|failure_rate|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|p50_ms|p95_ms|min_opp|ambiguous|failures_by_m|failures_by_boundary|failures_by_seed|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|6|0|0.000000|[0.0000,0.3903]|0|0|0|0|0|0|1|3816|107.981|98.077|108.125|129|0||||ablation|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|6|0|0.000000|[0.0000,0.3903]|0|0|0|0|0|0|1|3816|95.295|90.924|91.681|129|0||||ablation|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|6|0|0.000000|[0.0000,0.3903]|0|0|0|0|0|0|1|3816|89.805|90.533|91.192|129|0||||ablation|

## Paired vs r=0

|candidate|p|k|period|input|cases|baseline_failures|candidate_failures|candidate_rate|Wilson95|fixed_vs_r0|introduced_vs_r0|unchanged_pass|unchanged_fail|net_improvement|McNemar_p|avg_ms|p50_ms|p95_ms|failures_by_m|failures_by_boundary|failures_by_seed|migration_map|recommendation|
|---|---:|---:|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|---|
|exact_pre0d8|9|1|exact|controlled-zero-noise|6|0|0|0.000000|[0.0000,0.3903]|0|0|6|0|0|1.000000|107.981|98.077|108.125|||||keep_r0|
|exact_pre3d8|9|1|exact|controlled-zero-noise|6|0|0|0.000000|[0.0000,0.3903]|0|0|6|0|0|1.000000|95.295|90.924|91.681|||||keep_r0|
|exact_pre4d8|9|1|exact|controlled-zero-noise|6|0|0|0.000000|[0.0000,0.3903]|0|0|6|0|0|1.000000|89.805|90.533|91.192|||||keep_r0|

## Calibration Top

|tag|p|k|candidate|r|denominator|baseline_failures|candidate_failures|fixed|introduced|net|recommendation|
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
|calibration-top|9|1|exact_pre0d8|0|8|0|0|0|0|0|keep_r0|
