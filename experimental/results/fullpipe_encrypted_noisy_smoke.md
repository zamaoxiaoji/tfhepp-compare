# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|10|3|[0.1078,0.6032]|3|0|0|0|3|0|0|6360|97.707|129|0|256:3|ablation|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|10|2|[0.0567,0.5098]|2|0|0|0|2|1|4|6360|92.064|129|0|256:2|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|10|3|[0.1078,0.6032]|3|0|0|0|3|0|6|6360|94.508|129|0|256:3|ablation|
