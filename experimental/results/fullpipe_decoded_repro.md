# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|decoded-ablation|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|6|[0.0106,0.0494]|15|0|0|9|6|-5|68|165360|55.880|129|0|256:6|ablation|
|decoded-ablation|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|6|[0.0106,0.0494]|15|0|0|9|6|-5|0|165360|86.376|129|0|256:6|ablation|
|decoded-ablation|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|1|[0.0007,0.0215]|6|0|0|5|1|0|68|165360|55.702|129|0|255:1|ablation|
|decoded-ablation|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|1|[0.0007,0.0215]|6|0|0|5|1|0|0|165360|80.746|129|0|255:1|ablation|
|decoded-ablation|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|1|[0.0007,0.0215]|6|0|0|5|1|0|68|165360|60.288|129|0|255:1|ablation|
|decoded-ablation|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|1|[0.0007,0.0215]|6|0|0|5|1|0|0|165360|83.377|129|0|255:1|ablation|
