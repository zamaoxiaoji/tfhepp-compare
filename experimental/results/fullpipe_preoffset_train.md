# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|6|[0.0106,0.0494]|15|0|0|9|6|-4|68|165360|77.394|129|0|256:6|ablation|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|6|[0.0106,0.0494]|15|0|0|9|6|-4|0|165360|108.397|129|0|256:6|ablation|
|conversion-pbs|exact_pre2d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|3|[0.0039,0.0334]|8|0|0|5|3|-1|68|165360|81.381|129|0|255:1;256:2|ablation|
|conversion-pbs|exact_pre2d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|3|[0.0039,0.0334]|8|0|0|5|3|-1|0|165360|107.643|129|0|255:1;256:2|ablation|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|1|[0.0007,0.0215]|6|0|0|5|1|1|68|165360|76.252|129|0|255:1|probabilistic-candidate|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|1|[0.0007,0.0215]|6|0|0|5|1|1|0|165360|107.328|129|0|255:1|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|1|[0.0007,0.0215]|6|0|0|5|1|1|68|165360|81.601|129|0|255:1|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|1|[0.0007,0.0215]|6|0|0|5|1|1|0|165360|103.510|129|0|255:1|probabilistic-candidate|
|conversion-pbs|exact_pre5d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|2|[0.0021,0.0276]|7|0|0|5|2|0|68|165360|81.830|129|0|255:2|ablation|
|conversion-pbs|exact_pre5d8|exact|pbs|tfhepp_v10_poly|9|1|full|260|2|[0.0021,0.0276]|7|0|0|5|2|0|0|165360|108.208|129|0|255:2|ablation|
