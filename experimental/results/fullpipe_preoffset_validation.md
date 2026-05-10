# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|2600|86|[0.0269,0.0407]|178|0|0|92|86|-23|595|1653600|94.747|129|0|255:4;256:80;257:2|ablation|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|2600|49|[0.0143,0.0248]|113|0|0|64|49|14|595|1653600|94.539|129|0|255:22;256:27|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|2600|52|[0.0153,0.0261]|112|0|0|60|52|11|595|1653600|93.451|129|0|255:38;256:14|probabilistic-candidate|
|conversion-pbs|exact_pre5d8|exact|pbs|tfhepp_v10_poly|9|1|exact|2600|63|[0.0190,0.0309]|125|0|0|62|63|0|595|1653600|93.747|129|0|255:54;256:9|ablation|
