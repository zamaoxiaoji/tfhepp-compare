# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|conversion-pbs|exact_pre-1d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|7|[0.0131,0.0545]|18|0|0|11|7|-3|68|165360|92.308|129|0|256:7|ablation|
|conversion-pbs|exact_pre-2d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|8|[0.0157,0.0595]|20|0|0|12|8|-4|68|165360|88.141|129|0|256:8|ablation|
|conversion-pbs|exact_pre-3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|10|[0.0210,0.0693]|26|0|0|16|10|-6|68|165360|92.347|129|0|256:9;257:1|ablation|
|conversion-pbs|exact_pre-4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|15|[0.0353,0.0930]|33|0|0|18|15|-11|68|165360|92.112|129|0|256:13;257:2|ablation|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|6|[0.0106,0.0494]|15|0|0|9|6|-2|68|165360|93.670|129|0|256:6|ablation|
|conversion-pbs|exact_pre1d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|4|[0.0060,0.0389]|11|0|0|7|4|0|68|165360|93.344|129|0|256:4|ablation|
|conversion-pbs|exact_pre2d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|3|[0.0039,0.0334]|8|0|0|5|3|1|68|165360|90.083|129|0|255:1;256:2|probabilistic-candidate|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|1|[0.0007,0.0215]|6|0|0|5|1|3|68|165360|93.817|129|0|255:1|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|1|[0.0007,0.0215]|6|0|0|5|1|3|68|165360|93.989|129|0|255:1|probabilistic-candidate|
|conversion-pbs|exact_pre5d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|2|[0.0021,0.0276]|7|0|0|5|2|2|68|165360|90.465|129|0|255:2|probabilistic-candidate|
|conversion-pbs|exact_pre6d8|exact|pbs|tfhepp_v10_poly|9|1|exact|260|4|[0.0060,0.0389]|10|0|0|6|4|0|68|165360|92.997|129|0|255:4|ablation|
