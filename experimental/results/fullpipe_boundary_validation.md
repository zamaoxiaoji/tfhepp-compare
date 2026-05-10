# weighted_bitapprox summary

## Candidate results

|pipeline|candidate|final_model|source|policy|p|k|period|cases|failures|Wilson95|bit_failures|conversion_failures|final_failures_given_correct_guard|benign_bit_errors|fatal_bit_errors|failure_reduction_vs_exact|skipped_cmux|cmux_total|avg_ms|min_opp|ambiguous|failures_by_m|recommendation|
|---|---|---|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|exact|6800|91|[0.0109,0.0164]|187|0|0|96|91|-39|1625|4324800|96.791|129|0|0:42;1:3;255:2;256:43;257:1|ablation|
|conversion-pbs|exact_pre0d8|exact|pbs|tfhepp_v10_poly|9|1|full|6800|91|[0.0109,0.0164]|187|0|0|96|91|-39|0|4324800|131.236|129|0|0:42;1:3;255:2;256:43;257:1|ablation|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|exact|6800|51|[0.0057,0.0098]|112|0|0|61|51|1|1625|4324800|98.245|129|0|0:13;1:1;255:8;256:16;511:13|probabilistic-candidate|
|conversion-pbs|exact_pre3d8|exact|pbs|tfhepp_v10_poly|9|1|full|6800|51|[0.0057,0.0098]|112|0|0|61|51|1|0|4324800|123.986|129|0|0:13;1:1;255:8;256:16;511:13|probabilistic-candidate|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|exact|6800|52|[0.0058,0.0100]|105|0|0|53|52|0|1625|4324800|93.017|129|0|0:10;255:15;256:8;511:19|ablation|
|conversion-pbs|exact_pre4d8|exact|pbs|tfhepp_v10_poly|9|1|full|6800|52|[0.0058,0.0100]|105|0|0|53|52|0|0|4324800|119.854|129|0|0:10;255:15;256:8;511:19|ablation|
