#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${BUILD_DIR:-build-ethmsb-release}"
REPORT_DIR="${REPORT_DIR:-reports}"
mkdir -p "${REPORT_DIR}"

cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" \
  --target my_ethmsb_h3compat_targeted_tests \
           my_ethmsb_pruned_failure_sweep \
           my_ethmsb_pruned_failure_decomp \
           my_ethmsb_pruned_gapmsb_bench_smoke \
           my_ethmsb_pruned_tex_claim_audit \
           my_ethmsb_direct_mask_periodicity_audit \
           my_ethmsb_fused_lut_truth_table_audit \
           my_ethmsb_fused_final_tests \
           my_ethmsb_fused_vs_conversion_ab \
           my_ethmsb_pruned_fused_final_bench_smoke \
  -j"$(nproc)"

"${BUILD_DIR}/my_ethmsb_h3compat_targeted_tests" --adversarial --kappa 5

"${BUILD_DIR}/my_ethmsb_pruned_failure_sweep" \
  --impl ethmsb_h3compat_pruned_gap_single_round \
  --bits all \
  --k-range 1:11 \
  --trials-random 100 \
  --trials-boundary 5 \
  --seed 0 \
  --target-failure-rate 1e-4 \
  --csv "${REPORT_DIR}/pruned_failure_sweep.csv" \
  --json "${REPORT_DIR}/pruned_failure_sweep_summary.json"

"${BUILD_DIR}/my_ethmsb_pruned_tex_claim_audit" \
  --out "${REPORT_DIR}/pruned_tex_claim_audit.csv"

"${BUILD_DIR}/my_ethmsb_direct_mask_periodicity_audit" \
  > "${REPORT_DIR}/direct_mask_periodicity_audit.csv"

"${BUILD_DIR}/my_ethmsb_fused_lut_truth_table_audit" \
  > "${REPORT_DIR}/fused_lut_truth_table_audit.csv"

"${BUILD_DIR}/my_ethmsb_fused_final_tests" \
  > "${REPORT_DIR}/fused_final_pure_tests.txt"

selected_bits="$(
python3 - "${REPORT_DIR}/pruned_failure_sweep_summary.json" <<'PY'
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as f:
    data = json.load(f)
bits = [str(row["bits"]) for row in data["selected_candidates"]
        if row.get("selected_k") is not None]
print(",".join(bits))
PY
)"

unsupported_bits="$(
python3 - "${REPORT_DIR}/pruned_failure_sweep_summary.json" <<'PY'
import json, sys
with open(sys.argv[1], "r", encoding="utf-8") as f:
    data = json.load(f)
bits = [str(row["bits"]) for row in data["selected_candidates"]
        if row.get("selected_k") is None]
print(",".join(bits))
PY
)"

if [[ -n "${unsupported_bits}" ]]; then
  echo "unsupported_bits=${unsupported_bits}"
fi

if [[ -n "${selected_bits}" ]]; then
  "${BUILD_DIR}/my_ethmsb_pruned_gapmsb_bench_smoke" \
    --impl ethmsb_h3compat_pruned_gap_single_round \
    --bits "${selected_bits}" \
    --use-selected-k "${REPORT_DIR}/pruned_failure_sweep_summary.json" \
    --ops lt,gt,eq \
    --trials 3 \
    --warmup 1 \
    --csv "${REPORT_DIR}/pruned_gap_single_round_smoke.csv"
else
  echo "smoke_bench_skipped=no_selected_candidates"
fi
