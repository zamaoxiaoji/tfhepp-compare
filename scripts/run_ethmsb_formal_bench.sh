#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [[ -z "${ETHMSB_REPORT_DIR:-}" ]]; then
  if [[ -f "$ROOT/reports/latest_ethmsb_report_dir.txt" ]]; then
    ETHMSB_REPORT_DIR="$(cat "$ROOT/reports/latest_ethmsb_report_dir.txt")"
  else
    echo "warning: no prior correctness report found; running formal cert first"
    "$ROOT/scripts/run_ethmsb_formal_cert.sh"
    ETHMSB_REPORT_DIR="$(cat "$ROOT/reports/latest_ethmsb_report_dir.txt")"
  fi
fi

REPORT_DIR="$ETHMSB_REPORT_DIR"
BENCH_DIR="$REPORT_DIR/benchmark"
BUILD_DIR="$ROOT/build-ethmsb-release"
mkdir -p "$BENCH_DIR" "$REPORT_DIR/build_logs"

if [[ ! -f "$REPORT_DIR/correctness/PASS" ]]; then
  echo "BEGIN_NEED_INFO"
  echo "stage: bench_gate"
  echo "what_failed: correctness gate PASS marker is missing"
  echo "commands_run:"
  echo "  scripts/run_ethmsb_formal_bench.sh"
  echo "END_NEED_INFO"
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

TASKSET_CORE="${TASKSET_CORE:-0}"
TASKSET_CMD=()
if command -v taskset >/dev/null 2>&1; then
  TASKSET_CMD=(taskset -c "$TASKSET_CORE")
  export ETHMSB_TASKSET_CORE="$TASKSET_CORE"
else
  echo "warning: taskset not available; running without CPU pinning" | tee "$BENCH_DIR/taskset_warning.txt"
  export ETHMSB_TASKSET_CORE="unavailable"
fi

run_log() {
  local log="$1"
  shift
  {
    printf 'command:'
    printf ' %q' "$@"
    printf '\n'
    "$@"
  } >"$log" 2>&1
}

echo "command: rm -rf $BUILD_DIR"
rm -rf "$BUILD_DIR"
CMAKE_ARGS=(-S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release)
if [[ "${USE_CONCRETE:-OFF}" == "ON" ]]; then CMAKE_ARGS+=(-DUSE_CONCRETE=ON); fi
if [[ "${USE_FFTW3:-OFF}" == "ON" ]]; then CMAKE_ARGS+=(-DUSE_FFTW3=ON); fi
if [[ "${USE_MKL:-OFF}" == "ON" ]]; then CMAKE_ARGS+=(-DUSE_MKL=ON); fi

run_log "$REPORT_DIR/build_logs/bench_cmake_configure.log" cmake "${CMAKE_ARGS[@]}"
run_log "$REPORT_DIR/build_logs/bench_cmake_build.log" cmake --build "$BUILD_DIR" -j"$(nproc)" --target \
  my_ethmsb_bench_he3db_style \
  my_ethmsb_he3db_fairness_audit

run_log "$BENCH_DIR/he3db_fairness_audit.log" "$BUILD_DIR/my_ethmsb_he3db_fairness_audit"
HE3DB_FAIRNESS="$(grep -E '^baseline_fairness=' "$BENCH_DIR/he3db_fairness_audit.log" | tail -n1 | cut -d= -f2 || true)"
if [[ -z "$HE3DB_FAIRNESS" ]]; then HE3DB_FAIRNESS="not_run"; fi

MARGIN_LOG="$(ls "$REPORT_DIR"/correctness/margin_encrypted_t*.log 2>/dev/null | sort | tail -n1 || true)"
MIN_SAFE="$(grep -E 'MARGIN_CERT_RESULT' "$MARGIN_LOG" 2>/dev/null | sed -E 's/.*min_safe_margin=([^ ]+).*/\1/' | tail -n1)"
if [[ -z "$MIN_SAFE" ]]; then MIN_SAFE="unknown"; fi
SWITCH_WIDTH="${SWITCH_BAND_WIDTH:-0x0002000000000000}"

run_log "$BENCH_DIR/bench_ethmsb_h3compat_strict_t30.log" "${TASKSET_CMD[@]}" "$BUILD_DIR/my_ethmsb_bench_he3db_style" \
  --impl ethmsb_h3compat_strict \
  --bits all \
  --ops lt,gt,eq,le,ge,neq \
  --trials "${BENCH_TRIALS:-30}" \
  --warmup "${BENCH_WARMUP:-3}" \
  --seed 123 \
  --fixed-inputs random \
  --include-keygen-time no \
  --margin-cert-min-hex "$MIN_SAFE" \
  --switch-band-width-hex "$SWITCH_WIDTH" \
  --csv "$BENCH_DIR/bench_ethmsb_h3compat_strict_t${BENCH_TRIALS:-30}.csv"

if [[ "$HE3DB_FAIRNESS" == "passed" ]]; then
  echo "HE3DB fairness passed, but this script benchmarks only ethmsb_h3compat_strict unless HE3DB baseline implementation is explicitly enabled." > "$BENCH_DIR/he3db_baseline_status.txt"
else
  echo "HE3DB baseline marked non-comparable/logical-only; no speedup claim generated." > "$BENCH_DIR/he3db_baseline_status.txt"
fi

REPORT="$REPORT_DIR/ETHMSB_REPRO_REPORT.md"
BENCH_CSV="$BENCH_DIR/bench_ethmsb_h3compat_strict_t${BENCH_TRIALS:-30}.csv"
{
  echo "# ETHMSB Reproduction Report"
  echo
  echo "## A. Algorithm implemented"
  echo "This run certifies corrected ETHMSB strict / h3compat_l20_l02. It uses shift=kappa, guard_value=delta(k)*2^(k-kappa-1), unsigned comparison with k=t+1, and arithmetic PBS output 0/out_value. IdeGateBootstrapping is not used."
  echo
  echo "## B. Parameters"
  echo "- TFHEpp release code: see env_report.txt for git hash/tag."
  echo "- Custom h3compat params are in my_he3db_compat_params.hpp."
  echo "- P0: uint32_t, n=672, alpha=2^-16."
  echo "- P2: uint64_t, n=2048, alpha=2^-52."
  echo "- lvl20: t=2, basebit=10."
  echo "- lvl02: my_h3_lvl0param -> my_h3_lvl2param."
  echo "- kappa=5, BOOL_ONE=0x2000000000000000, delta(33)=0x0000000080000000."
  echo "- USE_CONCRETE does not alter h3compat params; this report records the actual macro state in CSV/env logs."
  echo
  echo "## C. Correctness certification"
  echo "- PBS valid matrix log: correctness/pbs_valid_matrix.log."
  echo "- Switch-band bool/guard CSVs: correctness/switchband_bool.csv and correctness/switchband_guard.csv."
  echo "- Margin trivial log: correctness/margin_trivial.log."
  echo "- Margin encrypted log: $(basename "$MARGIN_LOG")."
  echo "- Adversarial comparison: correctness/comparison_adversarial.csv."
  echo "- Random resumable: correctness/random_100.jsonl and correctness/random_1000.jsonl."
  echo "- Small-k exhaustive semantic plus representative PBS: correctness/small_k_exhaustive.log."
  echo "- Exact threshold phase+offset=Q/2-1 is treated as calibration artifact, not a valid encoding failure."
  echo
  echo "## D. Benchmark"
  echo "Key generation, encryption, and decryption are excluded from evaluation latency. Benchmark CSV: benchmark/$(basename "$BENCH_CSV")."
  echo
  echo '```csv'
  grep -v '^#' "$BENCH_CSV" | head -n 20
  echo '```'
  echo
  echo "eq/neq use lt+gt and therefore have 2x ETHMSB PBS count. le/ge use one strict comparison plus arithmetic complement."
  echo
  echo "## E. HE3DB comparison"
  echo "Fairness audit result: $HE3DB_FAIRNESS. If failed, HE3DB is logical-only/non-comparable and no speedup claim is made."
  echo
  echo "## F. Limitations"
  echo "- Random tests are statistical, not proof."
  echo "- Strict variant does not use IdeGateBootstrapping."
  echo "- Supported correctness scope is k<=33 / 32-bit unsigned comparison."
  echo "- The claim is TFHEpp release code plus HE3DB-compatible custom params, not default TFHEpp params."
  echo
  echo "## G. Reproduction commands"
  echo '```bash'
  echo "scripts/run_ethmsb_formal_cert.sh"
  echo "ETHMSB_REPORT_DIR=$REPORT_DIR scripts/run_ethmsb_formal_bench.sh"
  echo "$BUILD_DIR/my_ethmsb_random_resumable --impl h3compat_l20_l02 --bits all --kappa 5 --ops lt,gt,eq --trials-per-bit 1000 --seed 2 --max-seconds 7200 --checkpoint $REPORT_DIR/correctness/random_1000.jsonl --resume"
  echo '```'
} > "$REPORT"

echo "ETHMSB_FORMAL_BENCH_DONE report_dir=$REPORT_DIR bench_csv=$BENCH_CSV he3db_fairness=$HE3DB_FAIRNESS"
