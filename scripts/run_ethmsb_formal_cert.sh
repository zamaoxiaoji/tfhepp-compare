#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TS="${ETHMSB_REPORT_TS:-$(date +%Y%m%d_%H%M%S)}"
REPORT_DIR="${ETHMSB_REPORT_DIR:-$ROOT/reports/$TS}"
CORR_DIR="$REPORT_DIR/correctness"
BUILD_DIR="$ROOT/build-ethmsb-release"
mkdir -p "$CORR_DIR" "$REPORT_DIR/build_logs"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

TASKSET_CORE="${TASKSET_CORE:-0}"
TASKSET_CMD=()
if command -v taskset >/dev/null 2>&1; then
  TASKSET_CMD=(taskset -c "$TASKSET_CORE")
  export ETHMSB_TASKSET_CORE="$TASKSET_CORE"
else
  echo "warning: taskset not available; running without CPU pinning" | tee "$REPORT_DIR/taskset_warning.txt"
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

write_env_report() {
  local out="$1"
  {
    echo "command: pwd"; pwd
    echo "command: git status --short"; git status --short || true
    echo "command: git diff -- include/params/128bit.hpp include/params/concrete.hpp"; git diff -- include/params/128bit.hpp include/params/concrete.hpp || true
    echo "command: git rev-parse HEAD"; git rev-parse HEAD || true
    echo "command: git describe --tags --always"; git describe --tags --always || true
    echo "command: lscpu"; lscpu || true
    echo "command: grep -m1 flags /proc/cpuinfo"; grep -m1 flags /proc/cpuinfo || true
    echo "command: nproc"; nproc || true
    echo "command: uname -a"; uname -a
    echo "command: cmake --version"; cmake --version
    echo "command: c++ --version || g++ --version || clang++ --version"; c++ --version || g++ --version || clang++ --version
    echo "env: OMP_NUM_THREADS=$OMP_NUM_THREADS OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS MKL_NUM_THREADS=$MKL_NUM_THREADS"
    echo "env: USE_CONCRETE=${USE_CONCRETE:-OFF} USE_FFTW3=${USE_FFTW3:-OFF} USE_MKL=${USE_MKL:-OFF} USE_KEY_BUNDLE=${USE_KEY_BUNDLE:-OFF} USE_TERNARY_CMUX=${USE_TERNARY_CMUX:-OFF}"
    echo "command: grep -R \"TFHEpp::lvl[012]param\" -n my_ethmsb*.hpp my_ethmsb*.cpp my_he3db_compat_params.hpp"
    grep -R "TFHEpp::lvl[012]param" -n my_ethmsb*.hpp my_ethmsb*.cpp my_he3db_compat_params.hpp || true
    echo "note: h3compat strict targets are expected to use my_h3_lvl0param/my_h3_lvl2param; direct_v10_l22 and debug helpers may mention TFHEpp::lvl*param."
    echo "note: h3compat params are custom structs and do not change with USE_CONCRETE; USE_KEY_BUNDLE only changes Addends."
  } >"$out" 2>&1
}

write_env_report "$REPORT_DIR/env_report.txt"
cp "$REPORT_DIR/env_report.txt" "$ROOT/reports/env_report.txt"

echo "command: rm -rf $BUILD_DIR"
rm -rf "$BUILD_DIR"
CMAKE_ARGS=(-S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release)
if [[ "${USE_CONCRETE:-OFF}" == "ON" ]]; then CMAKE_ARGS+=(-DUSE_CONCRETE=ON); fi
if [[ "${USE_FFTW3:-OFF}" == "ON" ]]; then CMAKE_ARGS+=(-DUSE_FFTW3=ON); fi
if [[ "${USE_MKL:-OFF}" == "ON" ]]; then CMAKE_ARGS+=(-DUSE_MKL=ON); fi

run_log "$REPORT_DIR/build_logs/cmake_configure.log" cmake "${CMAKE_ARGS[@]}"
run_log "$REPORT_DIR/build_logs/cmake_build.log" cmake --build "$BUILD_DIR" -j"$(nproc)" --target \
  my_ethmsb_pbs_matrix_tests \
  my_ethmsb_switchband \
  my_ethmsb_margin_cert \
  my_ethmsb_h3compat_targeted_tests \
  my_ethmsb_random_resumable \
  my_ethmsb_static_tests \
  my_ethmsb_he3db_fairness_audit

BIN="$BUILD_DIR"
run_log "$CORR_DIR/params_audit.txt" "$BIN/my_ethmsb_he3db_fairness_audit"
run_log "$CORR_DIR/pbs_valid_matrix.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_pbs_matrix_tests" --impl h3compat_l20_l02 --valid-only
run_log "$CORR_DIR/switchband_bool.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_switchband" --impl h3compat_l20_l02 --all --out-value bool --csv "$CORR_DIR/switchband_bool.csv"
run_log "$CORR_DIR/switchband_guard.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_switchband" --impl h3compat_l20_l02 --all --out-value guard --csv "$CORR_DIR/switchband_guard.csv"

SWITCH_BAND_WIDTH="${SWITCH_BAND_WIDTH:-0x0002000000000000}"
MARGIN_TRIALS="${MARGIN_TRIALS:-20}"
run_log "$CORR_DIR/margin_trivial.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_margin_cert" --impl h3compat_l20_l02 --kappa 5 --max-k 33 --mode trivial --switch-band-width-hex "$SWITCH_BAND_WIDTH" --csv "$CORR_DIR/margin_trivial.csv"
run_log "$CORR_DIR/margin_encrypted_t${MARGIN_TRIALS}.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_margin_cert" --impl h3compat_l20_l02 --kappa 5 --max-k 33 --mode encrypted --trials-per-case "$MARGIN_TRIALS" --seed 0 --switch-band-width-hex "$SWITCH_BAND_WIDTH" --csv "$CORR_DIR/margin_encrypted_t${MARGIN_TRIALS}.csv"
run_log "$CORR_DIR/comparison_adversarial.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_h3compat_targeted_tests" --adversarial --kappa 5 --csv "$CORR_DIR/comparison_adversarial.csv"
RANDOM_SMOKE_TRIALS="${RANDOM_SMOKE_TRIALS:-100}"
RANDOM_STRONG_TRIALS="${RANDOM_STRONG_TRIALS:-1000}"
RANDOM_STRONG_MAX_SECONDS="${RANDOM_STRONG_MAX_SECONDS:-7200}"
run_log "$CORR_DIR/random_${RANDOM_SMOKE_TRIALS}.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_random_resumable" --impl h3compat_l20_l02 --bits all --kappa 5 --ops lt,gt,eq --trials-per-bit "$RANDOM_SMOKE_TRIALS" --seed 1 --max-seconds 1800 --checkpoint "$CORR_DIR/random_${RANDOM_SMOKE_TRIALS}.jsonl"
run_log "$CORR_DIR/random_${RANDOM_STRONG_TRIALS}.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_random_resumable" --impl h3compat_l20_l02 --bits all --kappa 5 --ops lt,gt,eq --trials-per-bit "$RANDOM_STRONG_TRIALS" --seed 2 --max-seconds "$RANDOM_STRONG_MAX_SECONDS" --checkpoint "$CORR_DIR/random_${RANDOM_STRONG_TRIALS}.jsonl"
run_log "$CORR_DIR/small_k_exhaustive.log" "${TASKSET_CMD[@]}" "$BIN/my_ethmsb_static_tests" --impl h3compat_l20_l02 --exhaustive-k 16 --kappa 5

touch "$CORR_DIR/PASS"
echo "$REPORT_DIR" > "$ROOT/reports/latest_ethmsb_report_dir.txt"
echo "ETHMSB_FORMAL_CERT_DONE report_dir=$REPORT_DIR"
