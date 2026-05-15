#!/bin/bash
# Resumable compliant Q14/Q3/Q5 runner.
# Each step writes a *.done sentinel when finished so reruns skip it.
# Logs go to bench_out/{q14,q3,q5}_compliant_16.{txt,status}.
# Status file is updated continuously so you can `cat bench_out/STATUS` anytime.

set -u
ROOT=/home/wangruoning/tfhepp-compare
OUT=$ROOT/bench_out
BIN=$ROOT/build/test
STATUS=$OUT/STATUS

ts() { date "+%Y-%m-%d %H:%M:%S"; }
echo "[$(ts)] runner started, pid=$$" > "$STATUS"

run_one() {
    local q=$1
    local sentinel=$OUT/${q}_compliant_16.done
    local log=$OUT/mine_${q}_compliant_16.txt
    if [ -f "$sentinel" ]; then
        echo "[$(ts)] $q SKIP (already done)" >> "$STATUS"
        return
    fi
    echo "[$(ts)] $q START" >> "$STATUS"
    cd "$BIN"
    "./tpch_${q}_3pbs" 16 > "$log" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        touch "$sentinel"
        echo "[$(ts)] $q DONE (rc=0)" >> "$STATUS"
    else
        echo "[$(ts)] $q FAILED (rc=$rc)" >> "$STATUS"
    fi
}

run_one q14
run_one q3
run_one q5

echo "[$(ts)] runner exit" >> "$STATUS"
