#!/usr/bin/env python3
"""Run standalone TPC-H 3PBS/CKKS benchmarks and parse timing/accuracy.

The runner is intentionally resumable. Each completed run writes one JSON
metric file; results.csv is rebuilt from those metric files after every run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "bench_out" / "tpch_standalone_grid"
LOG_DIR = OUT_DIR / "logs"
METRIC_DIR = OUT_DIR / "metrics"
STATUS_PATH = OUT_DIR / "STATUS.md"
CSV_PATH = OUT_DIR / "results.csv"

ROWS = [1024, 2048, 4096, 8192, 16384]
QUERIES = ["q14", "q3", "q5", "q6"]
STANDALONE_BINARIES = {
    "q14": ROOT / "build" / "test" / "tpch_q14_3pbs",
    "q3": ROOT / "build" / "test" / "tpch_q3_3pbs",
    "q5": ROOT / "build" / "test" / "tpch_q5_3pbs",
    "q6": ROOT / "build" / "test" / "tpch_q6_3pbs",
}
BUILD_TARGETS = [
    "tpch_q14_3pbs",
    "tpch_q3_3pbs",
    "tpch_q5_3pbs",
    "tpch_q6_3pbs",
]
HE3DB_BINARY = Path(
    "/home/wrn/AE_submit/SOTA/HE3DB/build/bin/relational_queries_test_16bits"
)

TOLERANCES = {
    "q14": 1e-2,
    "q3": 0.5,
    "q5": 0.5,
    "q6": 0.5,
}


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def ensure_dirs() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    METRIC_DIR.mkdir(parents=True, exist_ok=True)


def status(message: str) -> None:
    line = f"[{now()}] {message}"
    print(line, flush=True)
    with STATUS_PATH.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def metric_path(system: str, query: str, rows: int) -> Path:
    return METRIC_DIR / f"{system}_{query}_{rows}.json"


def log_path(system: str, query: str, rows: int | str) -> Path:
    return LOG_DIR / f"{system}_{query}_{rows}.log"


def run_logged(cmd: list[str], logfile: Path, cwd: Path = ROOT) -> tuple[int, float]:
    start = time.perf_counter()
    with logfile.open("w", encoding="utf-8", errors="replace") as out:
        out.write(f"$ {' '.join(cmd)}\n")
        out.write(f"# started: {now()}\n")
        out.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            errors="replace",
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            out.write(line)
            out.flush()
        rc = proc.wait()
        wall_ms = (time.perf_counter() - start) * 1000.0
        out.write(f"# finished: {now()}\n")
        out.write(f"# returncode: {rc}\n")
        out.write(f"# wall_ms: {wall_ms:.3f}\n")
    return rc, wall_ms


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def parse_first(pattern: str, text: str) -> float | None:
    m = re.search(pattern, text, flags=re.MULTILINE)
    if not m:
        return None
    return float(m.group(1))


def parse_query_rows(text: str) -> int | None:
    m = re.search(r"Records:\s*(\d+)", text)
    return int(m.group(1)) if m else None


def parse_table_section(text: str, heading: str, query: str) -> list[float]:
    start = text.find(heading)
    if start < 0:
        return []
    if heading == "Encrypted query result:":
        end = text.find("Plain query result:", start)
    else:
        end = len(text)
    block = text[start:end if end >= 0 else len(text)]
    values: list[float] = []

    if query in {"q3", "q5"}:
        for line in block.splitlines():
            if "|" not in line:
                continue
            last = line.rsplit("|", 1)[-1].strip()
            if re.match(r"^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$", last):
                values.append(float(last))
        return values

    for line in block.splitlines():
        stripped = line.strip()
        if not stripped or stripped in {
            heading,
            "revenue",
            "promo_revenue",
            "Plain query result:",
        }:
            continue
        if stripped.startswith("("):
            continue
        m = re.match(
            r"^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$",
            stripped,
        )
        if m:
            values.append(float(stripped))
            break
    return values


def parse_metric(
    *,
    system: str,
    query: str,
    rows: int,
    rc: int | None,
    wall_ms: float | None,
    log: Path,
    text: str,
    note: str = "",
) -> dict[str, object]:
    encrypted = parse_table_section(text, "Encrypted query result:", query)
    plain = parse_table_section(text, "Plain query result:", query)
    pair_count = min(len(encrypted), len(plain))
    diffs = [abs(encrypted[i] - plain[i]) for i in range(pair_count)]
    tolerance = TOLERANCES.get(query, 0.5)
    correct_count = sum(1 for d in diffs if d <= tolerance)
    accuracy = (100.0 * correct_count / pair_count) if pair_count else None
    max_abs_error = max(diffs) if diffs else None

    return {
        "timestamp": now(),
        "system": system,
        "query": query,
        "rows": rows,
        "returncode": rc,
        "wall_ms": wall_ms,
        "reported_total_ms": parse_first(r"Query Evaluation Time:\s*([0-9.eE+-]+)\s*ms", text),
        "filtering_ms": parse_first(r"Filtering Time:\s*([0-9.eE+-]+)\s*ms", text),
        "ckks_or_aggregation_ms": (
            parse_first(r"CKKS .* Time:\s*([0-9.eE+-]+)\s*ms", text)
            or parse_first(r"Aggregation Time:\s*([0-9.eE+-]+)\s*ms", text)
        ),
        "repack_avg_error": parse_first(r"Repack(?: mask)? average error\s*=\s*([0-9.eE+-]+)", text),
        "tolerance": tolerance,
        "result_count": pair_count,
        "correct_count": correct_count if pair_count else None,
        "accuracy_percent": accuracy,
        "max_abs_error": max_abs_error,
        "encrypted_values": encrypted,
        "plain_values": plain,
        "log": str(log.relative_to(ROOT)),
        "note": note,
    }


def write_metric(metric: dict[str, object]) -> None:
    path = metric_path(
        str(metric["system"]), str(metric["query"]), int(metric["rows"])
    )
    path.write_text(json.dumps(metric, indent=2, sort_keys=True), encoding="utf-8")


def load_metrics() -> list[dict[str, object]]:
    metrics = []
    for path in sorted(METRIC_DIR.glob("*.json")):
        metrics.append(json.loads(path.read_text(encoding="utf-8")))
    order = {q: i for i, q in enumerate(QUERIES)}
    metrics.sort(
        key=lambda m: (
            str(m.get("system", "")),
            order.get(str(m.get("query", "")), 99),
            int(m.get("rows", 0)),
        )
    )
    return metrics


def fmt_float(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return str(value)
        return f"{value:.6g}"
    return value


def write_csv() -> None:
    fields = [
        "system",
        "query",
        "rows",
        "returncode",
        "wall_ms",
        "reported_total_ms",
        "filtering_ms",
        "ckks_or_aggregation_ms",
        "repack_avg_error",
        "accuracy_percent",
        "max_abs_error",
        "tolerance",
        "result_count",
        "correct_count",
        "log",
        "note",
    ]
    with CSV_PATH.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for metric in load_metrics():
            writer.writerow({field: fmt_float(metric.get(field)) for field in fields})


def build_standalone(no_build: bool) -> None:
    if no_build:
        return
    logfile = LOG_DIR / "build_standalone.log"
    cmd = ["cmake", "--build", "build", "--target", *BUILD_TARGETS, "-j2"]
    status("building standalone benchmark binaries")
    rc, wall_ms = run_logged(cmd, logfile)
    status(f"build finished rc={rc} wall_ms={wall_ms:.0f} log={logfile.relative_to(ROOT)}")
    if rc != 0:
        raise SystemExit(rc)


def run_standalone(query: str, rows: int, force: bool) -> None:
    metric = metric_path("standalone_3pbs_ckks", query, rows)
    if metric.exists() and not force:
        status(f"skip existing standalone {query} rows={rows}")
        return

    exe = STANDALONE_BINARIES[query]
    logfile = log_path("standalone_3pbs_ckks", query, rows)
    if not exe.exists():
        status(f"missing executable for {query}: {exe}")
        write_metric(
            {
                "timestamp": now(),
                "system": "standalone_3pbs_ckks",
                "query": query,
                "rows": rows,
                "returncode": None,
                "wall_ms": None,
                "reported_total_ms": None,
                "filtering_ms": None,
                "ckks_or_aggregation_ms": None,
                "repack_avg_error": None,
                "tolerance": TOLERANCES.get(query, 0.5),
                "result_count": 0,
                "correct_count": None,
                "accuracy_percent": None,
                "max_abs_error": None,
                "encrypted_values": [],
                "plain_values": [],
                "log": "",
                "note": f"missing executable: {exe}",
            }
        )
        write_csv()
        return

    status(f"run standalone {query} rows={rows}")
    rc, wall_ms = run_logged([str(exe), str(rows)], logfile)
    text = read_text(logfile)
    parsed_rows = parse_query_rows(text) or rows
    metric_data = parse_metric(
        system="standalone_3pbs_ckks",
        query=query,
        rows=parsed_rows,
        rc=rc,
        wall_ms=wall_ms,
        log=logfile,
        text=text,
    )
    write_metric(metric_data)
    write_csv()
    status(
        "done standalone "
        f"{query} rows={parsed_rows} rc={rc} "
        f"time_ms={fmt_float(metric_data['reported_total_ms'])} "
        f"acc={fmt_float(metric_data['accuracy_percent'])}% "
        f"max_err={fmt_float(metric_data['max_abs_error'])}"
    )


def parse_he3db_blocks(text: str) -> Iterable[tuple[int, str]]:
    matches = list(re.finditer(r"Records:\s*(\d+)", text))
    for idx, match in enumerate(matches):
        start = match.start()
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(text)
        yield int(match.group(1)), text[start:end]


def run_he3db_q6(binary: Path, rows: list[int], force: bool) -> None:
    system = "he3db_original"
    query = "q6"
    already = all(metric_path(system, query, row).exists() for row in rows if row < 16384)
    missing_marker = metric_path(system, query, 16384).exists()
    if already and missing_marker and not force:
        status("skip existing HE3DB original q6")
        return
    if not binary.exists():
        status(f"missing HE3DB binary: {binary}")
        for row in rows:
            write_metric(
                {
                    "timestamp": now(),
                    "system": system,
                    "query": query,
                    "rows": row,
                    "returncode": None,
                    "wall_ms": None,
                    "reported_total_ms": None,
                    "filtering_ms": None,
                    "ckks_or_aggregation_ms": None,
                    "repack_avg_error": None,
                    "tolerance": TOLERANCES[query],
                    "result_count": 0,
                    "correct_count": None,
                    "accuracy_percent": None,
                    "max_abs_error": None,
                    "encrypted_values": [],
                    "plain_values": [],
                    "log": "",
                    "note": f"missing HE3DB binary: {binary}",
                }
            )
        write_csv()
        return

    logfile = log_path(system, query, "1024_8192_original_main")
    status("run HE3DB original q6 (original main emits 1024,2048,4096,8192)")
    rc, wall_ms = run_logged([str(binary)], logfile, cwd=binary.parent)
    text = read_text(logfile)
    emitted_rows = set()
    for row, block in parse_he3db_blocks(text):
        if row not in rows:
            continue
        emitted_rows.add(row)
        metric_data = parse_metric(
            system=system,
            query=query,
            rows=row,
            rc=rc,
            wall_ms=wall_ms,
            log=logfile,
            text=block,
            note="wall_ms is for the whole HE3DB original binary run",
        )
        write_metric(metric_data)

    for row in rows:
        if row in emitted_rows:
            continue
        note = (
            "not emitted by unmodified HE3DB original main "
            "(loop condition is i < 16384)"
            if row == 16384
            else "not emitted by HE3DB original binary"
        )
        write_metric(
            {
                "timestamp": now(),
                "system": system,
                "query": query,
                "rows": row,
                "returncode": rc,
                "wall_ms": None,
                "reported_total_ms": None,
                "filtering_ms": None,
                "ckks_or_aggregation_ms": None,
                "repack_avg_error": None,
                "tolerance": TOLERANCES[query],
                "result_count": 0,
                "correct_count": None,
                "accuracy_percent": None,
                "max_abs_error": None,
                "encrypted_values": [],
                "plain_values": [],
                "log": str(logfile.relative_to(ROOT)),
                "note": note,
            }
        )
    write_csv()
    status(f"done HE3DB original q6 rc={rc} wall_ms={wall_ms:.0f}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", nargs="+", type=int, default=ROWS)
    parser.add_argument("--queries", nargs="+", choices=QUERIES, default=QUERIES)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--skip-he3db", action="store_true")
    parser.add_argument("--he3db-binary", type=Path, default=HE3DB_BINARY)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_dirs()
    STATUS_PATH.write_text(f"# TPC-H benchmark status\nstarted: {now()}\n", encoding="utf-8")
    status(f"output directory: {OUT_DIR.relative_to(ROOT)}")
    status(f"pid={os.getpid()}")
    build_standalone(args.no_build)

    for rows in args.rows:
        for query in args.queries:
            run_standalone(query, rows, args.force)

    if not args.skip_he3db:
        run_he3db_q6(args.he3db_binary, args.rows, args.force)

    write_csv()
    status(f"all requested runs finished; csv={CSV_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
