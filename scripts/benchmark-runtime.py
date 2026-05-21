#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


FIELDNAMES = ["language", "case", "n", "repeat", "seconds"]


PYTHON_CASES = {
    "rddensity_fixed": lambda rddensity, rdbwdensity, x: rddensity(x, h=[0.6, 0.6], bino_flag=False),
    "rddensity_estimated": lambda rddensity, rdbwdensity, x: rddensity(x, bino_flag=False),
    "rdbwdensity_default": lambda rddensity, rdbwdensity, x: rdbwdensity(x),
    "rdbwdensity_plugin": lambda rddensity, rdbwdensity, x: rdbwdensity(x, vce="plugin"),
}


def benchmark_data(n: int):
    import numpy as np

    index = np.arange(n, dtype=float)
    return np.linspace(-1.5, 1.5, n) + 0.05 * np.sin(index)


def time_call(func, repeats: int, warmups: int, calls_per_repeat: int) -> list[float]:
    for _ in range(warmups):
        func()
    elapsed = []
    for _ in range(repeats):
        start = time.perf_counter()
        for _ in range(calls_per_repeat):
            func()
        elapsed.append((time.perf_counter() - start) / calls_per_repeat)
    return elapsed


def benchmark_python(sizes: list[int], repeats: int, warmups: int, calls_per_repeat: int) -> list[dict[str, object]]:
    sys.path.insert(0, str(ROOT / "Python" / "rddensity" / "src"))
    from rddensity import rdbwdensity, rddensity

    rows: list[dict[str, object]] = []
    for n in sizes:
        x = benchmark_data(n)
        for case, func in PYTHON_CASES.items():
            timings = time_call(lambda: func(rddensity, rdbwdensity, x), repeats, warmups, calls_per_repeat)
            rows.extend(
                {"language": "python", "case": case, "n": n, "repeat": i + 1, "seconds": seconds}
                for i, seconds in enumerate(timings)
            )
    return rows


def benchmark_r(rscript: str, sizes: list[int], repeats: int, warmups: int, calls_per_repeat: int) -> list[dict[str, object]]:
    code = f"""
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output <- args[[2]]
sizes <- as.integer(strsplit(args[[3]], ",")[[1]])
repeats <- as.integer(args[[4]])
warmups <- as.integer(args[[5]])
calls_per_repeat <- as.integer(args[[6]])

source(file.path(repo, "R/rddensity/R/rddensity_fun.R"))
source(file.path(repo, "R/rddensity/R/rdbwdensity.R"))
source(file.path(repo, "R/rddensity/R/rddensity.R"))

benchmark_data <- function(n) {{
  index <- seq(0, n - 1)
  seq(-1.5, 1.5, length.out = n) + 0.05 * sin(index)
}}

run_case <- function(case, x) {{
  if (case == "rddensity_fixed") {{
    rddensity(x, h = c(0.6, 0.6), bino = FALSE)
  }} else if (case == "rddensity_estimated") {{
    rddensity(x, bino = FALSE)
  }} else if (case == "rdbwdensity_default") {{
    rdbwdensity(x)
  }} else if (case == "rdbwdensity_plugin") {{
    rdbwdensity(x, vce = "plugin")
  }}
}}

cases <- c("rddensity_fixed", "rddensity_estimated", "rdbwdensity_default", "rdbwdensity_plugin")
rows <- data.frame(language=character(), case=character(), n=integer(), repeat_id=integer(), seconds=double())
for (n in sizes) {{
  x <- benchmark_data(n)
  for (case in cases) {{
    for (i in seq_len(warmups)) run_case(case, x)
    for (i in seq_len(repeats)) {{
      elapsed <- system.time(for (k in seq_len(calls_per_repeat)) run_case(case, x))[["elapsed"]]
      rows <- rbind(rows, data.frame(language="r", case=case, n=n, repeat_id=i, seconds=elapsed / calls_per_repeat))
    }}
  }}
}}
names(rows)[names(rows) == "repeat_id"] <- "repeat"
write.csv(rows, output, row.names=FALSE)
"""
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "benchmark-runtime.R"
        output = Path(tmp) / "benchmark-r.csv"
        script.write_text(code, encoding="utf-8")
        subprocess.run(
            [
                rscript,
                "--vanilla",
                str(script),
                str(ROOT),
                str(output),
                ",".join(map(str, sizes)),
                str(repeats),
                str(warmups),
                str(calls_per_repeat),
            ],
            check=True,
        )
        return read_rows(output)


def benchmark_stata(stata: str, sizes: list[int], repeats: int, warmups: int, calls_per_repeat: int) -> list[dict[str, object]]:
    code = r"""
version 16
clear all
set more off
set graphics off

args repo output sizes repeats warmups calls_per_repeat
adopath ++ "`repo'/stata"
local sizes : subinstr local sizes "|" " ", all

capture program drop _bench_make_data
program define _bench_make_data
    syntax, N(integer)
    clear
    set obs `n'
    gen double x = -1.5 + 3 * (_n - 1) / (`n' - 1) + 0.05 * sin(_n - 1)
end

capture program drop _bench_time
program define _bench_time, rclass
    syntax, CASE(string) N(integer) REPEAT(integer) CALLS(integer)
    _bench_make_data, n(`n')
    timer clear 1
    timer on 1
    forvalues call = 1/`calls' {
        if "`case'" == "rddensity_fixed" {
            quietly rddensity x, h(0.6 0.6) nobinomial
        }
        else if "`case'" == "rddensity_estimated" {
            quietly rddensity x, nobinomial
        }
        else if "`case'" == "rdbwdensity_default" {
            quietly rdbwdensity x
        }
        else if "`case'" == "rdbwdensity_plugin" {
            quietly rdbwdensity x, vce(plugin)
        }
    }
    timer off 1
    quietly timer list 1
    return scalar seconds = r(t1) / `calls'
end

tempfile results
tempname handle
postfile `handle' str16 language str32 case long n int repeat double seconds using "`results'", replace

foreach n of local sizes {
    foreach case in rddensity_fixed rddensity_estimated rdbwdensity_default rdbwdensity_plugin {
        forvalues i = 1/`warmups' {
            quietly _bench_time, case("`case'") n(`n') repeat(`i') calls(`calls_per_repeat')
        }
        forvalues i = 1/`repeats' {
            quietly _bench_time, case("`case'") n(`n') repeat(`i') calls(`calls_per_repeat')
            post `handle' ("stata") ("`case'") (`n') (`i') (r(seconds))
        }
    }
}

postclose `handle'
use "`results'", clear
export delimited using "`output'", replace
"""
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "benchmark-runtime.do"
        batch_log = Path(tmp) / "benchmark-runtime.log"
        output = Path(tmp) / "benchmark-stata.csv"
        script.write_text(code, encoding="utf-8")
        subprocess.run(
            [
                stata,
                "/e",
                "do",
                str(script),
                str(ROOT),
                str(output),
                "|".join(map(str, sizes)),
                str(repeats),
                str(warmups),
                str(calls_per_repeat),
            ],
            check=True,
            cwd=tmp,
        )
        deadline = time.monotonic() + 30
        while not output.exists() and time.monotonic() < deadline:
            time.sleep(0.25)
        if not output.exists():
            log_tail = ""
            if batch_log.exists():
                lines = batch_log.read_text(encoding="utf-8", errors="replace").splitlines()
                log_tail = "\n".join(lines[-80:])
            raise RuntimeError(f"Stata benchmark did not create {output}.\n{log_tail}")
        return read_rows(output)


def read_rows(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8-sig") as file:
        return [
            {
                "language": row["language"],
                "case": row["case"],
                "n": int(row["n"]),
                "repeat": int(row["repeat"]),
                "seconds": float(row["seconds"]),
            }
            for row in csv.DictReader(file)
        ]


def write_rows(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict[str, object]]) -> None:
    grouped: dict[tuple[str, str, int], list[float]] = {}
    for row in rows:
        key = (str(row["language"]), str(row["case"]), int(row["n"]))
        grouped.setdefault(key, []).append(float(row["seconds"]))

    print("language,case,n,repeats,median_seconds,min_seconds")
    for key in sorted(grouped):
        language, case, n = key
        values = grouped[key]
        print(f"{language},{case},{n},{len(values)},{statistics.median(values):.6g},{min(values):.6g}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark rddensity runtime across supported package implementations.")
    parser.add_argument("--languages", nargs="+", choices=["python", "r", "stata"], default=["python", "r"])
    parser.add_argument("--sizes", nargs="+", type=int, default=[200, 2000, 10000])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--calls-per-repeat", type=int, default=1, help="Run each case this many times per timing repeat and report seconds per call.")
    parser.add_argument("--rscript", default=shutil.which("Rscript"))
    parser.add_argument("--stata", default=None, help="Path to Stata executable, for example StataMP-64.exe.")
    parser.add_argument("--output", type=Path, help="Optional CSV path for raw timing rows.")
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    if "python" in args.languages:
        rows.extend(benchmark_python(args.sizes, args.repeats, args.warmups, args.calls_per_repeat))
    if "r" in args.languages:
        if not args.rscript:
            raise SystemExit("Rscript was requested but was not found. Pass --rscript or remove r from --languages.")
        rows.extend(benchmark_r(args.rscript, args.sizes, args.repeats, args.warmups, args.calls_per_repeat))
    if "stata" in args.languages:
        if not args.stata:
            raise SystemExit("Stata was requested but --stata was not provided.")
        rows.extend(benchmark_stata(args.stata, args.sizes, args.repeats, args.warmups, args.calls_per_repeat))

    print_summary(rows)
    if args.output:
        write_rows(rows, args.output)
        print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
