#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_value(value: str) -> float:
    if value.strip().lower() in {"nan", "."}:
        return math.nan
    return float(value)


def load_rows(path: Path, language: str | None = None, precision: str | None = None) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.DictReader(file))

    normalized = []
    for row in rows:
        normalized.append(
            {
                "language": row.get("language") or language,
                "precision": row.get("precision") or precision,
                "case": row["case"],
                "metric": row["metric"],
                "value": parse_value(row["value"]),
            }
        )
    return normalized


def as_map(rows: list[dict[str, object]]) -> dict[tuple[str, str], float]:
    return {
        (str(row["case"]), str(row["metric"])): float(row["value"])
        for row in rows
    }


def is_close(expected: float, actual: float, atol: float, rtol: float) -> bool:
    if math.isnan(expected) and math.isnan(actual):
        return True
    if math.isnan(expected) or math.isnan(actual):
        return False
    return math.isclose(expected, actual, abs_tol=atol, rel_tol=rtol)


def compare_source(
    reference: dict[tuple[str, str], float],
    rows: list[dict[str, object]],
    *,
    atol: float,
    rtol: float,
    tolerance_for=None,
) -> tuple[list[dict[str, object]], dict[str, float]]:
    source = as_map(rows)
    comparisons = []
    max_abs = 0.0
    max_rel = 0.0
    failures = 0

    for key in sorted(set(reference) & set(source)):
        expected = reference[key]
        actual = source[key]
        metric_atol, metric_rtol = tolerance_for(key[1], atol, rtol) if tolerance_for else (atol, rtol)
        abs_diff = abs(actual - expected) if not (math.isnan(actual) or math.isnan(expected)) else math.nan
        rel_diff = abs_diff / max(abs(expected), 1e-300) if not math.isnan(abs_diff) else math.nan
        passed = is_close(expected, actual, metric_atol, metric_rtol)
        failures += 0 if passed else 1
        if not math.isnan(abs_diff):
            max_abs = max(max_abs, abs_diff)
        if not math.isnan(rel_diff):
            max_rel = max(max_rel, rel_diff)
        comparisons.append(
            {
                "language": rows[0]["language"],
                "precision": rows[0]["precision"],
                "case": key[0],
                "metric": key[1],
                "reference": expected,
                "value": actual,
                "abs_diff": abs_diff,
                "rel_diff": rel_diff,
                "passed": passed,
            }
        )

    return comparisons, {"count": len(comparisons), "failures": failures, "max_abs": max_abs, "max_rel": max_rel}


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare rddensity numerical baselines across Python, R, and Stata.")
    parser.add_argument("--python", default=ROOT / "Python" / "rddensity" / "tests" / "fixtures" / "numerical-baseline.csv", type=Path)
    parser.add_argument("--r", default=ROOT / "R" / "rddensity" / "tests" / "testthat" / "fixtures" / "numerical-baseline.csv", type=Path)
    parser.add_argument("--stata", default=ROOT / "stata" / "tests" / "fixtures" / "numerical-baseline.csv", type=Path)
    parser.add_argument("--output", default=ROOT / "stata" / "tests" / "fixtures" / "numerical-comparison.csv", type=Path)
    parser.add_argument("--r-atol", default=1e-8, type=float)
    parser.add_argument("--r-rtol", default=1e-8, type=float)
    parser.add_argument("--stata-double-atol", default=1e-8, type=float)
    parser.add_argument("--stata-double-rtol", default=1e-8, type=float)
    parser.add_argument("--stata-single-atol", default=5e-6, type=float)
    parser.add_argument("--stata-single-rtol", default=5e-6, type=float)
    parser.add_argument("--stata-single-neff-atol", default=10, type=float)
    args = parser.parse_args()

    python_rows = load_rows(args.python, "python", "double")
    r_rows = load_rows(args.r, "r", "double")
    stata_rows = load_rows(args.stata)
    reference = as_map(python_rows)

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    grouped[("r", "double")].extend(r_rows)
    for row in stata_rows:
        grouped[(str(row["language"]), str(row["precision"]))].append(row)

    all_comparisons: list[dict[str, object]] = []
    any_failures = False
    for key, rows in sorted(grouped.items()):
        language, precision = key
        if language == "r":
            atol, rtol = args.r_atol, args.r_rtol
        elif precision == "single":
            atol, rtol = args.stata_single_atol, args.stata_single_rtol
        else:
            atol, rtol = args.stata_double_atol, args.stata_double_rtol

        def tolerance_for(metric: str, base_atol: float, base_rtol: float) -> tuple[float, float]:
            # Single-precision Stata can move observations exactly at bandwidth
            # boundaries, so effective local counts need a count-scale tolerance.
            if language == "stata" and precision == "single" and metric in {"n_eff_left", "n_eff_right"}:
                return max(base_atol, args.stata_single_neff_atol), base_rtol
            return base_atol, base_rtol

        comparisons, summary = compare_source(reference, rows, atol=atol, rtol=rtol, tolerance_for=tolerance_for)
        all_comparisons.extend(comparisons)
        any_failures = any_failures or summary["failures"] > 0
        print(
            f"{language}/{precision}: compared {summary['count']} common metrics; "
            f"failures={summary['failures']}; max_abs={summary['max_abs']:.3g}; max_rel={summary['max_rel']:.3g}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as file:
        fieldnames = ["language", "precision", "case", "metric", "reference", "value", "abs_diff", "rel_diff", "passed"]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_comparisons)
    print(f"Wrote {args.output.relative_to(ROOT)}")

    return 1 if any_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
