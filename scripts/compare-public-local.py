#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import textwrap
import time
from contextlib import nullcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

NUMERIC_FIELDNAMES = ["version", "language", "case", "metric", "value", "error"]
TIMING_FIELDNAMES = ["version", "language", "case", "n", "repeat", "seconds", "error"]


CASES = [
    "rddensity_fixed",
    "rddensity_estimated",
    "rddensity_uniform",
    "rddensity_restricted",
    "rddensity_all",
    "rddensity_epanechnikov_plugin",
    "rddensity_nonzero_cutoff",
    "rddensity_masspoints_adjusted",
    "rddensity_masspoints_unadjusted",
    "rdbwdensity_default",
    "rdbwdensity_plugin",
    "rdbwdensity_uniform_plugin",
    "rdbwdensity_nonzero_cutoff",
    "rdbwdensity_masspoints_adjusted",
    "rdbwdensity_masspoints_unadjusted",
    "illustration_rddensity_default",
    "illustration_rddensity_all",
    "illustration_rddensity_restricted_plugin",
    "illustration_rdbwdensity_default",
    "illustration_rddensity_h_10_hr",
    "illustration_rddensity_uniform",
    "illustration_rddensity_bwselect_diff",
    "illustration_rddensity_h_10_15",
    "illustration_rddensity_p2_q4",
    "illustration_rddensity_c5_all",
    "illustration_rdbwdensity_p3_restricted",
    "illustration_rdbwdensity_uniform_jackknife",
    "help_rddensity_h_10_20_plugin",
    "help_rdbwdensity_plugin",
    "replication_cjm_p1_each",
    "replication_cjm_p2_each",
    "replication_cjm_p3_each",
    "replication_cjm_p1_diff",
    "replication_cjm_p2_diff",
    "replication_cjm_p3_diff",
    "replication_cit_senate_default",
    "replication_cit_polecon_default",
    "replication_ckt_chemo",
    "replication_ckt_art",
]

BENCHMARK_CASES = [
    "rddensity_fixed",
    "rddensity_estimated",
    "rddensity_epanechnikov_plugin",
    "rddensity_nonzero_cutoff",
    "rddensity_masspoints_adjusted",
    "rddensity_masspoints_unadjusted",
    "rdbwdensity_default",
    "rdbwdensity_plugin",
    "rdbwdensity_uniform_plugin",
    "rdbwdensity_nonzero_cutoff",
    "rdbwdensity_masspoints_adjusted",
    "rdbwdensity_masspoints_unadjusted",
]


PYTHON_NUMERIC = r"""
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rddensity import rdbwdensity, rddensity

version, output, repo_arg = sys.argv[1], sys.argv[2], sys.argv[3]
REPO = Path(repo_arg)
REPLICATION_ROOT = Path(r"C:\Users\cattaneo\Dropbox\software\rdpackages-replication")


def baseline_data():
    index = np.arange(200, dtype=float)
    return np.linspace(-1.5, 1.5, 200) + 0.05 * np.sin(index)


def masspoint_data():
    index = np.arange(240, dtype=float)
    return np.round(np.linspace(-1.4, 1.6, 240) + 0.08 * np.sin(index / 2), 1)


def add_value(rows, case, metric, value):
    try:
        value = float(value)
    except Exception:
        value = math.nan
    rows.append(
        {
            "version": version,
            "language": "python",
            "case": case,
            "metric": metric,
            "value": value,
            "error": "",
        }
    )


def add_error(rows, case, error):
    rows.append(
        {
            "version": version,
            "language": "python",
            "case": case,
            "metric": "__error__",
            "value": math.nan,
            "error": str(error).replace("\n", " ")[:500],
        }
    )


def iter_items(values):
    return values.items() if hasattr(values, "items") else enumerate(values)


def summarize_rddensity(rows, case, fit):
    for prefix, attr in (
        ("hat", "hat"),
        ("sd_asy", "sd_asy"),
        ("sd_jk", "sd_jk"),
        ("test", "test"),
        ("h", "h"),
        ("n", "n"),
        ("hat_p", "hat_p"),
        ("sd_asy_p", "sd_asy_p"),
        ("sd_jk_p", "sd_jk_p"),
        ("test_p", "test_p"),
        ("bino", "bino"),
    ):
        values = getattr(fit, attr, None)
        if values is None:
            continue
        for name, value in iter_items(values):
            add_value(rows, case, f"{prefix}_{name}", value)


def summarize_rdbwdensity(rows, case, fit):
    h = fit.h
    for row_name in h.index:
        for col_name in h.columns:
            add_value(rows, case, f"h_{row_name}_{col_name}", h.loc[row_name, col_name])
    for name, value in iter_items(fit.n):
        add_value(rows, case, f"n_{name}", value)


def run_case(case, x, x_masspoints, source):
    if case == "rddensity_fixed":
        return "rddensity", rddensity(x, h=[0.6, 0.6], bino_flag=False)
    if case == "rddensity_estimated":
        return "rddensity", rddensity(x, bino_flag=False)
    if case == "rddensity_uniform":
        return "rddensity", rddensity(x, h=[0.7, 0.7], kernel="uniform", bino_flag=False)
    if case == "rddensity_restricted":
        return "rddensity", rddensity(x, h=[0.8, 0.8], fitselect="restricted", bino_flag=False)
    if case == "rddensity_all":
        return "rddensity", rddensity(x, h=[0.6, 0.6], useall=True, bino_flag=False)
    if case == "rddensity_epanechnikov_plugin":
        return "rddensity", rddensity(x, h=[0.75, 0.65], kernel="epanechnikov", vce="plugin", bino_flag=False)
    if case == "rddensity_nonzero_cutoff":
        return "rddensity", rddensity(x, c=0.15, h=[0.5, 0.7], bino_flag=False)
    if case == "rddensity_masspoints_adjusted":
        return "rddensity", rddensity(x_masspoints, h=[0.8, 0.8], bino_flag=False)
    if case == "rddensity_masspoints_unadjusted":
        return "rddensity", rddensity(x_masspoints, h=[0.8, 0.8], massPoints=False, bino_flag=False)
    if case == "rdbwdensity_default":
        return "rdbwdensity", rdbwdensity(x)
    if case == "rdbwdensity_plugin":
        return "rdbwdensity", rdbwdensity(x, vce="plugin")
    if case == "rdbwdensity_uniform_plugin":
        return "rdbwdensity", rdbwdensity(x, kernel="uniform", vce="plugin")
    if case == "rdbwdensity_nonzero_cutoff":
        return "rdbwdensity", rdbwdensity(x, c=0.15, p=3, kernel="epanechnikov")
    if case == "rdbwdensity_masspoints_adjusted":
        return "rdbwdensity", rdbwdensity(x_masspoints)
    if case == "rdbwdensity_masspoints_unadjusted":
        return "rdbwdensity", rdbwdensity(x_masspoints, massPoints=False)
    if case == "illustration_rddensity_default":
        return "rddensity", rddensity(source["senate"], bino_flag=False)
    if case == "illustration_rddensity_all":
        return "rddensity", rddensity(source["senate"], useall=True, bino_flag=False)
    if case == "illustration_rddensity_restricted_plugin":
        return "rddensity", rddensity(source["senate"], fitselect="restricted", vce="plugin", bino_flag=False)
    if case == "illustration_rdbwdensity_default":
        return "rdbwdensity", rdbwdensity(source["senate"])
    if case == "illustration_rddensity_h_10_hr":
        bw = rdbwdensity(source["senate"])
        return "rddensity", rddensity(source["senate"], h=[10, bw.h.loc["r", "bw"]], bino_flag=False)
    if case == "illustration_rddensity_uniform":
        return "rddensity", rddensity(source["senate"], kernel="uniform", bino_flag=False)
    if case == "illustration_rddensity_bwselect_diff":
        return "rddensity", rddensity(source["senate"], bwselect="diff", bino_flag=False)
    if case == "illustration_rddensity_h_10_15":
        return "rddensity", rddensity(source["senate"], h=[10, 15], bino_flag=False)
    if case == "illustration_rddensity_p2_q4":
        return "rddensity", rddensity(source["senate"], p=2, q=4, bino_flag=False)
    if case == "illustration_rddensity_c5_all":
        return "rddensity", rddensity(source["senate"], c=5, useall=True, bino_flag=False)
    if case == "illustration_rdbwdensity_p3_restricted":
        return "rdbwdensity", rdbwdensity(source["senate"], p=3, fitselect="restricted")
    if case == "illustration_rdbwdensity_uniform_jackknife":
        return "rdbwdensity", rdbwdensity(source["senate"], kernel="uniform", vce="jackknife")
    if case == "help_rddensity_h_10_20_plugin":
        return "rddensity", rddensity(source["senate"], h=[10, 20], vce="plugin", bino_flag=False)
    if case == "help_rdbwdensity_plugin":
        return "rdbwdensity", rdbwdensity(source["senate"], vce="plugin")
    if case == "replication_cjm_p1_each":
        return "rddensity", rddensity(source["cjm_headstart"], p=1, bwselect="each", bino_flag=False)
    if case == "replication_cjm_p2_each":
        return "rddensity", rddensity(source["cjm_headstart"], p=2, bwselect="each", bino_flag=False)
    if case == "replication_cjm_p3_each":
        return "rddensity", rddensity(source["cjm_headstart"], p=3, bwselect="each", bino_flag=False)
    if case == "replication_cjm_p1_diff":
        return "rddensity", rddensity(source["cjm_headstart"], p=1, bwselect="diff", bino_flag=False)
    if case == "replication_cjm_p2_diff":
        return "rddensity", rddensity(source["cjm_headstart"], p=2, bwselect="diff", bino_flag=False)
    if case == "replication_cjm_p3_diff":
        return "rddensity", rddensity(source["cjm_headstart"], p=3, bwselect="diff", bino_flag=False)
    if case == "replication_cit_senate_default":
        return "rddensity", rddensity(source["cit_senate"], bino_flag=False)
    if case == "replication_cit_polecon_default":
        return "rddensity", rddensity(source["cit_polecon"], bino_flag=False)
    if case == "replication_ckt_chemo":
        return "rddensity", rddensity(source["ckt_chemo"], c=25.5, vce="jackknife", h=[5, 5], q=2, binoWStep=1)
    if case == "replication_ckt_art":
        return "rddensity", rddensity(source["ckt_art"], c=350, vce="jackknife", bino_flag=False)
    raise ValueError(case)


rows = []
x = baseline_data()
x_masspoints = masspoint_data()
source = {
    "senate": pd.read_csv(REPO / "Python" / "rddensity_senate.csv")["margin"].dropna(),
    "cjm_headstart": (pd.read_csv(REPLICATION_ROOT / "CJM_2020_JASA" / "headstart.csv")["povrate60"] - 59.198).dropna(),
    "cit_senate": pd.read_csv(REPLICATION_ROOT / "CIT_2020_CUP" / "CIT_2020_CUP_senate.csv")["demmv"].dropna(),
    "cit_polecon": pd.read_csv(REPLICATION_ROOT / "CIT_2020_CUP" / "CIT_2020_CUP_polecon.csv")["X"].dropna(),
    "ckt_chemo": pd.read_stata(REPLICATION_ROOT / "CKT_2023_SIM" / "CKT_2023_Chemo.dta")["onc_score"].dropna(),
    "ckt_art": pd.read_stata(REPLICATION_ROOT / "CKT_2023_SIM" / "CKT_2023_SIM--ART.dta")["cd4"].dropna(),
}
for case in CASES:
    try:
        kind, fit = run_case(case, x, x_masspoints, source)
        if kind == "rddensity":
            summarize_rddensity(rows, case, fit)
        else:
            summarize_rdbwdensity(rows, case, fit)
    except Exception as exc:
        add_error(rows, case, exc)

with open(output, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=NUMERIC_FIELDNAMES)
    writer.writeheader()
    writer.writerows(rows)
"""


PYTHON_BENCHMARK = r"""
import csv
import sys
import time

import numpy as np
from rddensity import rdbwdensity, rddensity

version = sys.argv[1]
output = sys.argv[2]
sizes = [int(value) for value in sys.argv[3].split(",") if value]
repeats = int(sys.argv[4])
warmups = int(sys.argv[5])
calls_per_repeat = int(sys.argv[6])


def benchmark_data(n):
    index = np.arange(n, dtype=float)
    return np.linspace(-1.5, 1.5, n) + 0.05 * np.sin(index)


def benchmark_masspoint_data(n):
    index = np.arange(n, dtype=float)
    return np.round(np.linspace(-1.4, 1.6, n) + 0.08 * np.sin(index / 2), 1)


def run_case(case, data):
    if case == "rddensity_fixed":
        return rddensity(data["continuous"], h=[0.6, 0.6], bino_flag=False)
    if case == "rddensity_estimated":
        return rddensity(data["continuous"], bino_flag=False)
    if case == "rddensity_epanechnikov_plugin":
        return rddensity(data["continuous"], h=[0.75, 0.65], kernel="epanechnikov", vce="plugin", bino_flag=False)
    if case == "rddensity_nonzero_cutoff":
        return rddensity(data["continuous"], c=0.15, h=[0.5, 0.7], bino_flag=False)
    if case == "rddensity_masspoints_adjusted":
        return rddensity(data["masspoints"], h=[0.8, 0.8], bino_flag=False)
    if case == "rddensity_masspoints_unadjusted":
        return rddensity(data["masspoints"], h=[0.8, 0.8], massPoints=False, bino_flag=False)
    if case == "rdbwdensity_default":
        return rdbwdensity(data["continuous"])
    if case == "rdbwdensity_plugin":
        return rdbwdensity(data["continuous"], vce="plugin")
    if case == "rdbwdensity_uniform_plugin":
        return rdbwdensity(data["continuous"], kernel="uniform", vce="plugin")
    if case == "rdbwdensity_nonzero_cutoff":
        return rdbwdensity(data["continuous"], c=0.15, p=3, kernel="epanechnikov")
    if case == "rdbwdensity_masspoints_adjusted":
        return rdbwdensity(data["masspoints"])
    if case == "rdbwdensity_masspoints_unadjusted":
        return rdbwdensity(data["masspoints"], massPoints=False)
    raise ValueError(case)


rows = []
for n in sizes:
    data = {
        "continuous": benchmark_data(n),
        "masspoints": benchmark_masspoint_data(n),
    }
    for case in BENCHMARK_CASES:
        error = ""
        try:
            for _ in range(warmups):
                run_case(case, data)
            for repeat in range(1, repeats + 1):
                start = time.perf_counter()
                for _ in range(calls_per_repeat):
                    run_case(case, data)
                seconds = (time.perf_counter() - start) / calls_per_repeat
                rows.append(
                    {
                        "version": version,
                        "language": "python",
                        "case": case,
                        "n": n,
                        "repeat": repeat,
                        "seconds": seconds,
                        "error": "",
                    }
                )
        except Exception as exc:
            rows.append(
                {
                    "version": version,
                    "language": "python",
                    "case": case,
                    "n": n,
                    "repeat": 0,
                    "seconds": "",
                    "error": str(exc).replace("\n", " ")[:500],
                }
            )

with open(output, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=TIMING_FIELDNAMES)
    writer.writeheader()
    writer.writerows(rows)
"""


R_COMMON = r"""
baseline_data <- function() {
  index <- seq(0, 199)
  seq(-1.5, 1.5, length.out = 200) + 0.05 * sin(index)
}

masspoint_data <- function() {
  index <- seq(0, 239)
  round(seq(-1.4, 1.6, length.out = 240) + 0.08 * sin(index / 2), 1)
}

append_value <- function(rows, version, language, case_name, metric, value) {
  rbind(
    rows,
    data.frame(
      version = version,
      language = language,
      case = case_name,
      metric = metric,
      value = as.numeric(value),
      error = "",
      stringsAsFactors = FALSE
    )
  )
}

append_error <- function(rows, version, language, case_name, error) {
  rbind(
    rows,
    data.frame(
      version = version,
      language = language,
      case = case_name,
      metric = "__error__",
      value = NA_real_,
      error = substr(gsub("[\r\n]+", " ", conditionMessage(error)), 1, 500),
      stringsAsFactors = FALSE
    )
  )
}

append_series <- function(rows, version, language, case_name, prefix, values) {
  values <- unlist(values)
  for (name in names(values)) {
    rows <- append_value(rows, version, language, case_name, paste(prefix, name, sep = "_"), values[[name]])
  }
  rows
}

summarize_rddensity <- function(rows, version, language, case_name, fit) {
  for (prefix in c("hat", "sd_asy", "sd_jk", "test", "h", "N", "hat_p", "sd_asy_p", "sd_jk_p", "test_p", "bino")) {
    metric_prefix <- if (prefix == "N") "n" else prefix
    rows <- append_series(rows, version, language, case_name, metric_prefix, fit[[prefix]])
  }
  rows
}

summarize_rdbwdensity <- function(rows, version, language, case_name, fit) {
  h <- as.data.frame(fit$h)
  for (row_name in row.names(h)) {
    for (col_name in colnames(h)) {
      rows <- append_value(rows, version, language, case_name, paste("h", row_name, col_name, sep = "_"), h[row_name, col_name])
    }
  }
  rows <- append_series(rows, version, language, case_name, "n", fit$N)
  rows
}

run_numeric_case <- function(case_name, x, x_masspoints, source_data) {
  if (case_name == "rddensity_fixed") {
    list(kind = "rddensity", fit = rddensity(x, h = c(0.6, 0.6), bino = FALSE))
  } else if (case_name == "rddensity_estimated") {
    list(kind = "rddensity", fit = rddensity(x, bino = FALSE))
  } else if (case_name == "rddensity_uniform") {
    list(kind = "rddensity", fit = rddensity(x, h = c(0.7, 0.7), kernel = "uniform", bino = FALSE))
  } else if (case_name == "rddensity_restricted") {
    list(kind = "rddensity", fit = rddensity(x, h = c(0.8, 0.8), fitselect = "restricted", bino = FALSE))
  } else if (case_name == "rddensity_all") {
    list(kind = "rddensity", fit = rddensity(x, h = c(0.6, 0.6), all = TRUE, bino = FALSE))
  } else if (case_name == "rddensity_epanechnikov_plugin") {
    list(kind = "rddensity", fit = rddensity(x, h = c(0.75, 0.65), kernel = "epanechnikov", vce = "plugin", bino = FALSE))
  } else if (case_name == "rddensity_nonzero_cutoff") {
    list(kind = "rddensity", fit = rddensity(x, c = 0.15, h = c(0.5, 0.7), bino = FALSE))
  } else if (case_name == "rddensity_masspoints_adjusted") {
    list(kind = "rddensity", fit = rddensity(x_masspoints, h = c(0.8, 0.8), bino = FALSE))
  } else if (case_name == "rddensity_masspoints_unadjusted") {
    list(kind = "rddensity", fit = rddensity(x_masspoints, h = c(0.8, 0.8), massPoints = FALSE, bino = FALSE))
  } else if (case_name == "rdbwdensity_default") {
    list(kind = "rdbwdensity", fit = rdbwdensity(x))
  } else if (case_name == "rdbwdensity_plugin") {
    list(kind = "rdbwdensity", fit = rdbwdensity(x, vce = "plugin"))
  } else if (case_name == "rdbwdensity_uniform_plugin") {
    list(kind = "rdbwdensity", fit = rdbwdensity(x, kernel = "uniform", vce = "plugin"))
  } else if (case_name == "rdbwdensity_nonzero_cutoff") {
    list(kind = "rdbwdensity", fit = rdbwdensity(x, c = 0.15, p = 3, kernel = "epanechnikov"))
  } else if (case_name == "rdbwdensity_masspoints_adjusted") {
    list(kind = "rdbwdensity", fit = rdbwdensity(x_masspoints))
  } else if (case_name == "rdbwdensity_masspoints_unadjusted") {
    list(kind = "rdbwdensity", fit = rdbwdensity(x_masspoints, massPoints = FALSE))
  } else if (case_name == "illustration_rddensity_default") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, bino = FALSE))
  } else if (case_name == "illustration_rddensity_all") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, all = TRUE, bino = FALSE))
  } else if (case_name == "illustration_rddensity_restricted_plugin") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, fitselect = "restricted", vce = "plugin", bino = FALSE))
  } else if (case_name == "illustration_rdbwdensity_default") {
    list(kind = "rdbwdensity", fit = rdbwdensity(source_data$senate))
  } else if (case_name == "illustration_rddensity_h_10_hr") {
    bw <- rdbwdensity(source_data$senate)
    list(kind = "rddensity", fit = rddensity(source_data$senate, h = c(10, as.numeric(bw$h[2, 1])), bino = FALSE))
  } else if (case_name == "illustration_rddensity_uniform") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, kernel = "uniform", bino = FALSE))
  } else if (case_name == "illustration_rddensity_bwselect_diff") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, bwselect = "diff", bino = FALSE))
  } else if (case_name == "illustration_rddensity_h_10_15") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, h = c(10, 15), bino = FALSE))
  } else if (case_name == "illustration_rddensity_p2_q4") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, p = 2, q = 4, bino = FALSE))
  } else if (case_name == "illustration_rddensity_c5_all") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, c = 5, all = TRUE, bino = FALSE))
  } else if (case_name == "illustration_rdbwdensity_p3_restricted") {
    list(kind = "rdbwdensity", fit = rdbwdensity(source_data$senate, p = 3, fitselect = "restricted"))
  } else if (case_name == "illustration_rdbwdensity_uniform_jackknife") {
    list(kind = "rdbwdensity", fit = rdbwdensity(source_data$senate, kernel = "uniform", vce = "jackknife"))
  } else if (case_name == "help_rddensity_h_10_20_plugin") {
    list(kind = "rddensity", fit = rddensity(source_data$senate, h = c(10, 20), vce = "plugin", bino = FALSE))
  } else if (case_name == "help_rdbwdensity_plugin") {
    list(kind = "rdbwdensity", fit = rdbwdensity(source_data$senate, vce = "plugin"))
  } else if (case_name == "replication_cjm_p1_each") {
    list(kind = "rddensity", fit = rddensity(source_data$cjm_headstart, p = 1, bwselect = "each", bino = FALSE))
  } else if (case_name == "replication_cjm_p2_each") {
    list(kind = "rddensity", fit = rddensity(source_data$cjm_headstart, p = 2, bwselect = "each", bino = FALSE))
  } else if (case_name == "replication_cjm_p3_each") {
    list(kind = "rddensity", fit = rddensity(source_data$cjm_headstart, p = 3, bwselect = "each", bino = FALSE))
  } else if (case_name == "replication_cjm_p1_diff") {
    list(kind = "rddensity", fit = rddensity(source_data$cjm_headstart, p = 1, bwselect = "diff", bino = FALSE))
  } else if (case_name == "replication_cjm_p2_diff") {
    list(kind = "rddensity", fit = rddensity(source_data$cjm_headstart, p = 2, bwselect = "diff", bino = FALSE))
  } else if (case_name == "replication_cjm_p3_diff") {
    list(kind = "rddensity", fit = rddensity(source_data$cjm_headstart, p = 3, bwselect = "diff", bino = FALSE))
  } else if (case_name == "replication_cit_senate_default") {
    list(kind = "rddensity", fit = rddensity(source_data$cit_senate, bino = FALSE))
  } else if (case_name == "replication_cit_polecon_default") {
    list(kind = "rddensity", fit = rddensity(source_data$cit_polecon, bino = FALSE))
  } else if (case_name == "replication_ckt_chemo") {
    list(kind = "rddensity", fit = rddensity(source_data$ckt_chemo, c = 25.5, vce = "jackknife", h = c(5, 5), q = 2, binoWStep = 1))
  } else if (case_name == "replication_ckt_art") {
    list(kind = "rddensity", fit = rddensity(source_data$ckt_art, c = 350, vce = "jackknife", bino = FALSE))
  } else {
    stop(case_name)
  }
}
"""


R_NUMERIC = R_COMMON + r"""
args <- commandArgs(trailingOnly = TRUE)
version <- args[[1]]
repo <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output <- args[[3]]

if (version == "local") {
  source(file.path(repo, "R/rddensity/R/rddensity_fun.R"))
  source(file.path(repo, "R/rddensity/R/rdbwdensity.R"))
  source(file.path(repo, "R/rddensity/R/rddensity.R"))
} else {
  suppressPackageStartupMessages(library(rddensity))
}

cases <- c(
  "rddensity_fixed",
  "rddensity_estimated",
  "rddensity_uniform",
  "rddensity_restricted",
  "rddensity_all",
  "rddensity_epanechnikov_plugin",
  "rddensity_nonzero_cutoff",
  "rddensity_masspoints_adjusted",
  "rddensity_masspoints_unadjusted",
  "rdbwdensity_default",
  "rdbwdensity_plugin",
  "rdbwdensity_uniform_plugin",
  "rdbwdensity_nonzero_cutoff",
  "rdbwdensity_masspoints_adjusted",
  "rdbwdensity_masspoints_unadjusted",
  "illustration_rddensity_default",
  "illustration_rddensity_all",
  "illustration_rddensity_restricted_plugin",
  "illustration_rdbwdensity_default",
  "illustration_rddensity_h_10_hr",
  "illustration_rddensity_uniform",
  "illustration_rddensity_bwselect_diff",
  "illustration_rddensity_h_10_15",
  "illustration_rddensity_p2_q4",
  "illustration_rddensity_c5_all",
  "illustration_rdbwdensity_p3_restricted",
  "illustration_rdbwdensity_uniform_jackknife",
  "help_rddensity_h_10_20_plugin",
  "help_rdbwdensity_plugin",
  "replication_cjm_p1_each",
  "replication_cjm_p2_each",
  "replication_cjm_p3_each",
  "replication_cjm_p1_diff",
  "replication_cjm_p2_diff",
  "replication_cjm_p3_diff",
  "replication_cit_senate_default",
  "replication_cit_polecon_default",
  "replication_ckt_chemo",
  "replication_ckt_art"
)

rows <- data.frame(version=character(), language=character(), case=character(), metric=character(), value=double(), error=character())
x <- baseline_data()
x_masspoints <- masspoint_data()
replication_root <- normalizePath("C:/Users/cattaneo/Dropbox/software/rdpackages-replication", winslash = "/", mustWork = TRUE)
source_data <- list(
  senate = read.csv(file.path(repo, "R/rddensity_senate.csv"))$margin,
  cjm_headstart = na.omit(read.csv(file.path(replication_root, "CJM_2020_JASA/headstart.csv"))$povrate60 - 59.198),
  cit_senate = na.omit(read.csv(file.path(replication_root, "CIT_2020_CUP/CIT_2020_CUP_senate.csv"))$demmv),
  cit_polecon = na.omit(read.csv(file.path(replication_root, "CIT_2020_CUP/CIT_2020_CUP_polecon.csv"))$X),
  ckt_chemo = na.omit(foreign::read.dta(file.path(replication_root, "CKT_2023_SIM/CKT_2023_Chemo.dta"))$onc_score),
  ckt_art = na.omit(foreign::read.dta(file.path(replication_root, "CKT_2023_SIM/CKT_2023_SIM--ART.dta"))$cd4)
)
for (case_name in cases) {
  tryCatch({
    result <- run_numeric_case(case_name, x, x_masspoints, source_data)
    if (result$kind == "rddensity") {
      rows <- summarize_rddensity(rows, version, "r", case_name, result$fit)
    } else {
      rows <- summarize_rdbwdensity(rows, version, "r", case_name, result$fit)
    }
  }, error = function(e) {
    rows <<- append_error(rows, version, "r", case_name, e)
  })
}
write.csv(rows, output, row.names = FALSE, quote = TRUE, na = "NaN")
"""


R_BENCHMARK = r"""
args <- commandArgs(trailingOnly = TRUE)
version <- args[[1]]
repo <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output <- args[[3]]
sizes <- as.integer(strsplit(args[[4]], ",")[[1]])
repeats <- as.integer(args[[5]])
warmups <- as.integer(args[[6]])
calls_per_repeat <- as.integer(args[[7]])

if (version == "local") {
  source(file.path(repo, "R/rddensity/R/rddensity_fun.R"))
  source(file.path(repo, "R/rddensity/R/rdbwdensity.R"))
  source(file.path(repo, "R/rddensity/R/rddensity.R"))
} else {
  suppressPackageStartupMessages(library(rddensity))
}

benchmark_data <- function(n) {
  index <- seq(0, n - 1)
  seq(-1.5, 1.5, length.out = n) + 0.05 * sin(index)
}

benchmark_masspoint_data <- function(n) {
  index <- seq(0, n - 1)
  round(seq(-1.4, 1.6, length.out = n) + 0.08 * sin(index / 2), 1)
}

run_benchmark_case <- function(case_name, x, x_masspoints) {
  if (case_name == "rddensity_fixed") {
    rddensity(x, h = c(0.6, 0.6), bino = FALSE)
  } else if (case_name == "rddensity_estimated") {
    rddensity(x, bino = FALSE)
  } else if (case_name == "rddensity_epanechnikov_plugin") {
    rddensity(x, h = c(0.75, 0.65), kernel = "epanechnikov", vce = "plugin", bino = FALSE)
  } else if (case_name == "rddensity_nonzero_cutoff") {
    rddensity(x, c = 0.15, h = c(0.5, 0.7), bino = FALSE)
  } else if (case_name == "rddensity_masspoints_adjusted") {
    rddensity(x_masspoints, h = c(0.8, 0.8), bino = FALSE)
  } else if (case_name == "rddensity_masspoints_unadjusted") {
    rddensity(x_masspoints, h = c(0.8, 0.8), massPoints = FALSE, bino = FALSE)
  } else if (case_name == "rdbwdensity_default") {
    rdbwdensity(x)
  } else if (case_name == "rdbwdensity_plugin") {
    rdbwdensity(x, vce = "plugin")
  } else if (case_name == "rdbwdensity_uniform_plugin") {
    rdbwdensity(x, kernel = "uniform", vce = "plugin")
  } else if (case_name == "rdbwdensity_nonzero_cutoff") {
    rdbwdensity(x, c = 0.15, p = 3, kernel = "epanechnikov")
  } else if (case_name == "rdbwdensity_masspoints_adjusted") {
    rdbwdensity(x_masspoints)
  } else if (case_name == "rdbwdensity_masspoints_unadjusted") {
    rdbwdensity(x_masspoints, massPoints = FALSE)
  }
}

cases <- c(
  "rddensity_fixed",
  "rddensity_estimated",
  "rddensity_epanechnikov_plugin",
  "rddensity_nonzero_cutoff",
  "rddensity_masspoints_adjusted",
  "rddensity_masspoints_unadjusted",
  "rdbwdensity_default",
  "rdbwdensity_plugin",
  "rdbwdensity_uniform_plugin",
  "rdbwdensity_nonzero_cutoff",
  "rdbwdensity_masspoints_adjusted",
  "rdbwdensity_masspoints_unadjusted"
)

rows <- data.frame(version=character(), language=character(), case=character(), n=integer(), repeat_id=integer(), seconds=double(), error=character())
for (n in sizes) {
  x <- benchmark_data(n)
  x_masspoints <- benchmark_masspoint_data(n)
  for (case_name in cases) {
    tryCatch({
      for (i in seq_len(warmups)) run_benchmark_case(case_name, x, x_masspoints)
      for (i in seq_len(repeats)) {
        elapsed <- system.time(for (k in seq_len(calls_per_repeat)) run_benchmark_case(case_name, x, x_masspoints))[["elapsed"]]
        rows <- rbind(rows, data.frame(version=version, language="r", case=case_name, n=n, repeat_id=i, seconds=elapsed / calls_per_repeat, error=""))
      }
    }, error = function(e) {
      rows <<- rbind(rows, data.frame(version=version, language="r", case=case_name, n=n, repeat_id=0, seconds=NA_real_, error=substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1, 500)))
    })
  }
}
names(rows)[names(rows) == "repeat_id"] <- "repeat"
write.csv(rows, output, row.names = FALSE, quote = TRUE, na = "")
"""


STATA_NUMERIC = r"""
version 16
clear all
set more off
set graphics off

args version_label ado_path output repo_root
adopath + "`ado_path'"
local repo_root : subinstr local repo_root "\" "/", all
local replication_root "C:/Users/cattaneo/Dropbox/software/rdpackages-replication"
global RDD_REPO_ROOT "`repo_root'"
global RDD_REPLICATION_ROOT "`replication_root'"

tempfile results
tempname handle
postfile `handle' str16 version str16 language str64 case str32 metric double value str244 error using "`results'", replace
global RDD_COMPARE_HANDLE `handle'
global RDD_VERSION_LABEL "`version_label'"

capture program drop _rdd_post_error
program define _rdd_post_error
    syntax, CASE(string) RC(integer)
    local handle "$RDD_COMPARE_HANDLE"
    local version_label "$RDD_VERSION_LABEL"
    post `handle' ("`version_label'") ("stata") ("`case'") ("__error__") (.) ("rc=`rc'")
end

capture program drop _rdd_post_density
program define _rdd_post_density
    syntax, CASE(string) [ALL]
    local handle "$RDD_COMPARE_HANDLE"
    local version_label "$RDD_VERSION_LABEL"

    local vce_suffix "jk"
    if "`e(vce)'" == "plugin" local vce_suffix "asy"

    post `handle' ("`version_label'") ("stata") ("`case'") ("hat_left") (e(f_ql)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("hat_right") (e(f_qr)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("hat_diff") (e(f_qr) - e(f_ql)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("sd_`vce_suffix'_left") (e(se_ql)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("sd_`vce_suffix'_right") (e(se_qr)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("sd_`vce_suffix'_diff") (e(se_q)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("test_t_`vce_suffix'") (e(T_q)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("test_p_`vce_suffix'") (e(pv_q)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("h_left") (e(h_l)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("h_right") (e(h_r)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_full") (e(N_l) + e(N_r)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_left") (e(N_l)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_right") (e(N_r)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_eff_left") (e(N_h_l)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_eff_right") (e(N_h_r)) ("")

    if "`all'" != "" {
        post `handle' ("`version_label'") ("stata") ("`case'") ("hat_p_left") (e(f_pl)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("hat_p_right") (e(f_pr)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("hat_p_diff") (e(f_pr) - e(f_pl)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("sd_`vce_suffix'_p_left") (e(se_pl)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("sd_`vce_suffix'_p_right") (e(se_pr)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("sd_`vce_suffix'_p_diff") (e(se_p)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("test_p_t_`vce_suffix'") (e(T_p)) ("")
        post `handle' ("`version_label'") ("stata") ("`case'") ("test_p_p_`vce_suffix'") (e(pv_p)) ("")
    }
end

capture program drop _rdd_post_bandwidth
program define _rdd_post_bandwidth
    syntax, CASE(string)
    local handle "$RDD_COMPARE_HANDLE"
    local version_label "$RDD_VERSION_LABEL"

    matrix H = e(h)
    local rownames : rownames H
    local colnames : colnames H
    forvalues i = 1/`=rowsof(H)' {
        local row : word `i' of `rownames'
        if "`row'" == "f_left" local rowkey "l"
        if "`row'" == "f_right" local rowkey "r"
        if "`row'" == "f_diff" local rowkey "diff"
        if "`row'" == "f_sum" local rowkey "sum"
        forvalues j = 1/`=colsof(H)' {
            local col : word `j' of `colnames'
            if "`col'" == "bandwidth" local colkey "bw"
            if "`col'" == "var" local colkey "variance"
            if "`col'" == "bias2" local colkey "biassq"
            post `handle' ("`version_label'") ("stata") ("`case'") ("h_`rowkey'_`colkey'") (H[`i', `j']) ("")
        }
    }
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_full") (e(N_l) + e(N_r)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_left") (e(N_l)) ("")
    post `handle' ("`version_label'") ("stata") ("`case'") ("n_right") (e(N_r)) ("")
end

capture program drop _rdd_make_data
program define _rdd_make_data
    clear
    set obs 200
    gen double x = -1.5 + 3 * (_n - 1) / 199 + 0.05 * sin(_n - 1)
end

capture program drop _rdd_make_masspoint_data
program define _rdd_make_masspoint_data
    clear
    set obs 240
    gen double x = round(-1.4 + 3 * (_n - 1) / 239 + 0.08 * sin((_n - 1) / 2), .1)
end

capture program drop _rdd_run_case
program define _rdd_run_case
    syntax, CASE(string)
    local repo_root "$RDD_REPO_ROOT"
    local replication_root "$RDD_REPLICATION_ROOT"
    if "`case'" == "rddensity_fixed" {
        _rdd_make_data
        rddensity x, h(0.6 0.6) nobinomial
    }
    else if "`case'" == "rddensity_estimated" {
        _rdd_make_data
        rddensity x, nobinomial
    }
    else if "`case'" == "rddensity_uniform" {
        _rdd_make_data
        rddensity x, h(0.7 0.7) kernel(uniform) nobinomial
    }
    else if "`case'" == "rddensity_restricted" {
        _rdd_make_data
        rddensity x, h(0.8 0.8) fitselect(restricted) nobinomial
    }
    else if "`case'" == "rddensity_all" {
        _rdd_make_data
        rddensity x, h(0.6 0.6) all nobinomial
    }
    else if "`case'" == "rddensity_epanechnikov_plugin" {
        _rdd_make_data
        rddensity x, h(0.75 0.65) kernel(epanechnikov) vce(plugin) nobinomial
    }
    else if "`case'" == "rddensity_nonzero_cutoff" {
        _rdd_make_data
        rddensity x, c(0.15) h(0.5 0.7) nobinomial
    }
    else if "`case'" == "rddensity_masspoints_adjusted" {
        _rdd_make_masspoint_data
        rddensity x, h(0.8 0.8) nobinomial
    }
    else if "`case'" == "rddensity_masspoints_unadjusted" {
        _rdd_make_masspoint_data
        rddensity x, h(0.8 0.8) nomasspoints nobinomial
    }
    else if "`case'" == "rdbwdensity_default" {
        _rdd_make_data
        rdbwdensity x
    }
    else if "`case'" == "rdbwdensity_plugin" {
        _rdd_make_data
        rdbwdensity x, vce(plugin)
    }
    else if "`case'" == "rdbwdensity_uniform_plugin" {
        _rdd_make_data
        rdbwdensity x, kernel(uniform) vce(plugin)
    }
    else if "`case'" == "rdbwdensity_nonzero_cutoff" {
        _rdd_make_data
        rdbwdensity x, c(0.15) p(3) kernel(epanechnikov)
    }
    else if "`case'" == "rdbwdensity_masspoints_adjusted" {
        _rdd_make_masspoint_data
        rdbwdensity x
    }
    else if "`case'" == "rdbwdensity_masspoints_unadjusted" {
        _rdd_make_masspoint_data
        rdbwdensity x, nomasspoints
    }
    else if "`case'" == "illustration_rddensity_default" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, nobinomial
    }
    else if "`case'" == "illustration_rddensity_all" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, all nobinomial
    }
    else if "`case'" == "illustration_rddensity_restricted_plugin" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, fitselect(restricted) vce(plugin) nobinomial
    }
    else if "`case'" == "illustration_rdbwdensity_default" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rdbwdensity margin
    }
    else if "`case'" == "illustration_rddensity_h_10_hr" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        quietly rdbwdensity margin
        matrix H = e(h)
        local hr = H[2,1]
        rddensity margin, h(10 `hr') nobinomial
    }
    else if "`case'" == "illustration_rddensity_uniform" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, kernel(uniform) nobinomial
    }
    else if "`case'" == "illustration_rddensity_bwselect_diff" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, bwselect(diff) nobinomial
    }
    else if "`case'" == "illustration_rddensity_h_10_15" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, h(10 15) nobinomial
    }
    else if "`case'" == "illustration_rddensity_p2_q4" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, p(2) q(4) nobinomial
    }
    else if "`case'" == "illustration_rddensity_c5_all" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, c(5) all nobinomial
    }
    else if "`case'" == "illustration_rdbwdensity_p3_restricted" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rdbwdensity margin, p(3) fitselect(restricted)
    }
    else if "`case'" == "illustration_rdbwdensity_uniform_jackknife" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rdbwdensity margin, kernel(uniform) vce(jackknife)
    }
    else if "`case'" == "help_rddensity_h_10_20_plugin" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rddensity margin, h(10 20) vce(plugin) nobinomial
    }
    else if "`case'" == "help_rdbwdensity_plugin" {
        use "`repo_root'/stata/rddensity_senate.dta", clear
        rdbwdensity margin, vce(plugin)
    }
    else if "`case'" == "replication_cjm_p1_each" {
        use "`replication_root'/CJM_2020_JASA/headstart.dta", clear
        gen double x = povrate60 - 59.198
        rddensity x, p(1) bwselect(each) nobinomial
    }
    else if "`case'" == "replication_cjm_p2_each" {
        use "`replication_root'/CJM_2020_JASA/headstart.dta", clear
        gen double x = povrate60 - 59.198
        rddensity x, p(2) bwselect(each) nobinomial
    }
    else if "`case'" == "replication_cjm_p3_each" {
        use "`replication_root'/CJM_2020_JASA/headstart.dta", clear
        gen double x = povrate60 - 59.198
        rddensity x, p(3) bwselect(each) nobinomial
    }
    else if "`case'" == "replication_cjm_p1_diff" {
        use "`replication_root'/CJM_2020_JASA/headstart.dta", clear
        gen double x = povrate60 - 59.198
        rddensity x, p(1) bwselect(diff) nobinomial
    }
    else if "`case'" == "replication_cjm_p2_diff" {
        use "`replication_root'/CJM_2020_JASA/headstart.dta", clear
        gen double x = povrate60 - 59.198
        rddensity x, p(2) bwselect(diff) nobinomial
    }
    else if "`case'" == "replication_cjm_p3_diff" {
        use "`replication_root'/CJM_2020_JASA/headstart.dta", clear
        gen double x = povrate60 - 59.198
        rddensity x, p(3) bwselect(diff) nobinomial
    }
    else if "`case'" == "replication_cit_senate_default" {
        use "`replication_root'/CIT_2020_CUP/CIT_2020_CUP_senate.dta", clear
        rename demmv x
        rddensity x, nobinomial
    }
    else if "`case'" == "replication_cit_polecon_default" {
        use "`replication_root'/CIT_2020_CUP/CIT_2020_CUP_polecon.dta", clear
        rddensity X, nobinomial
    }
    else if "`case'" == "replication_ckt_chemo" {
        use "`replication_root'/CKT_2023_SIM/CKT_2023_Chemo.dta", clear
        rddensity onc_score, c(25.5) vce(jackknife) h(5 5) bino_wstep(1 1) q(2)
    }
    else if "`case'" == "replication_ckt_art" {
        use "`replication_root'/CKT_2023_SIM/CKT_2023_SIM--ART.dta", clear
        rddensity cd4, c(350) vce(jackknife) nobinomial
    }
end

local cases rddensity_fixed rddensity_estimated rddensity_uniform rddensity_restricted ///
    rddensity_all rddensity_epanechnikov_plugin rddensity_nonzero_cutoff ///
    rddensity_masspoints_adjusted rddensity_masspoints_unadjusted ///
    rdbwdensity_default rdbwdensity_plugin rdbwdensity_uniform_plugin rdbwdensity_nonzero_cutoff ///
    rdbwdensity_masspoints_adjusted rdbwdensity_masspoints_unadjusted ///
    illustration_rddensity_default illustration_rddensity_all illustration_rddensity_restricted_plugin ///
    illustration_rdbwdensity_default illustration_rddensity_h_10_hr illustration_rddensity_uniform ///
    illustration_rddensity_bwselect_diff illustration_rddensity_h_10_15 illustration_rddensity_p2_q4 ///
    illustration_rddensity_c5_all illustration_rdbwdensity_p3_restricted illustration_rdbwdensity_uniform_jackknife ///
    help_rddensity_h_10_20_plugin help_rdbwdensity_plugin ///
    replication_cjm_p1_each replication_cjm_p2_each replication_cjm_p3_each ///
    replication_cjm_p1_diff replication_cjm_p2_diff replication_cjm_p3_diff ///
    replication_cit_senate_default replication_cit_polecon_default replication_ckt_chemo replication_ckt_art

foreach case of local cases {
    capture noisily _rdd_run_case, case("`case'")
    if _rc {
        _rdd_post_error, case("`case'") rc(`=_rc')
    }
    else if strpos("`case'", "rddensity") > 0 & strpos("`case'", "rdbwdensity") == 0 {
        if inlist("`case'", "rddensity_all", "illustration_rddensity_all", "illustration_rddensity_c5_all") _rdd_post_density, case("`case'") all
        else _rdd_post_density, case("`case'")
    }
    else {
        _rdd_post_bandwidth, case("`case'")
    }
}

postclose `handle'
use "`results'", clear
format value %21.16g
export delimited using "`output'", replace
"""


STATA_BENCHMARK = r"""
version 16
clear all
set more off
set graphics off

args version_label ado_path output sizes repeats warmups calls_per_repeat
adopath + "`ado_path'"
local sizes : subinstr local sizes "|" " ", all

capture program drop _bench_make_data
program define _bench_make_data
    syntax, N(integer)
    clear
    set obs `n'
    gen double x = -1.5 + 3 * (_n - 1) / (`n' - 1) + 0.05 * sin(_n - 1)
end

capture program drop _bench_make_masspoint_data
program define _bench_make_masspoint_data
    syntax, N(integer)
    clear
    set obs `n'
    gen double x = round(-1.4 + 3 * (_n - 1) / (`n' - 1) + 0.08 * sin((_n - 1) / 2), .1)
end

capture program drop _bench_time
program define _bench_time, rclass
    syntax, CASE(string) N(integer) CALLS(integer)
    if inlist("`case'", "rddensity_masspoints_adjusted", "rddensity_masspoints_unadjusted", "rdbwdensity_masspoints_adjusted", "rdbwdensity_masspoints_unadjusted") {
        _bench_make_masspoint_data, n(`n')
    }
    else {
        _bench_make_data, n(`n')
    }
    timer clear 1
    timer on 1
    forvalues call = 1/`calls' {
        if "`case'" == "rddensity_fixed" {
            quietly rddensity x, h(0.6 0.6) nobinomial
        }
        else if "`case'" == "rddensity_estimated" {
            quietly rddensity x, nobinomial
        }
        else if "`case'" == "rddensity_epanechnikov_plugin" {
            quietly rddensity x, h(0.75 0.65) kernel(epanechnikov) vce(plugin) nobinomial
        }
        else if "`case'" == "rddensity_nonzero_cutoff" {
            quietly rddensity x, c(0.15) h(0.5 0.7) nobinomial
        }
        else if "`case'" == "rddensity_masspoints_adjusted" {
            quietly rddensity x, h(0.8 0.8) nobinomial
        }
        else if "`case'" == "rddensity_masspoints_unadjusted" {
            quietly rddensity x, h(0.8 0.8) nomasspoints nobinomial
        }
        else if "`case'" == "rdbwdensity_default" {
            quietly rdbwdensity x
        }
        else if "`case'" == "rdbwdensity_plugin" {
            quietly rdbwdensity x, vce(plugin)
        }
        else if "`case'" == "rdbwdensity_uniform_plugin" {
            quietly rdbwdensity x, kernel(uniform) vce(plugin)
        }
        else if "`case'" == "rdbwdensity_nonzero_cutoff" {
            quietly rdbwdensity x, c(0.15) p(3) kernel(epanechnikov)
        }
        else if "`case'" == "rdbwdensity_masspoints_adjusted" {
            quietly rdbwdensity x
        }
        else if "`case'" == "rdbwdensity_masspoints_unadjusted" {
            quietly rdbwdensity x, nomasspoints
        }
    }
    timer off 1
    quietly timer list 1
    return scalar seconds = r(t1) / `calls'
end

tempfile results
tempname handle
postfile `handle' str16 version str16 language str64 case long n int repeat double seconds str244 error using "`results'", replace

local cases rddensity_fixed rddensity_estimated rddensity_epanechnikov_plugin ///
    rddensity_nonzero_cutoff rddensity_masspoints_adjusted ///
    rddensity_masspoints_unadjusted rdbwdensity_default rdbwdensity_plugin ///
    rdbwdensity_uniform_plugin rdbwdensity_nonzero_cutoff ///
    rdbwdensity_masspoints_adjusted rdbwdensity_masspoints_unadjusted

foreach n of local sizes {
    foreach case of local cases {
        forvalues i = 1/`warmups' {
            capture quietly _bench_time, case("`case'") n(`n') calls(`calls_per_repeat')
        }
        forvalues i = 1/`repeats' {
            capture quietly _bench_time, case("`case'") n(`n') calls(`calls_per_repeat')
            if _rc {
                post `handle' ("`version_label'") ("stata") ("`case'") (`n') (`i') (.) ("rc=`=_rc'")
            }
            else {
                post `handle' ("`version_label'") ("stata") ("`case'") (`n') (`i') (r(seconds)) ("")
            }
        }
    }
}

postclose `handle'
use "`results'", clear
export delimited using "`output'", replace
"""


STATA_INSTALL_PUBLIC = r"""
version 16
clear all
set more off

args install_dir from_url
sysdir set PLUS "`install_dir'"
net install rddensity, from("`from_url'") replace
"""


def run(command: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None, timeout: int | None = None) -> None:
    subprocess.run(command, cwd=cwd, env=env, check=True, timeout=timeout)


def write_script(path: Path, code: str) -> Path:
    path.write_text(
        code.replace("BENCHMARK_CASES", repr(BENCHMARK_CASES)).replace("CASES", repr(CASES)).replace(
            "NUMERIC_FIELDNAMES", repr(NUMERIC_FIELDNAMES)
        ).replace("TIMING_FIELDNAMES", repr(TIMING_FIELDNAMES)),
        encoding="utf-8",
    )
    return path


def python_public_path(tmp: Path, version: str) -> Path:
    target = tmp / "python-public"
    if not target.exists():
        run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--target",
                str(target),
                f"rddensity=={version}",
                "lpdensity<3",
                "pandas<3",
            ],
            cwd=tmp,
            timeout=300,
        )
    return target


def run_python(script: Path, version: str, output: Path, *, public_path: Path | None = None) -> None:
    env = os.environ.copy()
    if public_path is None:
        env["PYTHONPATH"] = str(ROOT / "Python" / "rddensity" / "src")
    else:
        env["PYTHONPATH"] = str(public_path)
    run([sys.executable, str(script), version, str(output), str(ROOT)], cwd=script.parent, env=env, timeout=300)


def run_python_benchmark(
    script: Path,
    version: str,
    output: Path,
    sizes: list[int],
    repeats: int,
    warmups: int,
    calls_per_repeat: int,
    *,
    public_path: Path | None = None,
) -> None:
    env = os.environ.copy()
    if public_path is None:
        env["PYTHONPATH"] = str(ROOT / "Python" / "rddensity" / "src")
    else:
        env["PYTHONPATH"] = str(public_path)
    run(
        [
            sys.executable,
            str(script),
            version,
            str(output),
            ",".join(map(str, sizes)),
            str(repeats),
            str(warmups),
            str(calls_per_repeat),
        ],
        cwd=script.parent,
        env=env,
        timeout=600,
    )


def run_r(rscript: str, script: Path, version: str, output: Path) -> None:
    run([rscript, "--vanilla", str(script), version, str(ROOT), str(output)], cwd=script.parent, timeout=300)


def run_r_benchmark(
    rscript: str,
    script: Path,
    version: str,
    output: Path,
    sizes: list[int],
    repeats: int,
    warmups: int,
    calls_per_repeat: int,
) -> None:
    run(
        [
            rscript,
            "--vanilla",
            str(script),
            version,
            str(ROOT),
            str(output),
            ",".join(map(str, sizes)),
            str(repeats),
            str(warmups),
            str(calls_per_repeat),
        ],
        cwd=script.parent,
        timeout=900,
    )


def stata_install_public(stata: str, script: Path, install_dir: Path, from_url: str) -> None:
    install_dir.mkdir(parents=True, exist_ok=True)
    run([stata, "/e", "do", str(script), str(install_dir), from_url], cwd=script.parent, timeout=300)


def run_stata(stata: str, script: Path, version: str, ado_path: Path, output: Path, timeout: int = 300) -> None:
    run([stata, "/e", "do", str(script), version, str(ado_path), str(output), str(ROOT)], cwd=script.parent, timeout=timeout)
    wait_for(output, script.with_suffix(".log"))


def run_stata_benchmark(
    stata: str,
    script: Path,
    version: str,
    ado_path: Path,
    output: Path,
    sizes: list[int],
    repeats: int,
    warmups: int,
    calls_per_repeat: int,
) -> None:
    run(
        [
            stata,
            "/e",
            "do",
            str(script),
            version,
            str(ado_path),
            str(output),
            "|".join(map(str, sizes)),
            str(repeats),
            str(warmups),
            str(calls_per_repeat),
        ],
        cwd=script.parent,
        timeout=1200,
    )
    wait_for(output, script.with_suffix(".log"))


def wait_for(path: Path, log_path: Path, seconds: int = 30) -> None:
    deadline = time.monotonic() + seconds
    while not path.exists() and time.monotonic() < deadline:
        time.sleep(0.25)
    if path.exists():
        return
    tail = ""
    log_paths = [log_path]
    log_paths.extend(sorted(log_path.parent.glob("*.log")))
    seen = set()
    for candidate in log_paths:
        if candidate in seen or not candidate.exists():
            continue
        seen.add(candidate)
        lines = candidate.read_text(encoding="utf-8", errors="replace").splitlines()
        tail += f"\n--- {candidate.name} ---\n" + "\n".join(lines[-80:])
    raise RuntimeError(f"Expected output was not created: {path}\n{tail}")


def read_numeric(paths: list[Path]) -> list[dict[str, object]]:
    rows = []
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as file:
            for row in csv.DictReader(file):
                value = row.get("value", "")
                try:
                    parsed = float(value) if value not in {"", "."} else math.nan
                except ValueError:
                    parsed = math.nan
                rows.append(
                    {
                        "version": row["version"],
                        "language": row["language"],
                        "case": row["case"],
                        "metric": row["metric"],
                        "value": parsed,
                        "error": row.get("error", ""),
                    }
                )
    return rows


def read_timings(paths: list[Path]) -> list[dict[str, object]]:
    rows = []
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as file:
            for row in csv.DictReader(file):
                seconds = row.get("seconds", "")
                rows.append(
                    {
                        "version": row["version"],
                        "language": row["language"],
                        "case": row["case"],
                        "n": int(row["n"]),
                        "repeat": int(row["repeat"]),
                        "seconds": float(seconds) if seconds not in {"", "."} else math.nan,
                        "error": row.get("error", ""),
                    }
                )
    return rows


def is_close(left: float, right: float, atol: float, rtol: float) -> bool:
    if math.isnan(left) and math.isnan(right):
        return True
    if math.isnan(left) or math.isnan(right):
        return False
    return math.isclose(left, right, abs_tol=atol, rel_tol=rtol)


def compare_numeric(rows: list[dict[str, object]], atol: float, rtol: float) -> list[dict[str, object]]:
    output = []
    for language in sorted({str(row["language"]) for row in rows}):
        public = {
            (str(row["case"]), str(row["metric"])): float(row["value"])
            for row in rows
            if row["language"] == language and row["version"] == "public" and row["metric"] != "__error__"
        }
        local = {
            (str(row["case"]), str(row["metric"])): float(row["value"])
            for row in rows
            if row["language"] == language and row["version"] == "local" and row["metric"] != "__error__"
        }
        for case, metric in sorted(set(public) & set(local)):
            old = public[(case, metric)]
            new = local[(case, metric)]
            abs_diff = abs(new - old) if not (math.isnan(old) or math.isnan(new)) else math.nan
            rel_diff = abs_diff / max(abs(old), 1e-300) if not math.isnan(abs_diff) else math.nan
            output.append(
                {
                    "language": language,
                    "case": case,
                    "metric": metric,
                    "public": old,
                    "local": new,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_diff,
                    "passed": is_close(old, new, atol, rtol),
                }
            )
    return output


def summarize_numeric(rows: list[dict[str, object]], comparisons: list[dict[str, object]]) -> None:
    print("\nNUMERICAL: public vs local")
    for language in sorted({str(row["language"]) for row in rows}):
        errors = [
            row
            for row in rows
            if row["language"] == language and row["metric"] == "__error__" and str(row.get("error", "")).strip()
        ]
        language_rows = [row for row in comparisons if row["language"] == language]
        failures = [row for row in language_rows if not row["passed"]]
        max_abs = max((row["abs_diff"] for row in language_rows if not math.isnan(row["abs_diff"])), default=0.0)
        max_rel = max((row["rel_diff"] for row in language_rows if not math.isnan(row["rel_diff"])), default=0.0)
        print(
            f"{language}: compared {len(language_rows)} common metrics; "
            f"failures={len(failures)}; errors={len(errors)}; max_abs={max_abs:.3g}; max_rel={max_rel:.3g}"
        )
        by_case: dict[str, list[dict[str, object]]] = {}
        for failure in failures:
            by_case.setdefault(str(failure["case"]), []).append(failure)
        for case in sorted(by_case):
            case_rows = by_case[case]
            case_abs = max(row["abs_diff"] for row in case_rows if not math.isnan(row["abs_diff"]))
            print(f"  {case}: {len(case_rows)} differing metrics, max_abs={case_abs:.3g}")
        for error in errors[:8]:
            print(f"  ERROR {error['version']} {error['case']}: {error['error']}")

    largest = sorted(
        [row for row in comparisons if not math.isnan(row["abs_diff"])],
        key=lambda row: row["abs_diff"],
        reverse=True,
    )[:12]
    if largest:
        print("\nLargest public/local absolute differences")
        print("language,case,metric,public,local,abs_diff")
        for row in largest:
            print(
                f"{row['language']},{row['case']},{row['metric']},"
                f"{row['public']:.12g},{row['local']:.12g},{row['abs_diff']:.3g}"
            )


def summarize_speed(rows: list[dict[str, object]]) -> None:
    print("\nTIMING: public vs local median seconds")
    grouped: dict[tuple[str, str, str, int], list[float]] = {}
    for row in rows:
        if row.get("error") or math.isnan(float(row["seconds"])):
            continue
        key = (str(row["language"]), str(row["version"]), str(row["case"]), int(row["n"]))
        grouped.setdefault(key, []).append(float(row["seconds"]))

    print("language,case,n,public,local,speedup_public_over_local")
    ratios = []
    for language in sorted({key[0] for key in grouped}):
        for case in sorted({key[2] for key in grouped if key[0] == language}):
            for n in sorted({key[3] for key in grouped if key[0] == language and key[2] == case}):
                public_values = grouped.get((language, "public", case, n), [])
                local_values = grouped.get((language, "local", case, n), [])
                if not public_values or not local_values:
                    continue
                public_median = statistics.median(public_values)
                local_median = statistics.median(local_values)
                ratio = public_median / local_median if local_median > 0 else math.inf
                ratios.append((language, case, n, ratio))
                print(f"{language},{case},{n},{public_median:.6g},{local_median:.6g},{ratio:.3g}")

    if ratios:
        print("\nSpeedup rollup by language")
        for language in sorted({row[0] for row in ratios}):
            values = [ratio for lang, _, _, ratio in ratios if lang == language and math.isfinite(ratio)]
            print(f"{language}: median={statistics.median(values):.3g}x; best={max(values):.3g}x; worst={min(values):.3g}x")

    errors = [row for row in rows if row.get("error")]
    for error in errors[:12]:
        print(f"TIMING ERROR {error['language']} {error['version']} {error['case']} n={error['n']}: {error['error']}")


def copy_outputs(rows: list[dict[str, object]], output: Path, fieldnames: list[str]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare public rddensity releases against the local checkout.")
    parser.add_argument("--mode", choices=["numeric", "speed", "all"], default="all")
    parser.add_argument("--python-public-version", default="2.4.6")
    parser.add_argument("--rscript", default=shutil.which("Rscript"))
    parser.add_argument("--stata", default=r"C:\Program Files\StataNow19\StataMP-64.exe")
    parser.add_argument("--stata-public-from", default="https://raw.githubusercontent.com/rdpackages/rddensity/main/stata")
    parser.add_argument("--languages", nargs="+", choices=["python", "r", "stata"], default=["python", "r", "stata"])
    parser.add_argument("--sizes", nargs="+", type=int, default=[200, 2000, 10000])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--calls-per-repeat", type=int, default=3)
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--work-dir", type=Path, help="Optional directory to keep temporary scripts/logs.")
    args = parser.parse_args()

    if args.rscript is None:
        raise SystemExit("Rscript was not found.")
    if args.stata and not Path(args.stata).exists():
        raise SystemExit(f"Stata executable was not found: {args.stata}")

    tmp_context = nullcontext(args.work_dir) if args.work_dir else tempfile.TemporaryDirectory(prefix="rddensity-public-local-")
    with tmp_context as tmp_name:
        tmp = Path(tmp_name)
        tmp.mkdir(parents=True, exist_ok=True)
        py_numeric = write_script(tmp / "numeric_python.py", PYTHON_NUMERIC)
        py_benchmark = write_script(tmp / "benchmark_python.py", PYTHON_BENCHMARK)
        r_numeric = write_script(tmp / "numeric_r.R", R_NUMERIC)
        r_benchmark = write_script(tmp / "benchmark_r.R", R_BENCHMARK)
        stata_numeric = write_script(tmp / "numeric_stata.do", STATA_NUMERIC)
        stata_benchmark = write_script(tmp / "benchmark_stata.do", STATA_BENCHMARK)
        stata_install = write_script(tmp / "install_public_stata.do", STATA_INSTALL_PUBLIC)

        public_python = python_public_path(tmp, args.python_public_version) if "python" in args.languages else None
        public_stata = tmp / "stata-public"
        if "stata" in args.languages and args.stata:
            stata_install_public(args.stata, stata_install, public_stata, args.stata_public_from)
        public_stata_ado = public_stata / "r" if (public_stata / "r").exists() else public_stata

        numeric_paths = []
        if args.mode in {"numeric", "all"}:
            if "python" in args.languages:
                run_python(py_numeric, "public", tmp / "python-public-numeric.csv", public_path=public_python)
                run_python(py_numeric, "local", tmp / "python-local-numeric.csv")
            if "r" in args.languages:
                run_r(args.rscript, r_numeric, "public", tmp / "r-public-numeric.csv")
                run_r(args.rscript, r_numeric, "local", tmp / "r-local-numeric.csv")
            if "stata" in args.languages and args.stata:
                run_stata(args.stata, stata_numeric, "public", public_stata_ado, tmp / "stata-public-numeric.csv")
                run_stata(args.stata, stata_numeric, "local", ROOT / "stata", tmp / "stata-local-numeric.csv")
            numeric_paths = sorted(tmp.glob("*-numeric.csv"))
            numeric_rows = read_numeric(numeric_paths)
            comparisons = compare_numeric(numeric_rows, args.atol, args.rtol)
            summarize_numeric(numeric_rows, comparisons)
            if args.output_dir:
                copy_outputs(numeric_rows, args.output_dir / "public-local-numeric-raw.csv", NUMERIC_FIELDNAMES)
                copy_outputs(comparisons, args.output_dir / "public-local-numeric-comparison.csv", list(comparisons[0].keys()) if comparisons else [])

        if args.mode in {"speed", "all"}:
            timing_paths = []
            if "python" in args.languages:
                run_python_benchmark(
                    py_benchmark,
                    "public",
                    tmp / "python-public-timing.csv",
                    args.sizes,
                    args.repeats,
                    args.warmups,
                    args.calls_per_repeat,
                    public_path=public_python,
                )
                run_python_benchmark(
                    py_benchmark,
                    "local",
                    tmp / "python-local-timing.csv",
                    args.sizes,
                    args.repeats,
                    args.warmups,
                    args.calls_per_repeat,
                )
            if "r" in args.languages:
                run_r_benchmark(args.rscript, r_benchmark, "public", tmp / "r-public-timing.csv", args.sizes, args.repeats, args.warmups, args.calls_per_repeat)
                run_r_benchmark(args.rscript, r_benchmark, "local", tmp / "r-local-timing.csv", args.sizes, args.repeats, args.warmups, args.calls_per_repeat)
            if "stata" in args.languages and args.stata:
                run_stata_benchmark(args.stata, stata_benchmark, "public", public_stata_ado, tmp / "stata-public-timing.csv", args.sizes, args.repeats, args.warmups, args.calls_per_repeat)
                run_stata_benchmark(args.stata, stata_benchmark, "local", ROOT / "stata", tmp / "stata-local-timing.csv", args.sizes, args.repeats, args.warmups, args.calls_per_repeat)
            timing_paths = sorted(tmp.glob("*-timing.csv"))
            timing_rows = read_timings(timing_paths)
            summarize_speed(timing_rows)
            if args.output_dir:
                copy_outputs(timing_rows, args.output_dir / "public-local-timing-raw.csv", TIMING_FIELDNAMES)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
