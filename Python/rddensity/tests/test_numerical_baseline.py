from pathlib import Path

import numpy as np
import pandas as pd

from rddensity import rdbwdensity, rddensity


FIXTURES = Path(__file__).resolve().parent / "fixtures"


def baseline_data() -> np.ndarray:
    index = np.arange(200, dtype=float)
    return np.linspace(-1.5, 1.5, 200) + 0.05 * np.sin(index)


def summarize_rddensity(case: str, fit) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for prefix, values in (
        ("hat", fit.hat),
        ("sd_asy", fit.sd_asy),
        ("sd_jk", fit.sd_jk),
        ("test", fit.test),
        ("h", fit.h),
        ("n", fit.n),
        ("hat_p", fit.hat_p),
        ("sd_asy_p", fit.sd_asy_p),
        ("sd_jk_p", fit.sd_jk_p),
        ("test_p", fit.test_p),
    ):
        for name, value in values.items():
            rows.append({"case": case, "metric": f"{prefix}_{name}", "value": value})
    return rows


def summarize_rdbwdensity(case: str, fit) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for row in fit.h.index:
        for col in fit.h.columns:
            rows.append({"case": case, "metric": f"h_{row}_{col}", "value": fit.h.loc[row, col]})
    for name, value in fit.n.items():
        rows.append({"case": case, "metric": f"n_{name}", "value": value})
    return rows


def compute_baseline() -> pd.DataFrame:
    x = baseline_data()
    rows = []
    rows.extend(summarize_rddensity("rddensity_fixed", rddensity(x, h=[0.6, 0.6], bino_flag=False)))
    rows.extend(summarize_rddensity("rddensity_uniform", rddensity(x, h=[0.7, 0.7], kernel="uniform", bino_flag=False)))
    rows.extend(
        summarize_rddensity(
            "rddensity_restricted",
            rddensity(x, h=[0.8, 0.8], fitselect="restricted", bino_flag=False),
        )
    )
    rows.extend(summarize_rddensity("rddensity_all", rddensity(x, h=[0.6, 0.6], useall=True, bino_flag=False)))
    rows.extend(summarize_rdbwdensity("rdbwdensity_default", rdbwdensity(x)))
    rows.extend(summarize_rdbwdensity("rdbwdensity_uniform_plugin", rdbwdensity(x, kernel="uniform", vce="plugin")))
    return pd.DataFrame(rows)


def test_numerical_baseline_matches_current_contract():
    expected = pd.read_csv(FIXTURES / "numerical-baseline.csv")
    observed = compute_baseline()

    assert observed[["case", "metric"]].equals(expected[["case", "metric"]])
    np.testing.assert_allclose(
        observed["value"].to_numpy(dtype=float),
        expected["value"].to_numpy(dtype=float),
        rtol=1e-10,
        atol=1e-10,
        equal_nan=True,
    )


def test_estimated_bandwidth_path_works_with_current_pandas():
    fit = rddensity(baseline_data(), bino_flag=False)

    assert np.isfinite(fit.h["left"])
    assert np.isfinite(fit.h["right"])
