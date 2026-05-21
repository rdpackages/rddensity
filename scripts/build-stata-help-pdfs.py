#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


WINDOWS_CANDIDATES = [
    r"C:\Program Files\StataNow19\StataMP-64.exe",
    r"C:\Program Files\StataNow19\StataSE-64.exe",
    r"C:\Program Files\StataNow19\StataBE-64.exe",
    r"C:\Program Files\Stata19\StataMP-64.exe",
    r"C:\Program Files\Stata19\StataSE-64.exe",
    r"C:\Program Files\Stata19\StataBE-64.exe",
]

HELP_FILES = (
    "rddensity",
    "rdbwdensity",
)


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for path in (start, *start.parents):
        if (path / "stata" / "rddensity.pkg").exists():
            return path
    raise SystemExit("Could not find rddensity repository root.")


def find_stata(user_value: str | None) -> Path:
    candidates: list[str] = []
    if user_value:
        candidates.append(user_value)
    if env_value := os.environ.get("STATA_EXE"):
        candidates.append(env_value)

    path_names = [
        "StataMP-64.exe",
        "StataSE-64.exe",
        "StataBE-64.exe",
        "stata-mp",
        "stata-se",
        "stata",
    ]
    candidates.extend(found for name in path_names if (found := shutil.which(name)))
    if os.name == "nt":
        candidates.extend(WINDOWS_CANDIDATES)

    seen: set[Path] = set()
    for candidate in candidates:
        path = Path(candidate)
        try:
            resolved = path.resolve()
        except OSError:
            resolved = path
        if resolved in seen:
            continue
        seen.add(resolved)
        if path.exists():
            return path

    raise SystemExit("Could not find Stata. Pass --stata or set STATA_EXE.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Stata help PDFs from .sthlp files.")
    parser.add_argument("--stata", help="Path to the Stata executable.")
    args = parser.parse_args()

    repo_root = find_repo_root(Path.cwd())
    stata = find_stata(args.stata)
    do_file = repo_root / "scripts" / "build-stata-help-pdfs.do"
    stata_dir = repo_root / "stata"

    command = [str(stata), "/e", "do", str(do_file), str(repo_root)]
    print(f"Repository: {repo_root}")
    print(f"Stata:      {stata}")
    subprocess.run(command, cwd=repo_root, check=True)

    missing_or_empty = []
    for stem in HELP_FILES:
        pdf = stata_dir / f"{stem}.pdf"
        if not pdf.exists() or pdf.stat().st_size == 0:
            missing_or_empty.append(pdf.name)
        else:
            print(f"Wrote {pdf.relative_to(repo_root)} ({pdf.stat().st_size:,} bytes)")

    if missing_or_empty:
        raise SystemExit("Missing or empty PDF(s): " + ", ".join(missing_or_empty))

    print("Stata help PDF generation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
