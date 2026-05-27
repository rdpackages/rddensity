#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


DEFAULT_STATA16 = r"C:\Program Files\Stata16\StataMP-64.exe"
WINDOWS_CANDIDATES = [
    DEFAULT_STATA16,
]


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

    if os.name == "nt":
        candidates.extend(WINDOWS_CANDIDATES)

    path_names = [
        "StataMP-64.exe",
        "StataSE-64.exe",
        "StataBE-64.exe",
        "stata-mp",
        "stata-se",
        "stata",
    ]
    candidates.extend(found for name in path_names if (found := shutil.which(name)))

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

    raise SystemExit("Could not find Stata 16. Pass --stata or set STATA_EXE.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Stata Mata library for rddensity.")
    parser.add_argument("--stata", help="Path to the Stata 16 executable.")
    args = parser.parse_args()

    repo_root = find_repo_root(Path.cwd())
    stata = find_stata(args.stata)
    do_file = repo_root / "scripts" / "build-stata-mlib.do"
    library = repo_root / "stata" / "lrddensity.mlib"

    command = [str(stata), "/e", "do", str(do_file), str(repo_root)]
    print(f"Repository: {repo_root}")
    print(f"Stata:      {stata}")
    subprocess.run(command, cwd=repo_root, check=True)

    if not library.exists() or library.stat().st_size == 0:
        raise SystemExit("Missing or empty Stata Mata library: " + str(library))

    for object_file in (repo_root / "stata").glob("*.mo"):
        object_file.unlink()

    print(f"Wrote {library.relative_to(repo_root)} ({library.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
