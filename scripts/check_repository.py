"""Repository-level checks for the rdhte multi-language package repo."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def require_path(relative_path: str) -> Path:
    path = ROOT / relative_path
    if not path.exists():
        fail(f"Missing required path: {relative_path}")
    return path


def read_text(relative_path: str) -> str:
    return require_path(relative_path).read_text(encoding="utf-8")


def check_required_layout() -> None:
    required_paths = [
        "README.md",
        "LICENSE.md",
        "CHANGELOG.md",
        "R/rdhte/DESCRIPTION",
        "R/rdhte/NAMESPACE",
        "R/rdhte_illustration.R",
        "R/rdhte_dataset.csv",
        "R/rdhte/data/rdhte_dataset.rda",
        "stata/stata.toc",
        "stata/rdhte.pkg",
        "scripts/build-stata-help-pdfs.py",
    ]

    for relative_path in required_paths:
        require_path(relative_path)


def read_pyproject_string(pyproject: str, key: str) -> str:
    match = re.search(rf"^{re.escape(key)}\s*=\s*\"([^\"]+)\"\s*$", pyproject, re.MULTILINE)
    if not match:
        fail(f"Python/rdhte/pyproject.toml is missing project {key!r}.")
    return match.group(1).strip()


def check_python_metadata() -> str | None:
    pyproject = ROOT / "Python" / "rdhte" / "pyproject.toml"
    if not pyproject.exists():
        return None

    require_path("Python/rdhte/pyproject.toml")
    require_path("Python/rdhte/src/rdhte/__init__.py")
    require_path("Python/rdhte/tests")

    pyproject_text = pyproject.read_text(encoding="utf-8")
    name = read_pyproject_string(pyproject_text, "name")
    version = read_pyproject_string(pyproject_text, "version")

    if name != "rdhte":
        fail(f"Unexpected Python package name: {name!r}")
    if not version:
        fail("Python package version is empty.")

    return version


def check_r_metadata() -> str:
    description = read_text("R/rdhte/DESCRIPTION")
    package_match = re.search(r"^Package:\s*(\S+)\s*$", description, re.MULTILINE)
    version_match = re.search(r"^Version:\s*(\S+)\s*$", description, re.MULTILINE)

    if not package_match:
        fail("R/rdhte/DESCRIPTION is missing Package.")
    if package_match.group(1) != "rdhte":
        fail(f"Unexpected R package name: {package_match.group(1)!r}")
    if not version_match:
        fail("R/rdhte/DESCRIPTION is missing Version.")

    return version_match.group(1)


def check_stata_manifest() -> str:
    package_text = read_text("stata/rdhte.pkg")
    distribution_match = re.search(r"^d Distribution-Date:\s*(\d{8})\s*$", package_text, re.MULTILINE)

    if not distribution_match:
        fail("stata/rdhte.pkg is missing Distribution-Date.")

    missing_files = []
    for line in package_text.splitlines():
        if not line.startswith("f "):
            continue
        relative_file = line[2:].strip()
        if not (ROOT / "stata" / relative_file).is_file():
            missing_files.append(relative_file)

    if missing_files:
        missing = ", ".join(missing_files)
        fail(f"Files listed in stata/rdhte.pkg are missing: {missing}")

    return distribution_match.group(1)


def main() -> None:
    check_required_layout()
    r_version = check_r_metadata()
    python_version = check_python_metadata()
    stata_distribution_date = check_stata_manifest()

    print("Repository checks passed.")
    print(f"R package version: {r_version}")
    if python_version is None:
        print("Python package version: not present yet")
    else:
        print(f"Python package version: {python_version}")
    print(f"Stata distribution date: {stata_distribution_date}")


if __name__ == "__main__":
    main()
