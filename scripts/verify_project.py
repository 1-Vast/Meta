"""Run the maintained MetaSieve test suite."""
from __future__ import annotations

import sys

import pytest


def main() -> int:
    args = sys.argv[1:]
    if any(value in {"-h", "--help"} for value in args):
        print("usage: python main.py verify tests [pytest arguments]")
        print("\nRuns pytest with -q when no pytest arguments are supplied.")
        return 0
    return int(pytest.main(args or ["-q"]))


if __name__ == "__main__":
    raise SystemExit(main())
