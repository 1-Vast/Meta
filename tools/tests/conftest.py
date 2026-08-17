"""Shared pytest configuration: the regression-suite tiers.

The default `pytest tools/tests -q` run must stay a *contract* suite — structural,
algebraic, gradient, data and CLI gates that a change can break silently. Two
groups of tests do not belong in that budget:

``research_gate``
    Synthetic training gates for retained experimental families. Closed-family
    implementations and their dedicated tests are no longer kept in the active
    tree; their verdicts remain in the R-series reports and Git history.

``slow``
    Full-corpus subprocess smokes.

Skipping a research gate is a cost decision about a **settled** question. It
is not a lowering of the admission bar: a new or reopened family must run its
own Stage 1 gates in the research tier before any real-data training, and
`scripts/run_stage.py` does exactly that.

Enable with::

    RUN_RESEARCH_GATES=1 pytest tools/tests -q  # retained training gates
    RUN_SLOW=1 pytest tools/tests -q            # full-corpus subprocess smokes
"""
from __future__ import annotations

import os

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "research_gate: synthetic training gate for an already-decided model "
        "family; opt in with RUN_RESEARCH_GATES=1")
    config.addinivalue_line(
        "markers",
        "slow: loads the full corpus or spawns a subprocess; opt in with "
        "RUN_SLOW=1")


def pytest_collection_modifyitems(config: pytest.Config,
                                  items: list[pytest.Item]) -> None:
    if os.environ.get("RUN_RESEARCH_GATES") == "1":
        return
    skip = pytest.mark.skip(
        reason="research gate for a decided family: set RUN_RESEARCH_GATES=1 "
               "to run (verdict retained in report/meta_fewshot/)")
    for item in items:
        if "research_gate" in item.keywords:
            item.add_marker(skip)
