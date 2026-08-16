"""Shared pytest configuration: the regression-suite tiers.

The default `pytest tests -q` run must stay a *contract* suite — structural,
algebraic, gradient, data and CLI gates that a change can break silently. Two
groups of tests do not belong in that budget:

``research_gate``
    Synthetic **training** gates belonging to model families whose verdict is
    already recorded as immutable evidence. They train small synthetic models
    on CPU and dominate the wall time (~314 s of a 410 s run). Their measured
    outcomes live in `report/meta_fewshot/stageR8_stronger_shape_20260816/`
    and `report/meta_fewshot/stageR13_shape_direct_20260816/RESULT.json`, so
    the scientific record does not depend on re-running them.

``slow``
    Full-corpus subprocess smokes.

Skipping a research gate is a cost decision about a **settled** question. It
is not a lowering of the admission bar: a new or reopened family must run its
own Stage 1 gates in the research tier before any real-data training, and
`scripts/run_stage.py` does exactly that.

Enable with::

    RUN_RESEARCH_GATES=1 pytest tests -q      # the closed families' training gates
    RUN_SLOW=1           pytest tests -q      # full-corpus subprocess smokes
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
