"""Leak-free checkpoint selection: a held-out partition of `meta_train`.

Until 2026-08-18 the maintained trainer selected its checkpoint on `meta_val`
and then reported `meta_val`. Stage B measured what that costs — an identical
arm selected on `meta_val` scored 2.1246 at k=0 against the leak-free 2.7425,
resolved at k=2/3/5 — so **every** `meta_val` number produced by the standard
trainer, including the recorded incumbent band, is optimistic by about
0.62 pK^2 at k=0. That is 5.6x the largest mechanism effect measured in the
whole cycle, which is why no mechanism verdict could be read cleanly off it.

This module is the promotion of Stage B's rule
(`tools/research/stageB_complementary/train_stageb.py`) into the maintained
path, with its constants unchanged so that Stage B's recorded artifacts remain
reproducible:

* the partition is **by homology component, never by target** — two targets of
  one component are not independent evidence about an unseen component, so a
  target-level split would leak;
* training episodes are drawn only from the fit components;
* checkpoint selection reads only the withheld internal-validation components;
* `meta_val` is read once, after the candidate is frozen.

The constants below were frozen before any Stage B arm trained and must not be
retuned: changing them silently redefines what "leak-free" means for every
artifact that cites this module.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import EpisodeSpec, QPSMPData, stable_seed

# Frozen 2026-08-18, before any Stage B arm trained.
INTERNAL_VAL_FRACTION = 0.12
PARTITION_SEED = 20260818
INTERNAL_BANK_SEED = 20260819
SUPPORT_SIZES = (0, 1, 2, 3, 5)

# Stage B's bank seeds were namespaced "stageb-internal". The namespace is kept
# verbatim so a rebuilt bank is bit-identical to the recorded one.
BANK_NAMESPACE = "stageb-internal"


def partition_components(data: QPSMPData) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split meta_train components into fit / internal-validation, by component."""
    components = sorted(data.components["meta_train"])
    order = np.random.default_rng(PARTITION_SEED).permutation(len(components))
    held = max(2, int(round(INTERNAL_VAL_FRACTION * len(components))))
    internal = tuple(sorted(components[int(i)] for i in order[:held]))
    fit = tuple(sorted(components[int(i)] for i in order[held:]))
    return fit, internal


def eligible_targets(data: QPSMPData, components: tuple[str, ...],
                     needed: int) -> dict[str, tuple[str, ...]]:
    out: dict[str, tuple[str, ...]] = {}
    for component in components:
        targets = tuple(
            target for target in data.components["meta_train"][component]
            if data._unique_ligand_count(data.tasks["meta_train"][target]) >= needed)
        if targets:
            out[component] = targets
    return out


def draw_fit_episode(data: QPSMPData, components: tuple[str, ...],
                     support_size: int, query_size: int,
                     rng: np.random.Generator) -> EpisodeSpec:
    """`draw_episode` restricted to the fit components.

    Component-uniform then target-uniform, matching the incumbent sampler, so
    the only difference from production sampling is the withheld internal-
    validation components.
    """
    pool = eligible_targets(data, components, support_size + 1)
    if not pool:
        raise ValueError("no fit task can provide the requested episode")
    keys = sorted(pool)
    component = keys[int(rng.integers(len(keys)))]
    targets = pool[component]
    target = targets[int(rng.integers(len(targets)))]
    order = data._unique_ligand_order(data.tasks["meta_train"][target], rng)
    support = order[:support_size]
    query = order[support_size:support_size + min(query_size,
                                                  len(order) - support_size)]
    donors = [key for key in keys if key != component]
    donor_component = donors[int(rng.integers(len(donors)))]
    donor = pool[donor_component][int(rng.integers(len(pool[donor_component])))]
    return EpisodeSpec("meta_train", component, target, tuple(map(int, support)),
                       tuple(map(int, query)), donor)


def internal_validation_specs(data: QPSMPData, components: tuple[str, ...],
                              query_size: int = 16, draws: int = 1,
                              support_sizes: tuple[int, ...] = SUPPORT_SIZES,
                              ) -> dict[int, tuple[EpisodeSpec, ...]]:
    """A fixed nested episode bank over the withheld meta_train components."""
    max_support = max(support_sizes)
    pool = eligible_targets(data, components, max_support + 1)
    keys = sorted(pool)
    if len(keys) < 2:
        raise ValueError("internal validation needs at least two components")
    banks: dict[int, list[EpisodeSpec]] = {size: [] for size in support_sizes}
    for index, component in enumerate(keys):
        donor_component = keys[(index + 1) % len(keys)]
        for target in pool[component]:
            for draw in range(draws):
                rng = np.random.default_rng(stable_seed(
                    BANK_NAMESPACE, INTERNAL_BANK_SEED, target, draw))
                order = data._unique_ligand_order(
                    data.tasks["meta_train"][target], rng)
                query = tuple(map(int, order[
                    max_support:max_support + min(
                        query_size, len(order) - max_support)]))
                if len(query) < 2:
                    continue
                donor_targets = pool[donor_component]
                donor = donor_targets[int(rng.integers(len(donor_targets)))]
                for size in support_sizes:
                    banks[size].append(EpisodeSpec(
                        "meta_train", component, target,
                        tuple(map(int, order[:size])), query, donor))
    return {size: tuple(specs) for size, specs in banks.items()}


def selection_record(fit: tuple[str, ...], internal: tuple[str, ...],
                     rule: str) -> dict:
    """The provenance block every artifact must carry about its selection rule.

    `rule` is `"internal"` (leak-free) or `"meta_val"` (the disclosed
    diagnostic). An artifact that omits this block cannot be compared with one
    that carries it, because the two are not measuring the same thing.
    """
    leak_free = rule == "internal"
    return {
        "rule": rule,
        "leak_free": leak_free,
        "module": "scripts/internal_validation.py",
        "partition_seed": PARTITION_SEED,
        "internal_validation_fraction": INTERNAL_VAL_FRACTION,
        "fit_components": len(fit),
        "internal_validation_components": len(internal),
        "meta_val_read_during_training": not leak_free,
        "disclosure": (
            "checkpoint selected on withheld meta_train components; meta_val "
            "was not read during training"
            if leak_free else
            "DISCLOSED LEAK: checkpoint selected on meta_val, the population "
            "this run also reports. Stage B measured this rule at about "
            "+0.62 pK^2 optimism at k=0; the figures are development "
            "estimates, not held-out ones."),
    }
