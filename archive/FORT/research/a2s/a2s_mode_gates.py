"""A2S-MODE Gates A0-A4: is there a discrete target state to infer at all?

These gates are **measurements on the existing substrate**, not the mechanism.
They exist to kill A2S-MODE cheaply, before any adaptation module is built, in
the same way the Q1 stratum gate had to run before A2S-TRACE.

* **A0 - is there anything to select?**  A per-target head on a compact label-free
  ligand basis, fitted on that target's abundant labels and evaluated on held-out
  queries, reported **per relation stratum**.  If its gain lives only where the
  transport class already wins, the premise of A2S-MODE is wrong.
* **A1 - mode sufficiency.**  Cluster source-target heads into ``M`` modes.  On
  unseen probe targets, does the best *mode* beat the single global head?  This
  is the confound that killed the group-transfer route: a "mode" that is really
  just a better global ligand model.  The ceiling is measured with **split-half
  query selection**, never with hindsight on the scored queries.
* **A2 - k-shot identifiability (decisive).**  Can ``k`` support residuals pick
  the mode the oracle would pick, above chance ``1/M``?  This is the question the
  ``tau/sigma`` arithmetic cannot settle, because discrete separation ``D_k`` is
  an object training can shape while ``tau`` is not.
* **A3 - complementarity.**  Does any gain survive in the ``t < 0.35`` strata
  where the whole transport class measures zero?  Reported together with the P0
  readout: gain as a function of support diversity vs support-query Tanimoto.
* **A4 - synthetic positive control.**  Inject a world with true discrete modes
  and confirm the pipeline recovers them.  No null may be reported without it.

Only ``fit`` and ``probe`` roles are opened.  ``locked`` and the recipient roster
are never requested.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from research.a2s.a2s_information_gate import FEATURES, PROTEINS, REGISTRY, canonical, sha256_file
from research.a2s.a2s_trace import (
    BATCH_EPISODES,
    DEVICE,
    EVALUATION_EPISODE_SEED,
    TRAIN_POLICIES,
    Substrate,
    build_episodes,
    gather_batch,
    group_by_k,
    load_substrate,
    make_batch,
    tanimoto,
)
from research.a2s.a2s_trace_stratum import (
    DEFAULT_LOCK,
    DEFAULT_OOF,
    MIN_STRATUM_QUERY,
    STRATUM_NAMES,
    metric_loss,
    paired_bootstrap,
    stratum_of,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_mode_gates_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_mode_gates_records_2026-08-02.parquet"

SEED = 20260802
PCA_COMPONENTS = 16
BASIS_DIM = PCA_COMPONENTS + 10
HEAD_RIDGE = 1.0
MIN_ROWS_FOR_DICTIONARY = 40
MIN_ROWS_FOR_TARGET_ORACLE = 24
MODE_COUNTS = (2, 3, 4, 6)
HEADLINE_M = 4
MIN_SPLIT_QUERY = 12
MDE = 0.005
BOOTSTRAP_DRAWS = 2000
ADMITTED_STRATUM = "t55_100"
NULL_STRATA = ("t00_20", "t20_35")


# --------------------------------------------------------------------------- #
# Compact label-free ligand basis
# --------------------------------------------------------------------------- #


def build_basis(substrate: Substrate) -> tuple[np.ndarray, dict[str, object]]:
    """``g(x)``: 10 standardised descriptors + the top Morgan principal components.

    Every statistic is computed on ``fit``-role rows only and no label is read,
    so the basis carries no probe information.
    """

    fit_mask = torch.as_tensor(
        (substrate.labeled.role == "fit").to_numpy(), device=DEVICE
    )
    bits = substrate.bits
    centre = bits[fit_mask].mean(dim=0, keepdim=True)
    centred_fit = bits[fit_mask] - centre
    # ``torch.svd_lowrank`` draws a random projection and takes no generator, so
    # it returns a different subspace on every call.  Every gate built on this
    # basis was therefore irreproducible.  The fit-role Gram matrix is only
    # 1024x1024, so an exact symmetric eigendecomposition is affordable and
    # removes the randomness entirely.  Eigenvector *signs* are still arbitrary,
    # so they are fixed by a deterministic convention.
    gram = (centred_fit.T @ centred_fit).double()
    eigenvalues, eigenvectors = torch.linalg.eigh(gram)
    order = torch.argsort(eigenvalues, descending=True)[:PCA_COMPONENTS]
    projection = eigenvectors[:, order]
    leading = projection.abs().argmax(dim=0)
    signs = torch.sign(projection[leading, torch.arange(projection.shape[1], device=projection.device)])
    projection = (projection * torch.where(signs == 0, torch.ones_like(signs), signs)).to(bits.dtype)
    scores = (bits - centre) @ projection
    score_mean = scores[fit_mask].mean(dim=0, keepdim=True)
    score_scale = scores[fit_mask].std(dim=0, keepdim=True).clamp(min=1e-6)
    scores = (scores - score_mean) / score_scale
    basis = torch.cat((substrate.desc, scores), dim=1).cpu().numpy().astype(np.float64)
    stats = {
        "dimension": int(basis.shape[1]),
        "pca_components": PCA_COMPONENTS,
        "descriptors": 10,
        "fitted_on": "fit-role rows only, label-free",
        "decomposition": "exact symmetric eigendecomposition, sign-fixed (deterministic)",
        "explained_variance_ratio": float(
            (eigenvalues[order].sum() / eigenvalues.clamp(min=0).sum()).item()
        ),
    }
    return basis, stats


def fit_head(
    design: np.ndarray, residual: np.ndarray, ridge: float = HEAD_RIDGE
) -> tuple[np.ndarray, float]:
    """Ridge head on the compact basis; the intercept is not penalised."""

    centre = design.mean(axis=0)
    centred = design - centre
    gram = centred.T @ centred + ridge * np.eye(design.shape[1])
    weight = np.linalg.solve(gram, centred.T @ (residual - residual.mean()))
    intercept = float(residual.mean() - centre @ weight)
    return weight, intercept


# --------------------------------------------------------------------------- #
# Mode dictionary
# --------------------------------------------------------------------------- #


@dataclass
class Dictionary:
    modes: np.ndarray            # (M, d) mode weight vectors, mode 0 excluded
    global_head: np.ndarray      # (d,) the M=1 control
    sigma: float                 # within-target residual scale after a head fit
    level_sd: float              # across-target SD of the head intercept
    n_targets: int
    assignment: dict[str, int]

    @property
    def n_modes(self) -> int:
        return int(self.modes.shape[0])


def source_heads(
    substrate: Substrate, basis: np.ndarray, *, role: str = "fit"
) -> tuple[np.ndarray, list[str], float, float]:
    frame = substrate.labeled
    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    heads: list[np.ndarray] = []
    names: list[str] = []
    intercepts: list[float] = []
    dispersion: list[float] = []
    for target, group in frame.loc[frame.role == role].groupby("target", sort=True):
        rows = group.index.to_numpy()
        if len(rows) < MIN_ROWS_FOR_DICTIONARY:
            continue
        weight, intercept = fit_head(basis[rows], residual[rows])
        heads.append(weight)
        names.append(str(target))
        intercepts.append(intercept)
        fitted = basis[rows] @ weight + intercept
        dispersion.append(float(np.std(residual[rows] - fitted)))
    if not heads:
        raise RuntimeError("no source target had enough rows for a head")
    return (
        np.vstack(heads),
        names,
        float(np.mean(dispersion)),
        float(np.std(intercepts)),
    )


def build_dictionary(
    heads: np.ndarray, names: list[str], sigma: float, level_sd: float, n_modes: int
) -> Dictionary:
    from sklearn.cluster import KMeans

    model = KMeans(n_clusters=n_modes, n_init=10, random_state=SEED)
    labels = model.fit_predict(heads)
    return Dictionary(
        modes=np.asarray(model.cluster_centers_, dtype=np.float64),
        global_head=heads.mean(axis=0),
        sigma=sigma,
        level_sd=level_sd,
        n_targets=int(len(names)),
        assignment={name: int(label) for name, label in zip(names, labels)},
    )


def mode_matrix(dictionary: Dictionary) -> np.ndarray:
    """Candidate set including the null mode as row 0."""

    return np.vstack([np.zeros(dictionary.modes.shape[1]), dictionary.modes])


# --------------------------------------------------------------------------- #
# k-shot mode evidence
# --------------------------------------------------------------------------- #


def mode_log_evidence(
    support_basis: np.ndarray,
    support_residual: np.ndarray,
    modes: np.ndarray,
    sigma: float,
    level_sd: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Empirical-Bayes log evidence of each mode given k support residuals.

    The level is shrunk, not free.  At k=1 a shrunk level still absorbs almost
    all of the single observation, which is why the rank channel is structurally
    near-silent there; that silence is a registered property, not an accident.
    """

    k = len(support_residual)
    predictions = support_basis @ modes.T                     # (k, M+1)
    gap = support_residual[:, None] - predictions
    shrink = level_sd**2 / (level_sd**2 + sigma**2 / max(k, 1))
    level = shrink * gap.mean(axis=0)                          # (M+1,)
    error = gap - level[None, :]
    evidence = -0.5 * (error**2).sum(axis=0) / sigma**2 - 0.5 * level**2 / level_sd**2
    return evidence, level


# --------------------------------------------------------------------------- #
# Episode evaluation
# --------------------------------------------------------------------------- #


def evaluate_gates(
    substrate: Substrate,
    basis: np.ndarray,
    dictionary: Dictionary,
    episodes: list,
    *,
    synthetic: np.ndarray | None = None,
) -> pd.DataFrame:
    """One row per (episode, stratum, scope).

    ``synthetic`` supplies a per-target true mode index for Gate A4; when given,
    the query labels and support residuals are replaced by that world.
    """

    frame = substrate.labeled
    residual_all = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    affinity_all = substrate.affinity.cpu().numpy().astype(np.float64)
    base_all = substrate.base.cpu().numpy().astype(np.float64)
    rows_by_target = {
        str(target): group.index.to_numpy()
        for target, group in frame.groupby("target", sort=True)
    }
    modes = mode_matrix(dictionary)
    records: list[dict[str, object]] = []

    by_k = group_by_k(episodes, substrate)
    for k, (block, support_index, query_index) in sorted(by_k.items()):
        for start in range(0, len(block), BATCH_EPISODES):
            selection = np.arange(start, min(start + BATCH_EPISODES, len(block)))
            batch = make_batch(support_index, query_index, selection)
            data = gather_batch(substrate, batch)
            nearest_all = tanimoto(data["query_bits"], data["support_bits"]).amax(-1).cpu().numpy()
            # `make_batch` trims the padded query block to this batch's widest
            # episode, so the row ids must come from the batch, not from the
            # globally padded index array.
            support_rows = batch.support.cpu().numpy()
            query_rows = batch.query.cpu().numpy()
            mask_all = data["mask"].cpu().numpy()

            for offset, position in enumerate(selection):
                episode = block[int(position)]
                queries = query_rows[offset][mask_all[offset]]
                supports = support_rows[offset]
                if len(queries) == 0:
                    continue
                query_basis = basis[queries]
                support_basis = basis[supports]
                base_query = base_all[queries]

                if synthetic is None:
                    label_query = affinity_all[queries]
                    support_residual = residual_all[supports]
                else:
                    true_mode = int(synthetic[int(rows_by_target[episode.target][0])])
                    generator = np.random.default_rng(
                        int(sha256(f"{SEED}:{episode.target}:{episode.k}:{episode.draw}".encode()).hexdigest()[:8], 16)
                    )
                    truth = modes[true_mode]
                    level = generator.normal(0.0, dictionary.level_sd)
                    label_query = (
                        base_query + level + query_basis @ truth
                        + generator.normal(0.0, dictionary.sigma, len(queries))
                    )
                    support_residual = (
                        level + support_basis @ truth
                        + generator.normal(0.0, dictionary.sigma, len(supports))
                    )

                # --- candidate corrections -------------------------------- #
                mode_delta = query_basis @ modes.T                     # (Q, M+1)
                global_delta = query_basis @ dictionary.global_head    # (Q,)
                evidence, _ = mode_log_evidence(
                    support_basis, support_residual, modes, dictionary.sigma, dictionary.level_sd
                )
                chosen = int(np.argmax(evidence))

                # --- A0: per-target oracle head on held-out queries -------- #
                target_rows = rows_by_target[episode.target]
                train_rows = np.setdiff1d(target_rows, queries, assume_unique=False)
                if len(train_rows) >= MIN_ROWS_FOR_TARGET_ORACLE:
                    weight, intercept = fit_head(basis[train_rows], residual_all[train_rows])
                    target_delta = query_basis @ weight + intercept
                else:
                    target_delta = None

                nearest = nearest_all[offset][mask_all[offset]]
                strata = stratum_of(nearest)
                parity = np.arange(len(queries)) % 2

                for stratum in (*STRATUM_NAMES, "all"):
                    active = np.ones(len(queries), dtype=bool)
                    if stratum != "all":
                        active &= strata == stratum
                    if int(active.sum()) < MIN_STRATUM_QUERY:
                        continue
                    truth_values = label_query[active]
                    if float(np.std(truth_values)) < 1e-9:
                        continue

                    scopes: list[tuple[str, np.ndarray, np.ndarray | None]] = [
                        ("full", active, None)
                    ]
                    if int(active.sum()) >= MIN_SPLIT_QUERY:
                        first = active & (parity == 0)
                        second = active & (parity == 1)
                        if first.sum() >= MIN_STRATUM_QUERY // 2 and second.sum() >= MIN_STRATUM_QUERY // 2:
                            scopes.append(("heldout", second, first))

                    for scope, evaluated, selector in scopes:
                        labels = label_query[evaluated]
                        if float(np.std(labels)) < 1e-9:
                            continue
                        entry: dict[str, object] = {
                            "policy": episode.policy,
                            "k": episode.k,
                            "target": episode.target,
                            "component": episode.component,
                            "stratum": stratum,
                            "scope": scope,
                            "n_query": int(evaluated.sum()),
                            "nearest_tanimoto_mean": float(nearest[evaluated].mean()),
                            "support_spread": float(
                                np.mean(np.std(support_basis @ modes.T, axis=0))
                            ),
                            "kshot_mode": chosen,
                        }
                        candidates = {
                            "base": base_query[evaluated],
                            "global_head": base_query[evaluated] + global_delta[evaluated],
                            "kshot_mode": base_query[evaluated] + mode_delta[evaluated, chosen],
                        }
                        if target_delta is not None:
                            candidates["target_oracle"] = base_query[evaluated] + target_delta[evaluated]
                        # Hindsight mode on the scored queries: inflated, flagged.
                        hindsight = max(
                            range(modes.shape[0]),
                            key=lambda index: metric_loss(labels, base_query[evaluated] + mode_delta[evaluated, index])["ci"],
                        )
                        candidates["mode_hindsight"] = base_query[evaluated] + mode_delta[evaluated, hindsight]
                        entry["hindsight_mode"] = int(hindsight)
                        if selector is not None:
                            selector_labels = label_query[selector]
                            honest = max(
                                range(modes.shape[0]),
                                key=lambda index: metric_loss(
                                    selector_labels, base_query[selector] + mode_delta[selector, index]
                                )["ci"],
                            )
                            candidates["mode_splithalf"] = base_query[evaluated] + mode_delta[evaluated, honest]
                            entry["splithalf_mode"] = int(honest)
                            entry["kshot_matches_splithalf"] = int(chosen == honest)
                        for name, prediction in candidates.items():
                            for metric, value in metric_loss(labels, prediction).items():
                                entry[f"{name}__{metric}"] = float(value)
                        records.append(entry)
    return pd.DataFrame.from_records(records)


# --------------------------------------------------------------------------- #
# Summaries and gate verdicts
# --------------------------------------------------------------------------- #


CONTRASTS = {
    "a0_target_oracle_minus_base": ("target_oracle", "base"),
    "a1_splithalf_minus_global": ("mode_splithalf", "global_head"),
    "a1_splithalf_minus_base": ("mode_splithalf", "base"),
    "a1_global_minus_base": ("global_head", "base"),
    "a2_kshot_minus_global": ("kshot_mode", "global_head"),
    "a2_kshot_minus_base": ("kshot_mode", "base"),
    "hindsight_minus_global": ("mode_hindsight", "global_head"),
}
METRICS = ("ci", "ndcg10", "rmse")


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for scope in sorted(records.scope.unique()):
        for k in sorted(records.k.unique()):
            for stratum in sorted(records.stratum.unique()):
                frame = records.loc[
                    (records.scope == scope) & (records.k == k) & (records.stratum == stratum)
                ].copy()
                if frame.empty:
                    continue
                cell: dict[str, object] = {
                    "episodes": int(len(frame)),
                    "targets": int(frame.target.nunique()),
                    "components": int(frame.component.nunique()),
                    "nearest_tanimoto_mean": float(frame.nearest_tanimoto_mean.mean()),
                    "absolute_ci": {},
                    "contrasts": {},
                }
                for name in ("base", "global_head", "kshot_mode", "mode_splithalf", "mode_hindsight", "target_oracle"):
                    column = f"{name}__ci"
                    if column in frame:
                        cell["absolute_ci"][name] = float(np.nanmean(frame[column].to_numpy()))
                if "kshot_matches_splithalf" in frame:
                    accuracy = frame.kshot_matches_splithalf.dropna()
                    if len(accuracy):
                        cell["kshot_selection_accuracy"] = float(accuracy.mean())
                for label, (left, right) in CONTRASTS.items():
                    for metric in METRICS:
                        left_column, right_column = f"{left}__{metric}", f"{right}__{metric}"
                        if left_column not in frame or right_column not in frame:
                            continue
                        sign = -1.0 if metric == "rmse" else 1.0
                        frame["value"] = sign * (frame[left_column] - frame[right_column])
                        cell["contrasts"].setdefault(label, {})[metric] = paired_bootstrap(
                            frame, "value", draws=BOOTSTRAP_DRAWS
                        )
                summary.setdefault(scope, {}).setdefault(f"k{int(k)}", {})[stratum] = cell
    return summary


def decide(summary: dict[str, object], n_modes: int) -> dict[str, object]:
    def cell(scope: str, k: int, stratum: str) -> dict[str, object] | None:
        return summary.get(scope, {}).get(f"k{k}", {}).get(stratum)

    gates: dict[str, object] = {}

    # A0: per-target headroom outside the transport-admitted stratum.
    a0: list[dict[str, object]] = []
    for stratum in (*NULL_STRATA, "t35_55", ADMITTED_STRATUM, "all"):
        entry = cell("full", 5, stratum)
        if entry and "a0_target_oracle_minus_base" in entry["contrasts"]:
            interval = entry["contrasts"]["a0_target_oracle_minus_base"]["ci"]
            a0.append({"stratum": stratum, **interval})
    gates["A0"] = {
        "records": a0,
        "pass": any(
            item["stratum"] in NULL_STRATA and item["lower95"] > MDE for item in a0
        ),
        "criterion": "per-target head beats base with LCB > MDE in a stratum where transport is null",
    }

    # A1: mode sufficiency, split-half selection only.
    a1: list[dict[str, object]] = []
    for k in (3, 5):
        for stratum in (*NULL_STRATA, ADMITTED_STRATUM, "all"):
            entry = cell("heldout", k, stratum)
            if entry and "a1_splithalf_minus_global" in entry["contrasts"]:
                a1.append(
                    {"k": k, "stratum": stratum, **entry["contrasts"]["a1_splithalf_minus_global"]["ci"]}
                )
    gates["A1"] = {
        "records": a1,
        "pass": any(item["lower95"] > MDE for item in a1),
        "criterion": "split-half-selected mode beats the single global head with LCB > MDE",
    }

    # A2: k-shot identifiability.
    a2: list[dict[str, object]] = []
    for k in (1, 3, 5):
        for stratum in (*NULL_STRATA, ADMITTED_STRATUM, "all"):
            entry = cell("heldout", k, stratum)
            if entry is None:
                continue
            record = {"k": k, "stratum": stratum, "chance": 1.0 / (n_modes + 1)}
            if "kshot_selection_accuracy" in entry:
                record["selection_accuracy"] = entry["kshot_selection_accuracy"]
            if "a2_kshot_minus_global" in entry["contrasts"]:
                record.update(entry["contrasts"]["a2_kshot_minus_global"]["ci"])
            a2.append(record)
    gates["A2"] = {
        "records": a2,
        "pass": any(
            item.get("lower95", -1.0) > MDE
            and item.get("selection_accuracy", 0.0) > item["chance"]
            for item in a2
            if item["k"] >= 3
        ),
        "criterion": "k-shot selection beats chance and beats the global head with LCB > MDE at k>=3",
    }

    # A3: complementarity in the strata where transport is null.
    a3 = [item for item in a2 if item["stratum"] in NULL_STRATA and item["k"] >= 3]
    gates["A3"] = {
        "records": a3,
        "pass": any(item.get("lower95", -1.0) > MDE for item in a3),
        "criterion": "k-shot mode gain survives where the transport class measures zero",
    }

    verdict = "MODE_ROUTE_ADMITTED" if (gates["A1"]["pass"] and gates["A2"]["pass"]) else "MODE_ROUTE_NOT_ADMITTED"
    return {"gates": gates, "verdict": verdict}


# --------------------------------------------------------------------------- #
# Runner
# --------------------------------------------------------------------------- #


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    if DEVICE.type != "cuda":
        raise RuntimeError("run these gates with D:\\anaconda\\envs\\drug\\python.exe")
    substrate, context = load_substrate(lock_path, oof_cache)
    basis, basis_stats = build_basis(substrate)
    heads, names, sigma, level_sd = source_heads(substrate, basis)

    probe = [
        episode
        for episode in build_episodes(substrate.labeled, "probe")
        if episode.policy in TRAIN_POLICIES and episode.seed == EVALUATION_EPISODE_SEED
    ]

    sweep: dict[str, object] = {}
    headline_records: pd.DataFrame | None = None
    headline_dictionary: Dictionary | None = None
    for n_modes in MODE_COUNTS:
        dictionary = build_dictionary(heads, names, sigma, level_sd, n_modes)
        records = evaluate_gates(substrate, basis, dictionary, probe)
        summary = summarise(records)
        sweep[f"M{n_modes}"] = {
            "dictionary": {
                "source_targets": dictionary.n_targets,
                "sigma": dictionary.sigma,
                "level_sd": dictionary.level_sd,
                "cluster_sizes": sorted(
                    pd.Series(list(dictionary.assignment.values())).value_counts().to_dict().items()
                ),
                "mode_norms": [float(np.linalg.norm(row)) for row in dictionary.modes],
                "global_head_norm": float(np.linalg.norm(dictionary.global_head)),
            },
            "decision": decide(summary, n_modes),
        }
        if n_modes == HEADLINE_M:
            headline_records = records
            headline_dictionary = dictionary
            sweep[f"M{n_modes}"]["summary"] = summary

    assert headline_records is not None and headline_dictionary is not None
    records_path.parent.mkdir(parents=True, exist_ok=True)
    headline_records.to_parquet(records_path, index=False)

    # A4: synthetic positive control with the headline dictionary.
    generator = np.random.default_rng(SEED)
    true_mode = np.zeros(len(substrate.labeled), dtype=np.int64)
    for target, group in substrate.labeled.groupby("target", sort=True):
        true_mode[group.index.to_numpy()] = int(generator.integers(1, headline_dictionary.n_modes + 1))
    control_records = evaluate_gates(
        substrate, basis, headline_dictionary, probe, synthetic=true_mode
    )
    control_summary = summarise(control_records)
    control_decision = decide(control_summary, headline_dictionary.n_modes)

    result: dict[str, object] = {
        "schema": "a2s-mode-gates-v1",
        "status": "SOURCE_ONLY_DIAGNOSTIC",
        "question": "A0-A4: does a small, separable, k-shot-identifiable discrete target state exist?",
        "labels": {
            "opened_roles": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
        },
        "device": torch.cuda.get_device_name(0),
        "lock": {"path": str(lock_path), "content_sha256": context["lock"]["content_sha256"]},
        "inputs": {
            "registry_sha256": sha256_file(REGISTRY),
            "features_sha256": sha256_file(FEATURES),
            "proteins_sha256": sha256_file(PROTEINS),
        },
        "protocol": {
            "basis": basis_stats,
            "head_ridge": HEAD_RIDGE,
            "min_rows_for_dictionary": MIN_ROWS_FOR_DICTIONARY,
            "min_rows_for_target_oracle": MIN_ROWS_FOR_TARGET_ORACLE,
            "mode_counts": list(MODE_COUNTS),
            "headline_M": HEADLINE_M,
            "policies": list(TRAIN_POLICIES),
            "episode_seed": EVALUATION_EPISODE_SEED,
            "mde": MDE,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "selection": "split-half query selection; hindsight-on-scored-queries reported but never gated on",
            "level": "empirical-Bayes shrunk; k=1 rank silence retained by design",
        },
        "episodes": {"probe": len(probe)},
        "sweep": sweep,
        "synthetic_positive_control": {
            "design": "each probe target assigned a true mode from the headline dictionary; labels regenerated",
            "decision": control_decision,
            "summary": {
                scope: {
                    k: {
                        stratum: {
                            "absolute_ci": cell["absolute_ci"],
                            "kshot_selection_accuracy": cell.get("kshot_selection_accuracy"),
                            "a2_kshot_minus_global": cell["contrasts"].get("a2_kshot_minus_global", {}).get("ci"),
                        }
                        for stratum, cell in strata.items()
                    }
                    for k, strata in by_k.items()
                }
                for scope, by_k in control_summary.items()
            },
        },
        "decision": sweep[f"M{HEADLINE_M}"]["decision"],
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, default=float) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    args = parser.parse_args()
    result = run(
        args.lock.resolve(), args.out.resolve(), args.records.resolve(), args.oof_cache.resolve()
    )
    print(
        json.dumps(
            {
                "verdict": result["decision"]["verdict"],
                "gates": {
                    name: gate["pass"] for name, gate in result["decision"]["gates"].items()
                },
                "control_verdict": result["synthetic_positive_control"]["decision"]["verdict"],
            },
            indent=2,
            sort_keys=True,
            default=float,
        )
    )


if __name__ == "__main__":
    main()
