"""Support-conditioned adaptation operators on the PSEP substrate.

The question is no longer whether a target-specific object exists -- M4 settled
that (+0.0190 over a target-agnostic head of identical capacity, -0.0062 for a
wrong-target head).  The question is whether a *learned operator*

    Delta_t(x_q) = A_phi(S_t, x_q)

can convert k <= 5 support measurements into query-dependent corrections better
than the closed-form adaptations that need k ~ 20 to clear the MDE.

Why an operator rather than a latent state.  M2/M3 measured the task family to be
near-isotropic: 553 heads in 266 dimensions with participation ratio 115, rank-16
retention 29.6 %.  A low-dimensional `z_t` therefore cannot carry the object by
construction.  An operator that maps the support *set* to a *function* is not
bound by that: it never has to name a coordinate system for the task.

Structural predictions this design lets us test rather than assume:

  * Any operator of pure transport form `Delta = sum_i alpha_i(x_q) e_i` is a
    **structural no-op at k=1**: with one support point the convex weights are 1,
    the correction is a constant, and within-document concordance is invariant to
    constants.  `attention` must therefore score exactly 0 at k=1.
  * A decoder that reads `x_q` alongside the pooled support -- `cnp` -- is *not*
    so bound and may act at k=1.  The gap between them at k=1 is a clean read on
    whether anything beyond transport is being learned.

Arms
  base        Delta = 0                                     (support-free)
  intercept   Delta = mean(e_S)                             (calibration)
  krr         Delta = k_q^T (K+lambda I)^{-1} e_S, Tanimoto (kernel adaptation)
  ridge       closed-form ridge in the fixed 266-d basis     (what M2 measured)
  r2d2        meta-learned encoder + closed-form ridge to a meta-learned prior
              (the differentiable ANIL/MetaOptNet limit; also the Stage-1 lever)
  attention   alpha_i(x_q) = softmax(<u(x_q), u(x_i)>)       (learned transport)
  cnp         z_i = h(x_i, e_i); r_q = xattn(x_q, {z_i}); Delta = g(x_q, r_q)

Controls, run for every arm at every k
  correct | random_target | protein_hard | chemistry_matched | label_permuted

Components are split inside `discover` only: 60 % meta-train, 20 % meta-val
(early stopping), 20 % meta-test (reported).  `validate` and `confirm` are never
opened.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits
from research.psep.psep_m0 import rich_basis

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_operator_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_operator_records_2026-08-02.parquet"

CAPACITY = 266
K_SWEEP = (1, 3, 5)
EVAL_EPISODES = 24
TRAIN_EPISODES_PER_EPOCH = 4
MAX_EPOCHS = 40
PATIENCE = 6
LEARNING_RATE = 2e-3
WEIGHT_DECAY = 1e-5
BATCH_UNITS = 32
MAX_PAIRS = 1024
TAU = 0.5
EMB = 64
HIDDEN = 128
BOOTSTRAP_DRAWS = 2000
MIN_SUPPORT_POOL = 20
CONTROLS = ("correct", "random_target", "protein_hard", "chemistry_matched", "label_permuted")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# --------------------------------------------------------------------------- #
# Episodes
# --------------------------------------------------------------------------- #


@dataclass
class Unit:
    unit: str
    component: str
    endpoint: str
    accession: str
    role: str
    support_pool: np.ndarray      # training rows (documents disjoint from query)
    query: np.ndarray
    pair_left: np.ndarray         # within-document query pairs
    pair_right: np.ndarray


def within_document_pairs(documents: np.ndarray, label: np.ndarray, limit: int) -> tuple[np.ndarray, np.ndarray]:
    left, right = np.triu_indices(len(documents), k=1)
    keep = (documents[left] == documents[right]) & (label[left] != label[right])
    left, right = left[keep], right[keep]
    if len(left) > limit:
        rng = np.random.default_rng(SEED)
        pick = rng.choice(len(left), size=limit, replace=False)
        left, right = left[pick], right[pick]
    return left, right


def build_units(substrate, splits) -> list[Unit]:
    documents = substrate.rows.docs.astype(str).to_numpy()
    label = substrate.affinity
    units: list[Unit] = []
    for split in splits:
        if split.regime != "separated" or len(split.train) < MIN_SUPPORT_POOL:
            continue
        left, right = within_document_pairs(documents[split.evaluation], label[split.evaluation], MAX_PAIRS)
        if len(left) < 8:
            continue
        row = substrate.rows.iloc[split.evaluation[0]]
        units.append(Unit(
            split.unit, split.component, split.endpoint, str(row.accession), "",
            split.train, split.evaluation, left, right,
        ))
    return units


def assign_meta_role(component: str) -> str:
    draw = (int(sha256(f"{SEED}:meta:{component}".encode()).hexdigest()[:8], 16) % 10_000) / 10_000.0
    if draw < 0.60:
        return "meta_train"
    if draw < 0.80:
        return "meta_val"
    return "meta_test"


# --------------------------------------------------------------------------- #
# Operators
# --------------------------------------------------------------------------- #


class Attention(nn.Module):
    """Delta(x_q) = s * sum_i softmax_i(<u(x_q), u(x_i)>/T) * e_i.

    Pure transport: the correction is a convex combination of support residuals,
    so it is a structural no-op at k=1.
    """

    def __init__(self, dimension: int):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(dimension, HIDDEN), nn.GELU(), nn.Linear(HIDDEN, EMB)
        )
        self.log_temperature = nn.Parameter(torch.zeros(1))
        self.log_scale = nn.Parameter(torch.zeros(1))

    def forward(self, query: torch.Tensor, support: torch.Tensor, residual: torch.Tensor) -> torch.Tensor:
        u_query = F.normalize(self.encoder(query), dim=-1)
        u_support = F.normalize(self.encoder(support), dim=-1)
        weights = ((u_query @ u_support.T) * self.log_temperature.exp().clamp(max=30.0)).softmax(dim=-1)
        return self.log_scale.exp() * (weights @ residual)


class ConditionalNeuralProcess(nn.Module):
    """z_i = h([x_i, e_i]); r_q = crossattn(x_q, {x_i}, {z_i}); Delta = g([x_q, r_q]).

    The decoder reads the query alongside the pooled support, so the correction is
    not restricted to the span of the support residuals and can act at k=1.
    """

    def __init__(self, dimension: int):
        super().__init__()
        self.context = nn.Sequential(
            nn.Linear(dimension + 1, HIDDEN), nn.GELU(), nn.Linear(HIDDEN, EMB)
        )
        self.query_projection = nn.Linear(dimension, EMB)
        self.key_projection = nn.Linear(dimension, EMB)
        self.decoder = nn.Sequential(
            nn.Linear(dimension + EMB, HIDDEN), nn.GELU(), nn.Linear(HIDDEN, 1)
        )
        self.log_scale = nn.Parameter(torch.zeros(1))

    def forward(self, query: torch.Tensor, support: torch.Tensor, residual: torch.Tensor) -> torch.Tensor:
        context = self.context(torch.cat([support, residual[:, None]], dim=-1))
        weights = (
            (self.query_projection(query) @ self.key_projection(support).T) / np.sqrt(EMB)
        ).softmax(dim=-1)
        pooled = weights @ context
        return self.log_scale.exp() * self.decoder(torch.cat([query, pooled], dim=-1)).squeeze(-1)


class MetaRidge(nn.Module):
    """Meta-learned encoder + closed-form ridge to a meta-learned prior.

    w_t = w0 + Phi_S^T (Phi_S Phi_S^T + lambda I)^{-1} (e_S - Phi_S w0)

    This is the linear-head limit of ANIL and the R2D2/MetaOptNet construction:
    the inner loop is exact and differentiable, so the encoder, the prior `w0` and
    the regulariser are all meta-learned.  It is the direct test of whether a
    *learned representation* concentrates the head where the fixed basis does not.
    """

    def __init__(self, dimension: int):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(dimension, HIDDEN), nn.GELU(), nn.Linear(HIDDEN, EMB)
        )
        self.prior = nn.Parameter(torch.zeros(EMB))
        self.log_lambda = nn.Parameter(torch.zeros(1))
        self.log_scale = nn.Parameter(torch.zeros(1))

    def forward(self, query: torch.Tensor, support: torch.Tensor, residual: torch.Tensor) -> torch.Tensor:
        phi_support = self.encoder(support)
        phi_query = self.encoder(query)
        centre = phi_support.mean(dim=0, keepdim=True)
        phi_support = phi_support - centre
        phi_query = phi_query - centre
        target = residual - residual.mean() - phi_support @ self.prior
        gram = phi_support @ phi_support.T
        ridge = self.log_lambda.exp().clamp(1e-3, 1e6)
        identity = torch.eye(gram.shape[0], device=gram.device, dtype=gram.dtype)
        dual = torch.linalg.solve(gram + ridge * identity, target)
        weight = self.prior + phi_support.T @ dual
        return self.log_scale.exp() * (phi_query @ weight)


MODELS = {"attention": Attention, "cnp": ConditionalNeuralProcess, "r2d2": MetaRidge}


# --------------------------------------------------------------------------- #
# Closed-form baselines
# --------------------------------------------------------------------------- #


def tanimoto(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    intersection = left @ right.T
    sizes_left = left.sum(axis=1, keepdims=True)
    sizes_right = right.sum(axis=1, keepdims=True).T
    return intersection / np.maximum(sizes_left + sizes_right - intersection, 1e-9)


def baseline_delta(
    name: str, query_basis: np.ndarray, support_basis: np.ndarray, residual: np.ndarray,
    query_bits: np.ndarray, support_bits: np.ndarray, ridge: float, scale: float,
) -> np.ndarray:
    if name == "base":
        return np.zeros(len(query_basis))
    if name == "intercept":
        return np.full(len(query_basis), residual.mean())
    if name == "krr":
        gram = tanimoto(support_bits, support_bits)
        cross = tanimoto(query_bits, support_bits)
        centred = residual - residual.mean()
        alpha = np.linalg.solve(gram + ridge * np.eye(len(gram)), centred)
        return residual.mean() + scale * (cross @ alpha)
    if name == "ridge":
        centre = support_basis.mean(axis=0)
        centred = support_basis - centre
        gram = centred @ centred.T
        alpha = np.linalg.solve(gram + ridge * np.eye(len(gram)), residual - residual.mean())
        return residual.mean() + scale * ((query_basis - centre) @ (centred.T @ alpha))
    raise KeyError(name)


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #


def concordance(prediction: np.ndarray, label: np.ndarray, left: np.ndarray, right: np.ndarray) -> float:
    difference = np.sign(prediction[left] - prediction[right])
    truth = np.sign(label[left] - label[right])
    return float((difference == truth).mean() + 0.5 * (difference == 0).mean())


def ndcg_at_k(prediction: np.ndarray, label: np.ndarray, k: int = 10) -> float:
    if len(label) < 2:
        return float("nan")
    relevance = label - label.min()
    order = np.argsort(-prediction)[:k]
    ideal = np.argsort(-relevance)[:k]
    discount = 1.0 / np.log2(np.arange(2, min(k, len(label)) + 2))
    gain = float((relevance[order] * discount[: len(order)]).sum())
    best = float((relevance[ideal] * discount[: len(ideal)]).sum())
    return gain / best if best > 0 else float("nan")


def unit_metrics(prediction: np.ndarray, label: np.ndarray, unit: Unit) -> dict[str, float]:
    from scipy.stats import pearsonr, spearmanr

    metrics = {
        "ci_within": concordance(prediction, label, unit.pair_left, unit.pair_right),
        "rmse": float(np.sqrt(np.mean((prediction - label) ** 2))),
        "mae": float(np.mean(np.abs(prediction - label))),
        "ndcg10": ndcg_at_k(prediction, label),
    }
    if len(label) > 2 and np.std(prediction) > 1e-9:
        metrics["pearson"] = float(pearsonr(prediction, label)[0])
        metrics["spearman"] = float(spearmanr(prediction, label)[0])
    else:
        metrics["pearson"] = float("nan")
        metrics["spearman"] = float("nan")
    return metrics


def component_bootstrap(frame: pd.DataFrame, column: str, draws: int = BOOTSTRAP_DRAWS) -> dict[str, float]:
    usable = frame[["component", "unit", column]].dropna()
    if usable.empty:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    per_unit = usable.groupby(["component", "unit"], sort=True)[column].mean().reset_index()
    per_component = per_unit.groupby("component", sort=True)[column].mean().to_numpy(dtype=np.float64)
    per_component = per_component[np.isfinite(per_component)]
    if len(per_component) == 0:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    rng = np.random.default_rng(SEED)
    sample = rng.integers(0, len(per_component), size=(draws, len(per_component)))
    means = per_component[sample].mean(axis=1)
    return {
        "components": int(len(per_component)),
        "units": int(len(per_unit)),
        "mean": float(per_component.mean()),
        "lower95": float(np.quantile(means, 0.025)),
        "upper95": float(np.quantile(means, 0.975)),
    }


# --------------------------------------------------------------------------- #
# Donor maps for the wrong-support controls
# --------------------------------------------------------------------------- #


def sequence_similarity(accessions: list[str]) -> np.ndarray:
    """4-mer Jaccard between unit accessions -- the same measure that defines
    homology components, reused to pick the *hardest* wrong protein."""

    from sklearn.feature_extraction.text import CountVectorizer

    proteins = pd.read_csv(
        ROOT / "dataset" / "public" / "papyrus_05_7" / "raw" / "05.7_combined_set_protein_targets.tsv.xz",
        sep="\t", low_memory=False, usecols=["target_id", "Sequence"],
    )
    lookup: dict[str, str] = {}
    for target_id, sequence in zip(proteins.target_id, proteins.Sequence):
        if isinstance(sequence, str):
            lookup.setdefault(str(target_id).split("_")[0], sequence)
    sequences = [lookup.get(accession, "X") for accession in accessions]
    vectoriser = CountVectorizer(analyzer="char", ngram_range=(4, 4), binary=True, lowercase=False)
    matrix = vectoriser.fit_transform(sequences).astype(np.float32)
    matrix.data[:] = 1.0
    intersection = (matrix @ matrix.T).toarray()
    sizes = np.asarray(matrix.sum(axis=1)).ravel()
    return intersection / np.maximum(sizes[:, None] + sizes[None, :] - intersection, 1.0)


def donor_maps(units: list[Unit], bits: np.ndarray) -> dict[str, np.ndarray]:
    """For each unit, a donor index per control.  Donors always sit in a
    different homology component, so no control can leak the correct target."""

    count = len(units)
    components = np.asarray([unit.component for unit in units])
    different = components[:, None] != components[None, :]

    identity = sequence_similarity([unit.accession for unit in units])
    fingerprint = np.stack([bits[unit.support_pool].mean(axis=0) for unit in units])
    fingerprint = fingerprint / np.maximum(np.linalg.norm(fingerprint, axis=1, keepdims=True), 1e-9)
    chemistry = fingerprint @ fingerprint.T

    rng = np.random.default_rng(SEED)
    random_donor = np.zeros(count, dtype=np.int64)
    for position in range(count):
        allowed = np.flatnonzero(different[position])
        random_donor[position] = rng.choice(allowed) if len(allowed) else position

    def hardest(matrix: np.ndarray) -> np.ndarray:
        masked = np.where(different, matrix, -np.inf)
        return masked.argmax(axis=1)

    return {
        "random_target": random_donor,
        "protein_hard": hardest(identity),
        "chemistry_matched": hardest(chemistry),
    }


def draw_support(
    unit: Unit, units: list[Unit], donors: dict[str, np.ndarray], position: int,
    control: str, k: int, rng: np.random.Generator,
) -> tuple[np.ndarray, bool]:
    """Support rows and whether the residual labels should be permuted."""

    if control in ("correct", "label_permuted"):
        pool = unit.support_pool
    else:
        pool = units[int(donors[control][position])].support_pool
    if len(pool) < k:
        return np.array([], dtype=np.int64), False
    return rng.choice(pool, size=k, replace=False), control == "label_permuted"


# --------------------------------------------------------------------------- #
# Training
# --------------------------------------------------------------------------- #


def smoothed_concordance_loss(
    prediction: torch.Tensor, label: torch.Tensor, left: torch.Tensor, right: torch.Tensor
) -> torch.Tensor:
    """Bounded surrogate for 1 - CI.  A convex margin loss is the wrong choice on
    this substrate: the frozen base orders at chance with large spread, so its
    optimum sits far from the concordance optimum."""

    gap = prediction[left] - prediction[right]
    sign = torch.sign(label[left] - label[right])
    return torch.sigmoid(-gap * sign / TAU).mean()


def train_operator(
    name: str, basis: torch.Tensor, residual: torch.Tensor, affinity: torch.Tensor,
    base: torch.Tensor, units: list[Unit], roles: dict[str, str], seed: int,
) -> tuple[nn.Module, dict[str, object]]:
    torch.manual_seed(seed)
    model = MODELS[name](basis.shape[1]).to(DEVICE)
    optimiser = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    train_units = [u for u in units if roles[u.component] == "meta_train"]
    val_units = [u for u in units if roles[u.component] == "meta_val"]
    rng = np.random.default_rng(seed)

    history, best_score, best_state, stale = [], -np.inf, None, 0
    for epoch in range(MAX_EPOCHS):
        model.train()
        order = rng.permutation(len(train_units))
        total, steps = 0.0, 0
        for start in range(0, len(order), BATCH_UNITS):
            optimiser.zero_grad()
            losses = []
            for index in order[start : start + BATCH_UNITS]:
                unit = train_units[index]
                for _ in range(TRAIN_EPISODES_PER_EPOCH):
                    k = int(rng.choice(K_SWEEP))
                    if len(unit.support_pool) < k:
                        continue
                    support = rng.choice(unit.support_pool, size=k, replace=False)
                    delta = model(basis[unit.query], basis[support], residual[support])
                    losses.append(smoothed_concordance_loss(
                        base[unit.query] + delta, affinity[unit.query],
                        torch.as_tensor(unit.pair_left, device=DEVICE),
                        torch.as_tensor(unit.pair_right, device=DEVICE),
                    ))
            if not losses:
                continue
            loss = torch.stack(losses).mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimiser.step()
            total += float(loss.item())
            steps += 1

        score = validation_gain(model, basis, residual, affinity, base, val_units, seed + epoch)
        history.append({"epoch": epoch, "loss": total / max(steps, 1), "val_gain": score})
        if score > best_score + 1e-5:
            best_score, stale = score, 0
            best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
        else:
            stale += 1
            if stale >= PATIENCE:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"history": history, "best_val_gain": best_score, "epochs_run": len(history)}


@torch.no_grad()
def validation_gain(
    model: nn.Module, basis: torch.Tensor, residual: torch.Tensor, affinity: torch.Tensor,
    base: torch.Tensor, units: list[Unit], seed: int,
) -> float:
    model.eval()
    rng = np.random.default_rng(seed)
    per_component: dict[str, list[float]] = {}
    for unit in units:
        label = affinity[unit.query].cpu().numpy()
        reference = concordance(base[unit.query].cpu().numpy(), label, unit.pair_left, unit.pair_right)
        gains = []
        for k in K_SWEEP:
            if len(unit.support_pool) < k:
                continue
            for _ in range(4):
                support = rng.choice(unit.support_pool, size=k, replace=False)
                delta = model(basis[unit.query], basis[support], residual[support])
                prediction = (base[unit.query] + delta).cpu().numpy()
                gains.append(concordance(prediction, label, unit.pair_left, unit.pair_right) - reference)
        if gains:
            per_component.setdefault(unit.component, []).append(float(np.mean(gains)))
    if not per_component:
        return float("nan")
    return float(np.mean([np.mean(values) for values in per_component.values()]))


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def tune_baseline(
    name: str, basis_np: np.ndarray, bits: np.ndarray, residual_np: np.ndarray,
    affinity_np: np.ndarray, base_np: np.ndarray, units: list[Unit], seed: int,
) -> dict[int, tuple[float, float]]:
    """Grid-select ridge and a global transport scale per k on meta-val.

    Every baseline is granted its own scale: a gain measured against an unscaled
    kernel is a scale artefact, not a mechanism.
    """

    grid = [(ridge, scale) for ridge in (0.01, 0.1, 1.0, 10.0) for scale in (0.25, 0.5, 1.0)]
    chosen: dict[int, tuple[float, float]] = {}
    for k in K_SWEEP:
        best, best_score = (1.0, 1.0), -np.inf
        for ridge, scale in grid:
            rng = np.random.default_rng(seed)
            per_component: dict[str, list[float]] = {}
            for unit in units:
                if len(unit.support_pool) < k:
                    continue
                label = affinity_np[unit.query]
                reference = concordance(base_np[unit.query], label, unit.pair_left, unit.pair_right)
                gains = []
                for _ in range(4):
                    support = rng.choice(unit.support_pool, size=k, replace=False)
                    delta = baseline_delta(
                        name, basis_np[unit.query], basis_np[support], residual_np[support],
                        bits[unit.query], bits[support], ridge, scale,
                    )
                    gains.append(concordance(base_np[unit.query] + delta, label,
                                             unit.pair_left, unit.pair_right) - reference)
                if gains:
                    per_component.setdefault(unit.component, []).append(float(np.mean(gains)))
            if per_component:
                score = float(np.mean([np.mean(v) for v in per_component.values()]))
                if score > best_score:
                    best, best_score = (ridge, scale), score
        chosen[k] = best
    return chosen


@torch.no_grad()
def evaluate_all(
    models: dict[str, nn.Module], tuned: dict[str, dict[int, tuple[float, float]]],
    basis: torch.Tensor, basis_np: np.ndarray, bits: np.ndarray,
    residual_np: np.ndarray, affinity_np: np.ndarray,
    base_np: np.ndarray, units: list[Unit], all_units: list[Unit],
    donors: dict[str, np.ndarray], index_of: dict[str, int], seed: int,
) -> pd.DataFrame:
    for model in models.values():
        model.eval()
    records: list[dict[str, object]] = []
    for number, unit in enumerate(units):
        position = index_of[unit.unit]
        label = affinity_np[unit.query]
        reference = unit_metrics(base_np[unit.query], label, unit)
        for k in K_SWEEP:
            for control in CONTROLS:
                rng = np.random.default_rng(
                    int(sha256(f"{seed}:{unit.unit}:{k}:{control}".encode()).hexdigest()[:8], 16)
                )
                accumulated: dict[str, list[dict[str, float]]] = {}
                for _ in range(EVAL_EPISODES):
                    support, permute = draw_support(unit, all_units, donors, position, control, k, rng)
                    if len(support) == 0:
                        continue
                    values = residual_np[support].copy()
                    if permute:
                        values = values[rng.permutation(len(values))]
                    tensor_values = torch.as_tensor(values, dtype=basis.dtype, device=DEVICE)
                    for name in ("intercept", "krr", "ridge"):
                        ridge, scale = tuned[name][k]
                        delta = baseline_delta(
                            name, basis_np[unit.query], basis_np[support], values,
                            bits[unit.query], bits[support], ridge, scale,
                        )
                        accumulated.setdefault(name, []).append(
                            unit_metrics(base_np[unit.query] + delta, label, unit))
                    for name, model in models.items():
                        delta = model(basis[unit.query], basis[support], tensor_values)
                        prediction = base_np[unit.query] + delta.cpu().numpy().astype(np.float64)
                        accumulated.setdefault(name, []).append(unit_metrics(prediction, label, unit))
                for name, values in accumulated.items():
                    row: dict[str, object] = {
                        "unit": unit.unit, "component": unit.component, "endpoint": unit.endpoint,
                        "arm": name, "control": control, "k": k,
                        "n_support_pool": int(len(unit.support_pool)), "n_query": int(len(unit.query)),
                    }
                    for metric in ("ci_within", "rmse", "mae", "ndcg10", "pearson", "spearman"):
                        scores = [item[metric] for item in values if np.isfinite(item[metric])]
                        row[metric] = float(np.mean(scores)) if scores else float("nan")
                        row[f"base_{metric}"] = reference[metric]
                        row[f"gain_{metric}"] = (
                            row[metric] - reference[metric] if np.isfinite(row[metric]) else float("nan")
                        )
                    records.append(row)
        if number % 25 == 0:
            print(f"  eval {number}/{len(units)}", flush=True)
    return pd.DataFrame.from_records(records)


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for (arm, control, k), frame in records.groupby(["arm", "control", "k"]):
        cell = {metric: component_bootstrap(frame, f"gain_{metric}")
                for metric in ("ci_within", "rmse", "spearman", "ndcg10")}
        summary.setdefault(str(arm), {}).setdefault(str(control), {})[f"k{k}"] = cell
    return summary


def gate(summary: dict[str, object], records: pd.DataFrame) -> dict[str, object]:
    """The admission gate, evaluated arm by arm on within-document concordance."""

    def value(arm: str, control: str, k: int, field: str = "mean") -> float:
        try:
            return summary[arm][control][f"k{k}"]["ci_within"][field]
        except KeyError:
            return float("nan")

    verdicts: dict[str, object] = {}
    for arm in sorted(set(records.arm)):
        per_k = {}
        for k in K_SWEEP:
            correct = value(arm, "correct", k)
            lower = value(arm, "correct", k, "lower95")
            wrong = max(
                (value(arm, control, k) for control in CONTROLS if control != "correct"),
                default=float("nan"),
            )
            beats_base = bool(np.isfinite(lower) and lower > 0.005)
            beats_wrong = bool(np.isfinite(correct) and np.isfinite(wrong) and correct - wrong > 0.005)
            beats_calibration = bool(correct - value("intercept", "correct", k) > 0.005)
            beats_kernel = bool(
                correct - max(value("krr", "correct", k), value("ridge", "correct", k)) > 0.005
            )
            per_k[f"k{k}"] = {
                "gain": correct, "lower95": lower, "best_wrong_support": wrong,
                "beats_support_free": beats_base, "beats_wrong_support": beats_wrong,
                "beats_calibration": beats_calibration, "beats_kernel_and_ridge": beats_kernel,
                "admitted": bool(beats_base and beats_wrong and beats_calibration and beats_kernel),
            }
        verdicts[arm] = per_k
    admitted = sorted({
        arm for arm, cells in verdicts.items()
        if any(cell["admitted"] for cell in cells.values()) and arm in MODELS
    })
    return {
        "per_arm": verdicts,
        "admitted_mechanisms": admitted,
        "verdict": "MECHANISM_ADMITTED" if admitted else "NO_OPERATOR_PASSES_ADMISSION_GATE",
    }


def run(substrate_dir: Path, role: str, output: Path, records_path: Path, seed: int) -> dict[str, object]:
    started = time.time()
    from scipy.sparse import load_npz

    basis_np, substrate, stats = rich_basis(substrate_dir, role)
    basis_np = np.ascontiguousarray(basis_np[:, :CAPACITY])
    splits = build_splits(substrate.rows)
    units = build_units(substrate, splits)
    roles = {unit.component: assign_meta_role(unit.component) for unit in units}
    print(f"{len(units)} units, {len(set(u.component for u in units))} components; device {DEVICE}", flush=True)

    bits_all = load_npz(substrate_dir / "morgan.npz").tocsr()
    bits = np.asarray(bits_all[substrate.rows.structure_row.to_numpy()].todense(), dtype=np.float32)

    scale = np.maximum(basis_np.std(axis=0), 1e-6)
    basis_np = (basis_np - basis_np.mean(axis=0)) / scale
    basis = torch.as_tensor(basis_np, dtype=torch.float32, device=DEVICE)
    residual_np = np.asarray(substrate.residual, dtype=np.float64)
    residual = torch.as_tensor(residual_np, dtype=torch.float32, device=DEVICE)
    affinity_np = np.asarray(substrate.affinity, dtype=np.float64)
    base_np = np.asarray(substrate.base, dtype=np.float64)
    affinity = torch.as_tensor(affinity_np, dtype=torch.float32, device=DEVICE)
    base = torch.as_tensor(base_np, dtype=torch.float32, device=DEVICE)

    donors = donor_maps(units, bits)
    index_of = {unit.unit: position for position, unit in enumerate(units)}
    val_units = [u for u in units if roles[u.component] == "meta_val"]
    test_units = [u for u in units if roles[u.component] == "meta_test"]
    print(f"meta_train {sum(1 for u in units if roles[u.component]=='meta_train')} "
          f"meta_val {len(val_units)} meta_test {len(test_units)}", flush=True)

    tuned = {
        name: tune_baseline(name, basis_np, bits, residual_np, affinity_np, base_np, val_units, seed)
        for name in ("intercept", "krr", "ridge")
    }
    print("tuned:", tuned, flush=True)

    models, training = {}, {}
    for name in MODELS:
        print(f"training {name} ...", flush=True)
        model, info = train_operator(name, basis, residual, affinity, base, units, roles, seed)
        models[name], training[name] = model, info
        print(f"  {name}: best_val_gain={info['best_val_gain']:+.4f} epochs={info['epochs_run']}", flush=True)

    records = evaluate_all(models, tuned, basis, basis_np, bits, residual_np,
                           affinity_np, base_np, test_units, units, donors, index_of, seed)
    summary = summarise(records)
    decision = gate(summary, records)

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "device": str(DEVICE),
        "substrate": stats,
        "firewall": {
            "evaluated_role": role, "validate_read": False, "confirm_read": False,
            "meta_split": "components split 60/20/20 inside discover",
            "reported_on": "meta_test components only",
        },
        "protocol": {
            "seed": seed, "capacity": CAPACITY, "k_sweep": list(K_SWEEP),
            "eval_episodes": EVAL_EPISODES, "controls": list(CONTROLS),
            "loss": "bounded smoothed-concordance surrogate on within-document query pairs",
            "primary_metric": "within-document pair concordance gain over the support-free base",
            "baseline_tuning": "ridge and global transport scale grid-selected per k on meta_val",
        },
        "counts": {
            "units": len(units), "meta_test_units": len(test_units),
            "meta_test_components": len(set(u.component for u in test_units)),
        },
        "training": training,
        "baseline_hyperparameters": {n: {str(k): list(v) for k, v in t.items()} for n, t in tuned.items()},
        "summary": summary,
        **decision,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Support-conditioned adaptation operators")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--seed", type=int, default=SEED)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records, arguments.seed)
    print(json.dumps({"verdict": payload["verdict"],
                      "admitted": payload["admitted_mechanisms"]}, indent=2))
    header = "{:<11s}{:>3s}  {:>20s}  {:>11s}".format("arm", "k", "correct gain", "best wrong")
    print("\n" + header)
    for arm in sorted(payload["summary"]):
        for k in K_SWEEP:
            cell = payload["per_arm"][arm]["k{}".format(k)]
            flag = "ADMITTED" if cell["admitted"] else ""
            print("{:<11s}{:>3d}  {:+.4f} [{:+.4f}]  {:+.4f}  {}".format(
                arm, k, cell["gain"], cell["lower95"], cell["best_wrong_support"], flag))


if __name__ == "__main__":
    main()
