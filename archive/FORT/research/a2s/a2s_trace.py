"""A2S-TRACE Q2: amortised, label-free per-pair residual transport with abstention.

Q1 (`research/a2s_trace_stratum.py`) measured where correctly assigned support
labels carry transferable ranking information: only when the query is chemically
close to a support compound.  Inside that admitted stratum a zero-parameter
Tanimoto KRR already captures most of the available gain.

TRACE asks the next question.  A fixed isotropic similarity does not know *when*
short range lies: Tanimoto 0.7 is sometimes a bioisosteric swap that preserves
potency and sometimes an activity cliff.  TRACE meta-learns one label-free
function that says how much *more or less* than isotropic similarity a measured
residual on support compound ``i`` should transport to query compound ``q``, plus
a per-query gate that lets the transport abstain::

    delta_q = clip( c * alpha_q * sum_i m_qi * k_qi * r~_i , +/- max_i |r_i| )

Structural properties, provable from the functional form and asserted in
``tests/test_a2s_trace.py``:

* ``m`` and ``alpha`` never read a label or a residual, so a residual derangement
  leaves every weight bit-identical and the correct-minus-deranged contrast
  isolates compound-to-evidence assignment exactly;
* ``r_S == 0`` implies ``delta == 0`` bitwise;
* ``|delta_q| <= max_i |r_i|``, an observed quantity, so no unbounded learned
  scale can multiply a residual aggregate;
* switching the learned parts off recovers fixed Tanimoto KRR exactly, and
  switching the kernel row to its normalised form recovers Nadaraya-Watson
  exactly.

Zero target-specific quantities are estimated at meta-test.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn

from research.a2s.a2s_information_gate import (
    DEVICE,
    MORGAN_BITS,
    FEATURES,
    PROTEINS,
    REGISTRY,
    canonical,
    load_design,
    load_labeled_fit_probe,
    load_metadata,
    row_index,
    sha256_file,
    verify_lock,
)
from research.a2s.a2s_trace_stratum import (
    ADMISSION_MDE,
    DEFAULT_LOCK,
    DEFAULT_OOF,
    EPISODE_SEEDS,
    MIN_STRATUM_QUERY,
    STRATUM_NAMES,
    DrawnEpisode,
    build_episodes,
    derangement,
    metric_loss,
    paired_bootstrap,
    resolve_base,
    stratum_of,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_trace_q2_mechanism_2026-08-01.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_trace_q2_records_2026-08-01.parquet"

SEED = 20260801
TRAIN_POLICIES = ("random_within_target", "scaffold_disjoint")
ADMITTED_STRATA = ("t55_100",)
ADMITTED_K = (3, 5)
VALIDATION_FOLD = 4
EVALUATION_EPISODE_SEED = EPISODE_SEEDS[0]
EMBED_DIM = 64
PROTEIN_DIM = 32
HIDDEN = 128
PAIR_MARGIN = 0.10
MAX_EPOCHS = 24
PATIENCE = 4
BATCH_EPISODES = 64
LEARNING_RATE = 3e-3          # learned heads (validated on the synthetic control)
SCALE_LEARNING_RATE = 1e-2    # the two analytic scalars (transport scale, ridge)
WEIGHT_DECAY = 1e-5
RANK_TAU = 1.0
DEVIATION_WEIGHT = 0.0
EVALUATIONS_PER_EPOCH = 3
WARMUP_EPOCHS = 3
KRR_RIDGE_GRID = (0.03, 0.1, 0.3, 1.0)
MIXTURE_GAMMAS = (0.5, 1.0, 2.0, 4.0)
GATE_BIAS_INIT = 6.0


# --------------------------------------------------------------------------- #
# Substrate
# --------------------------------------------------------------------------- #


@dataclass
class Substrate:
    labeled: pd.DataFrame
    bits: torch.Tensor          # (n, 1024) binary Morgan
    desc: torch.Tensor          # (n, 10) standardised descriptors
    base: torch.Tensor          # (n,) frozen out-of-fold base prediction
    affinity: torch.Tensor      # (n,) pKi
    protein: torch.Tensor       # (n_targets, 1280) pooled ESM-2
    target_row: torch.Tensor    # (n,) index into `protein`
    position: dict[int, int]


def load_substrate(lock_path: Path, oof_cache: Path) -> tuple[Substrate, dict[str, object]]:
    lock = verify_lock(lock_path)
    metadata = load_metadata(lock)
    labeled = load_labeled_fit_probe(metadata)
    values, raw_features, target_to_protein = load_design(metadata, labeled)
    base, oof_stats = resolve_base(labeled, values, metadata, oof_cache)

    source_rows = labeled.source_row.to_numpy(dtype=np.int64)
    bits = torch.as_tensor(
        raw_features[source_rows, :MORGAN_BITS].astype(np.float32), device=DEVICE
    )
    desc = torch.as_tensor(values[:, MORGAN_BITS : MORGAN_BITS + 10].astype(np.float32), device=DEVICE)
    protein_archive = np.load(PROTEINS, allow_pickle=False)
    pooled = np.asarray(protein_archive["pooled"], dtype=np.float32)
    protein = torch.as_tensor(pooled, device=DEVICE)
    target_row = torch.as_tensor(
        np.asarray([target_to_protein[str(target)] for target in labeled.target], dtype=np.int64),
        device=DEVICE,
    )
    substrate = Substrate(
        labeled=labeled,
        bits=bits,
        desc=desc,
        base=torch.as_tensor(base.astype(np.float32), device=DEVICE),
        affinity=torch.as_tensor(labeled.affinity.to_numpy(np.float32), device=DEVICE),
        protein=protein,
        target_row=target_row,
        position=row_index(labeled),
    )
    return substrate, {"lock": lock, "oof": oof_stats}


# --------------------------------------------------------------------------- #
# Episode tensors
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Batch:
    support: torch.Tensor       # (B, k) row positions
    query: torch.Tensor         # (B, Q) row positions
    query_mask: torch.Tensor    # (B, Q) bool
    k: int


def index_episodes(episodes: list[DrawnEpisode], substrate: Substrate) -> tuple[np.ndarray, np.ndarray]:
    position = substrate.position
    support = np.asarray(
        [[position[int(row)] for row in episode.support_rows] for episode in episodes], dtype=np.int64
    )
    width = max(len(episode.query_rows) for episode in episodes)
    query = np.full((len(episodes), width), -1, dtype=np.int64)
    for offset, episode in enumerate(episodes):
        indices = [position[int(row)] for row in episode.query_rows]
        query[offset, : len(indices)] = indices
    return support, query


def make_batch(support: np.ndarray, query: np.ndarray, rows: np.ndarray) -> Batch:
    support_t = torch.as_tensor(support[rows], device=DEVICE)
    query_block = query[rows]
    mask = query_block >= 0
    width = int(mask.sum(axis=1).max())
    query_block = query_block[:, :width]
    mask = mask[:, :width]
    return Batch(
        support=support_t,
        query=torch.as_tensor(np.where(mask, query_block, 0), device=DEVICE),
        query_mask=torch.as_tensor(mask, device=DEVICE),
        k=int(support.shape[1]),
    )


def tanimoto(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    """Binary Tanimoto between (B, A, D) and (B, C, D) fingerprint blocks."""

    intersection = left @ right.transpose(1, 2)
    left_count = left.sum(-1).unsqueeze(2)
    right_count = right.sum(-1).unsqueeze(1)
    union = left_count + right_count - intersection
    return intersection / union.clamp(min=1.0)


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TraceConfig:
    """A declared point in the transport family.

    ``modulation=False, gate=False, weights="krr", whiten=True`` is *exactly*
    bound-clipped fixed Tanimoto KRR; ``weights="nw", whiten=False`` with both
    learned parts off is *exactly* the Nadaraya-Watson smoother.  Everything the
    mechanism claims is therefore the measured delta from a named restriction.
    """

    weights: str = "krr"               # "krr" (unnormalised kernel row) | "nw"
    whiten: bool = True                # label-free support-Gram conditioning
    global_scale: bool = True          # one meta-learned scalar, target-independent
    modulation: bool = True            # learned per-pair transport reliability
    gate: bool = True                  # learned per-query transport gate
    protein: str = "learned"           # "learned" | "zero"
    pair_features: str = "full"        # "full" (adds a learned bilinear) | "scalar"
    ridge: float = 0.1
    embed_dim: int = EMBED_DIM
    protein_dim: int = PROTEIN_DIM
    hidden: int = HIDDEN


class Trace(nn.Module):
    """Learned reliability modulation of a fixed label-free transport kernel.

    ``delta_q = clamp( alpha_q * sum_i m_qi * w_qi * r~_i , +/- max_i|r_i| )``

    * ``w_qi`` and ``r~`` are fixed label-free analytic objects (Tanimoto kernel
      row and support-Gram-whitened residual);
    * ``m_qi = 2*sigmoid(g_psi(pair, protein))`` is the learned per-pair
      reliability, zero-initialised so ``m == 1`` before training;
    * ``alpha_q = sigmoid(h_psi(query, protein, k, similarity stats) + b)`` is the
      learned per-query gate, initialised open.

    Neither ``g_psi`` nor ``h_psi`` reads a label or a residual.
    """

    def __init__(self, config: TraceConfig) -> None:
        super().__init__()
        self.config = config
        self.fingerprint = nn.Linear(MORGAN_BITS, config.embed_dim, bias=False)
        self.protein_projection = nn.Linear(1280, config.protein_dim)
        pair_dim = 6 + 20 + config.protein_dim
        if config.pair_features == "full":
            pair_dim += 2 * config.embed_dim
        self.pair = nn.Sequential(
            nn.Linear(pair_dim, config.hidden),
            nn.SiLU(),
            nn.Linear(config.hidden, config.hidden // 2),
            nn.SiLU(),
            nn.Linear(config.hidden // 2, 1),
        )
        gate_dim = config.embed_dim + config.protein_dim + 3 + 3
        self.query_gate = nn.Sequential(
            nn.Linear(gate_dim, config.hidden // 2),
            nn.SiLU(),
            nn.Linear(config.hidden // 2, 1),
        )
        self.log_ridge = nn.Parameter(torch.tensor(math.log(math.expm1(config.ridge))))
        # One global, target-independent transport scale.  The frozen base carries
        # roughly twice the within-episode spread of the transport while ordering
        # at chance, so the relative precision of base and transport is itself a
        # free scalar.  It is an episode constant, therefore it belongs to the
        # baseline and never to the mechanism's claim.
        self.log_scale = nn.Parameter(torch.tensor(math.log(math.expm1(1.0))))
        nn.init.zeros_(self.pair[-1].weight)
        nn.init.zeros_(self.pair[-1].bias)
        nn.init.zeros_(self.query_gate[-1].weight)
        # A large positive bias starts the gate fully open, so an untrained gated
        # model is numerically the ungated model.  Closing has to be learned.
        nn.init.constant_(self.query_gate[-1].bias, GATE_BIAS_INIT)

    # -- pieces ---------------------------------------------------------- #

    def protein_context(self, protein: torch.Tensor) -> torch.Tensor:
        if self.config.protein == "zero":
            return torch.zeros(
                protein.shape[0], self.config.protein_dim, device=protein.device, dtype=protein.dtype
            )
        return torch.tanh(self.protein_projection(protein))

    def embed(self, bits: torch.Tensor) -> torch.Tensor:
        return torch.tanh(self.fingerprint(bits))

    def whitened(self, gram: torch.Tensor, residual: torch.Tensor) -> torch.Tensor:
        if not self.config.whiten:
            return residual
        ridge = nn.functional.softplus(self.log_ridge)
        eye = torch.eye(gram.shape[-1], device=gram.device, dtype=gram.dtype)
        return torch.linalg.solve(gram + ridge * eye, residual.unsqueeze(-1)).squeeze(-1)

    def pair_features(
        self,
        query_bits: torch.Tensor,
        support_bits: torch.Tensor,
        query_desc: torch.Tensor,
        support_desc: torch.Tensor,
        protein: torch.Tensor,
        cross: torch.Tensor,
    ) -> torch.Tensor:
        batch, n_query, _ = query_bits.shape
        n_support = support_bits.shape[1]
        intersection = query_bits @ support_bits.transpose(1, 2)
        query_count = query_bits.sum(-1).unsqueeze(2)
        support_count = support_bits.sum(-1).unsqueeze(1)
        dice = 2.0 * intersection / (query_count + support_count).clamp(min=1.0)
        only_query = (query_count - intersection) / query_count.clamp(min=1.0)
        only_support = (support_count - intersection) / support_count.clamp(min=1.0)
        shared = intersection / torch.minimum(query_count, support_count).clamp(min=1.0)
        log_cross = torch.log(cross.clamp(min=1e-6))
        scalar = torch.stack((cross, log_cross, dice, only_query, only_support, shared), dim=-1)
        delta_desc = query_desc.unsqueeze(2) - support_desc.unsqueeze(1)
        context = self.protein_context(protein)
        context = context.unsqueeze(1).unsqueeze(1).expand(batch, n_query, n_support, -1)
        parts = [scalar, delta_desc, delta_desc.abs(), context]
        if self.config.pair_features == "full":
            embed_query = self.embed(query_bits)
            embed_support = self.embed(support_bits)
            parts.insert(3, (embed_query.unsqueeze(2) - embed_support.unsqueeze(1)).abs())
            parts.insert(3, embed_query.unsqueeze(2) * embed_support.unsqueeze(1))
        return torch.cat(parts, dim=-1)

    def gate_value(
        self, query_bits: torch.Tensor, protein: torch.Tensor, cross: torch.Tensor
    ) -> torch.Tensor:
        batch, n_query, _ = query_bits.shape
        n_support = cross.shape[-1]
        context = self.protein_context(protein)
        statistics = torch.stack(
            (cross.mean(-1), cross.amax(-1), cross.std(-1, unbiased=False)), dim=-1
        )
        budget = torch.zeros(batch, 3, device=cross.device, dtype=cross.dtype)
        budget[:, {1: 0, 3: 1, 5: 2}.get(n_support, 2)] = 1.0
        features = torch.cat(
            (
                self.embed(query_bits),
                context.unsqueeze(1).expand(batch, n_query, -1),
                budget.unsqueeze(1).expand(batch, n_query, -1),
                statistics,
            ),
            dim=-1,
        )
        return torch.sigmoid(self.query_gate(features).squeeze(-1))

    # -- forward ---------------------------------------------------------- #

    def forward(
        self,
        query_bits: torch.Tensor,
        support_bits: torch.Tensor,
        query_desc: torch.Tensor,
        support_desc: torch.Tensor,
        protein: torch.Tensor,
        residual: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        cross = tanimoto(query_bits, support_bits)
        gram = tanimoto(support_bits, support_bits)
        if self.config.weights == "krr":
            transport = cross
        elif self.config.weights == "nw":
            transport = cross / cross.sum(-1, keepdim=True).clamp(min=1e-9)
        else:
            raise ValueError(f"unknown transport weights {self.config.weights!r}")
        if self.config.modulation:
            features = self.pair_features(
                query_bits, support_bits, query_desc, support_desc, protein, cross
            )
            modulation = 2.0 * torch.sigmoid(self.pair(features).squeeze(-1))
        else:
            modulation = torch.ones_like(cross)
        gate = (
            self.gate_value(query_bits, protein, cross)
            if self.config.gate
            else torch.ones(cross.shape[0], cross.shape[1], device=cross.device, dtype=cross.dtype)
        )
        scale = (
            nn.functional.softplus(self.log_scale)
            if self.config.global_scale
            else torch.ones((), device=cross.device, dtype=cross.dtype)
        )
        transported = scale * gate * (
            (modulation * transport) * self.whitened(gram, residual).unsqueeze(1)
        ).sum(-1)
        bound = residual.abs().amax(-1, keepdim=True)
        delta = torch.clamp(transported, min=-bound, max=bound)
        # `bound == 0` already forces `delta == 0` through the clamp; the explicit
        # mask keeps the residual-null guarantee independent of clamp edge cases.
        delta = torch.where(bound > 0, delta, torch.zeros_like(delta))
        diagnostics = {
            "modulation": modulation,
            "gate": gate,
            "scale": scale,
            "clipped": (transported.abs() > bound).to(delta.dtype),
        }
        return delta, diagnostics


class StaticMixture(nn.Module):
    """Baseline: globally learned convex mixture of fixed Tanimoto KRR experts."""

    def __init__(self, gammas: tuple[float, ...] = MIXTURE_GAMMAS) -> None:
        super().__init__()
        self.gammas = gammas
        self.logits = nn.Parameter(torch.zeros(len(gammas)))
        self.log_ridge = nn.Parameter(torch.full((len(gammas),), math.log(math.expm1(0.1))))
        # Granted the same global transport scale as the analytic bar, so the
        # comparison is not decided by a magnitude the baseline cannot express.
        self.log_scale = nn.Parameter(torch.tensor(math.log(math.expm1(1.0))))

    def forward(
        self,
        query_bits: torch.Tensor,
        support_bits: torch.Tensor,
        query_desc: torch.Tensor,
        support_desc: torch.Tensor,
        protein: torch.Tensor,
        residual: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        cross = tanimoto(query_bits, support_bits)
        gram = tanimoto(support_bits, support_bits)
        ridge = nn.functional.softplus(self.log_ridge)
        weights = torch.softmax(self.logits, dim=0)
        eye = torch.eye(gram.shape[-1], device=gram.device, dtype=gram.dtype)
        total = torch.zeros(cross.shape[0], cross.shape[1], device=cross.device, dtype=cross.dtype)
        for index, gamma in enumerate(self.gammas):
            gram_m = gram.clamp(min=0.0) ** gamma
            cross_m = cross.clamp(min=0.0) ** gamma
            solved = torch.linalg.solve(gram_m + ridge[index] * eye, residual.unsqueeze(-1)).squeeze(-1)
            total = total + weights[index] * (cross_m @ solved.unsqueeze(-1)).squeeze(-1)
        total = nn.functional.softplus(self.log_scale) * total
        bound = residual.abs().amax(-1, keepdim=True)
        delta = torch.clamp(total, min=-bound, max=bound)
        delta = torch.where(bound > 0, delta, torch.zeros_like(delta))
        return delta, {"mixture": weights}


# --------------------------------------------------------------------------- #
# Analytic operators
# --------------------------------------------------------------------------- #


def analytic_delta(
    query_bits: torch.Tensor,
    support_bits: torch.Tensor,
    residual: torch.Tensor,
    *,
    estimator: str,
    ridge: float = 0.1,
) -> torch.Tensor:
    cross = tanimoto(query_bits, support_bits)
    if estimator == "level":
        return residual.mean(-1, keepdim=True).expand(-1, cross.shape[1])
    if estimator == "nw":
        weights = cross / cross.sum(-1, keepdim=True).clamp(min=1e-9)
        return (weights * residual.unsqueeze(1)).sum(-1)
    if estimator == "krr":
        gram = tanimoto(support_bits, support_bits)
        eye = torch.eye(gram.shape[-1], device=gram.device, dtype=gram.dtype)
        solved = torch.linalg.solve(gram + ridge * eye, residual.unsqueeze(-1)).squeeze(-1)
        return (cross * solved.unsqueeze(1)).sum(-1)
    if estimator == "cka_nnls":
        return cka_nnls_delta(query_bits, support_bits, residual, ridge=ridge)
    raise ValueError(f"unknown estimator {estimator!r}")


def cka_nnls_delta(
    query_bits: torch.Tensor,
    support_bits: torch.Tensor,
    residual: torch.Tensor,
    *,
    ridge: float = 0.1,
    gammas: tuple[float, ...] = MIXTURE_GAMMAS,
) -> torch.Tensor:
    """Centred-kernel-alignment NNLS mixture, solved independently per episode."""

    from scipy.optimize import nnls

    cross = tanimoto(query_bits, support_bits)
    gram = tanimoto(support_bits, support_bits)
    batch, n_support = residual.shape
    centre = torch.eye(n_support, device=gram.device) - 1.0 / n_support
    eye = torch.eye(n_support, device=gram.device)
    experts = []
    aligned = []
    for gamma in gammas:
        gram_m = gram.clamp(min=0.0) ** gamma
        cross_m = cross.clamp(min=0.0) ** gamma
        solved = torch.linalg.solve(gram_m + ridge * eye, residual.unsqueeze(-1)).squeeze(-1)
        experts.append((cross_m @ solved.unsqueeze(-1)).squeeze(-1))
        centred = centre @ gram_m @ centre
        aligned.append(centred.reshape(batch, -1))
    target = (centre @ (residual.unsqueeze(-1) @ residual.unsqueeze(1)) @ centre).reshape(batch, -1)
    design = torch.stack(aligned, dim=-1).cpu().numpy().astype(np.float64)
    response = target.cpu().numpy().astype(np.float64)
    stacked = torch.stack(experts, dim=-1)
    weights = np.zeros((batch, len(gammas)), dtype=np.float32)
    for index in range(batch):
        solution, _ = nnls(design[index], response[index])
        total = float(solution.sum())
        weights[index] = solution / total if total > 1e-12 else np.full(len(gammas), 1.0 / len(gammas))
    weight_tensor = torch.as_tensor(weights, device=stacked.device).unsqueeze(1)
    return (stacked * weight_tensor).sum(-1)


# --------------------------------------------------------------------------- #
# Residual arms
# --------------------------------------------------------------------------- #


def residual_arm(
    residual: torch.Tensor, arm: str, rng: np.random.Generator, donor: torch.Tensor | None = None
) -> torch.Tensor:
    if arm == "correct":
        return residual
    if arm == "null":
        return torch.zeros_like(residual)
    if arm == "signflip":
        return -residual
    if arm == "deranged":
        k = residual.shape[1]
        if k < 2:
            return residual
        order = np.stack([derangement(k, rng) for _ in range(residual.shape[0])]).astype(np.int64)
        index = torch.as_tensor(order, device=residual.device, dtype=torch.long)
        return torch.gather(residual, 1, index)
    if arm == "wrong_target":
        if donor is None:
            raise ValueError("the wrong-target arm needs donor residuals")
        scale = residual.norm(dim=1, keepdim=True) / donor.norm(dim=1, keepdim=True).clamp(min=1e-9)
        return donor * scale
    raise ValueError(f"unknown arm {arm!r}")


# --------------------------------------------------------------------------- #
# Synthetic positive control
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SyntheticSpec:
    """A world in which transport reliability really is pair-dependent.

    Support residuals are drawn i.i.d.  A query's synthetic label is built by
    transporting them with Tanimoto weights *restricted to reliable pairs*, where
    a pair is reliable exactly when the two molecules are close on one held-out
    descriptor.  An isotropic kernel cannot represent that restriction; a learned
    per-pair modulation can.  If the pipeline cannot recover this, no null may be
    reported from it.
    """

    amplitude: float = 1.5
    noise: float = 0.30
    descriptor: int = 1          # standardised molecular weight
    threshold: float = 0.35
    enabled: bool = True


def apply_synthetic(
    data: dict[str, torch.Tensor], spec: SyntheticSpec, seed: int
) -> dict[str, torch.Tensor]:
    generator = torch.Generator(device=DEVICE).manual_seed(int(seed) % (2**31))
    residual = torch.randn(
        data["residual"].shape, generator=generator, device=DEVICE, dtype=data["residual"].dtype
    )
    cross = tanimoto(data["query_bits"], data["support_bits"])
    delta_desc = (
        data["query_desc"][:, :, spec.descriptor].unsqueeze(2)
        - data["support_desc"][:, :, spec.descriptor].unsqueeze(1)
    ).abs()
    reliable = (delta_desc <= spec.threshold).to(cross.dtype)
    weight = cross * reliable
    weight = weight / weight.sum(-1, keepdim=True).clamp(min=1e-6)
    transported = (weight * residual.unsqueeze(1)).sum(-1)
    noise = spec.noise * torch.randn(
        transported.shape, generator=generator, device=DEVICE, dtype=transported.dtype
    )
    return {
        **data,
        "residual": residual,
        "label": data["base"] + spec.amplitude * transported + noise,
    }


def synthetic_oracle_delta(
    data: dict[str, torch.Tensor], spec: SyntheticSpec
) -> torch.Tensor:
    """The transport an oracle with the true reliability mask would apply."""

    cross = tanimoto(data["query_bits"], data["support_bits"])
    delta_desc = (
        data["query_desc"][:, :, spec.descriptor].unsqueeze(2)
        - data["support_desc"][:, :, spec.descriptor].unsqueeze(1)
    ).abs()
    weight = cross * (delta_desc <= spec.threshold).to(cross.dtype)
    weight = weight / weight.sum(-1, keepdim=True).clamp(min=1e-6)
    return spec.amplitude * (weight * data["residual"].unsqueeze(1)).sum(-1)


# --------------------------------------------------------------------------- #
# Training
# --------------------------------------------------------------------------- #


def gather_batch(substrate: Substrate, batch: Batch) -> dict[str, torch.Tensor]:
    support_bits = substrate.bits[batch.support]
    query_bits = substrate.bits[batch.query]
    support_desc = substrate.desc[batch.support]
    query_desc = substrate.desc[batch.query]
    protein = substrate.protein[substrate.target_row[batch.support[:, 0]]]
    residual = substrate.affinity[batch.support] - substrate.base[batch.support]
    return {
        "support_bits": support_bits,
        "query_bits": query_bits,
        "support_desc": support_desc,
        "query_desc": query_desc,
        "protein": protein,
        "residual": residual,
        "base": substrate.base[batch.query],
        "label": substrate.affinity[batch.query],
        "mask": batch.query_mask,
    }


def pairwise_loss(
    prediction: torch.Tensor,
    label: torch.Tensor,
    mask: torch.Tensor,
    *,
    surrogate: str = "smooth_ci",
    tau: float = RANK_TAU,
) -> torch.Tensor:
    """Within-episode pairwise ranking loss.

    ``smooth_ci`` is a bounded surrogate for the pairwise error rate, which is
    what the declared endpoint (CI) measures.  The convex ``logistic`` surrogate
    is retained for reference: because it keeps punishing confidently-wrong
    pairs without bound, its optimum sits at a much smaller transport magnitude
    than the CI optimum, and a model trained on it will not find the transport
    scale that CI rewards.  That mismatch is a measured property of this
    substrate, not a hyperparameter preference.
    """

    difference = label.unsqueeze(2) - label.unsqueeze(1)
    valid = (difference > PAIR_MARGIN) & mask.unsqueeze(2) & mask.unsqueeze(1)
    if not bool(valid.any()):
        return prediction.sum() * 0.0
    gap = prediction.unsqueeze(2) - prediction.unsqueeze(1)
    if surrogate == "smooth_ci":
        losses = torch.sigmoid(-gap / tau)
    elif surrogate == "logistic":
        losses = nn.functional.softplus(-gap)
    else:
        raise ValueError(f"unknown ranking surrogate {surrogate!r}")
    per_episode = (losses * valid).sum(dim=(1, 2)) / valid.sum(dim=(1, 2)).clamp(min=1)
    active = valid.sum(dim=(1, 2)) > 0
    return per_episode[active].mean()


def deviation_penalty(diagnostics: dict[str, torch.Tensor], mask: torch.Tensor) -> torch.Tensor:
    """Shrink the learned parts toward their analytic restriction.

    The modulation defaults to 1 and the gate to fully open, which is exactly the
    nested KRR restriction.  Penalising deviation makes the learned increment pay
    for itself rather than drift, and it keeps an untrained-looking model from
    being an unregularised free operator.
    """

    weight = mask.to(diagnostics["gate"].dtype)
    gate = (((1.0 - diagnostics["gate"]) ** 2) * weight).sum() / weight.sum().clamp(min=1)
    modulation = diagnostics["modulation"]
    pair_weight = weight.unsqueeze(-1).expand_as(modulation)
    pair = (((modulation - 1.0) ** 2) * pair_weight).sum() / pair_weight.sum().clamp(min=1)
    return gate + pair


def episode_ci(prediction: np.ndarray, label: np.ndarray) -> float:
    left, right = np.triu_indices(len(label), k=1)
    truth = np.sign(label[left] - label[right])
    active = truth != 0
    if not active.any():
        return float("nan")
    predicted = np.sign(prediction[left] - prediction[right])
    return float((predicted[active] == truth[active]).mean() + 0.5 * (predicted[active] == 0).mean())


def train_model(
    model: nn.Module,
    substrate: Substrate,
    train_episodes: list[DrawnEpisode],
    val_episodes: list[DrawnEpisode],
    *,
    seed: int,
    max_epochs: int = MAX_EPOCHS,
    synthetic: SyntheticSpec | None = None,
    head_lr: float = LEARNING_RATE,
    deviation_weight: float = DEVIATION_WEIGHT,
    warmup_epochs: int = WARMUP_EPOCHS,
) -> dict[str, object]:
    torch.manual_seed(seed)
    model.to(DEVICE)
    analytic = [
        parameter
        for name, parameter in model.named_parameters()
        if name in {"log_scale", "log_ridge", "logits"}
    ]
    analytic_ids = {id(parameter) for parameter in analytic}
    heads = [parameter for parameter in model.parameters() if id(parameter) not in analytic_ids]
    groups = [{"params": analytic, "lr": SCALE_LEARNING_RATE, "weight_decay": 0.0}]
    if heads:
        groups.append({"params": heads, "lr": head_lr, "weight_decay": WEIGHT_DECAY})
    optimiser = torch.optim.Adam(groups)
    by_k_train = group_by_k(train_episodes, substrate)
    by_k_val = group_by_k(val_episodes, substrate)
    rng = np.random.default_rng(seed)
    # The untrained model is an exact analytic restriction, so it is a legitimate
    # candidate for early stopping: training must beat *not training*.
    initial = validation_score(model, substrate, by_k_val, synthetic=synthetic)
    best = {
        "score": initial,
        "epoch": -1,
        "state": {key: value.detach().clone() for key, value in model.state_dict().items()},
    }
    history: list[dict[str, float]] = [{"epoch": -1, "train_loss": float("nan"), "val_ci": initial}]
    stale = 0
    for epoch in range(max_epochs):
        # Warm-up freezes the learned heads so the analytic scalars converge
        # first.  The ladder is then nested in optimisation as well as in
        # parameterisation: the learned rungs start from the best restriction
        # rather than having to rediscover it.
        if len(groups) > 1:
            groups[1]["lr"] = 0.0 if epoch < warmup_epochs else head_lr
        model.train()
        order: list[tuple[int, np.ndarray]] = []
        for k, (episodes, support, query) in by_k_train.items():
            permutation = rng.permutation(len(episodes))
            for start in range(0, len(permutation), BATCH_EPISODES):
                order.append((k, permutation[start : start + BATCH_EPISODES]))
        rng.shuffle(order)
        checkpoints = {
            int(len(order) * (index + 1) / EVALUATIONS_PER_EPOCH) - 1
            for index in range(EVALUATIONS_PER_EPOCH)
        }
        total = 0.0
        for step, (k, rows) in enumerate(order):
            _, support, query = by_k_train[k]
            batch = make_batch(support, query, rows)
            data = gather_batch(substrate, batch)
            if synthetic is not None:
                data = apply_synthetic(data, synthetic, seed * 1_000_003 + epoch * 9973 + step)
            delta, diagnostics = model(
                data["query_bits"],
                data["support_bits"],
                data["query_desc"],
                data["support_desc"],
                data["protein"],
                data["residual"],
            )
            prediction = data["base"] + delta
            loss = pairwise_loss(prediction, data["label"], data["mask"])
            if "modulation" in diagnostics and deviation_weight > 0.0:
                loss = loss + deviation_weight * deviation_penalty(diagnostics, data["mask"])
            optimiser.zero_grad(set_to_none=True)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimiser.step()
            total += float(loss.detach())
            if step in checkpoints:
                score = validation_score(model, substrate, by_k_val, synthetic=synthetic)
                history.append(
                    {"epoch": epoch + (step + 1) / len(order), "train_loss": total / (step + 1), "val_ci": score}
                )
                if score > float(best["score"]) + 1e-6:
                    best = {
                        "score": score,
                        "epoch": epoch,
                        "state": {key: value.detach().clone() for key, value in model.state_dict().items()},
                    }
                    stale = 0
                elif epoch >= warmup_epochs:
                    stale += 1
                model.train()
        if stale >= PATIENCE * EVALUATIONS_PER_EPOCH:
            break
    model.load_state_dict(best["state"])
    return {"history": history, "best_epoch": int(best["epoch"]), "best_val_ci": float(best["score"])}


def group_by_k(
    episodes: list[DrawnEpisode], substrate: Substrate
) -> dict[int, tuple[list[DrawnEpisode], np.ndarray, np.ndarray]]:
    grouped: dict[int, list[DrawnEpisode]] = {}
    for episode in episodes:
        grouped.setdefault(episode.k, []).append(episode)
    output: dict[int, tuple[list[DrawnEpisode], np.ndarray, np.ndarray]] = {}
    for k, block in grouped.items():
        support, query = index_episodes(block, substrate)
        output[k] = (block, support, query)
    return output


@torch.no_grad()
def validation_score(
    model: nn.Module,
    substrate: Substrate,
    by_k: dict[int, tuple[list[DrawnEpisode], np.ndarray, np.ndarray]],
    *,
    synthetic: SyntheticSpec | None = None,
) -> float:
    """Mean episode CI over admitted-stratum queries of the inner validation fold."""

    model.eval()
    scores: list[float] = []
    for k, (episodes, support, query) in by_k.items():
        if k not in ADMITTED_K:
            continue
        for start in range(0, len(episodes), BATCH_EPISODES):
            rows = np.arange(start, min(start + BATCH_EPISODES, len(episodes)))
            batch = make_batch(support, query, rows)
            data = gather_batch(substrate, batch)
            if synthetic is not None:
                data = apply_synthetic(data, synthetic, 5_000_011 + k * 101 + start)
            delta, _ = model(
                data["query_bits"],
                data["support_bits"],
                data["query_desc"],
                data["support_desc"],
                data["protein"],
                data["residual"],
            )
            prediction = (data["base"] + delta).cpu().numpy()
            label = data["label"].cpu().numpy()
            mask = data["mask"].cpu().numpy()
            nearest = tanimoto(data["query_bits"], data["support_bits"]).amax(-1).cpu().numpy()
            for row in range(len(rows)):
                active = mask[row] & (nearest[row] >= 0.55)
                if int(active.sum()) < MIN_STRATUM_QUERY:
                    continue
                value = episode_ci(prediction[row][active], label[row][active])
                if np.isfinite(value):
                    scores.append(value)
    return float(np.mean(scores)) if scores else float("nan")


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


@torch.no_grad()
def evaluate(
    substrate: Substrate,
    episodes: list[DrawnEpisode],
    methods: dict[str, object],
    *,
    arms: dict[str, tuple[str, ...]] | None = None,
    krr_ridge: float = 0.1,
    seed: int = SEED,
    synthetic: SyntheticSpec | None = None,
) -> pd.DataFrame:
    """Evaluate every method on `episodes`.

    ``arms`` maps a method name to the residual arms it is run under.  Only the
    headline model and the analytic bar need the control arms, so the default
    keeps the metric budget proportional to the contrasts actually reported.
    """

    arms = arms or {}
    by_k = group_by_k(episodes, substrate)
    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    for k, (block, support, query) in sorted(by_k.items()):
        for start in range(0, len(block), BATCH_EPISODES):
            rows = np.arange(start, min(start + BATCH_EPISODES, len(block)))
            batch = make_batch(support, query, rows)
            data = gather_batch(substrate, batch)
            if synthetic is not None:
                data = apply_synthetic(data, synthetic, 7_000_003 + k * 101 + start)
            nearest = tanimoto(data["query_bits"], data["support_bits"]).amax(-1).cpu().numpy()
            donor_rows = rng.permutation(len(block))[: len(rows)]
            donor_batch = make_batch(support, query, donor_rows)
            donor = substrate.affinity[donor_batch.support] - substrate.base[donor_batch.support]
            predictions: dict[str, np.ndarray] = {
                "base__correct": data["base"].cpu().numpy()
            }
            if synthetic is not None:
                predictions["oracle__correct"] = (
                    data["base"] + synthetic_oracle_delta(data, synthetic)
                ).cpu().numpy()
            arm_cache: dict[str, torch.Tensor] = {}
            for name, method in methods.items():
                for arm in arms.get(name, ("correct",)):
                    if arm == "deranged" and k < 2:
                        continue
                    if arm not in arm_cache:
                        arm_cache[arm] = residual_arm(data["residual"], arm, rng, donor=donor)
                    residual = arm_cache[arm]
                    if isinstance(method, str):
                        delta = analytic_delta(
                            data["query_bits"], data["support_bits"], residual,
                            estimator=method, ridge=krr_ridge,
                        )
                    else:
                        method.eval()
                        delta, _ = method(
                            data["query_bits"], data["support_bits"], data["query_desc"],
                            data["support_desc"], data["protein"], residual,
                        )
                    predictions[f"{name}__{arm}"] = (data["base"] + delta).cpu().numpy()
            label = data["label"].cpu().numpy()
            mask = data["mask"].cpu().numpy()
            for offset, row_index_value in enumerate(rows):
                episode = block[int(row_index_value)]
                strata = stratum_of(nearest[offset])
                for stratum in (*STRATUM_NAMES, "all"):
                    active = mask[offset].copy()
                    if stratum != "all":
                        active &= strata == stratum
                    if int(active.sum()) < MIN_STRATUM_QUERY:
                        continue
                    truth = label[offset][active]
                    if float(np.std(truth)) < 1e-9:
                        continue
                    record: dict[str, object] = {
                        "policy": episode.policy,
                        "seed": episode.seed,
                        "draw": episode.draw,
                        "k": episode.k,
                        "target": episode.target,
                        "component": episode.component,
                        "stratum": stratum,
                        "n_query": int(active.sum()),
                        "nearest_tanimoto_mean": float(nearest[offset][active].mean()),
                    }
                    for name, values in predictions.items():
                        for metric, value in metric_loss(truth, values[offset][active]).items():
                            record[f"{name}__{metric}"] = float(value)
                    records.append(record)
    return pd.DataFrame.from_records(records)


def contrast_table(
    records: pd.DataFrame, contrasts: dict[str, tuple[str, str]], metrics: tuple[str, ...]
) -> dict[str, object]:
    summary: dict[str, object] = {}
    for policy in sorted(records.policy.unique()):
        for k in sorted(records.k.unique()):
            for stratum in sorted(records.stratum.unique()):
                frame = records.loc[
                    (records.policy == policy) & (records.k == k) & (records.stratum == stratum)
                ].copy()
                if frame.empty:
                    continue
                cell: dict[str, object] = {
                    "episodes": int(len(frame)),
                    "targets": int(frame.target.nunique()),
                    "components": int(frame.component.nunique()),
                    "absolute": {},
                    "contrasts": {},
                }
                columns = {name.rsplit("__", 1)[0] for name in frame.columns if "__" in name and name.endswith("__ci")}
                for name in sorted(columns):
                    cell["absolute"][name] = {
                        metric: float(np.nanmean(frame[f"{name}__{metric}"].to_numpy()))
                        for metric in metrics
                        if f"{name}__{metric}" in frame
                    }
                for label, (left, right) in contrasts.items():
                    for metric in metrics:
                        left_column, right_column = f"{left}__{metric}", f"{right}__{metric}"
                        if left_column not in frame or right_column not in frame:
                            continue
                        sign = -1.0 if metric == "rmse" else 1.0
                        frame["value"] = sign * (frame[left_column] - frame[right_column])
                        cell["contrasts"].setdefault(label, {})[metric] = paired_bootstrap(frame, "value")
                summary.setdefault(policy, {}).setdefault(f"k{int(k)}", {})[stratum] = cell
    return summary


# --------------------------------------------------------------------------- #
# Runner
# --------------------------------------------------------------------------- #


def split_fit_episodes(
    episodes: list[DrawnEpisode], substrate: Substrate
) -> tuple[list[DrawnEpisode], list[DrawnEpisode]]:
    folds = (
        substrate.labeled[["component", "oof_fold"]]
        .dropna()
        .drop_duplicates()
        .set_index("component")
        .oof_fold.to_dict()
    )
    train = [e for e in episodes if int(folds.get(e.component, 0)) != VALIDATION_FOLD]
    validation = [e for e in episodes if int(folds.get(e.component, 0)) == VALIDATION_FOLD]
    return train, validation


def select_krr_ridge(substrate: Substrate, episodes: list[DrawnEpisode]) -> tuple[float, dict[str, float]]:
    """Pick the KRR ridge on `fit` episodes only, so the baseline is not handicapped."""

    scores: dict[str, float] = {}
    by_k = group_by_k(episodes, substrate)
    for ridge in KRR_RIDGE_GRID:
        values: list[float] = []
        for k, (block, support, query) in by_k.items():
            if k not in ADMITTED_K:
                continue
            for start in range(0, len(block), BATCH_EPISODES):
                rows = np.arange(start, min(start + BATCH_EPISODES, len(block)))
                batch = make_batch(support, query, rows)
                data = gather_batch(substrate, batch)
                delta = analytic_delta(
                    data["query_bits"], data["support_bits"], data["residual"],
                    estimator="krr", ridge=ridge,
                )
                prediction = (data["base"] + delta).cpu().numpy()
                label = data["label"].cpu().numpy()
                mask = data["mask"].cpu().numpy()
                nearest = tanimoto(data["query_bits"], data["support_bits"]).amax(-1).cpu().numpy()
                for row in range(len(rows)):
                    active = mask[row] & (nearest[row] >= 0.55)
                    if int(active.sum()) < MIN_STRATUM_QUERY:
                        continue
                    value = episode_ci(prediction[row][active], label[row][active])
                    if np.isfinite(value):
                        values.append(value)
        scores[str(ridge)] = float(np.mean(values)) if values else float("nan")
    best = max(KRR_RIDGE_GRID, key=lambda value: scores[str(value)])
    return float(best), scores


def ladder(ridge: float) -> dict[str, TraceConfig]:
    """One rung per claim.  Each rung differs from its predecessor by one part."""

    off = {"global_scale": False, "modulation": False, "gate": False, "ridge": ridge}
    return {
        # Exact restrictions, run through the learned code path with every learned
        # part switched off.  They must reproduce the analytic operators bit for bit.
        "R2_nw": TraceConfig(weights="nw", whiten=False, **off),
        "R2b_krr": TraceConfig(weights="krr", whiten=True, **off),
        # The strong analytic bar: KRR plus one meta-learned, target-independent
        # transport scale.  This is the baseline the mechanism must beat.
        "R2c_krr_scaled": TraceConfig(
            weights="krr", whiten=True, global_scale=True, modulation=False, gate=False, ridge=ridge
        ),
        # Learned rungs, each adding exactly one part.
        "R3_gate": TraceConfig(
            weights="krr", whiten=True, global_scale=True, modulation=False, gate=True, ridge=ridge
        ),
        "R4_trace": TraceConfig(
            weights="krr", whiten=True, global_scale=True, modulation=True, gate=True, ridge=ridge
        ),
        "R5_modulation_only": TraceConfig(
            weights="krr", whiten=True, global_scale=True, modulation=True, gate=False, ridge=ridge
        ),
        "R6_protein_zero": TraceConfig(
            weights="krr", whiten=True, global_scale=True, modulation=True, gate=True,
            protein="zero", ridge=ridge,
        ),
        # Low-capacity modulation: interpretable pair scalars only, no learned
        # fingerprint bilinear.  If the full model overfits, this one should not.
        "R7_scalar_pairs": TraceConfig(
            weights="krr", whiten=True, global_scale=True, modulation=True, gate=True,
            pair_features="scalar", ridge=ridge,
        ),
    }


HEADLINE = "R4_trace"
ANALYTIC_BAR = "R2c_krr_scaled"
UNTRAINED_RUNGS = ("R2_nw", "R2b_krr")


def run(
    lock_path: Path,
    output: Path,
    records_path: Path,
    oof_cache: Path,
    *,
    seeds: tuple[int, ...] = EPISODE_SEEDS,
    max_epochs: int = MAX_EPOCHS,
    synthetic_spec: SyntheticSpec = SyntheticSpec(),
) -> dict[str, object]:
    if DEVICE.type != "cuda":
        raise RuntimeError("run this mechanism with D:\\anaconda\\envs\\drug\\python.exe")
    substrate, context = load_substrate(lock_path, oof_cache)

    fit_all = [e for e in build_episodes(substrate.labeled, "fit") if e.policy in TRAIN_POLICIES]
    # The probe role is measured under one declared episode seed; the three model
    # seeds vary the learned object, not the episode draw, so the paired
    # component bootstrap keeps the same episodes on both sides of a contrast.
    probe_all = [
        episode
        for episode in build_episodes(substrate.labeled, "probe")
        if episode.policy in TRAIN_POLICIES and episode.seed == EVALUATION_EPISODE_SEED
    ]
    fit_train, fit_val = split_fit_episodes(fit_all, substrate)
    krr_ridge, ridge_scores = select_krr_ridge(substrate, fit_val)

    trained: dict[str, dict[str, object]] = {}
    per_seed_records: list[pd.DataFrame] = []
    for seed in seeds:
        methods: dict[str, object] = {
            "level": "level",
            "nw": "nw",
            "krr": "krr",
            "cka_nnls": "cka_nnls",
        }
        mixture = StaticMixture()
        mixture_stats = train_model(
            mixture, substrate, fit_train, fit_val, seed=seed, max_epochs=max_epochs
        )
        methods["static_mixture"] = mixture
        trained.setdefault("static_mixture", {})[str(seed)] = {
            "best_epoch": mixture_stats["best_epoch"],
            "best_val_ci": mixture_stats["best_val_ci"],
            "weights": torch.softmax(mixture.logits.detach(), 0).cpu().tolist(),
        }
        for name, config in ladder(krr_ridge).items():
            model = Trace(config).to(DEVICE)
            if name in UNTRAINED_RUNGS:
                # Exact analytic restrictions: no parameter of theirs is used, so
                # training them would only move `log_ridge` away from the value
                # the baseline was granted.
                methods[name] = model
                trained.setdefault(name, {})[str(seed)] = {"trained": False, "config": vars(config)}
                continue
            stats = train_model(model, substrate, fit_train, fit_val, seed=seed, max_epochs=max_epochs)
            methods[name] = model
            trained.setdefault(name, {})[str(seed)] = {
                "trained": True,
                "best_epoch": stats["best_epoch"],
                "best_val_ci": stats["best_val_ci"],
                "parameters": int(sum(p.numel() for p in model.parameters())),
                "history": stats["history"],
            }
        frame = evaluate(
            substrate,
            probe_all,
            methods,
            arms={
                HEADLINE: ("correct", "deranged", "null", "signflip", "wrong_target"),
                "krr": ("correct", "deranged"),
            },
            krr_ridge=krr_ridge,
            seed=seed,
        )
        frame["model_seed"] = seed
        per_seed_records.append(frame)

    records = pd.concat(per_seed_records, ignore_index=True)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)

    contrasts = {
        # Headline: the learned transport against the strong analytic bar it nests.
        "trace_minus_krr_scaled": (f"{HEADLINE}__correct", f"{ANALYTIC_BAR}__correct"),
        "trace_minus_krr": (f"{HEADLINE}__correct", "R2b_krr__correct"),
        "trace_minus_analytic_krr": (f"{HEADLINE}__correct", "krr__correct"),
        "trace_minus_nw": (f"{HEADLINE}__correct", "R2_nw__correct"),
        "trace_minus_base": (f"{HEADLINE}__correct", "base__correct"),
        "trace_minus_static_mixture": (f"{HEADLINE}__correct", "static_mixture__correct"),
        "trace_minus_cka_nnls": (f"{HEADLINE}__correct", "cka_nnls__correct"),
        # Rung-by-rung increments.
        "rung_global_scale_gain": (f"{ANALYTIC_BAR}__correct", "R2b_krr__correct"),
        "rung_gate_gain": ("R3_gate__correct", f"{ANALYTIC_BAR}__correct"),
        "rung_modulation_gain": (f"{HEADLINE}__correct", "R3_gate__correct"),
        "rung_modulation_only_gain": ("R5_modulation_only__correct", f"{ANALYTIC_BAR}__correct"),
        "rung_protein_gain": (f"{HEADLINE}__correct", "R6_protein_zero__correct"),
        "rung_scalar_pairs_gain": ("R7_scalar_pairs__correct", f"{ANALYTIC_BAR}__correct"),
        # Controls.
        "krr_minus_base": ("krr__correct", "base__correct"),
        "trace_correct_minus_deranged": (f"{HEADLINE}__correct", f"{HEADLINE}__deranged"),
        "trace_correct_minus_wrong_target": (f"{HEADLINE}__correct", f"{HEADLINE}__wrong_target"),
        "trace_correct_minus_signflip": (f"{HEADLINE}__correct", f"{HEADLINE}__signflip"),
        "trace_null_minus_base": (f"{HEADLINE}__null", "base__correct"),
        "krr_correct_minus_deranged": ("krr__correct", "krr__deranged"),
        "level_minus_base": ("level__correct", "base__correct"),
        # Restriction audit: these must be numerically zero.
        "nesting_krr_restriction": ("R2b_krr__correct", "krr__correct"),
        "nesting_nw_restriction": ("R2_nw__correct", "nw__correct"),
    }
    summary = contrast_table(records, contrasts, ("ci", "ndcg10", "spearman", "rmse"))
    control = synthetic_positive_control(
        substrate, fit_train, fit_val, probe_all, krr_ridge, synthetic_spec, max_epochs
    )
    verdict = decide(summary, control)
    result: dict[str, object] = {
        "schema": "a2s-trace-q2-mechanism-v1",
        "status": "SOURCE_ONLY_DEVELOPMENT",
        "question": "Q2: is per-pair residual transport reliability a learnable, transferable function inside the Q1-admitted stratum?",
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
            "train_policies": list(TRAIN_POLICIES),
            "admitted_strata": list(ADMITTED_STRATA),
            "admitted_k": list(ADMITTED_K),
            "seeds": list(seeds),
            "evaluation_episode_seed": EVALUATION_EPISODE_SEED,
            "inner_validation_fold": VALIDATION_FOLD,
            "ladder": {name: vars(config) for name, config in ladder(krr_ridge).items()},
            "krr_ridge_selected_on_fit": krr_ridge,
            "krr_ridge_scores": ridge_scores,
            "max_epochs": max_epochs,
            "batch_episodes": BATCH_EPISODES,
            "headline_model": HEADLINE,
            "mde": ADMISSION_MDE,
        },
        "episodes": {
            "fit_train": len(fit_train),
            "fit_validation": len(fit_val),
            "probe": len(probe_all),
        },
        "training": trained,
        "summary": summary,
        "synthetic_positive_control": control,
        "decision": verdict,
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, default=float) + "\n",
        encoding="utf-8",
    )
    return result


def synthetic_positive_control(
    substrate: Substrate,
    fit_train: list[DrawnEpisode],
    fit_val: list[DrawnEpisode],
    probe: list[DrawnEpisode],
    krr_ridge: float,
    spec: SyntheticSpec,
    max_epochs: int,
    *,
    seed: int = EPISODE_SEEDS[0],
) -> dict[str, object]:
    """Train and evaluate the same ladder in a world where reliability is real.

    Without this, a null on real data is uninterpretable: it could equally mean
    "no signal" or "no power".  The control reports how much of an injected
    oracle gap the identical learner recovers on unseen probe components.
    """

    configs = ladder(krr_ridge)
    models: dict[str, object] = {"krr": "krr"}
    for name in (ANALYTIC_BAR, HEADLINE):
        model = Trace(configs[name])
        train_model(
            model, substrate, fit_train, fit_val, seed=seed, max_epochs=max_epochs, synthetic=spec
        )
        models[name] = model
    records = evaluate(substrate, probe, models, krr_ridge=krr_ridge, seed=seed, synthetic=spec)
    output: dict[str, object] = {"spec": vars(spec), "cells": {}}
    for k in ADMITTED_K:
        for stratum in ADMITTED_STRATA:
            frame = records.loc[(records.k == k) & (records.stratum == stratum)].copy()
            if frame.empty:
                continue
            absolute = {
                name: float(np.nanmean(frame[f"{name}__correct__ci"].to_numpy()))
                for name in ("base", "oracle", "krr", ANALYTIC_BAR, HEADLINE)
                if f"{name}__correct__ci" in frame
            }
            frame["value"] = frame[f"{HEADLINE}__correct__ci"] - frame[f"{ANALYTIC_BAR}__correct__ci"]
            recovered = paired_bootstrap(frame, "value")
            oracle_gap = absolute.get("oracle", float("nan")) - absolute.get(ANALYTIC_BAR, float("nan"))
            output["cells"][f"k{k}_{stratum}"] = {
                "absolute_ci": absolute,
                "oracle_gap": float(oracle_gap),
                "trace_minus_bar": recovered,
                "recovered_fraction": float(recovered["mean"] / oracle_gap) if oracle_gap > 0 else None,
            }
    detected = [
        cell["trace_minus_bar"]["lower95"] > ADMISSION_MDE for cell in output["cells"].values()
    ]
    output["detects_injected_reliability"] = bool(any(detected))
    return output


def decide(summary: dict[str, object], control: dict[str, object] | None = None) -> dict[str, object]:
    """Apply the preregistered M1/M1b/M4/M5 gates in the admitted stratum."""

    gates: list[dict[str, object]] = []
    for policy, by_k in summary.items():
        for k_label, by_stratum in by_k.items():
            k = int(k_label[1:])
            if k not in ADMITTED_K:
                continue
            for stratum in ADMITTED_STRATA:
                cell = by_stratum.get(stratum)
                if cell is None:
                    continue
                contrasts = cell["contrasts"]
                entry = {
                    "policy": policy,
                    "k": k,
                    "stratum": stratum,
                    "components": contrasts["trace_minus_krr_scaled"]["ci"]["components"],
                    "m1_ci_vs_krr_lower95": contrasts["trace_minus_krr_scaled"]["ci"]["lower95"],
                    "m1_ci_vs_krr_mean": contrasts["trace_minus_krr_scaled"]["ci"]["mean"],
                    "m1b_ndcg_vs_krr_lower95": contrasts["trace_minus_krr_scaled"]["ndcg10"]["lower95"],
                    "m4_assignment_lower95": contrasts["trace_correct_minus_deranged"]["ci"]["lower95"],
                    "m5_wrong_target_lower95": contrasts["trace_correct_minus_wrong_target"]["ci"]["lower95"],
                }
                entry["m1_pass"] = bool(entry["m1_ci_vs_krr_lower95"] > ADMISSION_MDE)
                entry["m1b_pass"] = bool(entry["m1b_ndcg_vs_krr_lower95"] > ADMISSION_MDE)
                entry["m4_pass"] = bool(entry["m4_assignment_lower95"] > 0.0)
                entry["m5_pass"] = bool(entry["m5_wrong_target_lower95"] > 0.0)
                gates.append(entry)
    headline_pass = all(bool(entry["m1_pass"]) for entry in gates) and bool(gates)
    powered = bool(control and control.get("detects_injected_reliability"))
    if headline_pass:
        verdict = "TRACE_BEATS_ANALYTIC_BAR"
    elif powered:
        verdict = "POSITIVELY_CONTROLLED_NULL_LEARNED_TRANSPORT_NOT_ADMITTED"
    else:
        verdict = "UNDERPOWERED_NULL_NOT_INTERPRETABLE"
    return {
        "gates": gates,
        "m1_all_pass": headline_pass,
        "positive_control_detects_injected_reliability": powered,
        "verdict": verdict,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--seeds", type=int, nargs="*", default=list(EPISODE_SEEDS))
    parser.add_argument("--max-epochs", type=int, default=MAX_EPOCHS)
    args = parser.parse_args()
    result = run(
        args.lock.resolve(),
        args.out.resolve(),
        args.records.resolve(),
        args.oof_cache.resolve(),
        seeds=tuple(int(value) for value in args.seeds),
        max_epochs=int(args.max_epochs),
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.out.resolve()),
                "episodes": result["episodes"],
                "decision": result["decision"],
            },
            indent=2,
            sort_keys=True,
            default=float,
        )
    )


if __name__ == "__main__":
    main()
