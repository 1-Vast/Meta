"""Gate T0: is there a transferable target-adaptation object at all?

Seven mechanisms have now been falsified on this substrate (TRACE, MODE, IDA,
RIP, HOTSPOT, assay-coherence, MMP grammar, PIRS, CFES).  Every one of them
*assumed* two premises and tested neither:

``P1``  The per-target ranking headroom measured by Gate A0 is a **chemical**
        object -- a property of how the target discriminates compounds -- rather
        than measurement-context structure (per-document / per-assay offsets)
        that happens to be predictable from chemistry because a ChEMBL document
        reports a congeneric series.

``P2``  That object **transfers**: one target's fitted response head is useful
        for another target.  Gate G2 tested a weaker linear surrogate (does the
        head lie in a low-rank *subspace* of source heads) and refuted it, but
        no gate has ever applied source target ``j``'s head to recipient target
        ``t`` and measured the result.  Abundant-to-scarce transfer is
        meaningless if ``P2`` is false, and every mechanism so far has silently
        depended on it.

This gate measures both, plus the information budget that any support-label
estimator must live inside.  It fits closed-form ridge heads from labels -- for
each source target and for each recipient's own support -- but trains no
gradient model and learns no selector.  It is a measurement of the object, taken
before another estimator is built for it.

**Status: exploratory.** The source ``probe`` outcome was consumed once by PIRS,
so nothing measured here may be used for model selection.  Every arm is also
reported on same-document query pairs, because ``T0A`` finds a chemistry-free
per-document oracle that outscores the full chemical head.

``T0A``  nuisance decomposition -- split the own-head oracle gain into
         same-document query pairs (measurement offsets cancel; what survives is
         chemistry) and cross-document pairs, and compare against a
         document-offset oracle that uses labels but no chemistry.

``T0B``  transfer oracle -- apply every source head to every probe target and
         report the best, the median, the full-support-selected and the
         ``k``-shot-selected gain.  This is the learning curve of a *discrete*
         adaptation object, where Gate G4 measured it for a continuous one.

``T0C``  label-free proposal -- can protein embedding or library chemotype
         overlap rank a useful source head into a short list without labels?
         This decides how many bits the ``k`` support labels must supply.

``T0D``  information accounting -- measured ``tau``, ``sigma`` and the resulting
         Shannon bits available at ``k`` in {1,3,5}, reported beside the bits
         each candidate mechanism would need.

Only ``fit`` and ``probe`` roles are opened.  ``locked`` and the recipient
roster are never requested.
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

from research.a2s.a2s_information_gate import canonical, sha256_file
from research.a2s.a2s_trace import DEVICE, Substrate, load_substrate
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK, DEFAULT_OOF, paired_bootstrap
from research.a2s.a2s_mode_gates import build_basis, fit_head, source_heads
from research.a2s.a2s_mode_generalization import target_splits


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_transfer_object_gate_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_transfer_object_gate_records_2026-08-02.parquet"

SEED = 20260802
K_SWEEP = (1, 3, 5, 10, 20)
DRAWS = 16
BOOTSTRAP_DRAWS = 2000
MDE = 0.005
MIN_PAIRS = 8
PROPOSAL_SHORTLIST = 8
SELECT_RIDGE = 1.0

# Preregistered admission thresholds.  Set before the run; not revisable after.
T0A_MIN_WITHIN_DOCUMENT_GAIN = 0.005     # lower 95% bound on same-document CI gain
T0B_MIN_TRANSFER_RETENTION = 0.50        # best source head / own head, point estimate
T0B_MIN_FULL_SUPPORT_GAIN = 0.005        # lower 95% bound, head selected on train rows
T0C_MIN_PROPOSAL_LIFT = 0.0              # label-free shortlist minus random shortlist
T0D_MIN_K5_SELECTION_GAIN = 0.005        # lower 95% bound at k=5 over the frozen base
# Tightened after the first run, in the failing direction only.  The original
# T0D comparator was selection-minus-*random-selection*.  Beating a uniform draw
# from a library whose median member is harmful is not evidence of few-shot
# skill, so the binding comparators are now the frozen base and the single
# pooled head.  A stricter threshold cannot manufacture a positive result.
T0D_BINDING_COMPARATORS = ("minus_base", "minus_pooled")


# --------------------------------------------------------------------------- #
# Concordance, partitioned by an arbitrary pair grouping
# --------------------------------------------------------------------------- #


def pair_concordance(
    label: np.ndarray, prediction: np.ndarray, group: np.ndarray | None = None
) -> dict[str, float]:
    """Concordance index over all pairs, and over within-/between-group pairs.

    ``group`` labels the measurement context of each row (ChEMBL document).  A
    within-group pair shares its context, so any additive per-context offset
    cancels exactly and the comparison is a pure chemical one.
    """

    label = np.asarray(label, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if len(label) < 2:
        return {"ci": float("nan"), "ci_within": float("nan"), "ci_between": float("nan"),
                "pairs": 0, "pairs_within": 0, "pairs_between": 0}
    left, right = np.triu_indices(len(label), k=1)
    truth = np.sign(label[left] - label[right])
    pred = np.sign(prediction[left] - prediction[right])
    active = truth != 0
    correct = (pred == truth).astype(np.float64) + 0.5 * (pred == 0).astype(np.float64)

    def score(mask: np.ndarray) -> float:
        selected = mask & active
        if selected.sum() < MIN_PAIRS:
            return float("nan")
        return float(correct[selected].mean())

    everything = np.ones(len(left), dtype=bool)
    if group is None:
        same = np.zeros(len(left), dtype=bool)
    else:
        codes = pd.factorize(np.asarray(group, dtype=object))[0]
        same = codes[left] == codes[right]
    return {
        "ci": score(everything),
        "ci_within": score(same),
        "ci_between": score(~same),
        "pairs": int(active.sum()),
        "pairs_within": int((same & active).sum()),
        "pairs_between": int(((~same) & active).sum()),
    }


# --------------------------------------------------------------------------- #
# Head application and selection
# --------------------------------------------------------------------------- #


def head_scores(design: np.ndarray, heads: np.ndarray) -> np.ndarray:
    """(n, M) adaptation term of every source head on every row."""

    return design @ heads.T


def support_fit_loss(
    design: np.ndarray, residual: np.ndarray, heads: np.ndarray
) -> np.ndarray:
    """Centred squared-error of each source head against ``k`` support residuals.

    The level is removed on both sides, so this scores the *shape* of the head
    and never the target's mean offset -- the level channel is rank-null and an
    unshrunk anchor is measurably harmful.
    """

    if len(residual) < 2:
        return np.zeros(len(heads), dtype=np.float64)
    predicted = head_scores(design, heads)                    # (k, M)
    predicted = predicted - predicted.mean(axis=0, keepdims=True)
    observed = residual - residual.mean()
    return ((predicted - observed[:, None]) ** 2).mean(axis=0)


def document_offset_prediction(
    train_docs: np.ndarray, train_residual: np.ndarray, eval_docs: np.ndarray
) -> np.ndarray:
    """Label-using, chemistry-free oracle: the target's own per-document mean."""

    frame = pd.DataFrame({"doc": train_docs, "r": train_residual})
    means = frame.groupby("doc", sort=False)["r"].mean()
    return np.asarray([means.get(doc, 0.0) for doc in eval_docs], dtype=np.float64)


# --------------------------------------------------------------------------- #
# Label-free proposals
# --------------------------------------------------------------------------- #


def pairwise_tanimoto(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    """(n, m) Tanimoto between two unbatched binary fingerprint blocks."""

    intersection = left @ right.T
    union = left.sum(dim=1, keepdim=True) + right.sum(dim=1).unsqueeze(0) - intersection
    return intersection / union.clamp(min=1.0)


@dataclass
class Proposals:
    protein: np.ndarray          # (T, M) similarity, higher = proposed earlier
    chemotype: np.ndarray
    names: list[str]


def build_proposals(
    substrate: Substrate,
    splits: list,
    source_names: list[str],
) -> Proposals:
    """Two label-free orderings of the source heads for each probe target."""

    frame = substrate.labeled
    protein_row = {
        str(target): int(substrate.target_row[group.index.to_numpy()[0]].item())
        for target, group in frame.groupby("target", sort=False)
    }
    pooled = substrate.protein
    normalised = pooled / pooled.norm(dim=1, keepdim=True).clamp(min=1e-8)

    source_rows = {
        name: frame.index[(frame.role == "fit") & (frame.target == name)].to_numpy()
        for name in source_names
    }

    protein_similarity = np.zeros((len(splits), len(source_names)), dtype=np.float64)
    chemotype_similarity = np.zeros((len(splits), len(source_names)), dtype=np.float64)
    for index, split in enumerate(splits):
        query_protein = normalised[protein_row[split.target]]
        for column, name in enumerate(source_names):
            protein_similarity[index, column] = float(
                (query_protein * normalised[protein_row[name]]).sum().item()
            )
        train_bits = substrate.bits[torch.as_tensor(split.train_rows, device=DEVICE)]
        for column, name in enumerate(source_names):
            rows = source_rows[name]
            if len(rows) == 0:
                continue
            other = substrate.bits[torch.as_tensor(rows, device=DEVICE)]
            similarity = pairwise_tanimoto(train_bits, other)   # (n_train, n_source)
            chemotype_similarity[index, column] = float(similarity.max(dim=1).values.mean().item())
    return Proposals(protein_similarity, chemotype_similarity, list(source_names))


# --------------------------------------------------------------------------- #
# The gate
# --------------------------------------------------------------------------- #


def evaluate(
    substrate: Substrate,
    basis: np.ndarray,
    heads: np.ndarray,
    source_names: list[str],
    splits: list,
    proposals: Proposals,
    sigma: float,
) -> pd.DataFrame:
    residual_all = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    label_all = substrate.affinity.cpu().numpy().astype(np.float64)
    base_all = substrate.base.cpu().numpy().astype(np.float64)
    documents = substrate.labeled.docs.astype(str).to_numpy()
    rng_global = np.random.default_rng(SEED)
    records: list[dict[str, object]] = []

    for index, split in enumerate(splits):
        train_rows, eval_rows = split.train_rows, split.eval_rows
        train_design, eval_design = basis[train_rows], basis[eval_rows]
        train_residual = residual_all[train_rows]
        label, base = label_all[eval_rows], base_all[eval_rows]
        eval_docs = documents[eval_rows]

        row: dict[str, object] = {
            "target": split.target,
            "component": split.component,
            "split": split.split,
            "n_train": int(len(train_rows)),
            "n_eval": int(len(eval_rows)),
            "n_eval_documents": int(len(set(eval_docs.tolist()))),
        }

        # --- reference arms -------------------------------------------------
        base_metrics = pair_concordance(label, base, eval_docs)
        if not np.isfinite(base_metrics["ci"]):
            # Too few discordant pairs to score a ranking at all; a split that
            # cannot express a gain must not enter the denominator either.
            continue
        row["base__ci"] = base_metrics["ci"]
        row["base__ci_within"] = base_metrics["ci_within"]
        row["base__ci_between"] = base_metrics["ci_between"]
        row["pairs"] = base_metrics["pairs"]
        row["pairs_within"] = base_metrics["pairs_within"]
        row["pairs_between"] = base_metrics["pairs_between"]

        own_weight, _ = fit_head(train_design, train_residual)
        own = pair_concordance(label, base + eval_design @ own_weight, eval_docs)
        row["own__ci"] = own["ci"]
        row["own__ci_within"] = own["ci_within"]
        row["own__ci_between"] = own["ci_between"]

        offset = document_offset_prediction(documents[train_rows], train_residual, eval_docs)
        offset_metrics = pair_concordance(label, base + offset, eval_docs)
        row["docoffset__ci"] = offset_metrics["ci"]
        row["docoffset__ci_within"] = offset_metrics["ci_within"]

        # --- T0B: every source head ----------------------------------------
        scores = head_scores(eval_design, heads)              # (n_eval, M)
        gains = np.empty(len(heads), dtype=np.float64)
        gains_within = np.empty(len(heads), dtype=np.float64)
        for column in range(len(heads)):
            metrics = pair_concordance(label, base + scores[:, column], eval_docs)
            gains[column] = metrics["ci"]
            gains_within[column] = metrics["ci_within"]
        row["source_best__ci"] = float(np.nanmax(gains))
        row["source_median__ci"] = float(np.nanmedian(gains))
        row["source_mean__ci"] = float(np.nanmean(gains))
        row["source_best__ci_within"] = float(np.nanmax(gains_within))
        # Two definitions of "useful", reported side by side.  The half-of-oracle
        # rule is arbitrary and the oracle maximum is itself noisy, so the bit
        # account it feeds is definition-dependent and must never be quoted alone.
        row["source_positive_fraction"] = float(np.mean(gains > base_metrics["ci"]))
        row["source_useful_count_half_oracle"] = int(
            np.sum(gains >= base_metrics["ci"] + 0.5 * (np.nanmax(gains) - base_metrics["ci"]))
        )
        row["source_useful_count_beats_base"] = int(np.sum(gains > base_metrics["ci"]))
        row["source_useful_count_beats_mde"] = int(np.sum(gains > base_metrics["ci"] + MDE))
        pooled_head = heads.mean(axis=0)
        row["source_pooled__ci"] = pair_concordance(
            label, base + eval_design @ pooled_head, eval_docs
        )["ci"]

        # full-support selection: choose the head on the target's own train rows
        losses = support_fit_loss(train_design, train_residual, heads)
        chosen = int(np.argmin(losses))
        row["select_full__ci"] = gains[chosen]
        # T0A showed a per-document offset oracle beats the chemical head, so any
        # transfer claim must also be scored on same-document pairs, where those
        # offsets cancel.  Recording only the all-pair number would let document
        # structure masquerade as transferred chemistry.
        row["select_full__ci_within"] = gains_within[chosen]
        row["select_full__index"] = chosen
        row["select_full__is_oracle"] = bool(chosen == int(np.nanargmax(gains)))

        # --- T0C: label-free shortlists ------------------------------------
        for name, similarity in (
            ("protein", proposals.protein[index]),
            ("chemotype", proposals.chemotype[index]),
        ):
            order = np.argsort(-similarity, kind="stable")
            row[f"proposal_{name}_top1__ci"] = gains[order[0]]
            shortlist = order[:PROPOSAL_SHORTLIST]
            row[f"proposal_{name}_shortlist_best__ci"] = float(np.nanmax(gains[shortlist]))
            row[f"proposal_{name}_shortlist_mean__ci"] = float(np.nanmean(gains[shortlist]))
        random_shortlists = rng_global.integers(0, len(heads), size=(DRAWS, PROPOSAL_SHORTLIST))
        row["proposal_random_top1__ci"] = float(np.nanmean(gains[random_shortlists[:, 0]]))
        row["proposal_random_shortlist_best__ci"] = float(
            np.nanmean([np.nanmax(gains[draw]) for draw in random_shortlists])
        )

        # --- T0B/T0D: k-shot selection of a discrete head -------------------
        digest = int(sha256(f"{SEED}:{split.target}:{split.split}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(digest)
        for k in K_SWEEP:
            if len(train_rows) < k + 1:
                continue
            selected, random_pick, shortlist_selected = [], [], []
            selected_within, random_within = [], []
            for _ in range(DRAWS):
                support = rng.choice(len(train_rows), size=k, replace=False)
                loss = support_fit_loss(
                    train_design[support], train_residual[support], heads
                )
                pick = int(np.argmin(loss))
                draw = int(rng.integers(0, len(heads)))
                selected.append(gains[pick])
                selected_within.append(gains_within[pick])
                random_pick.append(gains[draw])
                random_within.append(gains_within[draw])
                order = np.argsort(-proposals.chemotype[index], kind="stable")[:PROPOSAL_SHORTLIST]
                shortlist_selected.append(gains[int(order[int(np.argmin(loss[order]))])])
            row[f"select_k{k}__ci"] = float(np.mean(selected))
            row[f"select_k{k}__ci_within"] = float(np.nanmean(selected_within))
            row[f"select_k{k}_random__ci"] = float(np.mean(random_pick))
            row[f"select_k{k}_random__ci_within"] = float(np.nanmean(random_within))
            row[f"select_k{k}_shortlist__ci"] = float(np.mean(shortlist_selected))

        records.append(row)
    return pd.DataFrame.from_records(records)


def provenance_audit(substrate: Substrate, splits: list) -> dict[str, object]:
    """How leaky is the 'leakage-free' split?

    The within-target split holds out Murcko scaffolds only.  It makes no claim
    about documents or assays, and this measures how far that falls short.
    """

    documents = substrate.labeled.docs.astype(str).to_numpy()
    assays = substrate.labeled.assays.astype(str).to_numpy()
    report: dict[str, object] = {}
    for name in sorted({split.split for split in splits}):
        chosen = [split for split in splits if split.split == name]
        document_reuse, assay_reuse, shared, saturated, sizes = [], [], 0, 0, []
        for split in chosen:
            train_documents = set(documents[split.train_rows].tolist())
            train_assays = set(assays[split.train_rows].tolist())
            fraction = float(np.mean([d in train_documents for d in documents[split.eval_rows]]))
            document_reuse.append(fraction)
            assay_reuse.append(
                float(np.mean([a in train_assays for a in assays[split.eval_rows]]))
            )
            shared += int(fraction > 0.0)
            saturated += int(fraction == 1.0)
            sizes.append(len(split.train_rows))
        report[name] = {
            "splits": len(chosen),
            "targets_sharing_any_document": shared,
            "targets_with_every_query_document_seen": saturated,
            "query_rows_reusing_a_support_document": float(np.mean(document_reuse)),
            "query_rows_reusing_a_support_assay": float(np.mean(assay_reuse)),
            "support_rows_mean": float(np.mean(sizes)),
            "support_rows_median": float(np.median(sizes)),
        }
    return report


def synthetic_control(
    basis: np.ndarray,
    heads: np.ndarray,
    splits: list,
    sigma: float,
) -> dict[str, object]:
    """Power check: a world in which the discrete transfer object exists exactly.

    Each synthetic target's residual is generated by one real source head plus
    measured noise.  If ``k``-shot selection cannot recover it here, the harness
    has no power and no negative result may be reported.
    """

    rng = np.random.default_rng(SEED + 1)
    outcome: dict[str, list[float]] = {f"k{k}": [] for k in K_SWEEP}
    recovery: dict[str, list[float]] = {f"k{k}": [] for k in K_SWEEP}
    for split in splits:
        truth = int(rng.integers(0, len(heads)))
        rows = np.concatenate([split.train_rows, split.eval_rows])
        signal = basis[rows] @ heads[truth]
        residual = signal + rng.normal(0.0, sigma, len(rows))
        label = residual                                      # base == 0 by construction
        n_train = len(split.train_rows)
        train_design, train_residual = basis[rows[:n_train]], residual[:n_train]
        eval_design, eval_label = basis[rows[n_train:]], label[n_train:]
        zero = np.zeros(len(eval_label))
        base_ci = pair_concordance(eval_label, zero)["ci"]
        scores = head_scores(eval_design, heads)
        gains = np.asarray(
            [pair_concordance(eval_label, scores[:, column])["ci"] for column in range(len(heads))]
        )
        for k in K_SWEEP:
            if n_train < k + 1:
                continue
            picked, correct = [], []
            for _ in range(DRAWS):
                support = rng.choice(n_train, size=k, replace=False)
                loss = support_fit_loss(train_design[support], train_residual[support], heads)
                choice = int(np.argmin(loss))
                picked.append(gains[choice] - base_ci)
                correct.append(float(choice == truth))
            outcome[f"k{k}"].append(float(np.mean(picked)))
            recovery[f"k{k}"].append(float(np.mean(correct)))
    summary = {
        "sigma": sigma,
        "modes": int(len(heads)),
        "chance_recovery": 1.0 / len(heads),
        "gain": {key: float(np.mean(value)) for key, value in outcome.items() if value},
        "recovery": {key: float(np.mean(value)) for key, value in recovery.items() if value},
    }
    summary["pass"] = bool(
        summary["gain"].get("k5", 0.0) > 0.02
        and summary["recovery"].get("k5", 0.0) > 5.0 * summary["chance_recovery"]
    )
    return summary


def information_budget(
    substrate: Substrate, basis: np.ndarray, heads: np.ndarray, splits: list
) -> dict[str, object]:
    """Measured ``tau``/``sigma`` and the Shannon bits they leave at each ``k``."""

    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    documents = substrate.labeled.docs.astype(str).to_numpy()
    assays = substrate.labeled.assays.astype(str).to_numpy()

    signal, noise, within_document, within_assay = [], [], [], []
    for split in splits:
        rows = np.concatenate([split.train_rows, split.eval_rows])
        weight, _ = fit_head(basis[split.train_rows], residual[split.train_rows])
        # ``tau`` must be a *held-out* signal scale.  The in-sample spread of a
        # 26-parameter ridge head fitted on ~60 rows is dominated by overfit and
        # would overstate the information the support labels can carry.  The
        # regression of held-out residual on the held-out fitted value recovers
        # the part of the head that is real.
        fitted = basis[split.eval_rows] @ weight
        centred = fitted - fitted.mean()
        held = residual[split.eval_rows] - residual[split.eval_rows].mean()
        variance = float(centred @ centred)
        if variance > 0:
            slope = float(centred @ held) / variance
            signal.append(abs(slope) * float(np.std(centred)))
        noise.append(float(np.std(residual[rows] - residual[rows].mean())))
        frame = pd.DataFrame(
            {"r": residual[rows], "doc": documents[rows], "assay": assays[rows]}
        )
        for column, sink in (("doc", within_document), ("assay", within_assay)):
            groups = frame.groupby(column)["r"]
            spread = groups.transform(lambda values: values - values.mean())
            counts = groups.transform("size")
            usable = spread[counts >= 3]
            if len(usable) >= 5:
                sink.append(float(np.std(usable)))

    tau = float(np.mean(signal))
    sigma = float(np.mean(noise))
    def bits(k: int, dispersion: float) -> float:
        if k < 2:
            return 0.0
        return float((k - 1) * 0.5 * np.log2(1.0 + tau**2 / dispersion**2))

    return {
        "status": "HEURISTIC_ORDER_OF_MAGNITUDE_ONLY",
        "caveats": [
            "assumes an independent Gaussian channel; support contrasts are not "
            "exchangeable channel uses because supports share chemical design",
            "tau is a held-out ridge projection scale, not the separation among "
            "candidate operators, which is the quantity a selector actually decodes",
            "sigma is total residual dispersion, not conditional observation noise",
            "no empirical mutual information or decoder-error bound was estimated",
            "bits_required depends on an arbitrary definition of a useful head; "
            "both definitions are reported and they differ by more than a factor of two",
        ],
        "tau_ranking_signal": tau,
        "sigma_residual": sigma,
        "sigma_within_document": float(np.mean(within_document)) if within_document else float("nan"),
        "sigma_within_assay": float(np.mean(within_assay)) if within_assay else float("nan"),
        "snr_per_observation": float(tau**2 / sigma**2),
        "bits_available": {f"k{k}": bits(k, sigma) for k in (1, 3, 5, 10, 20)},
        "bits_available_within_document": {
            f"k{k}": bits(k, float(np.mean(within_document))) if within_document else float("nan")
            for k in (1, 3, 5, 10, 20)
        },
        "bits_required_discrete_selection": float(np.log2(len(heads))),
        "bits_required_shortlist": float(np.log2(PROPOSAL_SHORTLIST)),
        "bits_required_dense_head": float(basis.shape[1]),
        "bits_required_note": (
            "log2(M / M_useful) under the two reported definitions of M_useful; "
            "the half-of-oracle rule and the beats-base rule differ by >2x, so no "
            "quantitative agreement with the measured break-even may be claimed"
        ),
    }


# --------------------------------------------------------------------------- #
# Summary and decision
# --------------------------------------------------------------------------- #


CONTRASTS = {
    "own_minus_base": ("own__ci", "base__ci"),
    "own_minus_base_within_document": ("own__ci_within", "base__ci_within"),
    "own_minus_base_between_document": ("own__ci_between", "base__ci_between"),
    "docoffset_minus_base": ("docoffset__ci", "base__ci"),
    "docoffset_minus_base_within_document": ("docoffset__ci_within", "base__ci_within"),
    "source_best_minus_base": ("source_best__ci", "base__ci"),
    "source_best_minus_base_within_document": ("source_best__ci_within", "base__ci_within"),
    "source_median_minus_base": ("source_median__ci", "base__ci"),
    "source_pooled_minus_base": ("source_pooled__ci", "base__ci"),
    "source_best_minus_own": ("source_best__ci", "own__ci"),
    "select_full_minus_base": ("select_full__ci", "base__ci"),
    "select_full_minus_base_within_document": ("select_full__ci_within", "base__ci_within"),
    "select_full_minus_pooled": ("select_full__ci", "source_pooled__ci"),
    "proposal_protein_top1_minus_random": ("proposal_protein_top1__ci", "proposal_random_top1__ci"),
    "proposal_chemotype_top1_minus_random": ("proposal_chemotype_top1__ci", "proposal_random_top1__ci"),
    "proposal_protein_shortlist_minus_random": (
        "proposal_protein_shortlist_best__ci", "proposal_random_shortlist_best__ci"),
    "proposal_chemotype_shortlist_minus_random": (
        "proposal_chemotype_shortlist_best__ci", "proposal_random_shortlist_best__ci"),
}


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for split in sorted(records.split.unique()):
        frame = records.loc[records.split == split].copy()
        cells: dict[str, object] = {}
        for name, (left, right) in CONTRASTS.items():
            if left not in frame.columns or right not in frame.columns:
                continue
            frame[name] = frame[left] - frame[right]
            cells[name] = paired_bootstrap(frame, name, draws=BOOTSTRAP_DRAWS)
        for k in K_SWEEP:
            column = f"select_k{k}__ci"
            if column not in frame.columns:
                continue
            for label, left, right in (
                ("minus_base", column, "base__ci"),
                ("minus_random", column, f"select_k{k}_random__ci"),
                ("minus_pooled", column, "source_pooled__ci"),
                ("within_document_minus_base", f"select_k{k}__ci_within", "base__ci_within"),
                (
                    "within_document_minus_random",
                    f"select_k{k}__ci_within",
                    f"select_k{k}_random__ci_within",
                ),
            ):
                if left not in frame.columns or right not in frame.columns:
                    continue
                name = f"select_k{k}_{label}"
                frame[name] = frame[left] - frame[right]
                cells[name] = paired_bootstrap(frame, name, draws=BOOTSTRAP_DRAWS)
            shortlist = f"select_k{k}_shortlist__ci"
            if shortlist in frame.columns:
                name = f"select_k{k}_shortlist_minus_full_library"
                frame[name] = frame[shortlist] - frame[column]
                cells[name] = paired_bootstrap(frame, name, draws=BOOTSTRAP_DRAWS)
        cells["descriptives"] = {
            "targets": int(frame.target.nunique()),
            "components": int(frame.component.nunique()),
            "mean_eval_rows": float(frame.n_eval.mean()),
            "mean_eval_documents": float(frame.n_eval_documents.mean()),
            "mean_pairs_within_fraction": float(
                (frame.pairs_within / frame.pairs.clip(lower=1)).mean()
            ),
            "mean_source_positive_fraction": float(frame.source_positive_fraction.mean()),
            "mean_source_useful_count_half_oracle": float(frame.source_useful_count_half_oracle.mean()),
            "mean_source_useful_count_beats_base": float(frame.source_useful_count_beats_base.mean()),
            "mean_source_useful_count_beats_mde": float(frame.source_useful_count_beats_mde.mean()),
            "full_support_selection_hits_oracle": float(frame.select_full__is_oracle.mean()),
        }
        summary[split] = cells
    return summary


def decide(
    summary: dict[str, object],
    control: dict[str, object],
    budget: dict[str, object],
) -> dict[str, object]:
    cells = summary.get("scaffold_disjoint", {})

    def cell(name: str) -> dict[str, float]:
        value = cells.get(name)
        return value if isinstance(value, dict) else {"mean": float("nan"), "lower95": float("nan")}

    gates: dict[str, object] = {}

    within = cell("own_minus_base_within_document")
    gates["T0A"] = {
        "question": "does the per-target headroom survive inside a single measurement context",
        "within_document": within,
        "between_document": cell("own_minus_base_between_document"),
        "document_offset_oracle": cell("docoffset_minus_base"),
        "pass": bool(within.get("lower95", float("nan")) > T0A_MIN_WITHIN_DOCUMENT_GAIN),
    }

    best = cell("source_best_minus_base")
    own = cell("own_minus_base")
    retention = (
        float(best.get("mean", float("nan")) / own.get("mean", float("nan")))
        if np.isfinite(own.get("mean", float("nan"))) and own.get("mean", 0.0) != 0.0
        else float("nan")
    )
    full = cell("select_full_minus_base")
    full_within = cell("select_full_minus_base_within_document")
    gates["T0B"] = {
        "question": "does any source target's head transfer to a recipient target",
        "oracle_best_source": best,
        "own_head": own,
        "retention": retention,
        "median_source": cell("source_median_minus_base"),
        "pooled_source": cell("source_pooled_minus_base"),
        "full_support_selected": full,
        # Binding.  An all-pair gain is confounded by the per-document offsets
        # that T0A measured as larger than the whole chemical head, so transfer
        # is only admitted if it survives on same-document pairs.
        "full_support_selected_within_document": full_within,
        "pass": bool(
            retention >= T0B_MIN_TRANSFER_RETENTION
            and full.get("lower95", float("nan")) > T0B_MIN_FULL_SUPPORT_GAIN
            and full_within.get("lower95", float("nan")) > T0B_MIN_FULL_SUPPORT_GAIN
        ),
    }

    protein = cell("proposal_protein_shortlist_minus_random")
    chemotype = cell("proposal_chemotype_shortlist_minus_random")
    gates["T0C"] = {
        "question": "can a label-free proposal shrink the hypothesis space",
        "protein_shortlist": protein,
        "chemotype_shortlist": chemotype,
        "pass": bool(
            max(protein.get("lower95", -1.0), chemotype.get("lower95", -1.0))
            > T0C_MIN_PROPOSAL_LIFT
        ),
    }

    curve = {
        f"k{k}": {
            comparator: cell(f"select_k{k}_{comparator}")
            for comparator in ("minus_base", "minus_pooled", "minus_random")
        }
        for k in K_SWEEP
    }
    binding = [cell(f"select_k5_{name}") for name in T0D_BINDING_COMPARATORS]
    crossing = next(
        (k for k in K_SWEEP if cell(f"select_k{k}_minus_base").get("mean", -1.0) > 0.0),
        None,
    )
    gates["T0D"] = {
        "question": "can k<=5 labels select the transferable object",
        "binding_comparators": list(T0D_BINDING_COMPARATORS),
        "learning_curve": curve,
        "full_support": cell("select_full_minus_base"),
        "break_even_k": crossing,
        "information_budget": budget,
        "pass": bool(
            all(band.get("lower95", float("nan")) > T0D_MIN_K5_SELECTION_GAIN for band in binding)
        ),
    }

    if not control.get("pass", False):
        verdict = "T0_HARNESS_INVALID_NO_POWER"
    elif not gates["T0A"]["pass"]:
        verdict = "NO_TRANSFERABLE_CHEMICAL_HEADROOM_OBJECT_IS_MEASUREMENT_CONTEXT"
    elif not gates["T0B"]["pass"]:
        verdict = "ADAPTATION_OBJECT_IS_NOT_TRANSFERABLE_ACROSS_TARGETS"
    elif not gates["T0D"]["pass"]:
        verdict = "TRANSFERABLE_AT_FULL_SUPPORT_BUT_NOT_IDENTIFIABLE_AT_K5"
    else:
        verdict = "DISCRETE_TRANSFER_ADMITTED_PROCEED_TO_MECHANISM"
    return {"gates": gates, "verdict": verdict, "synthetic_positive_control": control}


# --------------------------------------------------------------------------- #
# Runner
# --------------------------------------------------------------------------- #


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    started = time.time()
    substrate, provenance = load_substrate(lock_path, oof_cache)
    roles = set(substrate.labeled.role.unique())
    if roles - {"fit", "probe"}:
        raise RuntimeError(f"role firewall breached: {sorted(roles)}")

    basis, basis_stats = build_basis(substrate)
    heads, source_names, sigma, level_sd = source_heads(substrate, basis, role="fit")
    splits = target_splits(substrate, role="probe")
    if not splits:
        raise RuntimeError("no probe target produced a within-target split")

    proposals = build_proposals(substrate, splits, source_names)
    records = evaluate(substrate, basis, heads, source_names, splits, proposals, sigma)
    summary = summarise(records)
    budget = information_budget(substrate, basis, heads, splits)
    control = synthetic_control(basis, heads, splits, sigma)
    decision = decide(summary, control, budget)

    payload: dict[str, object] = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "device": str(DEVICE),
        "firewall": {
            "roles_opened": sorted(roles),
            "locked_requested": False,
            "recipient_requested": False,
            # Not "trains nothing": closed-form ridge heads ARE fitted from
            # labels, for source targets and for each recipient's own support.
            # What is absent is any gradient-trained model or learned selector.
            "fits_closed_form_ridge_heads_from_labels": True,
            "trains_no_gradient_model": True,
            "probe_role_status": (
                "probe outcomes were consumed once by PIRS; this gate is therefore "
                "EXPLORATORY evidence only and may not be used for model selection"
            ),
        },
        "provenance_audit": provenance_audit(substrate, splits),
        "protocol": {
            "seed": SEED,
            "draws": DRAWS,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "k_sweep": list(K_SWEEP),
            "shortlist": PROPOSAL_SHORTLIST,
            "mde": MDE,
            "thresholds": {
                "T0A_min_within_document_gain": T0A_MIN_WITHIN_DOCUMENT_GAIN,
                "T0B_min_transfer_retention": T0B_MIN_TRANSFER_RETENTION,
                "T0B_min_full_support_gain": T0B_MIN_FULL_SUPPORT_GAIN,
                "T0C_min_proposal_lift": T0C_MIN_PROPOSAL_LIFT,
                "T0D_min_k5_selection_gain": T0D_MIN_K5_SELECTION_GAIN,
            },
            "basis": basis_stats,
            "source_heads": {"count": int(len(heads)), "sigma": sigma, "level_sd": level_sd},
        },
        "data": {
            "rows": int(len(substrate.labeled)),
            "targets": int(substrate.labeled.target.nunique()),
            "probe_splits": int(len(splits)),
            "documents": int(substrate.labeled.docs.nunique()),
            "assays": int(substrate.labeled.assays.nunique()),
        },
        "information_budget": budget,
        "summary": summary,
        **decision,
    }
    # Write the records first so their digest can be embedded, then hash the
    # payload and write it exactly once.  The previous order hashed the file,
    # appended an artifact block and rewrote it, so the recorded digest never
    # matched the file on disk.
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    payload["artifacts"] = {
        "json": str(output),
        "records": str(records_path),
        "records_sha256": sha256_file(records_path),
        "lock": provenance["lock"].get("content_sha256") if isinstance(provenance.get("lock"), dict) else None,
        "json_sha256": "self-referential; verify with content_sha256 over the payload minus this block",
    }
    payload["content_sha256"] = sha256(canonical(payload).encode()).hexdigest()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="A2S Gate T0: transferable adaptation object")
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.lock, arguments.output, arguments.records, arguments.oof_cache)
    print(json.dumps({"verdict": payload["verdict"], "gates": {
        name: gate.get("pass") for name, gate in payload["gates"].items()
    }}, indent=2))


if __name__ == "__main__":
    main()
