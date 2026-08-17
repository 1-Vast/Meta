"""Distance-stratified wrong-protein donors, and the ligand/label controls.

The v1 probe used exactly one donor per target — the *nearest* cross-component
protein. That is the hardest donor, which is the right choice for a
specificity claim, but it means a null could be an artifact of donor choice:
the nearest donor might simply be too similar to move anything. Phase 1 draws
donors at fixed quantiles of the similarity distribution so that perturbation
magnitude is a measured variable rather than a confound.

Every rule here is frozen in `_frozen.py` and label-blind: only pooled protein
embeddings, homology components and ligand structures are read.
"""
from __future__ import annotations

import numpy as np

from ._frozen import DONOR_STRATA


def whitened_protein_space(data, split: str, donor_pool: str,
                           whitening_pool: str):
    """The R2 `esm_whitened` geometry, rebuilt exactly as `matched_donors` does.

    Whitening statistics come from `meta_train` only, always: letting the
    evaluation split contribute its own covariance would leak population
    structure into the control that is supposed to be independent of it.
    """
    whitening_targets = sorted(data.tasks[whitening_pool])
    donor_targets = sorted(data.tasks[donor_pool])
    needed = set(whitening_targets) | set(donor_targets) | set(data.tasks[split])
    pooled = {t: np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
              for t in needed}
    matrix = np.stack([pooled[t] for t in whitening_targets])
    center = matrix.mean(0, keepdims=True)
    deviation = matrix - center
    covariance = deviation.T @ deviation / max(len(matrix) - 1, 1)
    values, vectors = np.linalg.eigh(covariance.astype(np.float64))
    whiten = ((vectors / np.sqrt(np.maximum(values, 1e-3)))
              @ vectors.T).astype(np.float32)

    def transform(vector: np.ndarray) -> np.ndarray:
        out = (vector - center[0]) @ whiten.T
        return out / max(float(np.linalg.norm(out)), 1e-9)

    return pooled, transform, donor_targets


def stratified_donors(data, split: str, donor_pool: str,
                      whitening_pool: str) -> dict[str, dict[str, tuple[str, float]]]:
    """For every target, one donor per similarity stratum.

    Returns `{target: {stratum: (donor_target, whitened_cosine)}}`. Candidates
    are restricted to a different homology component, which is what makes the
    substitution a legal wrong protein rather than a near-duplicate.
    """
    pooled, transform, donor_targets = whitened_protein_space(
        data, split, donor_pool, whitening_pool)
    component = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    embedded = {t: transform(pooled[t]) for t in set(donor_targets) | set(data.tasks[split])}
    bank = np.stack([embedded[t] for t in donor_targets])

    out: dict[str, dict[str, tuple[str, float]]] = {}
    for target in sorted(data.tasks[split]):
        similarity = bank @ embedded[target]
        legal = [(float(similarity[i]), donor_targets[i])
                 for i in range(len(donor_targets))
                 if component[donor_targets[i]] != component[target]]
        if not legal:
            raise ValueError(f"no cross-component donor exists for {target}")
        legal.sort(key=lambda item: -item[0])          # most similar first
        last = len(legal) - 1
        picks = {
            "nearest": 0,
            "q25": int(round(0.25 * last)),
            "median": int(round(0.50 * last)),
            "q75": int(round(0.75 * last)),
            "farthest": last,
        }
        out[target] = {name: (legal[index][1], legal[index][0])
                       for name, index in picks.items()
                       if name in DONOR_STRATA}
    return out


def novelty_and_scaffold_strata(data, scaffolds: dict[str, str],
                                train_split: str = "meta_train"):
    """Per-ligand max Tanimoto to `meta_train`, and whether its scaffold is seen.

    The double-cold split already removes exact and scaffold overlap by
    construction; these strata exist to check that any measured effect is not
    concentrated in whatever residual similarity survived.
    """
    fingerprints = data.fingerprints
    train_ligands = sorted({
        data.cells[int(i)]["ligand_id"]
        for target in data.tasks[train_split]
        for i in data.tasks[train_split][target]})
    train_matrix = np.stack([fingerprints[key].numpy() for key in train_ligands])
    train_sum = train_matrix.sum(1)
    train_scaffolds = {scaffolds.get(key, "") for key in train_ligands}
    train_scaffolds.discard("")

    cache: dict[str, float] = {}

    def novelty(ligand: str) -> float:
        if ligand not in cache:
            vector = fingerprints[ligand].numpy()
            intersection = train_matrix @ vector
            union = train_sum + vector.sum() - intersection
            cache[ligand] = float(
                (intersection / np.maximum(union, 1e-9)).max())
        return cache[ligand]

    def scaffold_seen(ligand: str) -> bool:
        return scaffolds.get(ligand, "") in train_scaffolds

    return novelty, scaffold_seen
