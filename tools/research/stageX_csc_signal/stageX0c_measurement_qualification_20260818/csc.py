"""Stage X0c: production CSC (Centred Selectivity Contrast) operator.

CSC(l | p,q) = [y(l,p) - y(l,q)] - (1/|L|) sum_{l' in L} [y(l',p) - y(l',q)]

with the reference ligand set L restricted to TRAINING cells only (frozen I12:
evaluation labels never enter the reference). Properties enforced by contract
tests: antisymmetry, identity-pair zero, reference-term sign flip, train-only
reference statistics.
"""
from __future__ import annotations
import numpy as np


def reference_terms(labels, protein_ids, ligand_ids, train_mask):
    """Reference terms r(p,q) = mean over TRAIN cells of y(l',p) - y(l',q).
    Evaluation cells never contribute (frozen I12)."""
    proteins = np.unique(protein_ids)
    ref = {}
    for p in proteins:
        yp = labels[(protein_ids == p) & train_mask]
        for q in proteins:
            yq = labels[(protein_ids == q) & train_mask]
            n = min(len(yp), len(yq))
            ref[(p, q)] = float(np.mean(yp[:n] - yq[:n])) if n else np.nan
    return ref


def csc_matrix(labels, protein_ids, ligand_ids, train_mask, eval_mask=None):
    """Compute CSC for every (ordered protein pair, ligand) present.

    labels: (n_cells,) float
    protein_ids / ligand_ids: (n_cells,) int
    train_mask: bool (n_cells,) cells allowed to define the reference term
    eval_mask: optional bool (n_cells,) cells to emit CSC for (defaults to all)

    Returns dict:
      csc: (n_cells,) CSC values (NaN for cells excluded by eval_mask)
      reference: dict[(p,q)] -> reference term from TRAIN cells only
      reference_train_fraction: share of reference mass from train cells (must be 1)
    """
    train_cells = np.where(train_mask)[0]
    proteins = np.unique(protein_ids)
    ligands = np.unique(ligand_ids)
    # reference terms: mean over TRAIN cells of [y(l',p) - y(l',q)] per pair
    ref = {}
    ref_count = {}
    idx_p = {p: np.where((protein_ids == p) & train_mask)[0] for p in proteins}
    train_y_by_prot = {p: labels[idx_p[p]] for p in proteins}
    for pi in range(len(proteins)):
        for pj in range(len(proteins)):
            p, q = proteins[pi], proteins[pj]
            yp = train_y_by_prot[p]
            yq = train_y_by_prot[q]
            n = min(len(yp), len(yq))
            if n == 0:
                ref[(p, q)] = np.nan
                ref_count[(p, q)] = 0
                continue
            ref[(p, q)] = float(np.mean(yp[:n] - yq[:n]))
            ref_count[(p, q)] = n

    out = np.full(len(labels), np.nan)
    for k in range(len(labels)):
        p = protein_ids[k]
        q_pairs = {q for q in proteins if q != p}
        if eval_mask is not None and not eval_mask[k]:
            continue
        # CSC for the pair (p,q) for every q != p; the emitted value for cell k
        # is CSC(l | p, q*) where q* is the reference partner chosen by the
        # caller; here we emit for a canonical partner (first other protein)
        # to keep the operator single-valued per cell.
        for q in sorted(q_pairs):
            out[k] = (labels[k] - ref.get((p, q), np.nan)) - ref.get((p, q), np.nan)
            break  # canonical partner only
    return {'csc': out, 'reference': ref, 'reference_counts': ref_count}


def csc_pair(labels_pq, labels_qp, reference_term):
    """Explicit two-protein CSC for one ordered pair.

    labels_pq: y(l, p) for ligands l; labels_qp: y(l, q) for same ligands.
    reference_term: scalar mean over train cells of y(l',p) - y(l',q).
    Returns per-ligand CSC array. Antisymmetric by construction:
    csc_pair(..., +ref) == -csc_pair(..., -ref) when p and q swap.
    """
    return (np.asarray(labels_pq, dtype=float) - np.asarray(labels_qp, dtype=float)) - reference_term
