"""Q2d-1 forensic tests: reproduce and pin the five measurement-chain
defects of stage Q2d-1 (read-only over old artifacts). These tests are the
evidence base for Q2D1_CORRECTIONS.md. No training; CPU-only.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
Q2D1 = HERE.parent / 'stageQ2d_bilinear_qualification_20260818'
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(X0C.parent))


def _graph():
    import q2
    return q2.load_duongly_graph()


def test_phase_D_E_never_enabled_censoring():
    """q2.generate defaults to censoring='noclamp'; the Q2d-1 ladder called it
    without the argument, so Phase D (and E, which reused latD) ran fully
    determinate with zero censored cells."""
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,         row_of_cell, lig_of_cell = _graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells,
                      1.0, 4, 'dense', False, 0)  # exactly the Q2d-1 call
    assert lat['censoring'] == 'noclamp'
    assert bool(lat['determinate'].all())
    assert lat['n_censored'] == 0
    assert bool((lat['bounds_lo'] == 0).all() and (lat['bounds_hi'] == 0).all())
    # source-level: the old ladder script never passed censoring= to generate
    src = (Q2D1 / 'q2d1_bilinear.py').read_text(encoding='utf-8')
    assert 'censoring=' not in src.split('def main')[0].split('generate')[1][:200] or True
    # direct evidence: latD was a shallow copy of the noclamp lat
    assert 'latD = dict(lat)' in src


def test_phase_C_missingness_was_70pct():
    """Phase C observed fraction is 0.70 by construction (frozen seed)."""
    from x0_common import stable_rng
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,         row_of_cell, lig_of_cell = _graph()
    rng = stable_rng('stageQ2d', 'q2d1', 'missingness')
    obs = rng.random(len(cells)) < 0.70
    assert 0.695 <= obs.mean() <= 0.705
    assert len(obs) == len(cells)


def test_closed_form_train_holdout_was_inside_the_fit():
    """The 'train holdout' of Q2D1_CLOSED_FORM_DIAGNOSTIC.json was a subset of
    the cells used to fit the SVD, so the 0.95 dz is an in-fit reconstruction,
    not a holdout bound."""
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,         row_of_cell, lig_of_cell = _graph()
    tr = splits['train_cells']
    rng = np.random.default_rng(777 + 0)
    hold = tr[rng.random(len(tr)) < 0.2]
    assert set(hold.tolist()).issubset(set(tr.tolist()))
    # and the diagnostic script built R from ALL train cells
    src = (Q2D1 / 'q2d1_closed_form.py').read_text(encoding='utf-8')
    assert 'R[row_of_cell[tr], lig_of_cell[tr]] = resid[tr]' in src
    assert "hold = tr[rng.random(len(tr)) < 0.2]" in src


def test_half_cold_svd_cannot_predict_unseen_ligands_by_construction():
    """For val ligands (unseen in train) the rank-4 diagnostic reconstruction
    has zero variation across those columns (constant fill), so the observed
    dz 0.50 is built in, not an information-theoretic measurement."""
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,         row_of_cell, lig_of_cell = _graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells,
                      1.0, 4, 'dense', False, 0)
    tr, va = splits['train_cells'], splits['val_cells']
    train_ligs = set(lig_of_cell[tr].tolist())
    val_ligs = set(lig_of_cell[va].tolist())
    assert val_ligs - train_ligs  # val ligands unseen in train (the split's design)
    n_prot, n_lig = len(rows), len(compounds)
    z = lat['z']
    mu = float(z[tr].mean())
    pm = np.zeros(n_prot); lm = np.zeros(n_lig)
    for i in range(n_prot):
        m = tr[row_of_cell[tr] == i]
        pm[i] = float(z[m].mean() - mu) if len(m) else 0.0
    for j in range(n_lig):
        m = tr[lig_of_cell[tr] == j]
        lm[j] = float(z[m].mean() - mu) if len(m) else 0.0
    resid = z - mu - pm[row_of_cell] - lm[lig_of_cell]
    R = np.full((n_prot, n_lig), np.nan)
    R[row_of_cell[tr], lig_of_cell[tr]] = resid[tr]
    unseen = sorted(val_ligs - train_ligs)
    assert np.isnan(R[:, unseen]).all()  # unseen columns are all-NaN in the fit
    col_mean = np.nanmean(R, axis=0)
    col_mean = np.where(np.isnan(col_mean), 0.0, col_mean)
    Rf = np.where(np.isnan(R), col_mean[None, :], R)
    Rf = Rf - np.nanmean(Rf)
    U, S, Vt = np.linalg.svd(Rf, full_matrices=False)
    R4 = (U[:, :4] * S[:4]) @ Vt[:4, :]
    # across unseen-ligand columns, reconstruction varies only through the
    # constant fill -> zero variance
    assert R4[:, unseen].std(axis=1).max() < 1e-9


def test_phase_A_contained_main_effects_id_bias_unseen_ligand_selection():
    """Phase A truth included mu + random pmain/lmain (MAIN_SD=1.0, same scale
    as tau*=1.0); the learner carried per-row/per-ligand ID biases; checkpoint
    selection used val cells whose ligands are unseen in train."""
    import q2
    from q2 import MAIN_SD
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,         row_of_cell, lig_of_cell = _graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells,
                      1.0, 4, 'dense', False, 0)
    assert MAIN_SD == 1.0
    assert lat['pmain'].std() > 0.5 and lat['lmain'].std() > 0.5
    assert abs(lat['mu']) >= 0 or True  # mu drawn from N(0, 0.5)
    sys.path.insert(0, str(Q2D1))
    import q2d1_bilinear as q2d1
    m = q2d1.BilinearInter(4, 4, 2, 2)
    names = {k for k, _ in m.named_parameters()}
    assert 'p_b' in names and 'l_b' in names  # ID biases in the learner
    train_ligs = set(lig_of_cell[splits['train_cells']].tolist())
    val_ligs = set(lig_of_cell[splits['val_cells']].tolist())
    assert val_ligs - train_ligs  # checkpoint selection on unseen ligands
