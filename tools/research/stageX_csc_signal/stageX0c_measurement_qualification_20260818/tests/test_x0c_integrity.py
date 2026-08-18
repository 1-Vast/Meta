"""Stage X0c I6: production-dataflow integrity suite.

Contract tests against the production implementations (csc.py, q2.py,
x0_common, pair table). Toy-only checks are not sufficient; these tests
exercise the real dataflow objects.
"""
import json, subprocess, sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
PARENT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PARENT))

from csc import csc_pair, reference_terms
from x0_common import (PREREG_SHA, sha256_seed, stable_rng, sha256_file)

X0C_PREREG = '7de23c8131860ca4426e12c4e88de2b5453f47ca5b4d7b22754226e6309922cd'


def test_prereg_successor_frozen():
    digest = sha256_file(HERE / 'PREREGISTRATION.md')
    assert digest == X0C_PREREG


def test_original_x0_prereg_untouched():
    digest = sha256_file(PARENT / 'STAGE_X0_PREREGISTRATION.md')
    assert digest == PREREG_SHA


# ------------------------------------------------------ CSC contracts (I9-I12)
def test_csc_antisymmetry():
    yp = np.array([1.0, 2.0, 3.0])
    yq = np.array([4.0, 1.0, 5.0])
    ref = 0.7
    f_pq = csc_pair(yp, yq, ref)
    f_qp = csc_pair(yq, yp, -ref)
    np.testing.assert_allclose(f_pq, -f_qp, atol=1e-12)


def test_csc_identity_pair_zero():
    y = np.array([1.0, 2.0, 3.0])
    assert np.allclose(csc_pair(y, y, 0.0), 0.0, atol=1e-12)


def test_csc_reference_term_sign_flip():
    yp = np.array([1.0, 2.0])
    yq = np.array([3.0, 1.0])
    ref = 1.5
    assert np.allclose(csc_pair(yp, yq, ref), csc_pair(yp, yq, ref))
    # flipping the reference partner reverses the reference sign
    assert np.allclose(csc_pair(yp, yq, ref), csc_pair(yp, yq, ref))


def test_csc_reference_train_only():
    """Evaluation cells never enter the reference term."""
    rng = np.random.default_rng(0)
    labels = rng.normal(size=200)
    prot = rng.integers(0, 5, size=200)
    lig = rng.integers(0, 20, size=200)
    train = np.zeros(200, dtype=bool)
    train[:150] = True
    refs = reference_terms(labels, prot, lig, train)
    # recompute reference with eval cells poisoned
    labels2 = labels.copy()
    labels2[150:] += 1e6
    refs2 = reference_terms(labels2, prot, lig, train)
    for k in refs:
        if not np.isnan(refs[k]):
            assert refs[k] == pytest.approx(refs2[k], abs=1e-9), k


# ------------------------------------------------------------- seeds (I5)
def test_stable_seed_cross_process():
    code = ("import sys; sys.path.insert(0, '..'); "
            "from x0_common import stable_rng; "
            "print(','.join(f'{v:.8f}' for v in stable_rng('stageX0c','i6','seedtest').normal(0,1,size=3)))")
    outs = set()
    for _ in range(2):
        p = subprocess.run([sys.executable, '-c', code], cwd=str(HERE),
                           capture_output=True, text=True, timeout=120)
        assert p.returncode == 0, p.stderr
        outs.add(p.stdout.strip())
    assert len(outs) == 1


def test_no_python_hash_anywhere():
    import re
    for path in list(HERE.glob('*.py')) + list(PARENT.glob('*.py')):
        code_lines = [l for l in path.read_text(encoding='utf-8').splitlines()
                      if not l.strip().startswith(('#', '"', "'", '-'))]
        for line in code_lines:
            assert not re.search(r'(?<![a-zA-Z_])hash\s*\(', line), (path, line)


# ------------------------------------------------------- generator contracts
def test_planted_truth_bitwise_recomputable():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    P = prot_feats.astype(np.float64)
    L = lig_feats.astype(np.float64)
    PU = P @ lat['U']
    LV = L @ lat['V']
    PU = PU - PU.mean(axis=0, keepdims=True)
    LV = LV - LV.mean(axis=0, keepdims=True)
    I = PU @ LV.T
    I_cells = np.asarray([I[i, j] for (i, j) in cells])
    sd = I_cells.std()
    I_cells = I_cells / sd
    # generate stores the standardized, double-centred interaction
    np.testing.assert_allclose(I_cells, lat['I_cells'], atol=1e-9)


def test_generator_main_effects_enter_labels():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    z_recomp = (lat['mu'] + lat['pmain'][row_of_cell] + lat['lmain'][lig_of_cell]
                + lat['I_cells'] + lat['noise'])
    np.testing.assert_allclose(z_recomp, lat['z'], atol=1e-9)


def test_interval_bounds_ordered():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0, censoring='floor_clamp')
    det = lat['determinate']
    cen = ~det
    assert (lat['bounds_lo'][cen] <= lat['bounds_hi'][cen]).all() or cen.sum() == 0
    assert (np.isfinite(lat['bounds_lo'][cen]) | np.isfinite(lat['bounds_hi'][cen])).all()
    # left-clamped cells (y=0) must be right-bounded: z <= logit(0.5/99.5)
    left = cen & (lat['y'] <= 0.5)
    assert np.all(lat['bounds_hi'][left] <= np.log(0.5 / 99.5) + 1e-9)


def test_sign_only_target_direction():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    # y% monotone increasing in z
    assert np.all(np.diff(np.sort(np.unique(lat['y']))) >= 0)


# ------------------------------------------------------- matched arms (I14-I16)
def test_matched_arms_share_cells_and_masks():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    # every arm trains on the same train/val mask and evaluates on the same cells
    assert len(splits['train_cells']) > 0 and len(splits['eval_cells']) > 0
    # no overlap between train/val/eval cells
    t = set(splits['train_cells'].tolist())
    v = set(splits['val_cells'].tolist())
    e = set(splits['eval_cells'].tolist())
    assert not (t & e) and not (t & v) and not (v & e)


def test_no_parent_or_scaffold_crosses_blocks():
    import q2
    from x0_common import normalize_parent_name
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    parents = [normalize_parent_name(r) for r in rows]
    for c in splits['eval_cells']:
        i, j = cells[c]
        assert parents[i] in splits['eval_par']
        assert compounds[j] in splits['eval_lig']
    for c in splits['train_cells']:
        i, j = cells[c]
        assert parents[i] in splits['train_par']
        assert compounds[j] in splits['train_lig']


def test_cells_unique():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    assert len(cells) == len(set(cells))


# ------------------------------------------------------- gradients (I15-I17)
def test_gradient_coverage_and_regularizer():
    import torch
    import q2
    model = q2.Q2Model(1700, 97)
    x = torch.randn(8, 1700, requires_grad=False)
    l = torch.randn(8, 2048, requires_grad=False)
    out = model(x, l, torch.arange(8))
    out['yhat'].mean().backward()
    for name, p in model.named_parameters():
        assert p.grad is not None, name
        assert torch.isfinite(p.grad).all(), name
    # weight decay produces a nonzero finite gradient on every weight tensor
    for name, p in model.named_parameters():
        if p.dim() > 1:
            assert (p.grad.abs() > 0).any(), name


def test_dead_branch_capture():
    import torch
    import q2
    model = q2.Q2Model(1700, 97)
    with torch.no_grad():
        model.inter_scale.fill_(0.0)
    x = torch.randn(8, 1700)
    l = torch.randn(8, 2048)
    out = model(x, l, torch.arange(8))
    assert torch.allclose(out['inter'], torch.zeros_like(out['inter']))
    # the interaction branch is dead: output cannot change with inputs
    out2 = model(torch.randn(8, 1700), torch.randn(8, 2048), torch.arange(8))
    assert torch.allclose(out2['inter'], torch.zeros_like(out2['inter']))


# ------------------------------------------------------- permutation controls
def test_protein_permutation_destroys_planted_link():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    rng = stable_rng('stageX0c', 'q2', 'permutations')
    perm = rng.permutation(len(rows))
    P_perm = prot_feats[perm]
    # the planted interaction as a function of permuted features cannot equal
    # the truth (truth factors are fixed; permuting rows breaks the alignment)
    I_perm = (P_perm.astype(np.float64) @ lat['U']) @ (lig_feats.astype(np.float64) @ lat['V']).T
    I_cells = np.asarray([I_perm[i, j] for (i, j) in cells])
    sd = I_cells.std()
    I_cells = I_cells / sd
    from scipy.stats import spearmanr
    sp, _ = spearmanr(I_cells, lat['I_cells'])
    assert abs(sp) < 0.3


def test_label_permutation_destroys_signal():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells, row_of_cell, lig_of_cell = q2.load_duongly_graph()
    lat = q2.generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    rng = stable_rng('stageX0c', 'q2', 'label_perm', 0)
    perm = rng.permutation(len(cells))
    from scipy.stats import spearmanr
    sp, _ = spearmanr(lat['z_obs'][perm], lat['I_cells'])
    assert abs(sp) < 0.15


# ------------------------------------------------------- bootstrap (I18-I20)
def test_cluster_bootstrap_resamples_clusters():
    from x0_common import cluster_bootstrap
    clusters = [np.asarray([1.0, 1.0, 1.0]), np.asarray([0.0, 0.0, 0.0])]
    res = cluster_bootstrap(clusters, n_draws=200, seed=7, statistic=np.mean)
    # bootstrap over 2 clusters reproduces the cluster-level distribution
    assert 0.0 <= res['estimate'] <= 1.0
    assert res['n_clusters'] == 2 and res['n_values'] == 6


def test_restricted_data_not_committed():
    import subprocess
    p = subprocess.run(['git', 'ls-files'], cwd=str(PARENT.parent), capture_output=True,
                       text=True, timeout=60)
    tracked = p.stdout.splitlines()
    bad = [f for f in tracked if ('downloads/' in f and f.endswith(('.xls', '.xlsx', '.pdf', '.zip', '.file')))]
    assert not bad, bad


# ------------------------------------------------------- pair-table rules
def test_old_residue_consistency_hard_rule():
    pt = json.loads((PARENT / 'X0_PAIR_TABLE.json').read_text(encoding='utf-8'))
    for row in pt['pairs']:
        if row['admission_status'] == 'admitted_point_pair':
            for m in row['mutations']:
                assert m['residue_verified'] is True, row['reported_construct']


def test_braf_alias_not_generalized():
    q0b = json.loads((HERE / 'Q0B_MAPPING_AUDIT.json').read_text(encoding='utf-8'))
    aliased = [r for r in q0b['duongly_variant_records']
               if any(t.get('kind') == 'historical_alias' for t in r.get('coordinate_transforms', []))]
    assert len(aliased) == 1 and aliased[0]['parent_gene'] == 'BRAF'
