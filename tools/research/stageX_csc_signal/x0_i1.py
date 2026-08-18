"""Stage X0 round-2: corrected I1 synthetic planted-signal control.

Corrects every round-1 defect of the uncommitted x0_planted.py draft (kept as
negative evidence):

1. einsum operand failure  -> explicit low-rank bilinear construction.
2. main effects generated but unused -> protein_main and ligand_main enter the
   label additively and are absorbed by dedicated heads at training time.
3. planted signal added on top of real % remaining values -> labels are built
   ENTIRELY from the known generative process; the real matrix contributes only
   the observation graph (which cells exist), the parent groups, the ligand
   scaffold groups, the size and the censoring semantics (clamp to [0,100]).
4. train/eval rows identical -> strict train/val/eval blocks by parent
   component AND ligand scaffold; eval = eval_parents x eval_ligands cells.
5. raw endpoint scored against interaction truth -> the fitted INTERACTION
   COMPONENT (interaction head output) is compared to the planted interaction.
6. ligand-only was an untrained zero vector -> a real capacity-matched
   ligand-only arm (protein input zeroed, same parameters and budget).
7. split ignored protein/pocket/scaffold novelty -> enforced (see 4).
8. raw endpoint offset dominating sign accuracy -> sign accuracy is computed
   between fitted interaction and planted interaction (both zero-centred by
   construction), never on raw endpoints.

Frozen gate (tau = 0.8, on held-out eval cells):
  sign accuracy >= 0.70 AND Spearman >= 0.30 AND
  sign_accuracy(correct_protein) - sign_accuracy(ligand_only) >= 0.05.
Other taus are the sensitivity/power curve only; they never select the result.
"""
from __future__ import annotations
import json, re
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from x0_common import (PREREG_SHA, HERE, sha256_seed, stable_rng, write_artifact,
                       load_duongly, cluster_bootstrap)

TAUS = [0.0, 0.2, 0.4, 0.8, 1.6]
GATE_TAU = 0.8
SIGN_ACC_GATE = 0.70
SPEARMAN_GATE = 0.30
GAP_GATE = 0.05
NOISE_SD = 1.0
MAIN_SD = 2.0
INTER_RANK = 4
PROT_DIM = 85 * 20          # KLIFS pocket one-hot
LIG_DIM = 2048              # ECFP4
HID = 32
DEAD_ZONE = 0.25            # logit units
BOOT_DRAWS = 2000
BOOT_SEED = 20260820
EPOCHS = 300
PATIENCE = 30
BATCH = 512
LR = 1e-3
WD = 1e-4


# ------------------------------------------------------------ data loading --
def load_features():
    """Protein pocket one-hot (97, 1700) and ligand ECFP (183, 2048), scaffolds."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    from rdkit.Chem.Scaffolds import MurckoScaffold
    info, matrix, seqs = load_duongly()
    pair_table = json.loads((HERE / 'X0_PAIR_TABLE.json').read_text(encoding='utf-8'))
    klifs = json.loads((HERE / 'klifs' / 'klifs_kinase_lookup.json').read_text(encoding='utf-8'))
    from x0_i2 import (klifs_pocket_for_parent, align_pocket_to_sequence,
                       mutate_pocket, pocket_onehot, KLIFS_QUERY_NAMES)

    rows = [str(x).strip() for x in matrix.iloc[1:, 0].tolist()]
    compounds = [re.sub(r'\s+', ' ', str(c)).strip() for c in matrix.columns[1:]]

    # protein pocket per assay row
    pocket_cache = {}
    prot_feats = np.zeros((len(rows), PROT_DIM), dtype=np.float32)
    row_meta = []
    for i, row in enumerate(rows):
        parent = _parent_of(row)
        if parent not in pocket_cache:
            seq = _parent_seq(parent, seqs)
            pocket, kid, note = klifs_pocket_for_parent(parent, klifs)
            align = align_pocket_to_sequence(pocket, seq) if pocket else None
            pocket_cache[parent] = (pocket, align, kid, note)
        pocket, align, kid, note = pocket_cache[parent]
        if pocket is None:
            pocket_vec = np.zeros(PROT_DIM, dtype=np.float32)
            row_meta.append({'row': row, 'parent': parent, 'pocket': None,
                             'note': note or 'no KLIFS pocket'})
        else:
            mut = _mutation_for_row(row, pair_table)
            vec_pocket = pocket
            pidx = None
            if mut and align is not None:
                _mp, pidx, _note = mutate_pocket(pocket, mut['pos'], mut['old'],
                                                 mut['new'], _parent_seq(parent, seqs), align)
                if pidx is not None:
                    vec_pocket = _mp
            pocket_vec = pocket_onehot(vec_pocket)
            row_meta.append({'row': row, 'parent': parent, 'pocket': bool(pocket),
                             'note': note, 'mutated_pocket_index': pidx})
        prot_feats[i] = pocket_vec

    # ligand ECFP
    pc = json.loads((HERE / 'pubchem_compounds.json').read_text(encoding='utf-8'))
    pc_clean = {re.sub(r'\s+', ' ', k).strip(): v for k, v in pc.items()}
    lig_feats = np.zeros((len(compounds), LIG_DIM), dtype=np.float32)
    lig_smiles = {}
    lig_resolved = []
    for j, name in enumerate(compounds):
        entry = pc_clean.get(name) or pc.get(name)
        smi = (entry.get('CanonicalSMILES') or entry.get('SMILES')
               or entry.get('IsomericSMILES')) if entry else None
        if smi:
            try:
                mol = Chem.MolFromSmiles(smi)
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=LIG_DIM)
                arr = np.zeros(LIG_DIM, dtype=np.float32)
                DataStructs.ConvertToNumpyArray(fp, arr)
                lig_feats[j] = arr
                lig_smiles[name] = smi
                lig_resolved.append(j)
            except Exception:
                smi = None
        if not smi:
            # deterministic hash features for unresolved names (train-only by split rule)
            rng = stable_rng('stageX0', 'ligand_hash_fallback', name)
            lig_feats[j, :128] = (rng.random(128) > 0.5).astype(np.float32)

    # Murcko scaffolds for resolved compounds
    scaffolds = {}
    for j, name in enumerate(compounds):
        smi = lig_smiles.get(name)
        if smi:
            try:
                mol = Chem.MolFromSmiles(smi)
                sc = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
                scaffolds[name] = sc or f'nosca:{j}'
            except Exception:
                scaffolds[name] = f'nosca:{j}'
        else:
            scaffolds[name] = f'unresolved:{j}'
    return (rows, compounds, prot_feats, lig_feats, scaffolds, row_meta)


def _parent_of(row):
    from x0_common import normalize_parent_name
    return normalize_parent_name(row)


def _parent_seq(parent, seqs):
    acc = {'ABL1': 'P00519', 'ALK': 'Q9UM73', 'BRAF': 'P15056', 'BTK': 'Q06187',
           'KIT': 'P10721', 'MET': 'P08581', 'SRC': 'P12931', 'CHEK2': 'O96017',
           'EGFR': 'P00533', 'FGFR1': 'P11362', 'FGFR2': 'P21802', 'FGFR3': 'P22607',
           'FGFR4': 'P22455', 'FLT3': 'P36888', 'JAK2': 'O60674', 'LRRK2': 'Q5S007',
           'MAP2K1': 'Q02750', 'MAPK14': 'Q16539', 'PDGFRA': 'P16234', 'RET': 'P07949',
           'TEK': 'Q02763'}[parent]
    return seqs[acc]['sequence']


def _mutation_for_row(row, pair_table):
    for r in pair_table['pairs']:
        if r['assay_row_s2'] == row and r['admission_status'] == 'admitted_point_pair':
            return r['mutations'][0]
    return None


# ---------------------------------------------------------------- splitting --
def make_splits(rows, compounds, prot_feats, lig_feats, scaffolds, matrix, seed=20260818):
    rng = np.random.default_rng(seed)
    parents = sorted({_parent_of(r) for r in rows}, key=lambda p: (-sum(_parent_of(r) == p for r in rows), p))
    parent_order = list(rng.permutation(len(parents)))
    eval_par = {parents[i] for i in parent_order[:5]}
    val_par = {parents[i] for i in parent_order[5:8]}
    train_par = {parents[i] for i in parent_order[8:]}

    sca_by_ligand = {c: scaffolds[c] for c in compounds}
    sca_names = sorted(set(sca_by_ligand.values()))
    sca_order = list(rng.permutation(len(sca_names)))
    eval_lig, val_lig, train_lig = set(), set(), set()
    for i, si in enumerate(sca_order):
        sca = sca_names[si]
        members = [c for c in compounds if sca_by_ligand[c] == sca]
        if si % 4 == 0:
            eval_lig.update(members)
        elif si % 4 == 1:
            val_lig.update(members)
        else:
            train_lig.update(members)
    # unresolved ligands never enter eval/val (no transferable features)
    eval_lig = {c for c in eval_lig if c in scaffolds and not scaffolds[c].startswith('unresolved')}
    val_lig = {c for c in val_lig if c in scaffolds and not scaffolds[c].startswith('unresolved')}

    labels = matrix.iloc[1:, 1:].to_numpy(dtype=float)
    row_parent = [_parent_of(r) for r in rows]
    cells = []
    for i, r in enumerate(rows):
        for j, c in enumerate(compounds):
            v = labels[i, j]
            if not np.isnan(v):
                cells.append((i, j))
    cell_idx = {c: k for k, c in enumerate(cells)}

    def block(pi, lj):
        out = []
        for (i, j) in cells:
            if row_parent[i] in pi and compounds[j] in lj:
                out.append(cell_idx[(i, j)])
        return np.asarray(out, dtype=np.int64)

    train_cells = block(train_par, train_lig)
    val_cells = np.concatenate([block(train_par, val_lig), block(val_par, train_lig)])
    eval_cells = block(eval_par, eval_lig)
    return {'parents': parents, 'train_par': sorted(train_par), 'val_par': sorted(val_par),
            'eval_par': sorted(eval_par), 'train_lig': sorted(train_lig),
            'val_lig': sorted(val_lig), 'eval_lig': sorted(eval_lig),
            'train_cells': train_cells, 'val_cells': val_cells, 'eval_cells': eval_cells,
            'n_scaffolds': len(sca_order), 'seed': seed}


# --------------------------------------------------------------- generation --
def generate_latents(cells, rows, compounds, prot_feats, lig_feats, tau, seed,
                     censoring='noclamp'):
    """Full generative process on the observed graph only.

    z(p,l) = protein_main(p) + ligand_main(l) + tau*I(p,l) + noise(p,l)
    I = low-rank bilinear in (protein pocket one-hot, ligand ECFP), unit-sd
    across observed cells, then scaled by tau (logit units).

    censoring modes (both apply the observation transform AFTER full latent
    generation; the primary realization mirrors the real Duong-Ly panel):
      'noclamp'     : y% = 100*logistic(z), fully determinate. This IS the real
                      Duong-Ly pattern: the panel reports continuous % remaining
                      including off-scale values (observed range -12.5 .. 191.3,
                      4023/17710 cells > 100); there is no hard assay floor.
      'floor_clamp' : y% = round(clip(100*logistic(z), 0, 100)) on the integer %
                      scale; cells at 0 / 100 are interval-censored with bounds
                      z <= logit(0.5/99.5) / z >= logit(99.5/0.5). Emulated
                      assay-floor realization used ONLY for censoring-machinery
                      validation and the floor-imputation control.
    """
    rng = stable_rng('stageX0', 'i1', 'tau', tau, 'seed', seed)
    n_cells = len(cells)
    pmain = rng.normal(0, MAIN_SD, size=len(rows)).astype(np.float64)
    lmain = rng.normal(0, MAIN_SD, size=len(compounds)).astype(np.float64)
    noise = rng.normal(0, NOISE_SD, size=n_cells).astype(np.float64)

    U = rng.normal(0, 1, size=(PROT_DIM, INTER_RANK)).astype(np.float64)
    V = rng.normal(0, 1, size=(LIG_DIM, INTER_RANK)).astype(np.float64)
    P = prot_feats.astype(np.float64)          # (n_prot, 1700)
    L = lig_feats.astype(np.float64)           # (n_lig, 2048)
    I_full = (P @ U) @ (L @ V).T               # (n_prot, n_lig)
    I_cells = np.asarray([I_full[i, j] for (i, j) in cells], dtype=np.float64)
    sd = I_cells.std()
    I_cells = (I_cells / sd) * tau if sd > 0 else I_cells * 0.0
    I_full = (I_full / sd) * tau if sd > 0 else I_full * 0.0

    z_cells = pmain[[i for (i, j) in cells]] + lmain[[j for (i, j) in cells]] + I_cells + noise
    y_pct = 100.0 / (1.0 + np.exp(-z_cells))
    if censoring == 'noclamp':
        determinate = np.ones(n_cells, dtype=bool)
        z_obs = z_cells.copy()
        bounds_lo = np.zeros(n_cells)
        bounds_hi = np.zeros(n_cells)
    elif censoring == 'floor_clamp':
        y_int = np.round(np.clip(y_pct, 0.0, 100.0))
        lo_b = np.log(0.5 / 99.5)
        hi_b = np.log(99.5 / 0.5)
        determinate = (y_int > 0) & (y_int < 100)
        z_obs = np.where(determinate, np.log(y_int / (100.0 - y_int)), np.nan)
        bounds_lo = np.where(determinate, 0.0, np.where(y_int <= 0, -np.inf, hi_b))
        bounds_hi = np.where(determinate, 0.0, np.where(y_int >= 100, np.inf, lo_b))
    else:
        raise ValueError(censoring)
    pmain_cells = pmain[[i for (i, j) in cells]]
    lmain_cells = lmain[[j for (i, j) in cells]]
    def corr(a, b):
        if np.std(a) == 0 or np.std(b) == 0:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])
    return {'tau': tau, 'seed': seed, 'censoring': censoring,
            'pmain': pmain, 'lmain': lmain, 'noise': noise,
            'U': U, 'V': V, 'I_full': I_full, 'I_cells': I_cells, 'z': z_cells,
            'y_pct': y_pct, 'determinate': determinate, 'z_obs': z_obs,
            'bounds_lo': bounds_lo, 'bounds_hi': bounds_hi,
            'n_censored': int((~determinate).sum()),
            'y_min': float(y_pct.min()), 'y_max': float(y_pct.max()),
            'I_mean': float(I_cells.mean()), 'I_sd': float(I_cells.std()),
            'I_rank': int(np.linalg.matrix_rank(I_full)),
            'I_pos_fraction': float((I_cells > 0).mean()),
            'corr_I_pmain': corr(I_cells, pmain_cells),
            'corr_I_lmain': corr(I_cells, lmain_cells)}


# -------------------------------------------------------------------- model --
class I1Model(nn.Module):
    def __init__(self, inter_rank=INTER_RANK):
        super().__init__()
        self.prot_enc = nn.Sequential(nn.Linear(PROT_DIM, HID), nn.Tanh())
        self.lig_enc = nn.Sequential(nn.Linear(LIG_DIM, HID), nn.Tanh())
        self.p_head = nn.Linear(HID, 1)
        self.l_head = nn.Linear(HID, 1)
        self.Wp = nn.Parameter(torch.randn(HID, inter_rank) * 0.1)
        self.Wl = nn.Parameter(torch.randn(HID, inter_rank) * 0.1)
        self.mu = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, prot, lig):
        pe = self.prot_enc(prot)
        le = self.lig_enc(lig)
        pmain = self.p_head(pe).squeeze(-1)
        lmain = self.l_head(le).squeeze(-1)
        inter = ((pe @ self.Wp) * (le @ self.Wl)).sum(-1) * self.inter_scale
        yhat = self.mu + pmain + lmain + inter
        return {'yhat': yhat, 'pmain': pmain, 'lmain': lmain, 'inter': inter}


class FreeIDModel(nn.Module):
    """Non-transferable upper bound: per-construct learned ID embedding."""
    def __init__(self, n_prot, inter_rank=INTER_RANK):
        super().__init__()
        self.emb = nn.Parameter(torch.randn(n_prot, HID) * 0.1)
        self.lig_enc = nn.Sequential(nn.Linear(LIG_DIM, HID), nn.Tanh())
        self.p_head = nn.Linear(HID, 1)
        self.l_head = nn.Linear(HID, 1)
        self.Wp = nn.Parameter(torch.randn(HID, inter_rank) * 0.1)
        self.Wl = nn.Parameter(torch.randn(HID, inter_rank) * 0.1)
        self.mu = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, prot_idx, lig):
        pe = self.emb[prot_idx]
        le = self.lig_enc(lig)
        pmain = self.p_head(pe).squeeze(-1)
        lmain = self.l_head(le).squeeze(-1)
        inter = ((pe @ self.Wp) * (le @ self.Wl)).sum(-1) * self.inter_scale
        yhat = self.mu + pmain + lmain + inter
        return {'yhat': yhat, 'pmain': pmain, 'lmain': lmain, 'inter': inter}


def _one_sided_huber(yhat, bound, margin=1.0):
    d = (yhat - bound) / margin
    return torch.where(d.abs() <= 1.0, 0.5 * d.square(), d.abs() - 0.5).mean()


def censored_loss(out, z_obs, det, blo, bhi):
    yhat = out['yhat']
    det_mask = torch.from_numpy(det).to(yhat.device)
    zt = torch.from_numpy(np.nan_to_num(z_obs, nan=0.0)).float().to(yhat.device)
    lo_t = torch.from_numpy(np.nan_to_num(blo, nan=0.0)).float().to(yhat.device)
    hi_t = torch.from_numpy(np.nan_to_num(bhi, nan=0.0)).float().to(yhat.device)
    mse = ((yhat - zt).square() * det_mask.float()).sum() / det_mask.float().sum().clamp(min=1)
    n_cen = (~det_mask).sum()
    left = ~det_mask & torch.isfinite(lo_t)
    right = ~det_mask & torch.isfinite(hi_t)
    loss = mse
    if left.any():
        loss = loss + _one_sided_huber(yhat[left], lo_t[left])
    if right.any():
        loss = loss + _one_sided_huber(yhat[right], hi_t[right])
    return loss


def train_arm(prot_feats, lig_feats, cells, row_of_cell, lig_of_cell, mask,
              z_obs, det, blo, bhi, arm, model, seed, device, report_interval=None):
    torch.manual_seed(seed)
    np.random.seed(seed % (2**32))
    for p in model.parameters():
        if p.dim() > 1:
            nn.init.xavier_uniform_(p)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    P = torch.from_numpy(prot_feats).float().to(device)
    L = torch.from_numpy(lig_feats).float().to(device)
    rows = torch.from_numpy(row_of_cell).long()
    ligs = torch.from_numpy(lig_of_cell).long()
    idx = np.arange(len(cells))
    best = None
    best_state = None
    n_train = len(mask)
    for ep in range(EPOCHS):
        model.train()
        rng = np.random.default_rng(seed + ep)
        perm = rng.permutation(n_train)
        tot = 0.0
        n = 0
        for b0 in range(0, n_train, BATCH):
            b = perm[b0:b0 + BATCH]
            c = idx[mask][b]
            if arm == 'free_target_id':
                out = model(rows[c], L[ligs[c]])
            else:
                out = model(P[rows[c]], L[ligs[c]])
            loss = censored_loss(out, z_obs[c], det[c], blo[c], bhi[c])
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss) * len(b)
            n += len(b)
        model.eval()
        with torch.no_grad():
            c = idx[mask]
            if arm == 'free_target_id':
                out = model(rows[c], L[ligs[c]])
            else:
                out = model(P[rows[c]], L[ligs[c]])
            train_loss = float(censored_loss(out, z_obs[c], det[c], blo[c], bhi[c]))
        if best is None or train_loss < best - 1e-6:
            best = train_loss
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            best_epoch = ep
        if ep - best_epoch >= PATIENCE:
            break
    model.load_state_dict(best_state)
    model.eval()
    return {'best_epoch': best_epoch, 'train_loss': best}


def predict(arm, model, prot_feats, lig_feats, rows, ligs, device, freeid_fallback=None):
    P = torch.from_numpy(prot_feats).float().to(device)
    L = torch.from_numpy(lig_feats).float().to(device)
    with torch.no_grad():
        if arm == 'free_target_id':
            out = model(rows, L[ligs])
        else:
            out = model(P[rows], L[ligs])
        return {k: v.cpu().numpy() for k, v in out.items()}


def spearman(a, b):
    from scipy.stats import spearmanr, pearsonr
    return float(spearmanr(a, b).correlation), float(pearsonr(a, b)[0])


def eval_metrics(inter_hat, I_true, det_mask=None):
    if det_mask is not None:
        inter_hat, I_true = inter_hat[det_mask], I_true[det_mask]
    sp, pe = spearman(inter_hat, I_true)
    nz = I_true != 0
    sign_acc = float((np.sign(inter_hat) == np.sign(I_true))[nz].mean()) if nz.any() else float('nan')
    dz = np.abs(I_true) > DEAD_ZONE
    dz_acc = float((np.sign(inter_hat) == np.sign(I_true))[dz].mean()) if dz.any() else float('nan')
    mse = float(np.mean((inter_hat - I_true) ** 2))
    slope = float(np.cov(inter_hat, I_true)[0, 1] / np.var(I_true)) if np.var(I_true) > 0 else float('nan')
    intercept = float(np.mean(inter_hat - I_true))
    return {'spearman': sp, 'pearson': pe, 'sign_accuracy': sign_acc,
            'dead_zone_sign_accuracy': dz_acc, 'interaction_mse': mse,
            'slope': slope, 'calibration_intercept': intercept}


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print('device:', device)
    rows, compounds, prot_feats, lig_feats, scaffolds, row_meta = load_features()
    info, matrix, seqs = load_duongly()
    splits = make_splits(rows, compounds, prot_feats, lig_feats, scaffolds, matrix)
    cells = []
    labels = matrix.iloc[1:, 1:].to_numpy(dtype=float)
    for i in range(len(rows)):
        for j in range(len(compounds)):
            if not np.isnan(labels[i, j]):
                cells.append((i, j))
    row_of_cell = np.asarray([i for (i, j) in cells], dtype=np.int64)
    lig_of_cell = np.asarray([j for (i, j) in cells], dtype=np.int64)
    print('cells:', len(cells), 'train/val/eval:',
          len(splits['train_cells']), len(splits['val_cells']), len(splits['eval_cells']))

    results = {}
    for tau in TAUS:
        lat = generate_latents(cells, rows, compounds, prot_feats, lig_feats, tau, seed=20260818)
        print(f"tau={tau}: n_censored={lat['n_censored']} I_sd={lat['I_sd']:.3f} "
              f"rank={lat['I_rank']} corr(I,pmain)={lat['corr_I_pmain']:.4f} "
              f"corr(I,lmain)={lat['corr_I_lmain']:.4f}")

        # ---- control constructions (permuted labels/features) ----
        rng_perm = stable_rng('stageX0', 'i1', 'tau', tau, 'permutations')
        lab_perm = rng_perm.permutation(len(cells))
        prot_perm = rng_perm.permutation(len(rows))
        prot_feats_perm = prot_feats[prot_perm]
        rnd_prot = stable_rng('stageX0', 'i1', 'tau', tau, 'random_prot').normal(0, 1, size=prot_feats.shape).astype(np.float32)
        prot_feats_zero = np.zeros_like(prot_feats)

        arms = {}
        for arm in ('ligand_only', 'correct_protein', 'shuffled_protein',
                    'random_protein', 'no_interaction_head'):
            P_arm = {'ligand_only': prot_feats_zero,
                     'correct_protein': prot_feats,
                     'shuffled_protein': prot_feats_perm,
                     'random_protein': rnd_prot,
                     'no_interaction_head': prot_feats}[arm]
            model = I1Model()
            model.to(device)
            if arm == 'no_interaction_head':
                model.inter_scale = nn.Parameter(torch.zeros(1).to(device))
            seed = sha256_seed('stageX0', 'i1', 'tau', tau, 'arm', arm)
            hist = train_arm(P_arm, lig_feats, cells, row_of_cell, lig_of_cell,
                             splits['train_cells'], lat['z_obs'], lat['determinate'],
                             lat['bounds_lo'], lat['bounds_hi'], arm, model, seed, device)
            arms[arm] = {'model': model, 'history': hist}

        # free target-ID upper bound (non-transferable; train cells only)
        fid = FreeIDModel(len(rows)).to(device)
        seed = sha256_seed('stageX0', 'i1', 'tau', tau, 'arm', 'free_target_id')
        hist = train_arm(prot_feats, lig_feats, cells, row_of_cell, lig_of_cell,
                         splits['train_cells'], lat['z_obs'], lat['determinate'],
                         lat['bounds_lo'], lat['bounds_hi'], 'free_target_id', fid, seed, device)
        arms['free_target_id'] = {'model': fid, 'history': hist}

        # label-permutation negative control: same graph, permuted labels
        if tau > 0:
            lab_perm_model = I1Model().to(device)
            seed = sha256_seed('stageX0', 'i1', 'tau', tau, 'arm', 'label_permuted')
            hist = train_arm(prot_feats, lig_feats, cells, row_of_cell, lig_of_cell,
                             splits['train_cells'], lat['z_obs'][lab_perm],
                             lat['determinate'][lab_perm], lat['bounds_lo'][lab_perm],
                             lat['bounds_hi'][lab_perm], 'label_permuted',
                             lab_perm_model, seed, device)
            arms['label_permuted'] = {'model': lab_perm_model, 'history': hist}

        tau_res = {}
        for arm, a in arms.items():
            model = a['model']
            if arm == 'no_interaction_head':
                with torch.no_grad():
                    model.inter_scale.fill_(0.0)
            eval_c = splits['eval_cells']
            train_c = splits['train_cells']
            pe = predict(arm, model, prot_feats if arm not in ('ligand_only',) else prot_feats_zero,
                         lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)
            inter_hat_eval = pe['inter']
            I_eval = lat['I_cells'][eval_c]
            m_eval = eval_metrics(inter_hat_eval, I_eval)
            # interaction component on TRAIN cells for the non-transferable ID arm
            pe_tr = predict(arm, model, prot_feats if arm not in ('ligand_only',) else prot_feats_zero,
                            lig_feats, row_of_cell[train_c], lig_of_cell[train_c], device)
            m_train = eval_metrics(pe_tr['inter'], lat['I_cells'][train_c])
            tau_res[arm] = {'eval': m_eval, 'train_interaction_recovery': m_train,
                            'best_epoch': a['history']['best_epoch'],
                            'train_loss': a['history']['train_loss']}

        # cluster bootstrap over eval parents for headline metrics
        bs = {}
        def signacc_boot(arm):
            inter_hat = predict(arm, arms[arm]['model'],
                                prot_feats if arm != 'ligand_only' else prot_feats_zero,
                                lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
            I_eval = lat['I_cells'][eval_c]
            clusters = []
            for p in sorted(set(_parent_of(rows[row_of_cell[c]]) for c in eval_c)):
                vals = np.asarray([(np.sign(inter_hat[k]) == np.sign(I_eval[k])) and I_eval[k] != 0
                                   for k, c in enumerate(eval_c)
                                   if _parent_of(rows[row_of_cell[c]]) == p], dtype=float)
                clusters.append(vals)
            return cluster_bootstrap(clusters, n_draws=BOOT_DRAWS, seed=BOOT_SEED, statistic=np.mean)
        bs['correct_protein'] = {'sign_accuracy': signacc_boot('correct_protein')}
        bs['ligand_only'] = {'sign_accuracy': signacc_boot('ligand_only')}

        # proper spearman bootstrap: resample parents, concat (inter_hat, I) pairs
        def spearman_boot(arm):
            inter_hat = predict(arm, arms[arm]['model'],
                                prot_feats if arm != 'ligand_only' else prot_feats_zero,
                                lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
            I_eval = lat['I_cells'][eval_c]
            clusters = []
            for p in sorted(set(_parent_of(rows[row_of_cell[c]]) for c in eval_c)):
                vals = [(inter_hat[k], I_eval[k]) for k, c in enumerate(eval_c)
                        if _parent_of(rows[row_of_cell[c]]) == p]
                clusters.append(np.asarray(vals))
            rng = np.random.default_rng(BOOT_SEED)
            stats = []
            for _ in range(BOOT_DRAWS):
                pick = rng.integers(0, len(clusters), size=len(clusters))
                pooled = np.concatenate([clusters[i] for i in pick])
                stats.append(spearman(pooled[:, 0], pooled[:, 1])[0])
            return {'estimate': spearman(inter_hat, I_eval)[0],
                    'ci_lo': float(np.percentile(stats, 2.5)),
                    'ci_hi': float(np.percentile(stats, 97.5)),
                    'n_clusters': len(clusters), 'seed': BOOT_SEED, 'draws': BOOT_DRAWS}
        bs['correct_protein']['spearman'] = spearman_boot('correct_protein')
        bs['ligand_only']['spearman'] = spearman_boot('ligand_only')

        # gap bootstrap (correct - ligand_only sign accuracy, same resamples)
        def gap_boot():
            ih_c = predict('correct_protein', arms['correct_protein']['model'], prot_feats,
                           lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
            ih_l = predict('ligand_only', arms['ligand_only']['model'], prot_feats_zero,
                           lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
            I_eval = lat['I_cells'][eval_c]
            clusters = []
            for p in sorted(set(_parent_of(rows[row_of_cell[c]]) for c in eval_c)):
                vals = [(ih_c[k], ih_l[k], I_eval[k]) for k, c in enumerate(eval_c)
                        if _parent_of(rows[row_of_cell[c]]) == p]
                clusters.append(np.asarray(vals))
            rng = np.random.default_rng(BOOT_SEED)
            stats = []
            for _ in range(BOOT_DRAWS):
                pick = rng.integers(0, len(clusters), size=len(clusters))
                pooled = np.concatenate([clusters[i] for i in pick])
                nz = pooled[:, 2] != 0
                if nz.any():
                    stats.append(float(
                        np.mean(np.sign(pooled[nz, 0]) == np.sign(pooled[nz, 2])) -
                        np.mean(np.sign(pooled[nz, 1]) == np.sign(pooled[nz, 2]))))
            gap0 = float(np.mean(np.sign(ih_c[I_eval != 0]) == np.sign(I_eval[I_eval != 0])) -
                         np.mean(np.sign(ih_l[I_eval != 0]) == np.sign(I_eval[I_eval != 0])))
            return {'estimate': gap0, 'ci_lo': float(np.percentile(stats, 2.5)),
                    'ci_hi': float(np.percentile(stats, 97.5)),
                    'n_clusters': len(clusters), 'seed': BOOT_SEED, 'draws': BOOT_DRAWS}
        bs['gap'] = gap_boot()

        # floor-imputation negative-control diagnostic at this tau
        floored = np.where(lat['determinate'], lat['z_obs'],
                           np.where(lat['y_pct'] <= 0.5, np.log(0.5 / 99.5), np.log(99.5 / 0.5)))
        inter_hat = predict('correct_protein', arms['correct_protein']['model'], prot_feats,
                            lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
        m_floor = eval_metrics(inter_hat, lat['I_cells'][eval_c])

        gate = None
        if tau == GATE_TAU:
            c = tau_res['correct_protein']['eval']
            l = tau_res['ligand_only']['eval']
            gate = {'sign_accuracy': c['sign_accuracy'], 'gate': SIGN_ACC_GATE,
                    'spearman': c['spearman'], 'spearman_gate': SPEARMAN_GATE,
                    'gap': c['sign_accuracy'] - l['sign_accuracy'], 'gap_gate': GAP_GATE,
                    'pass': (c['sign_accuracy'] >= SIGN_ACC_GATE and
                             c['spearman'] >= SPEARMAN_GATE and
                             c['sign_accuracy'] - l['sign_accuracy'] >= GAP_GATE)}

        neg = {}
        if tau == 0.0:
            # main-effect-only process: gate quantities undefined; the control is
            # that the correct arm shows no interaction structure
            inter_c = predict('correct_protein', arms['correct_protein']['model'], prot_feats,
                              lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
            inter_l = predict('ligand_only', arms['ligand_only']['model'], prot_feats_zero,
                              lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c], device)['inter']
            neg['tau0_no_spurious_interaction'] = {
                'correct_mean_abs_interaction': float(np.mean(np.abs(inter_c))),
                'ligand_only_mean_abs_interaction': float(np.mean(np.abs(inter_l))),
                'pass': float(np.mean(np.abs(inter_c))) <= 0.1,
                'note': 'main-effect-only generative process must not be classified as an interaction PASS: gate undefined at tau=0 and the interaction head must stay inert'}
        if tau > 0:
            neg['label_permutation'] = {
                'sign_accuracy': tau_res.get('label_permuted', {}).get('eval', {}).get('sign_accuracy'),
                'spearman': tau_res.get('label_permuted', {}).get('eval', {}).get('spearman'),
                'pass': (tau_res.get('label_permuted', {}).get('eval', {}).get('sign_accuracy') or 1.0) < SIGN_ACC_GATE,
                'note': 'permuting labels across cells must destroy the planted signal'}
            neg['protein_permutation'] = {
                'arm': 'shuffled_protein',
                'sign_accuracy': tau_res['shuffled_protein']['eval']['sign_accuracy'],
                'spearman': tau_res['shuffled_protein']['eval']['spearman'],
                'pass': tau_res['shuffled_protein']['eval']['sign_accuracy'] < SIGN_ACC_GATE,
                'note': 'permuting protein feature rows must destroy protein-conditioned interaction recovery'}
            neg['no_interaction_head'] = {
                'sign_accuracy': tau_res['no_interaction_head']['eval']['sign_accuracy'],
                'spearman': tau_res['no_interaction_head']['eval']['spearman'],
                'pass': tau_res['no_interaction_head']['eval']['sign_accuracy'] < SIGN_ACC_GATE,
                'note': 'removing the interaction head must fail interaction recovery'}
            neg['floor_imputation'] = {
                'sign_accuracy': m_floor['sign_accuracy'],
                'spearman': m_floor['spearman'],
                'note': 'floor-imputed censoring negative control; bias shown vs the censored-loss estimate (correct_protein arm above)'}
        results[str(tau)] = {
            'latent': {k: v for k, v in lat.items() if k not in ('U', 'V', 'I_full',
                                                                 'pmain', 'lmain', 'noise')},
            'arms': tau_res,
            'bootstrap': bs,
            'floor_imputation_negative_control': m_floor,
            'negative_controls': neg,
            'gate': gate,
        }
        print(f"  tau={tau}: correct sign_acc={tau_res['correct_protein']['eval']['sign_accuracy']:.4f} "
              f"spearman={tau_res['correct_protein']['eval']['spearman']:.4f} | "
              f"ligand_only sign_acc={tau_res['ligand_only']['eval']['sign_accuracy']:.4f} | "
              f"gap={tau_res['correct_protein']['eval']['sign_accuracy'] - tau_res['ligand_only']['eval']['sign_accuracy']:+.4f}"
              + (f" | GATE {'PASS' if gate and gate['pass'] else 'FAIL'}" if gate else ''))

    out = {
        'schema': 'MetaSieve.StageX.X0I1.v2',
        'stage': 'stageX_csc_signal',
        'preregistration_sha256': PREREG_SHA,
        'design': {
            'graph': 'Duong-Ly 2016 S2 observation graph (97 constructs x 183 compounds), real missingness, real parent groups, real ligand scaffold groups (PubChem Murcko via RDKit), clamp censoring to [0,100] applied after full latent generation',
            'label': 'z(p,l) = protein_main(p) + ligand_main(l) + tau*I(p,l) + noise(p,l); I = low-rank-4 bilinear of (KLIFS pocket one-hot, ECFP4), unit-sd then scaled by tau; y% = clamp(100*sigmoid(z), 0, 100); determinateness from the clamp (logit invertible on (0.5, 99.5))',
            'tau_units': 'logit (log-odds of % remaining activity) units: tau is the SD of the planted interaction across observed cells',
            'noise_sd': NOISE_SD, 'main_effect_sd': MAIN_SD, 'interaction_rank': INTER_RANK,
            'protein_feature': 'KLIFS 85-position pocket one-hot (1700 dims), mutated at the pocket index when the mutation maps',
            'ligand_feature': 'ECFP4 2048-bit (PubChem SMILES; deterministic hash fallback for unresolved names, train-only)',
            'model': 'protein encoder 1700->32 Tanh; ligand encoder 2048->32 Tanh; protein main-effect head; ligand main-effect head; rank-4 interaction head; AdamW lr=1e-3 wd=1e-4; batch 512; up to 300 epochs; early stopping patience 30 on val loss only',
            'arms': {
                'ligand_only': 'protein input zeroed (identical capacity/budget/optimizer)',
                'correct_protein': 'true protein features',
                'shuffled_protein': 'protein feature rows permuted (fixed seed)',
                'random_protein': 'SHA-256-seeded Gaussian, same shape',
                'no_interaction_head': 'interaction head scale pinned to 0',
                'free_target_id': 'per-construct learned ID embedding; NON-TRANSFERABLE upper bound, never a candidate model',
            },
            'split': 'train/val/eval blocks by parent component AND ligand scaffold; eval = eval_parents x eval_ligands observed cells; no parent construct and no scaffold crosses blocks; early stopping uses train+val only',
            'negative_controls': ['no_interaction_head must fail', 'tau=0 must show no protein-conditioned gain',
                                  'label permutation destroys the planted signal',
                                  'protein permutation destroys protein-conditioned interaction',
                                  'main-effect-only data must not be classified as interaction PASS',
                                  'floor-imputed censoring must show its bias'],
        },
        'splits': {k: v for k, v in splits.items() if 'cells' not in k},
        'cell_counts': {'total': len(cells), 'train': len(splits['train_cells']),
                        'val': len(splits['val_cells']), 'eval': len(splits['eval_cells'])},
        'results': results,
        'frozen_gate': {'tau': GATE_TAU, 'sign_accuracy': SIGN_ACC_GATE,
                        'spearman': SPEARMAN_GATE, 'gap': GAP_GATE,
                        'note': 'all three must hold on held-out eval cells for the correct_protein arm'},
    }
    inputs = [HERE / 'X0_PAIR_TABLE.json', HERE / 'downloads/duongly_mmc2.xlsx',
              HERE / 'downloads/duongly_mmc3.xlsx', HERE / 'klifs/klifs_kinase_lookup.json',
              HERE / 'pubchem_compounds.json']
    write_artifact(HERE / 'X0_I1.json', out, inputs)
    print(json.dumps({'gate': (out['results'].get(str(GATE_TAU)) or {}).get('gate'),
                      'cell_counts': out['cell_counts']}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
