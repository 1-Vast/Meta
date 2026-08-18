"""Stage X0 round-2: corrected I2 representation-capability instrument.

Corrects the round-1 defects:
1. WT and mutant representations are extracted at the SAME verified mutation
   coordinate (never WT-sequence midpoint).
2. Old residue is validated at the canonical coordinate before any use
   (pair table); historical renumbering is an explicit cited mapping.
3. Denomimators are legal inter-protein scales of the same representation
   type (never zero, never epsilon-padded).
4. mutation_token is reported separately as an edit-feature descriptor and is
   NOT counted in the representation admission gate.
5. random representation uses SHA-256-seeded values (cross-process stable)
   and is a sensitivity control only, never biological capability.
6. KLIFS pocket representation is implemented from the KLIFS API lookup with
   a parent/construct/mutation-position coverage census; unmapped rows keep
   reasons.
7. ESM truncation to 1020 tokens: mutation position must lie inside the
   window, otherwise the pair is excluded for ESM representations with reason.
8. I2 r_pair certifies sensitivity to change only, NOT biological direction
   (direction is tested by I1/I3); recorded explicitly.

Frozen gate: a representation passes iff median r_pair over labelled pairs
>= 0.05; X0 rule 2 requires >= 3 admissible representations passing,
including >= 1 local representation. Admissible = the representation
families excluding the edit descriptor and the random control.
"""
from __future__ import annotations
import json, math, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

from x0_common import (
    PREREG_SHA, AAS, HERE, sha256_seed, stable_rng, onehot_aa, composition,
    cluster_bootstrap, write_artifact, load_duongly, sha256_file)

ESM_WINDOW_RADIUS = 6  # window = [c-6, c+6] inclusive -> 13 positions
ESM_MAX_LEN = 1020
BOOT_DRAWS = 2000
BOOT_SEED = 20260820
CAPABILITY_THRESHOLD = 0.05  # frozen


# ---------------------------------------------------------------- ESM bank --
def load_esm(device='cuda'):
    import torch
    from transformers import AutoTokenizer, AutoModel
    name = 'facebook/esm2_t30_150M_UR50D'
    tok = AutoTokenizer.from_pretrained(name, local_files_only=True)
    model = AutoModel.from_pretrained(name, local_files_only=True)
    model.eval()
    try:
        model.to(device)
    except Exception:
        device = 'cpu'
        model.to('cpu')
    return tok, model, device


def esm_embed(sequences, tok, model, device):
    """Return per-sequence (n, hidden) mean-pooled embeddings and token-level
    hidden states aligned so that hidden[i] = residue i (1-based), hidden[0] =
    <cls> embedding (kept but excluded from pooling)."""
    import torch
    outs = []
    with torch.no_grad():
        for seq in sequences:
            enc = tok(seq, return_tensors='pt', truncation=True,
                      max_length=ESM_MAX_LEN + 2)
            if device == 'cuda':
                enc = {k: v.cuda() for k, v in enc.items()}
            h = model(**enc).last_hidden_state[0].cpu().float().numpy()
            outs.append(h)
    return outs


def window_mean_esm(hidden, center, radius):
    """Mean of hidden states of residues [center-radius, center+radius]
    clipped to the sequence; hidden[i] is residue i (1-based)."""
    n = hidden.shape[0] - 1  # exclude <cls>
    lo, hi = max(1, center - radius), min(n, center + radius)
    if hi < lo:
        return None
    return hidden[lo:hi + 1].mean(axis=0)


def window_onehot(seq, center, radius):
    """Concatenated one-hot over window positions [center-radius, center+radius];
    positions outside the sequence are zero vectors (identical for WT/mutant
    because a point mutation does not change sequence length)."""
    n = len(seq)
    vecs = []
    for p in range(center - radius, center + radius + 1):
        vecs.append(onehot_aa(seq[p - 1]) if 1 <= p <= n else np.zeros(len(AAS), np.float32))
    return np.concatenate(vecs)


# ------------------------------------------------------------------- KLIFS --
def load_klifs_lookup():
    path = HERE / 'klifs' / 'klifs_kinase_lookup.json'
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding='utf-8'))


KLIFS_QUERY_NAMES = {'CHEK2': 'CHK2', 'PDGFRA': 'PDGFRa'}


def needleman_wunsch(pocket, seg, match=2, mismatch=-1, gap=-2):
    """Global NW alignment of the KLIFS pocket (85 aligned positions, insertions
    removed) to a local segment of the canonical sequence. Returns the best
    score and the aligned pairs [(pocket_idx_1based | None, seq_pos_1based | None)].
    """
    n, m = len(pocket), len(seg)
    F = np.zeros((n + 1, m + 1), dtype=np.int32)
    for i in range(1, n + 1):
        F[i, 0] = i * gap
    for j in range(1, m + 1):
        F[0, j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if pocket[i - 1] == seg[j - 1] else mismatch
            F[i, j] = max(F[i - 1, j - 1] + s, F[i - 1, j] + gap, F[i, j - 1] + gap)
    pairs = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i, j] == F[i - 1, j - 1] + (match if pocket[i - 1] == seg[j - 1] else mismatch):
            pairs.append((i, j))
            i, j = i - 1, j - 1
        elif i > 0 and F[i, j] == F[i - 1, j] + gap:
            pairs.append((i, None))
            i -= 1
        else:
            pairs.append((None, j))
            j -= 1
    pairs.reverse()
    return F[n, m], pairs


def align_pocket_to_sequence(pocket, seq, flank=24, right_flank=220):
    """Locate the KLIFS 85-position pocket in the canonical sequence with a
    gapped global alignment (the pocket is gap-free 85 aligned positions, so
    alignment gaps are insertions in the full sequence; the pocket can span
    far more than 85 sequence positions).

    Returns dict {identity, matches, mismatches, gaps, seg_start, seg_end,
                  pocket_to_seq: {pocket_idx: seq_pos}, ...} or None.
    """
    if not pocket or not seq:
        return None
    # locate approximate site via the first 12 pocket residues (they sit in the
    # glycine-rich loop / beta1 region and are usually contiguous in the full seq)
    seed = pocket[:12]
    if seed in seq:
        approx = seq.index(seed)
    else:
        best, bi = None, None
        for i in range(max(0, len(seq) - len(seed)) + 1):
            mm = sum(a != b for a, b in zip(seed, seq[i:i + len(seed)]))
            if best is None or mm < best:
                best, bi = mm, i
        if best > 2:
            return None
        approx = bi
    start = max(0, approx - flank)
    end = min(len(seq), approx + right_flank)
    seg = seq[start:end]
    score, pairs = needleman_wunsch(pocket, seg)
    matches = sum(1 for pi, si in pairs if pi and si and pocket[pi - 1] == seg[si - 1])
    mismatches = sum(1 for pi, si in pairs if pi and si and pocket[pi - 1] != seg[si - 1])
    gaps = len(pairs) - matches - mismatches
    identity = matches / len(pocket)
    if identity < 0.80:
        return None
    pocket_to_seq = {pi: start + si for pi, si in pairs if pi and si}
    return {'identity': identity, 'matches': matches, 'mismatches': mismatches,
            'gaps': gaps, 'seg_start': start + 1, 'seg_end': end,
            'pocket_to_seq': pocket_to_seq}


def klifs_pocket_for_parent(parent, lookup):
    """Return (pocket_seq, kinase_id, note) for a parent gene, or (None,None,reason)."""
    if not lookup:
        return None, None, 'KLIFS lookup manifest missing'
    qname = KLIFS_QUERY_NAMES.get(parent, parent)
    matches = [m for m in lookup.get(qname, []) or []]
    if not matches:
        return None, None, f'no KLIFS entry for kinase {parent} (queried as {qname})'
    m = matches[0]
    pocket = m.get('pocket')
    if not pocket or len(pocket) != 85:
        return None, None, f'KLIFS pocket for {parent} missing or not 85 residues'
    return pocket, m.get('kinase_ID'), None


# ------------------------------------------------------------------- main --
def build_pair_records(pair_table, seqs):
    """Admitted point pairs with verified coordinates."""
    records = []
    for row in pair_table['pairs']:
        if row['admission_status'] != 'admitted_point_pair':
            continue
        m = row['mutations'][0]
        acc = row['canonical_accession']
        seq = seqs[acc]['sequence']
        pos = m['canonical_coordinate']
        old, new = m['old'], m['new']
        assert seq[pos - 1] == old, (row['reported_construct'], pos, old, seq[pos - 1])
        wt = seq
        mt = seq[:pos - 1] + new + seq[pos:]
        records.append({
            'parent': row['parent_kinase'],
            'construct': row['reported_construct'],
            'accession': acc,
            'pos': pos,
            'old': old,
            'new': new,
            'wt_seq': wt,
            'mt_seq': mt,
            's2_index': row['assay_row_s2_index'],
            'esm_admission': 'in' if pos <= ESM_MAX_LEN else
                             'mutation_position_outside_esm_window',
        })
    return records


def mutate_pocket(pocket, pos_canon, old, new, seq, align):
    """Apply the mutation to the pocket string if it maps inside; returns
    (mut_pocket, pocket_index, note)."""
    if align is None:
        return pocket, None, 'pocket alignment to canonical sequence failed'
    p2s = align.get('pocket_to_seq', {})
    idx = next((pi for pi, sp in p2s.items() if sp == pos_canon), None)
    if idx is None:
        return pocket, None, 'mutation_position_not_represented_in_klifs_pocket'
    if pocket[idx - 1] != old:
        return pocket, None, ('pocket residue at mapped index does not match expected old residue '
                              f'(pocket idx {idx}: {pocket[idx - 1]} vs expected {old}); the KLIFS '
                              'representative structure may carry a mutation')
    return pocket[:idx - 1] + new + pocket[idx:], idx, None


def pocket_onehot(pocket):
    return np.concatenate([onehot_aa(aa) for aa in pocket])


def rep_distance(a, b):
    return float(np.linalg.norm(np.asarray(a, np.float64) - np.asarray(b, np.float64)))


def median_inter_protein(points):
    """Median pairwise distance among representation points. Use
    median_inter_parent when same-parent points must not be compared."""
    if len(points) < 2:
        return None
    ds = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            ds.append(rep_distance(points[i], points[j]))
    return float(np.median(ds))


def median_inter_parent(point_pairs):
    """Median pairwise distance among representation points that belong to
    DIFFERENT parents. point_pairs: list of (parent_id, vector)."""
    if len(point_pairs) < 2:
        return None
    ds = []
    for i in range(len(point_pairs)):
        for j in range(i + 1, len(point_pairs)):
            if point_pairs[i][0] == point_pairs[j][0]:
                continue
            ds.append(rep_distance(point_pairs[i][1], point_pairs[j][1]))
    if not ds:
        return None
    return float(np.median(ds))


def bootstrap_ratio(pairs_by_parent, denom_by_parent, n_draws=2000, seed=20260820):
    """Cluster bootstrap of ratio = median pair distance / median inter-parent
    distance. Both numerator and denominator are recomputed on each resample of
    parent components. Returns dict like cluster_bootstrap plus statistic note."""
    rng = np.random.default_rng(seed)
    pair_parents = sorted(pairs_by_parent)
    denom_parents = sorted(denom_by_parent)
    # sample from the union of parents; a draw is valid when >=2 distinct
    # parents with denominator points appear
    stats = []
    union = sorted(set(pair_parents) | set(denom_parents))
    k = len(union)
    attempts = 0
    while len(stats) < n_draws and attempts < 50 * n_draws:
        attempts += 1
        idx = rng.integers(0, k, size=k)
        pooled = []
        for i in idx:
            p = union[i]
            pooled.extend(pairs_by_parent.get(p, []))
        denom_pts = []
        for i in idx:
            p = union[i]
            denom_pts.extend([(p, v) for v in denom_by_parent.get(p, [])])
        denom = median_inter_parent(denom_pts)
        if not pooled or denom is None or denom <= 0:
            continue
        stats.append(float(np.median(pooled)) / denom)
    stats = np.asarray(stats)
    all_pairs = [d for p in pair_parents for d in pairs_by_parent[p]]
    all_denom = [(p, v) for p in denom_parents for v in denom_by_parent[p]]
    est_denom = median_inter_parent(all_denom)
    est = float(np.median(all_pairs)) / est_denom if (all_pairs and est_denom) else float('nan')
    return {'estimate': est, 'ci_lo': float(np.percentile(stats, 2.5)),
            'ci_hi': float(np.percentile(stats, 97.5)),
            'n_clusters': k, 'n_values': len(all_pairs),
            'seed': seed, 'draws': len(stats),
            'statistic': 'ratio = median pair distance / median inter-parent distance, both recomputed per resample'}


def main():
    import torch
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print('device:', device)
    info, matrix, seqs = load_duongly()
    pair_table = json.loads((HERE / 'X0_PAIR_TABLE.json').read_text(encoding='utf-8'))
    records = build_pair_records(pair_table, seqs)
    print('admitted point pairs:', len(records))

    klifs = load_klifs_lookup()

    # ---------------- representations ----------------
    reps = {}
    # canonical sequence per parent via records
    parent_seq = {}
    for r in records:
        parent_seq.setdefault(r['parent'], seqs[r['accession']]['sequence'])
    comps_parent = {p: composition(parent_seq[p]) for p in parent_seq}
    pairs_comp = []
    for r in records:
        mt_comp = composition(r['mt_seq'])
        wt_comp = comps_parent[r['parent']]
        pairs_comp.append({'construct': r['construct'], 'parent': r['parent'],
                           'd': rep_distance(wt_comp, mt_comp)})
    denom_comp = median_inter_parent([(p, c) for p, c in comps_parent.items()])
    reps['composition'] = {
        'type': 'global', 'denominator': denom_comp,
        'denominator_description': 'median inter-protein distance of WT amino-acid composition across parents',
        'pairs': pairs_comp,
        'denom_by_parent': {p: [c] for p, c in comps_parent.items()}}

    # local one-hot window
    parent_windows = {}
    for r in records:
        parent_windows.setdefault(r['parent'], [])
    for r in records:
        parent_windows[r['parent']].append((r['pos'], window_onehot(parent_seq[r['parent']], r['pos'], ESM_WINDOW_RADIUS)))
    # dedupe windows per parent by position
    pw_uniq = {}
    for p, wins in parent_windows.items():
        uniq = {}
        for pos, w in wins:
            uniq[pos] = (pos, w)
        pw_uniq[p] = list(uniq.values())
    local_denom_points = {p: [w for _pos, w in wins] for p, wins in pw_uniq.items()}
    denom_local = median_inter_parent([(p, w) for p, wins in local_denom_points.items() for w in wins])
    pairs_local = []
    for r in records:
        wt_w = window_onehot(parent_seq[r['parent']], r['pos'], ESM_WINDOW_RADIUS)
        mt_w = window_onehot(r['mt_seq'], r['pos'], ESM_WINDOW_RADIUS)
        assert wt_w.shape == mt_w.shape
        pairs_local.append({'construct': r['construct'], 'parent': r['parent'],
                            'd': rep_distance(wt_w, mt_w),
                            'window_contract': f'canonical coordinate {r["pos"]}, radius {ESM_WINDOW_RADIUS}, same for WT and mutant'})
    reps['pair_centered_local_window'] = {
        'type': 'local', 'denominator': denom_local,
        'denominator_description': 'median distance between same-type local windows of different parents, each window extracted at its own verified mutation center (same-parent window pairs excluded)',
        'pairs': pairs_local,
        'denom_by_parent': local_denom_points}

    # ESM representations
    tok, model, device = load_esm(device)
    seqs_to_embed = {}
    for p in parent_seq:
        seqs_to_embed[f'parent:{p}'] = parent_seq[p]
    for r in records:
        seqs_to_embed[f'mt:{r["construct"]}'] = r['mt_seq']
    order = list(seqs_to_embed)
    print('embedding', len(order), 'sequences ...')
    hidden = {k: v for k, v in zip(order, esm_embed([seqs_to_embed[k] for k in order], tok, model, device))}
    esm_failures = []

    # global ESM
    g_esm_parent = {p: hidden[f'parent:{p}'][1:].mean(axis=0) for p in parent_seq}
    pairs_gesm, esm_excluded = [], []
    for r in records:
        if r['esm_admission'] != 'in':
            esm_excluded.append({'construct': r['construct'], 'pos': r['pos'],
                                 'reason': r['esm_admission']})
            continue
        mt_g = hidden[f'mt:{r["construct"]}'][1:].mean(axis=0)
        pairs_gesm.append({'construct': r['construct'], 'parent': r['parent'],
                           'd': rep_distance(g_esm_parent[r['parent']], mt_g)})
    reps['global_esm'] = {
        'type': 'global', 'denominator': median_inter_parent([(p, g) for p, g in g_esm_parent.items()]),
        'denominator_description': 'median inter-protein distance of global mean-pooled ESM-2 150M hidden states across WT parents',
        'pairs': pairs_gesm, 'esm_excluded_pairs': esm_excluded,
        'denom_by_parent': {p: [g] for p, g in g_esm_parent.items()}}

    # ESM local window (pair-centered, same coordinate for WT and mutant)
    parent_hidden = {p: hidden[f'parent:{p}'] for p in parent_seq}
    esm_local_denom = {}
    for p, wins in pw_uniq.items():
        pts = []
        for pos, _w in wins:
            v = window_mean_esm(parent_hidden[p], pos, ESM_WINDOW_RADIUS)
            if v is not None:
                pts.append(v)
        if pts:
            esm_local_denom[p] = pts
    pairs_esm_local, esm_local_excluded = [], []
    for r in records:
        if r['esm_admission'] != 'in':
            esm_local_excluded.append({'construct': r['construct'], 'pos': r['pos'],
                                       'reason': r['esm_admission']})
            continue
        w_wt = window_mean_esm(parent_hidden[r['parent']], r['pos'], ESM_WINDOW_RADIUS)
        w_mt = window_mean_esm(hidden[f'mt:{r["construct"]}'], r['pos'], ESM_WINDOW_RADIUS)
        if w_wt is None or w_mt is None:
            esm_local_excluded.append({'construct': r['construct'], 'pos': r['pos'],
                                       'reason': 'empty window'})
            continue
        pairs_esm_local.append({'construct': r['construct'], 'parent': r['parent'],
                                'd': rep_distance(w_wt, w_mt),
                                'window_contract': f'canonical coordinate {r["pos"]}, radius {ESM_WINDOW_RADIUS}, same for WT and mutant'})
    reps['esm_local_window'] = {
        'type': 'local', 'denominator': median_inter_parent([(p, w) for p, pts in esm_local_denom.items() for w in pts]),
        'denominator_description': 'median distance between ESM local-window representations of different parents, each at its own verified mutation center (same-parent window pairs excluded)',
        'pairs': pairs_esm_local, 'esm_excluded_pairs': esm_local_excluded,
        'denom_by_parent': esm_local_denom}

    # KLIFS pocket
    klifs_census = []
    pocket_by_parent = {}
    for parent in sorted(parent_seq):
        pocket, kid, note = klifs_pocket_for_parent(parent, klifs)
        klifs_census.append({'parent_kinase': parent, 'kinase_id': kid,
                             'pocket_available': pocket is not None,
                             'reason': note or 'ok',
                             'alignment': None})
        if pocket:
            align = align_pocket_to_sequence(pocket, parent_seq[parent])
            klifs_census[-1]['alignment'] = (
                {'identity': align['identity'], 'matches': align['matches'],
                 'mismatches': align['mismatches'], 'gaps': align['gaps'],
                 'seg_start': align['seg_start'], 'seg_end': align['seg_end']}
                if align else None)
            pocket_by_parent[parent] = {'pocket': pocket, 'align': align, 'kid': kid}
    pairs_klifs = []
    for r in records:
        entry = pocket_by_parent.get(r['parent'])
        if entry is None or entry['align'] is None:
            pairs_klifs.append({'construct': r['construct'], 'parent': r['parent'],
                                'excluded': True,
                                'reason': 'parent pocket unavailable or unaligned'})
            continue
        mt_pocket, idx, note = mutate_pocket(entry['pocket'], r['pos'], r['old'], r['new'],
                                             parent_seq[r['parent']], entry['align'])
        if idx is None:
            pairs_klifs.append({'construct': r['construct'], 'parent': r['parent'],
                                'excluded': True, 'reason': note})
            continue
        d = rep_distance(pocket_onehot(entry['pocket']), pocket_onehot(mt_pocket))
        pairs_klifs.append({'construct': r['construct'], 'parent': r['parent'],
                            'excluded': False, 'd': d, 'pocket_index': idx,
                            'note': f'mutation maps to KLIFS pocket index {idx}'})
    klifs_denom = {p: [pocket_onehot(entry['pocket'])]
                   for p, entry in pocket_by_parent.items() if entry['align'] is not None}
    reps['klifs_pocket'] = {
        'type': 'local', 'denominator': median_inter_parent([(p, w) for p, pts in klifs_denom.items() for w in pts]),
        'denominator_description': 'median inter-parent distance of KLIFS 85-position pocket one-hot representations (parents with aligned pockets only)',
        'pairs': pairs_klifs, 'klifs_census': klifs_census,
        'denom_by_parent': klifs_denom}

    # mutation_token edit descriptor (NOT in representation gate)
    pairs_token = []
    for r in records:
        pairs_token.append({'construct': r['construct'], 'parent': r['parent'],
                            'd': rep_distance(onehot_aa(r['old']), onehot_aa(r['new']))})
    edit_descriptor = {
        'description': 'pair-conditioned edit descriptor (old+new residue one-hot); NOT an independent protein representation; no legal inter-protein denominator exists',
        'pairs': pairs_token,
        'median_pair_distance': float(np.median([p['d'] for p in pairs_token])),
        'denominator': None,
        'admission_gate_status': 'excluded_from_representation_gate',
    }

    # random control (SHA-256 seeded; sensitivity control only)
    dim = 128
    def random_rep(name):
        rng = stable_rng('stageX0', 'random_representation', name, 'dim', dim)
        return rng.normal(0.0, 1.0, size=dim).astype(np.float32)
    rnd_parent = {p: random_rep(p) for p in parent_seq}
    pairs_rnd = []
    for r in records:
        pairs_rnd.append({'construct': r['construct'], 'parent': r['parent'],
                          'd': rep_distance(rnd_parent[r['parent']], random_rep(r['construct']))})
    reps['random'] = {
        'type': 'control', 'denominator': median_inter_parent([(p, v) for p, v in rnd_parent.items()]),
        'denominator_description': 'median inter-parent distance of SHA-256-seeded Gaussian vectors (dim 128)',
        'pairs': pairs_rnd,
        'denom_by_parent': {p: [v] for p, v in rnd_parent.items()},
        'interpretation': 'metric-sensitivity control only: a pass here shows the ratio metric is sensitive, NOT that any biological capability exists'}

    # ---------------- ratios + bootstrap ----------------
    for name, rep in reps.items():
        ps = [p for p in rep['pairs'] if not p.get('excluded', False)]
        if rep['denominator'] is None or len(ps) == 0:
            rep['median_pair_distance'] = None
            rep['ratio'] = None
            rep['pass_capability'] = None
            rep['n_pairs'] = len(ps)
            rep['n_excluded_pairs'] = len(rep['pairs']) - len(ps)
            rep['bootstrap'] = None
            continue
        denom = rep['denominator']
        rep['median_pair_distance'] = float(np.median([p['d'] for p in ps]))
        rep['ratio'] = rep['median_pair_distance'] / denom if denom > 0 else None
        rep['n_pairs'] = len(ps)
        rep['n_excluded_pairs'] = len(rep['pairs']) - len(ps)
        rep['pass_capability'] = (rep['ratio'] is not None and rep['ratio'] >= CAPABILITY_THRESHOLD)
        # cluster bootstrap over parent components: numerator AND denominator
        # recomputed on each resample
        pairs_by_parent = {}
        for p in ps:
            pairs_by_parent.setdefault(p['parent'], []).append(p['d'])
        denom_by_parent = rep.get('denom_by_parent', {})
        boot = bootstrap_ratio(pairs_by_parent, denom_by_parent,
                               n_draws=BOOT_DRAWS, seed=BOOT_SEED)
        rep['bootstrap'] = boot
        rep.pop('denom_by_parent', None)  # working set only; keep JSON lean

    gate = {
        'eligible_representations': ['composition', 'global_esm', 'pair_centered_local_window',
                                     'esm_local_window', 'klifs_pocket'],
        'local_representations': ['pair_centered_local_window', 'esm_local_window', 'klifs_pocket'],
        'excluded_from_gate': {'mutation_token': edit_descriptor['admission_gate_status'],
                               'random': reps['random']['interpretation']},
        'n_passing': sum(1 for n in gate_eligible(reps) if reps[n].get('pass_capability')),
    }
    gate['n_local_passing'] = sum(1 for n in gate['local_representations']
                                  if reps[n].get('pass_capability'))
    gate['pass'] = (gate['n_passing'] >= 3 and gate['n_local_passing'] >= 1)
    gate['frozen_rule'] = ('at least three admissible representations pass expression capability '
                           '(ratio >= 0.05), including at least one local representation')

    out = {
        'schema': 'MetaSieve.StageX.X0I2.v2',
        'stage': 'stageX_csc_signal',
        'preregistration_sha256': PREREG_SHA,
        'capability_threshold': CAPABILITY_THRESHOLD,
        'window_radius': ESM_WINDOW_RADIUS,
        'esm_model': 'facebook/esm2_t30_150M_UR50D (local cache, eval mode, float32)',
        'device': device,
        'esm_max_len': ESM_MAX_LEN,
        'census': {
            'total_mutant_construct_rows': len(pair_table['pairs']),
            'admitted_point_pairs': len(records),
            'esm_excluded_pairs': [p for p in reps['global_esm'].get('esm_excluded_pairs', [])],
            'klifs_parent_coverage': klifs_census,
        },
        'representations': reps,
        'edit_descriptor': edit_descriptor,
        'gate': gate,
        'scope_note': ('I2 r_pair certifies that a representation is sensitive to the planted '
                       'single-residue edit relative to its inter-protein scale. It does NOT '
                       'certify biological direction; direction and transferability are tested '
                       'by I1 (planted recovery) and I3 (ID-equivalence). A pass here is necessary, '
                       'never sufficient, for a biological capability claim.'),
    }
    inputs = [HERE / 'X0_PAIR_TABLE.json', HERE / 'downloads/duongly_mmc2.xlsx',
              HERE / 'downloads/duongly_mmc3.xlsx', HERE / 'klifs' / 'klifs_kinase_lookup.json']
    inputs += sorted((HERE / 'uniprot').glob('*.fasta'))
    write_artifact(HERE / 'X0_I2.json', out, inputs)
    summary = {name: {'ratio': rep.get('ratio'), 'pass': rep.get('pass_capability'),
                      'n_pairs': rep.get('n_pairs')} for name, rep in reps.items()}
    summary['edit_descriptor'] = {'median_pair_distance': edit_descriptor['median_pair_distance'],
                                  'gate': edit_descriptor['admission_gate_status']}
    summary['gate'] = gate
    print(json.dumps(summary, indent=1))
    return 0


def gate_eligible(reps):
    return [n for n in ['composition', 'global_esm', 'pair_centered_local_window',
                        'esm_local_window', 'klifs_pocket'] if n in reps]


if __name__ == '__main__':
    raise SystemExit(main())
