"""Stage X0 round-2 shared infrastructure (corrected instruments).

Design rules implemented here:
- Every random seed derives from SHA-256 (never Python hash()).
- Frozen preregistration SHA is a module-level constant; nothing rewrites it.
- Data loaders read the local governed caches only (downloads/ + uniprot/).
- Bootstrapping is component-cluster bootstrapping by default.
"""
from __future__ import annotations
import hashlib, json, re, warnings
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PREREG_SHA = '03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'
PREREG_FILE = HERE / 'STAGE_X0_PREREGISTRATION.md'
AAS = 'ARNDCQEGHILKMFPSTWYV'

# Canonical parent gene symbol for each parent label used in Duong-Ly S2.
PARENT_ALIASES = {
    'P38α/MAPK14': 'MAPK14',
    'P38α': 'MAPK14',
    'TIE2/TEK': 'TEK',
    'TIE2': 'TEK',
    'C-KIT': 'KIT',
    'C-MET': 'MET',
    'C-SRC': 'SRC',
    'PDGFRα': 'PDGFRA',
    'MEK1': 'MAP2K1',
    'CHK2': 'CHEK2',
}


def sha256_seed(*parts) -> int:
    """Deterministic cross-process seed from SHA-256 of stringified parts."""
    text = '|'.join(str(p) for p in parts)
    digest = hashlib.sha256(text.encode('utf-8')).digest()
    return int.from_bytes(digest[:8], 'big')


def stable_rng(*parts) -> np.random.Generator:
    return np.random.default_rng(sha256_seed(*parts))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def verify_prereg_frozen() -> bool:
    return sha256_file(PREREG_FILE) == PREREG_SHA


def git_commit() -> str:
    import subprocess
    try:
        out = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=str(HERE.parents[3]),
                             capture_output=True, text=True, timeout=30)
        return out.stdout.strip() or 'unavailable'
    except Exception:
        return 'unavailable'


def load_duongly():
    """Return (s1_info, s2_matrix, sequences) from the local governed caches."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Unknown extension is not supported and will be removed',
            category=UserWarning,
        )
        info = pd.read_excel(HERE / 'downloads/duongly_mmc2.xlsx', sheet_name='Table S1')
        matrix = pd.read_excel(HERE / 'downloads/duongly_mmc3.xlsx', sheet_name='Table S2')
    seqs = {}
    for p in (HERE / 'uniprot').glob('*.fasta'):
        text = p.read_text(encoding='utf-8')
        seq = ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))
        header = next((line for line in text.splitlines() if line.startswith('>')), '')
        seqs[p.stem] = {'sequence': seq, 'header': header}
    return info, matrix, seqs


def normalize_construct_name(s2_label: str) -> str:
    """Normalize an S2 row label (e.g. 'P38α(T106M)') to the S1 construct key form."""
    s = s2_label.strip()
    s = s.replace(' ', '')
    for alias, canon in PARENT_ALIASES.items():
        s = s.replace(alias, canon)
    if '(' in s and not s.endswith(')'):
        s = s + ')'
    return s


def normalize_parent_name(s2_label: str) -> str:
    s = s2_label.strip().replace(' ', '')
    for alias, canon in PARENT_ALIASES.items():
        if s.startswith(alias):
            s = canon + s[len(alias):]
            break
    return s.split('(')[0]


def parse_point_notation(s: str):
    """Parse a single point-mutation notation like 'E255K' -> (old, pos, new)."""
    m = re.match(r'^([A-Za-z])([0-9]+)([A-Za-z])$', s.strip())
    if m:
        return (m.group(1).upper(), int(m.group(2)), m.group(3).upper())
    return None


def parse_mutation_list(reported: str):
    """Parse a mutation list into structured entries.

    Returns list of dicts:
      point: {'kind':'point','old':..,'pos':..,'new':..,'reported':..}
      deletion: {'kind':'deletion','start':..,'end':..,'reported':..}
      insertion: {'kind':'insertion','start':..,'end':..,'reported':..}
      unknown: {'kind':'unknown','reported':..}
    """
    entries = []
    s = str(reported).strip()
    if s in ('-', 'nan', 'None', ''):
        return entries
    for part in re.split(r'[,;]+', s):
        part = part.strip()
        if not part:
            continue
        pt = parse_point_notation(part)
        if pt:
            entries.append({'kind': 'point', 'old': pt[0], 'pos': pt[1], 'new': pt[2],
                            'reported': part})
            continue
        m = re.match(r'^[dD]([0-9]+)-([0-9]+)$', part)
        if m:
            entries.append({'kind': 'deletion', 'start': int(m.group(1)), 'end': int(m.group(2)),
                            'reported': part})
            continue
        m = re.match(r'^[dD]([0-9]+)-([0-9]+)/([A-Za-z])([0-9]+)([A-Za-z])$', part)
        if m:
            entries.append({'kind': 'deletion', 'start': int(m.group(1)), 'end': int(m.group(2)),
                            'reported': part})
            entries.append({'kind': 'point', 'old': m.group(3).upper(), 'pos': int(m.group(4)),
                            'new': m.group(5).upper(), 'reported': part})
            continue
        m2 = re.search(r'([0-9]+)[-\u2013]([0-9]+)', part)
        if any(w in part.lower() for w in ('tandem', 'insertion', 'duplication')):
            if m2:
                entries.append({'kind': 'insertion', 'start': int(m2.group(1)),
                                'end': int(m2.group(2)), 'reported': part})
                continue
        entries.append({'kind': 'unknown', 'reported': part})
    return entries


def parse_construct_range(clone: str, canon_len: int):
    """Parse the S1 Clone column into a construct range on the canonical sequence.

    Returns dict with 'kind' ('full','range','parts','unresolved'), 'start','end',
    'parts' (list of [a,b] for non-contiguous), 'note'.
    """
    s = str(clone).strip().lower()
    if not s or s == 'nan':
        return {'kind': 'unresolved', 'note': 'clone field empty'}
    if 'full-length' in s or s == 'full length':
        return {'kind': 'full', 'start': 1, 'end': canon_len, 'note': 'full-length'}
    # 'aa 27-end'
    m = re.search(r'aa\s*([0-9]+)\s*-\s*end', s)
    if m and 'end)' not in s:
        return {'kind': 'range', 'start': int(m.group(1)), 'end': canon_len,
                'note': 'aa X-end'}
    if 'end)' in s and 'aa' in s:
        m = re.search(r'aa\s*([0-9]+)\s*-\s*([0-9]+)', s)
        if m:
            return {'kind': 'range', 'start': int(m.group(1)), 'end': int(m.group(2)),
                    'note': 'aa X-Y (end)'}
    # non-contiguous: 'cytoplasmic domain [669-745, 751-1210(end)'
    parts = re.findall(r'([0-9]+)\s*-\s*([0-9]+)', s)
    if '[' in s and len(parts) >= 1:
        return {'kind': 'parts', 'parts': [[int(a), int(b)] for a, b in parts],
                'note': 'non-contiguous construct parts'}
    if len(parts) >= 1:
        a, b = int(parts[0][0]), int(parts[0][1])
        if a == 1 and b == canon_len:
            return {'kind': 'full', 'start': 1, 'end': canon_len, 'note': 'full-length'}
        return {'kind': 'range', 'start': a, 'end': b, 'note': 'aa X-Y'}
    return {'kind': 'unresolved', 'note': f'unparsed clone field: {clone!r}'}


def construct_sequence(canon_seq: str, crange: dict):
    """Return the construct sequence for a parsed construct range, or None."""
    if crange['kind'] in ('full', 'range'):
        a, b = crange['start'], crange['end']
        if b > len(canon_seq):
            return None  # construct extends beyond canonical sequence
        return canon_seq[a - 1:b]
    if crange['kind'] == 'parts':
        out = ''
        for a, b in crange['parts']:
            if b > len(canon_seq):
                return None
            out += canon_seq[a - 1:b]
        return out
    return None


def map_canonical_to_construct(pos_canon: int, crange: dict):
    """Map a canonical coordinate to a 1-based construct coordinate; None if unmappable."""
    if crange['kind'] in ('full', 'range'):
        a, b = crange['start'], crange['end']
        if a <= pos_canon <= b:
            return pos_canon - a + 1
        return None
    if crange['kind'] == 'parts':
        offset = 0
        for a, b in crange['parts']:
            if a <= pos_canon <= b:
                return offset + (pos_canon - a + 1)
            offset += (b - a + 1)
        return None
    return None


def onehot_aa(aa):
    v = np.zeros(len(AAS), dtype=np.float32)
    if aa in AAS:
        v[AAS.index(aa)] = 1.0
    return v


def composition(seq):
    v = np.zeros(len(AAS), dtype=np.float32)
    for aa in seq:
        if aa in AAS:
            v[AAS.index(aa)] += 1
    n = max(len(seq), 1)
    return v / n


def cluster_bootstrap(values_by_cluster, n_draws=2000, seed=20260820,
                      statistic=np.median):
    """Component-cluster bootstrap over groups of values.

    values_by_cluster: list of arrays (one array of values per cluster).
    statistic: callable array -> scalar. Returns dict with estimate, ci_lo, ci_hi,
    n_clusters, n_values, seed, draws.
    """
    rng = np.random.default_rng(seed)
    stats = np.empty(n_draws)
    k = len(values_by_cluster)
    if k == 0:
        return {'estimate': float('nan'), 'ci_lo': float('nan'), 'ci_hi': float('nan'),
                'n_clusters': 0, 'n_values': 0, 'seed': seed, 'draws': n_draws,
                'note': 'no clusters'}
    for d in range(n_draws):
        idx = rng.integers(0, k, size=k)
        pooled = np.concatenate([values_by_cluster[i] for i in idx])
        stats[d] = statistic(pooled)
    est = statistic(np.concatenate(values_by_cluster))
    lo, hi = float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))
    return {'estimate': float(est), 'ci_lo': lo, 'ci_hi': hi, 'n_clusters': k,
            'n_values': int(sum(len(v) for v in values_by_cluster)), 'seed': seed,
            'draws': n_draws}


def write_artifact(path: Path, obj: dict, input_paths=None):
    """Write an instrument JSON artifact with provenance fields."""
    out = dict(obj)
    out.setdefault('preregistration_sha256', PREREG_SHA)
    out.setdefault('code_commit', git_commit())
    out.setdefault('input_sha256', {str(p.relative_to(HERE)): sha256_file(p)
                                    for p in (input_paths or []) if p.exists()})
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    return out
