"""S2 — deterministic 3D mechanistic teacher from raw holo coordinates.

Frozen by TEACHER_CONTRACT.md.  No affinity label, no identity of any kind.
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np

CHANNELS = ["hbond_directional", "electrostatic_signed", "hydrophobic_burial",
            "aromatic_orientation", "steric_overlap", "pocket_burial"]

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "F": 1.47,
       "CL": 1.75, "BR": 1.85, "I": 1.98}
POLAR = {"N", "O"}
POS_ATOMS = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ"),
             ("HIS", "ND1"), ("HIS", "NE2")}
NEG_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
AROM_RES = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
            "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}
SKIP_LIG = {"HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "ACT", "CL",
            "NA", "K", "MG", "ZN", "CA", "MN", "FE", "CU", "NI", "CO", "CD", "HG",
            "IOD", "BR", "NO3", "TRS", "DMS", "FMT", "ACY", "IMD", "BME", "MES",
            "EPE", "TLA", "CIT", "SIN", "MLI", "AZI", "CAC", "PG4", "1PE", "P6G",
            "SCN", "NH4", "LI", "CS", "RB", "SR", "BA", "PB", "AU", "AG", "PT"}


def fcut(d, a, b):
    out = np.zeros_like(d)
    out[d <= a] = 1.0
    m = (d > a) & (d < b)
    out[m] = 0.5 * (1 + np.cos(np.pi * (d[m] - a) / (b - a)))
    return out


def _load(path):
    import gemmi
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
            doc = gemmi.cif.read_string(f.read())
    else:
        doc = gemmi.cif.read(path)
    st = gemmi.make_structure_from_block(doc.sole_block())
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    return st


def extract(path, min_lig_atoms=8, max_lig_atoms=120):
    """Return (protein atoms, ligand atoms) for the largest drug-like ligand."""
    import gemmi
    st = _load(path)
    model = st[0]
    prot, ligs = [], {}
    for ch in model:
        for res in ch:
            nm = res.name.upper()
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and info.is_amino_acid():
                for at in res:
                    if at.element.name.upper() != "H":
                        prot.append((nm, at.name.upper(), at.element.name.upper(),
                                     np.array([at.pos.x, at.pos.y, at.pos.z]),
                                     f"{ch.name}:{res.seqid.num}"))
            elif (nm not in SKIP_LIG and not res.is_water()
                  and not (info is not None and info.is_nucleic_acid())):
                key = (ch.name, res.seqid.num, nm)
                for at in res:
                    if at.element.name.upper() != "H":
                        ligs.setdefault(key, []).append(
                            (at.element.name.upper(),
                             np.array([at.pos.x, at.pos.y, at.pos.z])))
    if not prot or not ligs:
        return None, None, None
    key = max(ligs, key=lambda k: len(ligs[k]))
    lig = ligs[key]
    if not (min_lig_atoms <= len(lig) <= max_lig_atoms):
        return None, None, None
    return prot, lig, key[2]


def teacher(prot, lig):
    """Six named channels: pair-local contributions and complex aggregates."""
    P = np.stack([p[3] for p in prot])
    L = np.stack([a[1] for a in lig])
    pe = np.array([p[2] for p in prot])
    le = np.array([a[0] for a in lig])
    pres = np.array([p[0] for p in prot])
    pnam = np.array([p[1] for p in prot])
    pkey = np.array([p[4] for p in prot])

    # Restrict to whole protein residues with any atom within 10 A of the ligand.
    # Every channel is defined at <= 8 A and bonded neighbours are <= 1.8 A away,
    # so this is EXACTLY equivalent to using the full protein while making the
    # O(nP^2) neighbour computation tractable on large homo-oligomers.
    d0 = np.linalg.norm(L[:, None, :] - P[None, :, :], axis=2).min(axis=0)
    keep_res = {k for k in np.unique(pkey) if d0[pkey == k].min() < 10.0}
    if not keep_res:
        return None
    sel = np.array([k in keep_res for k in pkey])
    P, pe, pres, pnam, pkey = P[sel], pe[sel], pres[sel], pnam[sel], pkey[sel]

    D = np.linalg.norm(L[:, None, :] - P[None, :, :], axis=2)   # (nL, nP_site)
    near = D < 8.0
    if not near.any():
        return None

    # ---- 1 directional hydrogen bond -----------------------------------
    pol_p = np.isin(pe, list(POLAR))
    pol_l = np.isin(le, list(POLAR))
    # heavy-neighbour centroid of each protein polar atom, for the angle term
    nbr = np.zeros_like(P)
    has_nbr = np.zeros(len(P), bool)
    PD = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    for i in np.where(pol_p)[0]:
        j = np.where((PD[i] > 0.1) & (PD[i] < 1.8))[0]
        if len(j):
            nbr[i] = P[j].mean(0)
            has_nbr[i] = True
    v1 = P[None, :, :] - L[:, None, :]
    v2 = P[None, :, :] - nbr[None, :, :]
    n1 = np.linalg.norm(v1, axis=2) + 1e-9
    n2 = np.linalg.norm(v2, axis=2) + 1e-9
    cos = np.einsum("ijk,ijk->ij", v1, v2) / (n1 * n2)
    ang = np.where(has_nbr[None, :], cos ** 2, 0.0)
    hb = fcut(D, 3.2, 4.0) * (pol_l[:, None] & pol_p[None, :]) * ang

    # ---- 2 signed electrostatics ---------------------------------------
    qp = np.zeros(len(P))
    for i, (rn, an) in enumerate(zip(pres, pnam)):
        if (rn, an) in POS_ATOMS:
            qp[i] = 1.0
        elif (rn, an) in NEG_ATOMS or an == "OXT":
            qp[i] = -1.0
    ql = np.zeros(len(L))
    for i, e in enumerate(le):
        if e == "N":
            ql[i] = 1.0
        elif e == "O":
            d_on = np.linalg.norm(L[i] - L, axis=1)
            if np.sum((d_on > 0.1) & (d_on < 1.7)) <= 1:
                ql[i] = -1.0
    es = fcut(D, 4.0, 6.0) * (ql[:, None] * qp[None, :])

    # ---- 3 hydrophobic burial ------------------------------------------
    apolar_p = np.isin(pe, ["C", "S"])
    for i in np.where(apolar_p)[0]:
        j = np.where((PD[i] > 0.1) & (PD[i] < 1.8))[0]
        if len(j) and np.isin(pe[j], list(POLAR)).any():
            apolar_p[i] = False
    apolar_l = np.isin(le, ["C", "S"])
    LD = np.linalg.norm(L[:, None, :] - L[None, :, :], axis=2)
    for i in np.where(apolar_l)[0]:
        j = np.where((LD[i] > 0.1) & (LD[i] < 1.8))[0]
        if len(j) and np.isin(le[j], list(POLAR)).any():
            apolar_l[i] = False
    ph = fcut(D, 4.0, 5.0) * (apolar_l[:, None] & apolar_p[None, :])

    # ---- 4 aromatic orientation ----------------------------------------
    ar = np.zeros_like(D)
    rings = []
    for k in np.unique(pkey):
        m = pkey == k
        rn = pres[m][0]
        if rn in AROM_RES:
            sel = m & np.isin(pnam, AROM_RES[rn])
            if sel.sum() >= 5:
                c = P[sel].mean(0)
                u, s, vt = np.linalg.svd(P[sel] - c)
                rings.append((c, vt[2], np.where(sel)[0]))
    if rings and len(L) >= 6:
        cl = L.mean(0)
        u, s, vt = np.linalg.svd(L - cl)
        planar = s[2] / (s[0] + 1e-9) < 0.25
        if planar:
            nl = vt[2]
            for c, n, idx in rings:
                dr = np.linalg.norm(cl - c)
                w = float(fcut(np.array([dr]), 4.5, 6.5)[0])
                if w > 0:
                    a = abs(float(np.dot(nl, n))) ** 2
                    ar[:, idx] += w * a / max(len(idx), 1)

    # ---- 5 steric overlap ----------------------------------------------
    rp = np.array([VDW.get(e, 1.70) for e in pe])
    rl = np.array([VDW.get(e, 1.70) for e in le])
    ov = np.maximum(0.0, rl[:, None] + rp[None, :] - 0.4 - D)

    # ---- 6 pocket burial -----------------------------------------------
    bur = (D < 6.0).astype(float)

    stack = np.stack([hb, es, ph, ar, ov, bur])          # (6, nL, nP)
    stack = stack * near[None, :, :]
    agg = stack.sum(axis=(1, 2)) / max(len(L), 1)
    # pair-local: aggregate protein atoms to residues
    res_ids = np.unique(pkey)
    ridx = {r: i for i, r in enumerate(res_ids)}
    pl = np.zeros((6, len(L), len(res_ids)))
    for j, k in enumerate(pkey):
        pl[:, :, ridx[k]] += stack[:, :, j]
    return {"aggregate": agg, "pair_local": pl, "residue_ids": res_ids,
            "n_ligand_atoms": len(L), "n_protein_atoms": len(P),
            "ligand_elements": le, "residue_names":
                np.array([pres[pkey == r][0] for r in res_ids])}


def invariance_audit(path, seed=0):
    prot, lig, comp = extract(path)
    if prot is None:
        return None
    base = teacher(prot, lig)
    if base is None:
        return None
    rng = np.random.default_rng(seed)
    out = {}
    # rotation + translation
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    t = rng.normal(size=3) * 10
    prot_r = [(a, b, c, Q @ d + t, e) for a, b, c, d, e in prot]
    lig_r = [(e, Q @ p + t) for e, p in lig]
    r = teacher(prot_r, lig_r)
    out["rototranslation_max_abs_diff"] = float(
        np.max(np.abs(base["aggregate"] - r["aggregate"])))
    # atom permutation
    pi = rng.permutation(len(prot))
    li = rng.permutation(len(lig))
    p2 = teacher([prot[i] for i in pi], [lig[i] for i in li])
    out["permutation_max_abs_diff"] = float(
        np.max(np.abs(base["aggregate"] - p2["aggregate"])))
    # determinism
    d2 = teacher(prot, lig)
    out["determinism_max_abs_diff"] = float(
        np.max(np.abs(base["aggregate"] - d2["aggregate"])))
    out["aggregate"] = {c: float(v) for c, v in zip(CHANNELS, base["aggregate"])}
    return out
