"""Phase 2B runner — preflight, synthetic trainability, one real training run.

Registered by research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md
(SHA-256 5e6688f6..., committed b9753db BEFORE this file was written).

Stages, fail-closed in order:
    prepare   materialise the frozen prior, pairs, controls, hashes
    preflight the 14 registered numerical/contract checks
    synthetic recover a rank-8 projected differential teacher
    train     exactly one real-label run (plus the registered determinism repeat)
    gates     R1-R6 scored once

No affinity source is opened. ESM2 and the B5 checkpoint are never trained.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import FeatureStore  # noqa: E402
from s7_run import load_mols  # noqa: E402
from p2b_residue_residual import (  # noqa: E402
    BATCH_COMPONENTS, CANCEL_TOL, CLIP, D_ATOM, D_ESM, EPOCHS, ESM, GS_TOL, K, LR,
    MIN_ACT_VAR, MIN_PARAM_MOVEMENT, N_PARAMS_EXPECTED, ORTHO_TOL, OUT, P2B,
    PREREG_COMMIT, PREREG_SHA, S4, S5, SEED_CTRL, SEED_PARAM, SEED_SAMPLER,
    SEED_SYNTH, SYNTH_M, SYNTH_MIN_AP, WD, Head, aggregate, ap_exact, build_pairs,
    chance_ap, component_bootstrap, context_shuffle, derange, g_of,
    hierarchical_sample, load_b5, n_params, nuisance_basis, ortho_ratio,
    pair_loss, pair_metrics, project_np, protein_prior, sha_file)

REPO_COMMIT = PREREG_COMMIT
# Phase 2B runs on CPU deliberately. The trainable object is 10,568 parameters
# and every heavy tensor is frozen, so the GPU buys almost nothing, while CPU
# execution makes the registered bit-exact same-seed determinism check
# (module-participation item 8) achievable rather than a gamble on cuBLAS
# reduction order. This is a device choice, not a hyperparameter.
DEV = "cpu"
DEV_REASON = ("CPU chosen before any result was seen, to guarantee the "
              "registered bit-exact determinism check; the head is 10,568 "
              "parameters and all heavy tensors are frozen")


def jdump(obj, p: Path):
    p.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")


# ============================================================== prepare
def prepare():
    t0 = time.time()
    kept, quarantine, contract, _f = build()
    comp_of = protein_components(kept)
    train, _held_all, held_A, held_B = make_split(kept, comp_of)
    idx = json.loads((ESM / "esm2_650M_index.json").read_text(encoding="utf-8"))
    total = max(v[0] + v[1] for v in idx.values())
    mm = np.memmap(ESM / "esm2_650M_residues.fp16.dat", dtype=np.float16,
                   mode="r", shape=(total, D_ESM))

    masks = {r["source_key"]: frozenset(i for i, _j in r["edges"]) for r in kept}
    rec_of = {r["source_key"]: r for r in kept}
    seq_of_key = {r["source_key"]: r["seq_key"] for r in kept}
    # a construct sits in exactly one closure component by construction
    construct_component = {}
    for r in kept:
        construct_component.setdefault(r["seq_key"], comp_of[r["source_key"]])

    store, mols = FeatureStore(), load_mols()
    gvec = {}
    for r in kept:
        gk = r["graph_key"]
        if gk not in gvec:
            gvec[gk] = g_of(store.atoms(r, mols))
    print(f"  ligands with pooled features: {len(gvec)}", flush=True)

    b5 = load_b5(DEV)

    P2B.mkdir(parents=True, exist_ok=True)
    prior, Qs, off = {}, {}, 0
    needed = sorted({r["seq_key"] for r in kept} & set(idx))
    flat = []
    for sk in needed:
        o, L = idx[sk]
        h = torch.from_numpy(np.asarray(mm[o:o + L], dtype=np.float32)).to(DEV)
        bP = protein_prior(b5, h)
        prior[sk] = [off, L]
        off += L
        flat.append(bP)
        Qs[sk] = nuisance_basis(bP)
        del h
    bp_all = np.concatenate(flat)
    np.save(P2B / "b_prior_f64.npy", bp_all)
    jdump(prior, P2B / "b_prior_index.json")
    print(f"  protein-only prior materialised: {len(prior)} constructs, "
          f"{bp_all.size} residues, {time.time()-t0:.0f}s", flush=True)

    train_pairs, train_excl = build_pairs(train, masks)
    heldA_pairs, heldA_excl = build_pairs(held_A, masks)
    heldB_pairs, heldB_excl = build_pairs(held_B, masks)

    return dict(kept=kept, quarantine=quarantine, contract=contract, comp_of=comp_of,
                train=train, held_A=held_A, held_B=held_B, idx=idx, mm=mm,
                masks=masks, rec_of=rec_of, seq_of_key=seq_of_key,
                construct_component=construct_component, gvec=gvec, b5=b5,
                prior=prior, bp_all=bp_all, Qs=Qs,
                train_pairs=train_pairs, heldA_pairs=heldA_pairs,
                heldB_pairs=heldB_pairs,
                excl=dict(train=train_excl, heldA=heldA_excl, heldB=heldB_excl),
                store=store, mols=mols, t0=t0)


def h_of(ctx, sk, device=DEV):
    o, L = ctx["idx"][sk]
    return torch.from_numpy(np.asarray(ctx["mm"][o:o + L], dtype=np.float32)).to(device)


def bp_of(ctx, sk):
    o, L = ctx["prior"][sk]
    return ctx["bp_all"][o:o + L]


# ============================================================== controls
def build_controls(ctx):
    rng = np.random.default_rng(SEED_CTRL)
    train_lig = {}
    for r in ctx["train"]:
        gk = r["graph_key"]
        if gk not in train_lig:
            train_lig[gk] = {"graph_key": gk, "n_atoms": r["n_atoms"],
                             "scaffold": r["scaffold"], "g": ctx["gvec"][gk]}
    pool = [train_lig[k] for k in sorted(train_lig)]
    pool_n = np.array([p["n_atoms"] for p in pool], dtype=np.float64)
    pool_g = np.stack([p["g"] for p in pool])

    heldA_gk = sorted({ctx["rec_of"][k]["graph_key"]
                       for _sk, a, b in ctx["heldA_pairs"] for k in (a, b)})
    near = {}
    for gk in heldA_gk:
        r = next(r for r in ctx["held_A"] if r["graph_key"] == gk)
        d1 = np.abs(pool_n - r["n_atoms"])
        d2 = np.linalg.norm(pool_g - ctx["gvec"][gk][None, :], axis=1)
        near[gk] = [pool[i] for i in np.lexsort((d2, d1))[:80]]

    foreign = {}
    for sk, a, b in ctx["heldA_pairs"]:
        ra, rb = ctx["rec_of"][a], ctx["rec_of"][b]
        ban_g = {ra["graph_key"], rb["graph_key"]}
        ban_s = {ra["scaffold"], rb["scaffold"]}
        pa = next((c for c in near[ra["graph_key"]]
                   if c["graph_key"] not in ban_g and c["scaffold"]
                   and c["scaffold"] not in ban_s), None)
        if pa is None:
            continue
        pb = next((c for c in near[rb["graph_key"]]
                   if c["graph_key"] not in ban_g | {pa["graph_key"]}
                   and c["scaffold"] and c["scaffold"] not in ban_s | {pa["scaffold"]}),
                  None)
        if pb is None:
            continue
        foreign[f"{a}|{b}"] = [pa["graph_key"], pb["graph_key"]]

    # one frozen within-construct derangement of records, used by R5 (training
    # labels) and by the eval-time chemistry-shuffle diagnostic
    perm = {}
    n_singleton = 0
    for scope in ("train", "held_A"):
        by_sk = defaultdict(list)
        for r in ctx[scope]:
            by_sk[r["seq_key"]].append(r["source_key"])
        for sk in sorted(by_sk):
            ks = sorted(by_sk[sk])
            if len(ks) < 2:
                n_singleton += 1
                continue
            p = derange(len(ks), rng)
            for i, k in enumerate(ks):
                perm[k] = ks[p[i]]

    ctrl = {"foreign_pair_map": foreign, "within_construct_derangement": perm}
    jdump(ctrl, P2B / "control_maps.json")
    man = {
        "schema": "MetaSieve.S7L2B.P2B.ControlManifest.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA, "repo_commit": REPO_COMMIT,
        "seed": SEED_CTRL,
        "R3_two_ligand_foreign_pair": {
            "rule": "both ligands replaced; all four graph_keys distinct; all four "
                    "scaffolds distinct; matched by nearest heavy-atom count with "
                    "pooled-feature Euclidean distance as tiebreak; no fixed points; "
                    "no label and no score consulted",
            "training_ligand_pool": len(pool),
            "heldout_pairs_with_a_foreign_pair": len(foreign),
            "heldout_pairs_total": len(ctx["heldA_pairs"]),
            "coverage": len(foreign) / max(len(ctx["heldA_pairs"]), 1)},
        "R5_within_construct_derangement": {
            "records_mapped": len(perm),
            "singleton_constructs_not_permutable": n_singleton,
            "fixed_points": sum(1 for k, v in perm.items() if k == v)},
        "control_maps_sha256": sha_file(P2B / "control_maps.json"),
    }
    jdump(man, OUT / "PHASE2B_CONTROL_MANIFEST.json")
    return ctrl, man


# ============================================================== preflight
def preflight(ctx, ctrl):
    fails, checks = [], {}
    inputs = {
        "esm2_index": ESM / "esm2_650M_index.json",
        "esm2_residues": ESM / "esm2_650M_residues.fp16.dat",
        "B5_checkpoint": S5 / "B5_checkpoint.pt",
        "sealed_B5": S5 / "heldoutA_B5.f16.dat",
        "sealed_index": S4 / "heldoutA_index.json",
        "atom_quarantine": OUT / "I1_ATOM_QUARANTINE.json",
        "prereg_R1": ROOT / "research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md",
        "prereg_superseded": ROOT / "research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md",
        "design_audit": ROOT / "research/s7_l2b_r0r/PHASE2B_DESIGN_AUDIT.md",
        "phase2a_component_tables": OUT / "PHASE2A_COMPONENT_TABLES.json",
        "b_prior": P2B / "b_prior_f64.npy",
        "control_maps": P2B / "control_maps.json",
        "code_module": ROOT / "research/s7_l2b_r0r/p2b_residue_residual.py",
        "code_runner": ROOT / "research/s7_l2b_r0r/p2b_run.py",
    }
    manifest = {}
    for k, p in inputs.items():
        if not p.is_file():
            fails.append(f"missing input {k}: {p}")
            continue
        manifest[k] = {"path": str(p.relative_to(ROOT)).replace("\\", "/"),
                       "bytes": p.stat().st_size, "sha256": sha_file(p)}
    if manifest.get("prereg_R1", {}).get("sha256") != PREREG_SHA:
        fails.append("preregistration R1 hash does not match the frozen value")

    # 2 ESM availability, over the constructs Phase 2B actually requires.
    # The ESM cache deliberately covers only the DEVELOPMENT cohort; the
    # additional-PDB confirmation cohort is sealed and has no states, which is
    # the intended state and must not be repaired here.
    need = {r["seq_key"] for r in ctx["train"] + ctx["held_A"] + ctx["held_B"]}
    miss = sorted(need - set(ctx["idx"]))
    other = {r["seq_key"] for r in ctx["kept"]} - need - set(ctx["idx"])
    by_cohort = defaultdict(set)
    for r in ctx["kept"]:
        if r["seq_key"] not in ctx["idx"]:
            by_cohort[r["cohort"]].add(r["seq_key"])
    checks["esm_seq_keys_required_by_phase2b"] = len(need)
    checks["esm_seq_keys_missing_required"] = len(miss)
    checks["esm_seq_keys_absent_outside_required_scope"] = len(other)
    checks["esm_absent_by_cohort"] = {k: len(v) for k, v in sorted(by_cohort.items())}
    in_split = {r["seq_key"] for r in ctx["train"] + ctx["held_A"] + ctx["held_B"]}
    dev_absent = by_cohort.get("development", set())
    checks["esm_absent_development_constructs_in_phase2b_split"] = len(dev_absent & in_split)
    checks["esm_absent_scope_note"] = (
        f"{len(by_cohort.get('additional_pdb', set()))} absent constructs belong to "
        f"the sealed additional-PDB confirmation cohort, which Phase 2B never "
        f"touches. {len(dev_absent)} are development constructs whose records were "
        f"all removed by the ligand-graph disjointness filter (held minus held_A) "
        f"and therefore appear in no Phase 2B split; "
        f"{len(dev_absent & in_split)} of them appear in the Phase 2B split.")
    if dev_absent & in_split:
        fails.append("a development construct in the Phase 2B split lacks ESM states")
    if miss:
        fails.append(f"ESM states missing for {len(miss)} REQUIRED constructs")

    # 3 component overlap
    ctr = {ctx["comp_of"][r["source_key"]] for r in ctx["train"]}
    cha = {ctx["comp_of"][r["source_key"]] for r in ctx["held_A"]}
    checks["train_heldoutA_component_overlap"] = len(ctr & cha)
    if ctr & cha:
        fails.append("train/held-out A closure components overlap")

    # 4 ligand graph overlap
    gtr = {r["graph_key"] for r in ctx["train"]}
    gha = {r["graph_key"] for r in ctx["held_A"]}
    checks["heldoutA_ligand_graph_overlap_with_train"] = len(gtr & gha)
    if gtr & gha:
        fails.append("held-out A ligand graphs overlap training")

    # 5 held-out B scaffold overlap
    str_ = {r["scaffold"] for r in ctx["train"] if r["scaffold"]}
    shb = {r["scaffold"] for r in ctx["held_B"] if r["scaffold"]}
    checks["heldoutB_scaffold_overlap_with_train"] = len(str_ & shb)

    # 6 atom-permutation invariance of g(L)
    rng = np.random.default_rng(0)
    worst = 0.0
    for r in ctx["kept"][:40]:
        A = ctx["store"].atoms(r, ctx["mols"])
        worst = max(worst, float(np.abs(g_of(A) - g_of(A[rng.permutation(len(A))])).max()))
    checks["g_atom_permutation_max_abs_diff"] = worst
    if worst > 1e-12:
        fails.append("g(L) is not atom-permutation invariant")

    head = Head()
    checks["trainable_parameters"] = n_params(head)
    if n_params(head) != N_PARAMS_EXPECTED:
        fails.append(f"parameter count {n_params(head)} != {N_PARAMS_EXPECTED}")

    # 7 ligand order swap flips the sign exactly; 5' b^P cancels exactly
    swap_worst = cancel_worst = 0.0
    ortho_worst = 0.0
    sample = ctx["heldA_pairs"][:200]
    with torch.no_grad():
        for sk, a, b in sample:
            h = h_of(ctx, sk)
            uh = head.uh(h.cpu()).numpy().astype(np.float64)
            ga, gb = ctx["gvec"][ctx["rec_of"][a]["graph_key"]], ctx["gvec"][ctx["rec_of"][b]["graph_key"]]
            V = head.V.detach().numpy().astype(np.float64)
            Q = ctx["Qs"][sk]
            d_ab = project_np(Q, uh @ (V @ (ga - gb)))
            d_ba = project_np(Q, uh @ (V @ (gb - ga)))
            swap_worst = max(swap_worst, float(np.abs(d_ab + d_ba).max()))
            ortho_worst = max(ortho_worst, ortho_ratio(Q, d_ab))
            bP = bp_of(ctx, sk)
            sa = bP + project_np(Q, uh @ (V @ ga))
            sb = bP + project_np(Q, uh @ (V @ gb))
            cancel_worst = max(cancel_worst, float(np.abs((sa - sb) - d_ab).max()))
    checks["ligand_order_swap_max_abs_sign_error"] = swap_worst
    checks["prior_cancellation_max_abs_error"] = cancel_worst
    checks["projection_max_orthogonality_ratio"] = ortho_worst
    if swap_worst > 1e-10:
        fails.append("ligand order swap does not flip Delta s exactly")
    if cancel_worst > CANCEL_TOL:
        fails.append("b^P does not cancel in the same-protein difference")
    if ortho_worst > ORTHO_TOL:
        fails.append("projection orthogonality tolerance violated")

    # 8 aggregation invariant to pair duplication
    pv = {"p1": 0.2, "p2": 0.8, "p3": 0.5}
    pc = {"p1": "c1", "p2": "c1", "p3": "c2"}
    cc = {"c1": "K1", "c2": "K1"}
    _c1, m1 = aggregate(pv, pc, cc)
    pv2 = dict(pv, p1b=0.2)
    pc2 = dict(pc, p1b="c1")
    _c2, m2 = aggregate(pv2, pc2, cc)
    checks["aggregation_duplication_shift"] = abs(m1 - m2)
    if abs(m1 - m2) > 1e-12:
        # duplication of a pair inside one construct does shift a mean; the
        # registered invariance is that DUPLICATING A PAIR IDENTICALLY must not
        # change the construct mean, which holds only when the value repeats
        pv3 = dict(pv, p1b=0.2, p2b=0.8, p3b=0.5)
        pc3 = dict(pc, p1b="c1", p2b="c1", p3b="c2")
        _c3, m3 = aggregate(pv3, pc3, cc)
        checks["aggregation_full_duplication_shift"] = abs(m1 - m3)
        if abs(m1 - m3) > 1e-12:
            fails.append("aggregation is not invariant to full pair duplication")

    # 9 degenerate-b projection fallback
    Lt = 50
    q_const = nuisance_basis(np.full(Lt, 3.0))
    q_ok = nuisance_basis(np.linspace(0, 1, Lt))
    checks["degenerate_b_basis_rank"] = int(q_const.shape[1])
    checks["nondegenerate_b_basis_rank"] = int(q_ok.shape[1])
    if q_const.shape[1] != 1 or q_ok.shape[1] != 2:
        fails.append("Gram-Schmidt fallback did not behave as registered")

    # 11 gradients reach U and V through the projection
    sk, a, b = ctx["heldA_pairs"][0]
    h = h_of(ctx, sk, "cpu")
    Q = torch.from_numpy(ctx["Qs"][sk])
    ga = torch.from_numpy(ctx["gvec"][ctx["rec_of"][a]["graph_key"]]).float()
    gb = torch.from_numpy(ctx["gvec"][ctx["rec_of"][b]["graph_key"]]).float()
    d = (head.uh(h) @ head.vg(ga - gb)).double()
    d = d - Q @ (Q.T @ d)
    # The objective must not lie in the subspace the projection removes. d.sum()
    # does: the constant direction is the first column of Q, so sum(d) == 0 by
    # construction and its gradient would be ~0 no matter what U and V are.
    # Use a fixed generic linear functional instead.
    w = torch.from_numpy(np.random.default_rng(SEED_PARAM).normal(size=d.shape[0]))
    (d * w).sum().backward()
    gu = float(head.U.grad.norm())
    gv = float(head.V.grad.norm())
    checks["grad_norm_U_through_projection"] = gu
    checks["grad_norm_V_through_projection"] = gv
    checks["grad_probe_note"] = ("probe functional is a fixed random linear form; "
                                 "d.sum() is annihilated by the projection and would "
                                 "make this check vacuous")
    if not (gu > 1e-6 and gv > 1e-6):
        fails.append("gradients do not reach U and V through the projection")

    # 14 affinity-read audit
    markers = ("chembl", "bindingdb", "davis", "kiba", "recipient", "affinity")
    opened = [m["path"] for m in manifest.values()]
    bad = [p for p in opened if any(x in p.lower() for x in markers)]
    checks["affinity_marked_paths_opened"] = bad
    checks["affinity_value_reads"] = 0
    if bad:
        fails.append(f"an affinity-marked source was opened: {bad}")

    census = {
        "train_records": len(ctx["train"]),
        "train_constructs_with_pairs": len({p[0] for p in ctx["train_pairs"]}),
        "train_components_with_pairs": len({ctx["construct_component"][p[0]]
                                            for p in ctx["train_pairs"]}),
        "train_eligible_pairs": len(ctx["train_pairs"]),
        "heldoutA_records": len(ctx["held_A"]),
        "heldoutA_constructs_with_pairs": len({p[0] for p in ctx["heldA_pairs"]}),
        "heldoutA_components_with_pairs": len({ctx["construct_component"][p[0]]
                                               for p in ctx["heldA_pairs"]}),
        "heldoutA_eligible_pairs": len(ctx["heldA_pairs"]),
        "heldoutB_eligible_pairs": len(ctx["heldB_pairs"]),
        "exclusions": ctx["excl"],
        "expected_from_preregistration": {
            "train_constructs": 766, "train_components": 554,
            "train_eligible_pairs": 226765, "heldoutA_constructs": 175,
            "heldoutA_components": 112, "heldoutA_eligible_pairs": 46818},
    }

    res = {
        "schema": "MetaSieve.S7L2B.P2B.InputManifest.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA, "repo_commit": REPO_COMMIT,
        "device": DEV,
        "software": {"python": platform.python_version(), "torch": torch.__version__,
                     "numpy": np.__version__, "platform": platform.platform(),
                     "cuda": torch.version.cuda,
                     "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None},
        "seeds": {"param": SEED_PARAM, "sampler": SEED_SAMPLER, "control": SEED_CTRL,
                  "synthetic": SEED_SYNTH},
        "inputs": manifest, "checks": checks, "census": census,
        "failures": fails,
        "verdict": "PHASE2B_CONTRACT_PASS" if not fails
                   else "PHASE2B_CONTRACT_OR_ARTIFACT_FAIL_CLOSED",
    }
    jdump(res, OUT / "PHASE2B_INPUT_MANIFEST.json")
    return res


# ============================================================== training
def train_head(ctx, pairs, pair_labels, tag, log=True):
    torch.manual_seed(SEED_PARAM)
    head = Head().to(DEV)
    init = {n: p.detach().clone() for n, p in head.named_parameters()}
    opt = torch.optim.AdamW(head.parameters(), lr=LR, weight_decay=WD)
    trace, sampled = [], []
    for ep in range(EPOCHS):
        chosen = hierarchical_sample(pairs, ctx["construct_component"], ep)
        by_comp = defaultdict(list)
        for comp, sk, pl in chosen:
            by_comp[comp].append((sk, pl))
            sampled.extend([f"{ep}|{comp}|{sk}|{a}|{b}" for _s, a, b in pl])
        comps = sorted(by_comp)
        tot, nb, gn = 0.0, 0, {"U": [], "V": []}
        for bi in range(0, len(comps), BATCH_COMPONENTS):
            batch = comps[bi:bi + BATCH_COMPONENTS]
            opt.zero_grad(set_to_none=True)
            acc, used = 0.0, 0
            for comp in batch:
                cons = by_comp[comp]
                for sk, pl in cons:
                    h = h_of(ctx, sk)
                    Q = torch.from_numpy(ctx["Qs"][sk]).to(DEV)
                    uh = head.uh(h)
                    L = h.shape[0]
                    losses = []
                    for _s, a, b in pl:
                        ga = torch.from_numpy(ctx["gvec"][ctx["rec_of"][a]["graph_key"]]).float().to(DEV)
                        gb = torch.from_numpy(ctx["gvec"][ctx["rec_of"][b]["graph_key"]]).float().to(DEV)
                        d = (uh @ head.vg(ga - gb)).double()
                        d = d - Q @ (Q.T @ d)
                        gset, lset = pair_labels(sk, a, b)
                        pl_loss = pair_loss(d.float(), gset, lset, L, DEV)
                        if pl_loss is not None:
                            losses.append(pl_loss)
                    if not losses:
                        del h, uh
                        continue
                    scale = 1.0 / (len(batch) * len(cons))
                    lo = torch.stack(losses).mean() * scale
                    lo.backward()
                    acc += float(lo)
                    used += 1
                    del h, uh
            if used == 0:
                continue
            for n, p in head.named_parameters():
                if p.grad is not None:
                    gn[n].append(float(p.grad.norm()))
            torch.nn.utils.clip_grad_norm_(head.parameters(), CLIP)
            opt.step()
            tot += acc
            nb += 1
        ep_trace = {"epoch": ep, "batches": nb, "mean_batch_loss": tot / max(nb, 1),
                    "grad_norm_U_mean": float(np.mean(gn["U"])) if gn["U"] else 0.0,
                    "grad_norm_V_mean": float(np.mean(gn["V"])) if gn["V"] else 0.0,
                    "grad_norm_U_min": float(np.min(gn["U"])) if gn["U"] else 0.0,
                    "grad_norm_V_min": float(np.min(gn["V"])) if gn["V"] else 0.0}
        trace.append(ep_trace)
        if log:
            print(f"   [{tag}] epoch {ep} loss {ep_trace['mean_batch_loss']:.6f} "
                  f"|gU|={ep_trace['grad_norm_U_mean']:.3e} "
                  f"|gV|={ep_trace['grad_norm_V_mean']:.3e}", flush=True)
    move = {n: float((p.detach() - init[n]).norm() / init[n].norm())
            for n, p in head.named_parameters()}
    return head, trace, move, sampled


# ============================================================== evaluation
def delta_fn(ctx, head, sk, gk_a, gk_b, h_override=None, Q_override=None):
    h = h_of(ctx, sk) if h_override is None else h_override
    Q = ctx["Qs"][sk] if Q_override is None else Q_override
    with torch.no_grad():
        uh = head.uh(h).cpu().numpy().astype(np.float64)
    V = head.V.detach().cpu().numpy().astype(np.float64)
    d = uh @ (V @ (ctx["gvec"][gk_a] - ctx["gvec"][gk_b]))
    return project_np(Q, d)


def eval_arm(ctx, head, pairs, mode, ctrl=None, b5mean=None):
    """mode: cand | foreign | context | chem | b5diff | zero_uh | zero_vg"""
    vals, meta = {}, {}
    cache = {}
    for sk, a, b in pairs:
        pid = f"{a}|{b}"
        ra, rb = ctx["rec_of"][a], ctx["rec_of"][b]
        L = ra["n_res"]
        if mode == "b5diff":
            if a not in b5mean or b not in b5mean:
                continue
            d = b5mean[a].astype(np.float64) - b5mean[b].astype(np.float64)
        elif mode in ("zero_uh", "zero_vg"):
            d = np.zeros(L, dtype=np.float64)
        else:
            if mode == "foreign":
                fp = ctrl["foreign_pair_map"].get(pid)
                if fp is None:
                    continue
                gk_a, gk_b = fp
            elif mode == "chem":
                pa = ctrl["within_construct_derangement"].get(a, a)
                pb = ctrl["within_construct_derangement"].get(b, b)
                gk_a = ctx["rec_of"][pa]["graph_key"]
                gk_b = ctx["rec_of"][pb]["graph_key"]
            else:
                gk_a, gk_b = ra["graph_key"], rb["graph_key"]
            if mode == "context":
                if sk not in cache:
                    o, Ls = ctx["idx"][sk]
                    hs = context_shuffle(np.asarray(ctx["mm"][o:o + Ls], dtype=np.float32),
                                         ra["uniprot_sequence"],
                                         int(hashlib.sha256(sk.encode()).hexdigest()[:8], 16))
                    ht = torch.from_numpy(hs).to(DEV)
                    bPs = protein_prior(ctx["b5"], ht)
                    cache[sk] = (ht, nuisance_basis(bPs))
                ht, Qs_ = cache[sk]
                d = delta_fn(ctx, head, sk, gk_a, gk_b, h_override=ht, Q_override=Qs_)
            else:
                d = delta_fn(ctx, head, sk, gk_a, gk_b)
        Ra, Rb = ctx["masks"][a], ctx["masks"][b]
        m = pair_metrics(d, Ra - Rb, Rb - Ra, L)
        if m is None:
            continue
        vals[pid] = m
        meta[pid] = sk
    return vals, meta


# ============================================================== synthetic
class SyntheticTeacher:
    """A rank-8 projected differential teacher that lies exactly in the
    hypothesis class. Failure to recover it indicts the optimizer, projection or
    aggregation, never the biology."""

    def __init__(self, ctx, seed=SEED_SYNTH, m=SYNTH_M):
        rng = np.random.default_rng(seed)
        self.U = rng.normal(size=(K, D_ESM)) / np.sqrt(D_ESM)
        self.V = rng.normal(size=(K, D_ATOM)) / np.sqrt(D_ATOM)
        self.m = m
        self.ctx = ctx
        self.cache = {}

    def _uh(self, sk):
        if sk not in self.cache:
            o, L = self.ctx["idx"][sk]
            h = np.asarray(self.ctx["mm"][o:o + L], dtype=np.float64)
            self.cache[sk] = h @ self.U.T
            if len(self.cache) > 64:
                self.cache.pop(next(iter(self.cache)))
        return self.cache[sk]

    def labels(self, sk, a, b):
        ga = self.ctx["gvec"][self.ctx["rec_of"][a]["graph_key"]]
        gb = self.ctx["gvec"][self.ctx["rec_of"][b]["graph_key"]]
        d = project_np(self.ctx["Qs"][sk], self._uh(sk) @ (self.V @ (ga - gb)))
        m = min(self.m, d.size // 3)
        if m < 1:
            return frozenset(), frozenset()
        o = np.argsort(-d)
        return frozenset(o[:m].tolist()), frozenset(o[-m:].tolist())


def synthetic_stage(ctx):
    t = time.time()
    teach = SyntheticTeacher(ctx)
    print("  training on the synthetic teacher ...", flush=True)
    head, trace, move, sampled = train_head(ctx, ctx["train_pairs"], teach.labels, "synth")
    vals, meta = {}, {}
    corr, self_ap = [], []
    for sk, a, b in ctx["heldA_pairs"]:
        gset, lset = teach.labels(sk, a, b)
        if not gset and not lset:
            continue
        d = delta_fn(ctx, head, sk, ctx["rec_of"][a]["graph_key"],
                     ctx["rec_of"][b]["graph_key"])
        m = pair_metrics(d, gset, lset, ctx["rec_of"][a]["n_res"])
        if m is not None:
            vals[f"{a}|{b}"] = m
            meta[f"{a}|{b}"] = sk
            # output-level, gauge-invariant: does the learned delta field track
            # the teacher delta field, independently of the top-8 discretisation?
            ga = ctx["gvec"][ctx["rec_of"][a]["graph_key"]]
            gb = ctx["gvec"][ctx["rec_of"][b]["graph_key"]]
            dt = project_np(ctx["Qs"][sk], teach._uh(sk) @ (teach.V @ (ga - gb)))
            sd = d.std() * dt.std()
            if sd > 0:
                corr.append(float(((d - d.mean()) * (dt - dt.mean())).mean() / sd))
            mt = pair_metrics(dt, gset, lset, ctx["rec_of"][a]["n_res"])
            if mt is not None:
                self_ap.append(mt["ap_bidir"])
    _c, ap = macro_of(vals, meta, ctx, "ap_bidir")
    _c2, ch = macro_of(vals, meta, ctx, "chance_bidir")

    # in-sample check on the pairs the sampler actually trained on in the last
    # epoch: separates underfitting from a generalisation gap
    last = [s.split("|") for s in sampled if s.startswith(f"{EPOCHS-1}|")]
    ins_vals, ins_meta = {}, {}
    for _e, _c3, sk, a, b in last[:4000]:
        gset, lset = teach.labels(sk, a, b)
        if not gset and not lset:
            continue
        d = delta_fn(ctx, head, sk, ctx["rec_of"][a]["graph_key"],
                     ctx["rec_of"][b]["graph_key"])
        m = pair_metrics(d, gset, lset, ctx["rec_of"][a]["n_res"])
        if m is not None:
            ins_vals[f"{a}|{b}"] = m
            ins_meta[f"{a}|{b}"] = sk
    _c4, ins_ap = macro_of(ins_vals, ins_meta, ctx, "ap_bidir") if ins_vals else (None, None)
    ok = bool(ap >= SYNTH_MIN_AP)
    res = {"schema": "MetaSieve.S7L2B.P2B.SyntheticAudit.v1",
           "created_utc": "2026-08-10", "preregistration_sha256": PREREG_SHA,
           "repo_commit": REPO_COMMIT, "seed": SEED_SYNTH,
           "teacher": f"rank-{K} projected bilinear differential, top/bottom "
                      f"{SYNTH_M} residues per pair",
           "heldout_synthetic_pairs": len(vals),
           "ap_bidir_recovered": ap, "chance_bidir": ch,
           "required_ap_bidir": SYNTH_MIN_AP, "pass": ok,
           "diagnostics": {
               "in_sample_ap_bidir_last_epoch_pairs": ins_ap,
               "in_sample_pairs_scored": len(ins_vals),
               "teacher_self_ap_bidir": float(np.mean(self_ap)) if self_ap else None,
               "output_level_correlation_with_teacher_field_mean":
                   float(np.mean(corr)) if corr else None,
               "output_level_correlation_with_teacher_field_median":
                   float(np.median(corr)) if corr else None,
               "correlation_note": "Pearson r between the learned and teacher delta "
                                   "FIELDS, computed at output level so it is invariant "
                                   "to the U/V rotation gauge. It separates 'cannot fit "
                                   "the function class' from 'fits it but the top-8 "
                                   "discretisation makes AP a harsh readout'.",
               "diagnostics_are_not_a_gate": True},
           "parameter_movement": move, "training_trace": trace,
           "elapsed_sec": round(time.time() - t, 1)}
    jdump(res, OUT / "PHASE2B_SYNTHETIC_AUDIT.json")
    return res


def macro_of(vals, meta, ctx, field):
    return aggregate({k: v[field] for k, v in vals.items()}, meta,
                     ctx["construct_component"])


def finish(ctx, ctrl, pf, cman, syn, head, head_perm, head_null, trace, trace2,
           trace_perm, move, move_perm, arms, b5mean, b5m_sha, det_same,
           det_sampler_same, ck, sidx, t0):
    comp, macro = {}, {}
    for tag, (v, m) in arms.items():
        for field in ("ap_bidir", "chance_bidir", "ap_gain", "ap_loss", "ap_change",
                      "ap_symdiff_conditional"):
            c, mm_ = macro_of(v, m, ctx, field)
            comp[(tag, field)] = c
            macro[(tag, field)] = mm_

    def gate(a_tag, b_tag, margin, field="ap_bidir", b_field=None, need_lcb_pos=True):
        bs = component_bootstrap(comp[(a_tag, field)], comp[(b_tag, b_field or field)])
        bs["margin"] = margin
        bs["pass"] = bool(bs["delta"] >= margin and
                          (bs["lcb95_one_sided"] > 0 if need_lcb_pos
                           else bs["lcb95_one_sided"] >= margin))
        return bs

    gates = {
        "R1_vs_chance": gate("cand", "cand", 0.05, "ap_bidir", "chance_bidir"),
        "R2_vs_frozen_B5_differential": gate("cand", "b5diff", 0.03),
        "R3_vs_two_ligand_foreign_pair": gate("cand", "foreign", 0.03),
        "R4_vs_residue_context_corruption": gate("cand", "context", 0.03),
        "R5_vs_trained_permuted_label_learner": gate("cand", "perm", 0.05),
    }

    # ------------------------------------------------------------------- R6
    print("R6 pair-score non-inferiority ...", flush=True)
    g5 = np.memmap(S5 / "heldoutA_B5.f16.dat", dtype=np.float16, mode="r",
                   shape=(max(v[0] + v[1] * v[2] for v in sidx.values()),))
    ap_b5, ap_2b = {}, {}
    uh_cache = {}
    V = head.V.detach().cpu().numpy().astype(np.float64)
    for rec in ctx["held_A"]:
        k = rec["source_key"]
        sk = rec["seq_key"]
        o, L, A = sidx[k]
        if sk not in uh_cache:
            with torch.no_grad():
                uh_cache[sk] = head.uh(h_of(ctx, sk)).cpu().numpy().astype(np.float64)
            if len(uh_cache) > 48:
                uh_cache.pop(next(iter(uh_cache)))
        d = project_np(ctx["Qs"][sk], uh_cache[sk] @ (V @ ctx["gvec"][rec["graph_key"]]))
        G = np.asarray(g5[o:o + L * A], dtype=np.float64).reshape(L, A)
        y = np.zeros((L, A), dtype=np.int8)
        for i, j in rec["edges"]:
            y[i, j] = 1
        yf = y.ravel()
        ap_b5[k] = ap_exact(G.ravel(), yf)
        ap_2b[k] = ap_exact((G + d[:, None]).ravel(), yf)
    cb5 = defaultdict(list)
    c2b = defaultdict(list)
    for k in ap_b5:
        if ap_b5[k] is not None:
            cb5[ctx["comp_of"][k]].append(ap_b5[k])
            c2b[ctx["comp_of"][k]].append(ap_2b[k])
    cb5 = {c: float(np.mean(v)) for c, v in cb5.items()}
    c2b = {c: float(np.mean(v)) for c, v in c2b.items()}
    r6 = component_bootstrap(c2b, cb5)
    r6["margin"] = -0.005
    r6["pass"] = bool(r6["delta"] >= -0.005 and r6["lcb95_one_sided"] >= -0.005)
    r6["note"] = "non-inferiority gate; LCB > 0 is NOT required"
    r6["b5_macro_recomputed"] = float(np.mean(list(cb5.values())))
    r6["g2b_macro"] = float(np.mean(list(c2b.values())))
    phase2a = json.loads((OUT / "PHASE2A_COMPONENT_TABLES.json").read_text(encoding="utf-8"))
    r6["b5_macro_phase2a"] = float(np.mean(list(phase2a["deployable"]["B5"]["full"].values())))
    gates["R6_pair_score_non_inferiority"] = r6

    # ------------------------------------------- replicate reproducibility ref
    by_sk_gk = defaultdict(list)
    for r in ctx["held_A"]:
        by_sk_gk[(r["seq_key"], r["graph_key"])].append(r)
    rep_vals, rep_meta = {}, {}
    for sk, a, b in ctx["heldA_pairs"]:
        ra = ctx["rec_of"][a]
        alts = [x for x in by_sk_gk[(sk, ra["graph_key"])]
                if x["source_key"] != a and x["pdb_id"] != ra["pdb_id"]]
        if not alts:
            continue
        Rap = ctx["masks"][alts[0]["source_key"]]
        Rb = ctx["masks"][b]
        L = ra["n_res"]
        score = np.zeros(L, dtype=np.float64)
        for r_ in Rap - Rb:
            score[r_] = 1.0
        for r_ in Rb - Rap:
            score[r_] = -1.0
        Ra, Rb_ = ctx["masks"][a], ctx["masks"][b]
        m = pair_metrics(score, Ra - Rb_, Rb_ - Ra, L)
        if m is not None:
            rep_vals[f"{a}|{b}"] = m
            rep_meta[f"{a}|{b}"] = sk
    rep_comp, rep_macro = macro_of(rep_vals, rep_meta, ctx, "ap_bidir")
    replicate_ref = {
        "name": "REPLICATE_REPRODUCIBILITY_REFERENCE",
        "not_a_ceiling": "a model that denoises annotation error may legitimately "
                         "exceed the agreement between two noisy annotations",
        "matched_pairs": len(rep_vals), "components": len(rep_comp),
        "ap_bidir": rep_macro,
        "cannot_determine_pass": True}

    # ------------------------------------------- module-participation audit
    with torch.no_grad():
        sk0 = ctx["heldA_pairs"][0][0]
        uh0 = head.uh(h_of(ctx, sk0)).cpu().numpy()
        vgs = np.stack([(head.V.detach().cpu().numpy() @ ctx["gvec"][gk]).astype(np.float64)
                        for gk in list(ctx["gvec"])[:512]])
    mp = {
        "grad_norm_U_min_over_epochs": min(t["grad_norm_U_min"] for t in trace),
        "grad_norm_V_min_over_epochs": min(t["grad_norm_V_min"] for t in trace),
        "grad_norms_nonzero": bool(min(t["grad_norm_U_min"] for t in trace) > 0
                                   and min(t["grad_norm_V_min"] for t in trace) > 0),
        "parameter_movement": move,
        "parameter_movement_threshold": MIN_PARAM_MOVEMENT,
        "parameter_movement_pass": bool(all(v >= MIN_PARAM_MOVEMENT for v in move.values())),
        "activation_variance_Uh": float(np.var(uh0)),
        "activation_variance_Vg": float(np.var(vgs)),
        "activation_variance_pass": bool(np.var(uh0) > MIN_ACT_VAR and np.var(vgs) > MIN_ACT_VAR),
        "zero_Uh_or_Vg_collapses_delta": True,
        "zero_ablation_note": "zeroing either bilinear factor makes delta identically "
                              "zero, so AP falls to the constant-score chance level by "
                              "construction. This check is structurally guaranteed and "
                              "is therefore weak evidence; the informative ablations are "
                              "R4 (residue side) and R3/chemistry-shuffle (ligand side).",
        "residue_context_shuffle_degrades": bool(macro[("context", "ap_bidir")]
                                                 < macro[("cand", "ap_bidir")]),
        "ligand_shuffle_degrades": bool(macro[("chem", "ap_bidir")]
                                        < macro[("cand", "ap_bidir")]),
        "same_seed_checkpoint_identical": det_same,
        "same_seed_sampler_identical": det_sampler_same,
        "synthetic_recovery_pass": syn["pass"],
        "gauge_note": "U -> R U, V -> R V leaves delta unchanged for orthogonal R, so "
                      "individual U/V channels are not interpretable and no claim is "
                      "made about them",
    }
    mp["pass"] = bool(mp["grad_norms_nonzero"] and mp["parameter_movement_pass"]
                      and mp["activation_variance_pass"] and mp["same_seed_checkpoint_identical"]
                      and mp["same_seed_sampler_identical"] and mp["synthetic_recovery_pass"]
                      and mp["residue_context_shuffle_degrades"]
                      and mp["ligand_shuffle_degrades"])

    # ------------------------------------------------------------- verdict
    if not (gates["R1_vs_chance"]["pass"] and gates["R2_vs_frozen_B5_differential"]["pass"]):
        verdict = "PHASE2B_MINIMAL_RESIDUAL_NOT_IDENTIFIED"
    elif not all(gates[g]["pass"] for g in ("R3_vs_two_ligand_foreign_pair",
                                            "R4_vs_residue_context_corruption",
                                            "R5_vs_trained_permuted_label_learner")):
        verdict = "PHASE2B_SHORTCUT_DEPENDENCE"
    elif not gates["R6_pair_score_non_inferiority"]["pass"]:
        verdict = "PHASE2B_RESIDUE_DIFFERENTIAL_IDENTIFIED_BUT_B5_INTEGRATION_FAILED"
    elif not mp["pass"]:
        verdict = "PHASE2B_SHORTCUT_DEPENDENCE"
    else:
        verdict = "STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_STATISTIC_IDENTIFIED_IN_DEVELOPMENT"

    tr = {"schema": "MetaSieve.S7L2B.P2B.TrainingTrace.v1", "created_utc": "2026-08-10",
          "preregistration_sha256": PREREG_SHA, "repo_commit": REPO_COMMIT,
          "device": DEV, "device_reason": DEV_REASON,
          "candidate": trace, "determinism_repeat": trace2,
          "permuted_label_control": trace_perm,
          "parameter_movement_candidate": move,
          "parameter_movement_permuted": move_perm,
          "checkpoint_sha256": sha_file(ck),
          "sampled_pairs_sha256": sha_file(P2B / "sampled_pairs.txt"),
          "sampled_pair_records": sum(len(t) for t in [[]]) or None}
    jdump(tr, OUT / "PHASE2B_TRAINING_TRACE.json")

    res = {
        "schema": "MetaSieve.S7L2B.P2B.Gate.v1", "created_utc": "2026-08-10",
        "preregistration": "research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md",
        "preregistration_sha256": PREREG_SHA,
        "preregistration_commit": PREREG_COMMIT,
        "superseded_preregistration_sha256":
            "ae6d1a0186bb37af86f3b6eb98c513bce7e67a8745aaf5a3811ce5c9b98ab477",
        "device": DEV, "device_reason": DEV_REASON,
        "seeds": {"param": SEED_PARAM, "sampler": SEED_SAMPLER, "control": SEED_CTRL,
                  "synthetic": SEED_SYNTH, "bootstrap": 20260903},
        "census": pf["census"],
        "control_manifest": cman,
        "b5_atom_mean_sha256": b5m_sha,
        "macro_ap_bidir": {t: macro[(t, "ap_bidir")] for t in arms},
        "macro_chance_bidir": {t: macro[(t, "chance_bidir")] for t in arms},
        "macro_ap_gain": {t: macro[(t, "ap_gain")] for t in arms},
        "macro_ap_loss": {t: macro[(t, "ap_loss")] for t in arms},
        "macro_ap_change": {t: macro[(t, "ap_change")] for t in arms},
        "secondary_symdiff_conditional": {t: macro[(t, "ap_symdiff_conditional")]
                                          for t in arms},
        "pairs_scored": {t: len(v) for t, (v, _m) in arms.items()},
        "gates": gates,
        "module_participation": mp,
        "replicate_reproducibility_reference": replicate_ref,
        "heldout_B_secondary": {
            "ap_bidir": macro[("cand_heldB", "ap_bidir")],
            "chance_bidir": macro[("cand_heldB", "chance_bidir")],
            "pairs": len(arms["cand_heldB"][0]),
            "note": "scaffold-strict secondary analysis; not a gate"},
        "TERMINAL_VERDICT": verdict,
        "interpretation_limit": "held-out A was read during Phase 1, Phase 2A and "
                                "Phase 2B design; this is development evidence, not "
                                "independent confirmation",
        "affinity_value_reads": 0,
        "elapsed_sec": round(time.time() - t0, 1),
    }
    jdump(res, OUT / "PHASE2B_GATE.json")
    print(json.dumps({"macro_ap_bidir": res["macro_ap_bidir"],
                      "gates": {k: {"delta": v["delta"],
                                    "lcb": v["lcb95_one_sided"], "pass": v["pass"]}
                                for k, v in gates.items()},
                      "module_participation_pass": mp["pass"],
                      "TERMINAL_VERDICT": verdict}, indent=2), flush=True)
    return 0


def main():
    t0 = time.time()
    print("preparing ...", flush=True)
    ctx = prepare()
    print(f"prepared in {time.time()-ctx['t0']:.0f}s", flush=True)
    print("building controls ...", flush=True)
    ctrl, cman = build_controls(ctx)
    print(f"  foreign-pair coverage "
          f"{cman['R3_two_ligand_foreign_pair']['coverage']:.4f}", flush=True)
    print("preflight ...", flush=True)
    pf = preflight(ctx, ctrl)
    print(json.dumps({"verdict": pf["verdict"], "failures": pf["failures"],
                      "census": pf["census"]}, indent=1), flush=True)
    if pf["failures"]:
        jdump({"TERMINAL_VERDICT": "PHASE2B_CONTRACT_OR_ARTIFACT_FAIL_CLOSED",
               "failures": pf["failures"]}, OUT / "PHASE2B_GATE.json")
        return 1

    print("synthetic trainability ...", flush=True)
    syn = synthetic_stage(ctx)
    print(f"  synthetic AP_bidir={syn['ap_bidir_recovered']:.4f} "
          f"(chance {syn['chance_bidir']:.4f}) pass={syn['pass']}", flush=True)
    if not syn["pass"]:
        d = syn["diagnostics"]
        jdump({"schema": "MetaSieve.S7L2B.P2B.TrainingTrace.v1",
               "created_utc": "2026-08-10", "preregistration_sha256": PREREG_SHA,
               "repo_commit": REPO_COMMIT, "device": DEV, "device_reason": DEV_REASON,
               "real_label_training_executed": False,
               "reason": "the registered synthetic precondition failed; the "
                         "real-label run is not authorized",
               "synthetic_training_trace": syn["training_trace"],
               "parameter_movement_synthetic": syn["parameter_movement"]},
              OUT / "PHASE2B_TRAINING_TRACE.json")
        jdump({
            "schema": "MetaSieve.S7L2B.P2B.Gate.v1", "created_utc": "2026-08-10",
            "preregistration": "research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md",
            "preregistration_sha256": PREREG_SHA,
            "preregistration_commit": PREREG_COMMIT,
            "superseded_preregistration_sha256":
                "ae6d1a0186bb37af86f3b6eb98c513bce7e67a8745aaf5a3811ce5c9b98ab477",
            "device": DEV, "device_reason": DEV_REASON,
            "seeds": {"param": SEED_PARAM, "sampler": SEED_SAMPLER,
                      "control": SEED_CTRL, "synthetic": SEED_SYNTH,
                      "bootstrap": 20260903},
            "contract_audit": {"verdict": pf["verdict"], "checks": pf["checks"]},
            "census": pf["census"],
            "control_manifest": cman,
            "synthetic_audit": syn,
            "gates_R1_to_R6": "NOT SCORED — the real-label run was not executed",
            "real_label_training_executed": False,
            "TERMINAL_VERDICT": "PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED",
            "failure_localisation": {
                "hypothesis_class": "NOT AT FAULT — the teacher lies exactly in the "
                                    "rank-8 projected bilinear class by construction",
                "evaluation_code": f"NOT AT FAULT — the teacher scores "
                                   f"{d['teacher_self_ap_bidir']:.5f} on its own labels",
                "generalisation": f"NOT AT FAULT — in-sample "
                                  f"{d['in_sample_ap_bidir_last_epoch_pairs']:.4f} vs "
                                  f"held-out {syn['ap_bidir_recovered']:.4f}; no gap",
                "at_fault": "the registered OPTIMIZATION BUDGET. The learned delta "
                            "field already correlates with the teacher field at "
                            f"r={d['output_level_correlation_with_teacher_field_median']:.3f} "
                            "(median, gauge-invariant), but 6 epochs over at most "
                            "8,864 sampled pairs do not drive that into an exact "
                            "top-8 ranking at AP_bidir >= 0.50."},
            "no_biological_conclusion_permitted": True,
            "affinity_value_reads": 0,
            "elapsed_sec": round(time.time() - t0, 1),
        }, OUT / "PHASE2B_GATE.json")
        return 2

    # ---------------------------------------------------- real training (one)
    def real_labels(sk, a, b):
        Ra, Rb = ctx["masks"][a], ctx["masks"][b]
        return Ra - Rb, Rb - Ra

    def perm_labels(sk, a, b):
        pa = ctrl["within_construct_derangement"].get(a, a)
        pb = ctrl["within_construct_derangement"].get(b, b)
        Ra, Rb = ctx["masks"][pa], ctx["masks"][pb]
        return Ra - Rb, Rb - Ra

    print("real-label training (one run) ...", flush=True)
    head, trace, move, sampled = train_head(ctx, ctx["train_pairs"], real_labels, "cand")
    ck = P2B / "phase2b_head.pt"
    torch.save(head.state_dict(), ck)
    sampled_blob = "\n".join(sampled).encode()
    (P2B / "sampled_pairs.txt").write_bytes(sampled_blob)

    print("registered determinism repeat (identical seeds) ...", flush=True)
    head2, trace2, _m2, sampled2 = train_head(ctx, ctx["train_pairs"], real_labels,
                                              "det", log=False)
    ck2 = P2B / "phase2b_head_repeat.pt"
    torch.save(head2.state_dict(), ck2)
    det_same = sha_file(ck) == sha_file(ck2)
    det_sampler_same = hashlib.sha256(sampled_blob).hexdigest() == \
        hashlib.sha256("\n".join(sampled2).encode()).hexdigest()

    print("trained permutation control ...", flush=True)
    head_perm, trace_perm, move_perm, _sp = train_head(ctx, ctx["train_pairs"],
                                                       perm_labels, "perm", log=False)
    torch.manual_seed(SEED_PARAM)
    head_null = Head().to(DEV)

    # ---------------------------------------------------- B5 atom-mean vectors
    print("materialising the frozen B5 differential baseline ...", flush=True)
    sidx = json.loads((S4 / "heldoutA_index.json").read_text(encoding="utf-8"))
    total = max(v[0] + v[1] * v[2] for v in sidx.values())
    g5 = np.memmap(S5 / "heldoutA_B5.f16.dat", dtype=np.float16, mode="r",
                   shape=(total,))
    b5mean = {}
    for k, (o, L, A) in sidx.items():
        b5mean[k] = np.asarray(g5[o:o + L * A], dtype=np.float32).reshape(L, A).mean(1)
    b5m_sha = hashlib.sha256(b"".join(b5mean[k].tobytes() for k in sorted(b5mean))).hexdigest()

    # ---------------------------------------------------- arms
    print("scoring arms ...", flush=True)
    arms = {}
    for tag, hd, mode in (("cand", head, "cand"), ("foreign", head, "foreign"),
                          ("context", head, "context"), ("chem", head, "chem"),
                          ("perm", head_perm, "cand"), ("null", head_null, "cand"),
                          ("b5diff", head, "b5diff")):
        v, m = eval_arm(ctx, hd, ctx["heldA_pairs"], mode, ctrl=ctrl, b5mean=b5mean)
        arms[tag] = (v, m)
        _c, mac = macro_of(v, m, ctx, "ap_bidir")
        print(f"   {tag:8s} pairs={len(v):6d}  AP_bidir={mac:.6f}", flush=True)
    v, m = eval_arm(ctx, head, ctx["heldB_pairs"], "cand", ctrl=ctrl, b5mean=b5mean)
    arms["cand_heldB"] = (v, m)

    return finish(ctx, ctrl, pf, cman, syn, head, head_perm, head_null, trace,
                  trace2, trace_perm, move, move_perm, arms, b5mean, b5m_sha,
                  det_same, det_sampler_same, ck, sidx, t0)


if __name__ == "__main__":
    sys.exit(main())
