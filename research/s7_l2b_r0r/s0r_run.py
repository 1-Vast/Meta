"""Complete-panel, label-blind replay for Phase 2B synthetic control.

This runner is governed by PREREG_PHASE2B_S0R_COMPLETE_PANEL_REPLAY.md.
It never imports p2b_run or opens MONN residue-edge labels.  Its pair universe
is the frozen metadata-only census committed before this file was written.
"""
from __future__ import annotations

import hashlib
import json
import math
import pickle
import platform
import subprocess
import sys
import time
from collections import OrderedDict, defaultdict
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "s7_l2b_r0r"
sys.path.insert(0, str(HERE))

from p2b_residue_residual import (  # noqa: E402
    CLIP, D_ATOM, D_ESM, K, LR, WD, Head, g_of, nuisance_basis, pair_loss,
    project_np, sha_file,
)
from s0_synth import (  # noqa: E402
    DEV_SEED, HCache, Teacher, build_panels, build_stream, hier_macro, pair_ap,
    pair_bce_np, panel_hash, stream_hash, stream_stats,
)
from s7_dataset import atom_features  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S0R_COMPLETE_PANEL_REPLAY.md"
PREREG_SHA = "6e30ead522a1217cf8442eb7185a459ffdc65690cad5ea1c73f85705cfd20400"
META_DIR = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "phase2b_s0r"
EXEC_DIR = META_DIR / "execution"
OUT = ROOT / "report" / "s7_l2b_r0r"
ESM = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "esm2_650M"
P2B = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "phase2b"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"

INPUT_HASHES = {
    "metadata_only_records.jsonl": "dcc2a0cd0cd958640f6548ed7cd3e5076e6c0ae9f6e455e35ad3e82f534dc856",
    "train_pairs.jsonl": "847a0770856e7e26903e56bc40d7249bb4bc21082392aa7edf1e3166014c1195",
    "heldoutA_pairs.jsonl": "29925ed5139b1d054aa3b2a26d5a0b281336e717fc801dcfae553b5e6cc340ae",
}
EXPECTED = {
    "train_universe": 228845,
    "train_panel": 14333,
    "train_panel_components": 298,
    "heldout": 44746,
    "heldout_components": 81,
}
CHECKPOINTS = (0, 1, 10, 100, 210)
N_UPDATES = 210
N_BOOT = 2000
BOOT_SEED = 20260903
AP_DROP = 0.05
DEVICE = "cpu"


class ContractError(RuntimeError):
    pass


def load_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def write_jsonl(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def git_head():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()


def state_hash(head: Head) -> str:
    h = hashlib.sha256()
    for name, tensor in sorted(head.state_dict().items()):
        h.update(name.encode())
        h.update(tensor.detach().cpu().contiguous().numpy().tobytes())
    return h.hexdigest()


def file_manifest(paths):
    return {str(p.relative_to(ROOT)).replace("\\", "/"): sha_file(p)
            for p in paths if p.exists() and p.is_file()}


def stamp():
    return {
        "stage": "P1R2B-PHASE2B-S0R_COMPLETE_PANEL_REPLAY",
        "preregistration_sha256": PREREG_SHA,
        "execution_commit": git_head(),
        "created_utc": "2026-08-10",
        "device": DEVICE,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "torch": torch.__version__,
            "cuda_available_but_unused": bool(torch.cuda.is_available()),
        },
        "synthetic_only": True,
        "real_structural_edge_label_reads": 0,
        "affinity_value_reads": 0,
    }


def load_molecules():
    out = {}
    for name in ("mol_dict", "independent_dataset_mol_dict"):
        with (MONN / name).open("rb") as f:
            d = pickle.load(f, encoding="bytes")
        for k, v in d.items():
            key = k.decode() if isinstance(k, bytes) else str(k)
            out.setdefault(key, v)
    return out


def prepare_label_blind_context():
    """Build the minimum context from an explicit metadata whitelist only."""
    if sha_file(PREREG) != PREREG_SHA:
        raise ContractError("preregistration hash mismatch")
    for name, expected in INPUT_HASHES.items():
        observed = sha_file(META_DIR / name)
        if observed != expected:
            raise ContractError(f"input hash mismatch: {name}: {observed}")

    records = load_jsonl(META_DIR / "metadata_only_records.jsonl")
    train_rows = load_jsonl(META_DIR / "train_pairs.jsonl")
    held_rows = load_jsonl(META_DIR / "heldoutA_pairs.jsonl")
    forbidden = {"edges", "positive_binary_edges", "mask", "gain_set",
                 "loss_set", "symmetric_difference"}
    for row in records:
        if forbidden & set(row):
            raise ContractError("forbidden structural label field reached runner")

    rec_of = {r["source_key"]: r for r in records}
    pairs = {}
    for name, rows in (("train", train_rows), ("heldout", held_rows)):
        parsed = []
        for row in rows:
            sk, a, b = row["seq_key"], row["a"], row["b"]
            if a not in rec_of or b not in rec_of:
                raise ContractError(f"unknown source key in {name} pair")
            ra, rb = rec_of[a], rec_of[b]
            if not (ra["seq_key"] == rb["seq_key"] == sk):
                raise ContractError("pair crosses exact sequence")
            if ra["graph_key"] == rb["graph_key"]:
                raise ContractError("pair has identical ligand graph")
            if not (ra["scaffold"] and rb["scaffold"]
                    and ra["scaffold"] != rb["scaffold"]):
                raise ContractError("pair violates scaffold-distinct contract")
            parsed.append((sk, a, b))
        pairs[name] = parsed

    if len(pairs["train"]) != EXPECTED["train_universe"]:
        raise ContractError("train-universe count mismatch")
    if len(pairs["heldout"]) != EXPECTED["heldout"]:
        raise ContractError("heldout count mismatch")

    idx_raw = json.loads((ESM / "esm2_650M_index.json").read_text(encoding="utf-8"))
    idx = {k: (int(v["offset"]), int(v["length"])) if isinstance(v, dict)
           else (int(v[0]), int(v[1])) for k, v in idx_raw.items()}
    total = max(o + length for o, length in idx.values())
    mm = np.memmap(ESM / "esm2_650M_residues.fp16.dat", dtype=np.float16,
                   mode="r", shape=(total, D_ESM))

    prior_idx = json.loads((P2B / "b_prior_index.json").read_text(encoding="utf-8"))
    prior = np.load(P2B / "b_prior_f64.npy", mmap_mode="r")
    needed_sk = sorted({sk for rows in pairs.values() for sk, _a, _b in rows})
    Qs = {}
    construct_component = {}
    for sk in needed_sk:
        if sk not in idx or sk not in prior_idx:
            raise ContractError(f"missing frozen state for construct {sk}")
        p = prior_idx[sk]
        off, length = (int(p["offset"]), int(p["length"])) if isinstance(p, dict) \
            else (int(p[0]), int(p[1]))
        if length != idx[sk][1]:
            raise ContractError(f"prior/ESM length mismatch for {sk}")
        Qs[sk] = nuisance_basis(np.asarray(prior[off:off + length], dtype=np.float64))
        comps = {r["component"] for r in records if r["seq_key"] == sk}
        if len(comps) != 1:
            raise ContractError(f"non-unique component for {sk}")
        construct_component[sk] = next(iter(comps))

    molecules = load_molecules()
    used_sources = {x for rows in pairs.values() for _sk, a, b in rows for x in (a, b)}
    needed_graphs = {rec_of[s]["graph_key"]: rec_of[s]["ligand_ccd"]
                     for s in used_sources}
    gvec = {}
    for graph, ccd in sorted(needed_graphs.items()):
        if ccd not in molecules:
            raise ContractError(f"missing molecule {ccd}")
        vec = g_of(atom_features(molecules[ccd]))
        if vec.shape != (D_ATOM,) or not np.isfinite(vec).all():
            raise ContractError(f"invalid atom feature vector for {ccd}")
        old = gvec.get(graph)
        if old is not None and not np.array_equal(old, vec):
            raise ContractError(f"graph key maps to inconsistent features: {graph}")
        gvec[graph] = vec

    ctx = {
        "records": records,
        "rec_of": rec_of,
        "train_pairs": pairs["train"],
        "heldA_pairs": pairs["heldout"],
        "construct_component": construct_component,
        "idx": idx,
        "mm": mm,
        "Qs": Qs,
        "gvec": gvec,
    }
    return ctx


def serialize_stream(path: Path, stream):
    rows = []
    for update, (epoch, batch_index, batch) in enumerate(stream, start=1):
        rows.append({
            "update": update,
            "epoch": epoch,
            "batch_index": batch_index,
            "batch": [{
                "component": comp,
                "constructs": [{"seq_key": sk,
                                "pairs": [[a, b] for _same_sk, a, b in pl]}
                               for sk, pl in constructs],
            } for comp, constructs in batch],
        })
    write_jsonl(path, rows)


def deserialize_stream(path: Path):
    stream = []
    for row in load_jsonl(path):
        batch = []
        for c in row["batch"]:
            constructs = []
            for x in c["constructs"]:
                sk = x["seq_key"]
                constructs.append((sk, [(sk, a, b) for a, b in x["pairs"]]))
            batch.append((c["component"], constructs))
        stream.append((int(row["epoch"]), int(row["batch_index"]), batch))
    return stream


def freeze_stream(ctx):
    generated = build_stream(ctx, ctx["train_pairs"], e_max=6)[:N_UPDATES]
    if len(generated) != N_UPDATES:
        raise ContractError("insufficient deterministic updates")
    path = EXEC_DIR / "frozen_stream.jsonl"
    serialize_stream(path, generated)
    disk = deserialize_stream(path)
    if stream_hash(disk) != stream_hash(generated):
        raise ContractError("serialized stream semantic hash mismatch")
    known = {(sk, a, b) for sk, a, b in ctx["train_pairs"]}
    for _ep, _bi, batch in disk:
        for comp, constructs in batch:
            for sk, pl in constructs:
                if ctx["construct_component"][sk] != comp:
                    raise ContractError("stream component mismatch")
                if any(p not in known for p in pl):
                    raise ContractError("stream contains unknown pair")
    return disk, {
        "file": str(path.relative_to(ROOT)).replace("\\", "/"),
        "file_sha256": sha_file(path),
        "semantic_sha256": stream_hash(disk),
        **stream_stats(disk),
    }


def head_arrays(head: Head):
    return (head.U.detach().cpu().numpy().astype(np.float64),
            head.V.detach().cpu().numpy().astype(np.float64))


def score_panel(ctx, head, rows, hcache, artifact_prefix: Path | None = None):
    """Score every row and optionally persist raw residue scores plus labels."""
    U, V = head_arrays(head)
    by_sk = defaultdict(list)
    for row in rows:
        by_sk[row["sk"]].append(row)
    aps, bces, ordered_rows = [], [], []
    flat_scores = []
    index_rows = []
    offset = 0
    for sk in sorted(by_sk):
        h = hcache.get(sk).astype(np.float64)
        uh = h @ U.T
        Q = ctx["Qs"][sk]
        for row in by_sk[sk]:
            gd = ctx["gvec"][row["gk_a"]] - ctx["gvec"][row["gk_b"]]
            d = project_np(Q, uh @ (V @ gd))
            ap = pair_ap(d, row)[0]
            bce = pair_bce_np(d, row)
            aps.append(ap)
            bces.append(bce)
            ordered_rows.append(row)
            if artifact_prefix is not None:
                flat_scores.append(d)
                index_rows.append({
                    "sk": sk, "a": row["a"], "b": row["b"],
                    "component": ctx["construct_component"][sk],
                    "offset": offset, "length": int(d.size),
                    "gain_indices": row["gi"].tolist(),
                    "loss_indices": row["li"].tolist(),
                    "ap_bidir": ap, "bce": bce,
                })
                offset += int(d.size)
    comp_ap, macro_ap = hier_macro(aps, ordered_rows, ctx)
    comp_bce, macro_bce = hier_macro(bces, ordered_rows, ctx)
    by_construct_ap = defaultdict(list)
    by_construct_bce = defaultdict(list)
    for ap, bce, row in zip(aps, bces, ordered_rows):
        if ap is not None:
            by_construct_ap[row["sk"]].append(ap)
        if bce is not None:
            by_construct_bce[row["sk"]].append(bce)
    construct_ap = {k: float(np.mean(v)) for k, v in by_construct_ap.items()}
    construct_bce = {k: float(np.mean(v)) for k, v in by_construct_bce.items()}
    result = {
        "pairs": len(rows),
        "constructs": len({r["sk"] for r in rows}),
        "components": len(comp_ap),
        "component_macro_ap_bidir": macro_ap,
        "component_macro_bce": macro_bce,
        "component_ap": comp_ap,
        "component_bce": comp_bce,
        "construct_ap": construct_ap,
        "construct_bce": construct_bce,
    }
    artifacts = {}
    if artifact_prefix is not None:
        scores_path = artifact_prefix.with_name(artifact_prefix.name + "_scores_f64.npy")
        index_path = artifact_prefix.with_name(artifact_prefix.name + "_index.jsonl")
        np.save(scores_path, np.concatenate(flat_scores).astype(np.float64, copy=False))
        write_jsonl(index_path, index_rows)
        artifacts = file_manifest([scores_path, index_path])
    return result, artifacts


def replay_prediction_artifacts(scores_path: Path, index_path: Path,
                                construct_component):
    scores = np.load(scores_path, mmap_mode="r")
    aps, bces, rows = [], [], []
    for item in load_jsonl(index_path):
        start, length = int(item["offset"]), int(item["length"])
        d = np.asarray(scores[start:start + length], dtype=np.float64)
        yg = np.zeros(length, dtype=np.int8)
        yl = np.zeros(length, dtype=np.int8)
        yg[np.asarray(item["gain_indices"], dtype=np.int64)] = 1
        yl[np.asarray(item["loss_indices"], dtype=np.int64)] = 1
        row = {
            "sk": item["sk"], "a": item["a"], "b": item["b"],
            "L": length, "gi": np.flatnonzero(yg), "li": np.flatnonzero(yl),
            "yg": yg, "yl": yl,
        }
        aps.append(pair_ap(d, row)[0])
        bces.append(pair_bce_np(d, row))
        rows.append(row)
    tiny_ctx = {"construct_component": construct_component}
    comp_ap, macro_ap = hier_macro(aps, rows, tiny_ctx)
    _comp_bce, macro_bce = hier_macro(bces, rows, tiny_ctx)
    return {"pairs": len(rows), "components": len(comp_ap),
            "component_macro_ap_bidir": macro_ap,
            "component_macro_bce": macro_bce}


def candidate_path_witness(ctx, teacher, panels, hcache):
    head = Head()
    with torch.no_grad():
        head.U.copy_(torch.from_numpy(teacher.U).float())
        head.V.copy_(torch.from_numpy(teacher.V).float())
    U, V = head_arrays(head)
    rel_max = 0.0
    teacher_ap, candidate_ap, rows = [], [], []
    anti_max = 0.0
    identical_max = 0.0
    for row in panels["heldout"]:
        t = teacher.field(row["sk"], row["gk_a"], row["gk_b"])
        h = hcache.get(row["sk"]).astype(np.float64)
        gd = ctx["gvec"][row["gk_a"]] - ctx["gvec"][row["gk_b"]]
        d = project_np(ctx["Qs"][row["sk"]], (h @ U.T) @ (V @ gd))
        db = project_np(ctx["Qs"][row["sk"]], (h @ U.T) @ (V @ -gd))
        dz = project_np(ctx["Qs"][row["sk"]], (h @ U.T) @ (V @ (gd - gd)))
        rel_max = max(rel_max, float(np.linalg.norm(d - t)
                                     / (1e-30 + np.linalg.norm(t))))
        anti_max = max(anti_max, float(np.max(np.abs(d + db))))
        identical_max = max(identical_max, float(np.max(np.abs(dz))))
        teacher_ap.append(pair_ap(t, row)[0])
        candidate_ap.append(pair_ap(d, row)[0])
        rows.append(row)
    _ct, tap = hier_macro(teacher_ap, rows, ctx)
    _cc, cap = hier_macro(candidate_ap, rows, ctx)
    wc = U.T @ V
    wt = teacher.W
    w_rel = float(np.linalg.norm(wc - wt) / np.linalg.norm(wt))
    failures = []
    if rel_max > 1e-4:
        failures.append("candidate field relative error exceeds 1e-4")
    if abs(cap - tap) > 1e-3:
        failures.append("candidate/teacher AP difference exceeds 1e-3")
    if w_rel > 1e-4:
        failures.append("candidate/teacher product error exceeds 1e-4")
    if anti_max > 1e-12 or identical_max > 1e-15:
        failures.append("antisymmetry or identical-ligand invariant failed")
    return {
        **stamp(), "complete_heldout_pairs_checked": len(rows),
        "complete_heldout_components": len({ctx["construct_component"][r["sk"]]
                                              for r in rows}),
        "max_relative_field_error": rel_max,
        "teacher_ap_bidir": tap, "candidate_ap_bidir": cap,
        "ap_abs_difference": abs(cap - tap),
        "relative_product_error": w_rel,
        "max_antisymmetry_error": anti_max,
        "max_identical_ligand_error": identical_max,
        "failures": failures,
        "verdict": "PASS" if not failures else "SYNTHETIC_CONTRACT_INVALID",
    }


def estimate_ray_scale(ctx, teacher, train_rows):
    fields = [(row, teacher.field(row["sk"], row["gk_a"], row["gk_b"]))
              for row in train_rows]

    def objective(scale):
        vals = [pair_bce_np(scale * d, row) for row, d in fields]
        return hier_macro(vals, [row for row, _d in fields], ctx)[1]

    grid = np.logspace(-3.0, 3.0, 61)
    values = [objective(float(a)) for a in grid]
    best = int(np.argmin(values))
    lo = float(grid[max(0, best - 1)])
    hi = float(grid[min(len(grid) - 1, best + 1)])
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(40):
        x1 = hi - phi * (hi - lo)
        x2 = lo + phi * (hi - lo)
        if objective(x1) <= objective(x2):
            hi = x2
        else:
            lo = x1
    scale = (lo + hi) / 2.0
    return scale, {
        "estimated_on": "complete label-blind hash-stratified train panel",
        "train_pairs": len(train_rows),
        "train_components": len({ctx["construct_component"][r["sk"]]
                                  for r in train_rows}),
        "grid": grid.tolist(), "grid_bce": values,
        "a_star": scale, "bce_at_a_star": objective(scale),
    }


def factor_balance(U, V):
    a = U @ U.T
    b = V @ V.T
    return float(np.linalg.norm(a - b) / (1e-30 + np.linalg.norm(a) + np.linalg.norm(b)))


def train_one_update(ctx, head, optimizer, batch, teacher, hcache):
    optimizer.zero_grad(set_to_none=True)
    batch_loss = torch.zeros((), dtype=torch.float32)
    for _comp, constructs in batch:
        comp_loss = torch.zeros((), dtype=torch.float32)
        for sk, pair_list in constructs:
            h_t = torch.from_numpy(hcache.get(sk).copy()).float()
            uh = head.uh(h_t)
            Q = ctx["Qs"][sk]
            losses = []
            for _same_sk, a, b in pair_list:
                ga = ctx["gvec"][ctx["rec_of"][a]["graph_key"]]
                gb = ctx["gvec"][ctx["rec_of"][b]["graph_key"]]
                raw = uh @ head.vg(torch.from_numpy((ga - gb).copy()).float())
                q = torch.from_numpy(Q).to(dtype=raw.dtype)
                d = raw - q @ (q.T @ raw) if q.shape[1] else raw
                gain, loss = teacher.labels(sk, a, b)
                losses.append(pair_loss(d, gain, loss, len(d), DEVICE))
            if losses:
                comp_loss = comp_loss + torch.stack(losses).mean() / len(constructs)
        batch_loss = batch_loss + comp_loss / len(batch)
    batch_loss.backward()
    torch.nn.utils.clip_grad_norm_(head.parameters(), CLIP)
    optimizer.step()
    if not all(torch.isfinite(p).all() for p in head.parameters()):
        raise ContractError("non-finite model parameter")
    return float(batch_loss.detach())


def run_trajectory(tag, ctx, panels, stream, teacher, hcache, init_U, init_V):
    run_dir = EXEC_DIR / tag
    run_dir.mkdir(parents=True, exist_ok=True)
    head = Head()
    with torch.no_grad():
        head.U.copy_(torch.from_numpy(init_U).float())
        head.V.copy_(torch.from_numpy(init_V).float())
    optimizer = torch.optim.AdamW(head.parameters(), lr=LR, weight_decay=WD)
    trajectory = []
    reload_rows = []
    artifact_paths = []

    def capture(update, last_loss):
        prefix = run_dir / f"u{update:04d}_heldout"
        train_metrics, _ = score_panel(ctx, head, panels["train"], hcache)
        held_metrics, pred_artifacts = score_panel(
            ctx, head, panels["heldout"], hcache, prefix
        )
        checkpoint = run_dir / f"u{update:04d}_checkpoint.pt"
        torch.save({
            "update": update,
            "model": head.state_dict(),
            "optimizer": optimizer.state_dict(),
            "last_batch_loss": last_loss,
        }, checkpoint)
        artifact_paths.extend([
            checkpoint,
            prefix.with_name(prefix.name + "_scores_f64.npy"),
            prefix.with_name(prefix.name + "_index.jsonl"),
        ])
        model_sha = state_hash(head)
        loaded = torch.load(checkpoint, map_location="cpu", weights_only=False)
        reloaded_head = Head()
        reloaded_head.load_state_dict(loaded["model"])
        checkpoint_match = state_hash(reloaded_head) == model_sha
        replay = replay_prediction_artifacts(
            prefix.with_name(prefix.name + "_scores_f64.npy"),
            prefix.with_name(prefix.name + "_index.jsonl"),
            ctx["construct_component"],
        )
        metric_match = (
            replay["pairs"] == held_metrics["pairs"]
            and replay["components"] == held_metrics["components"]
            and abs(replay["component_macro_ap_bidir"]
                    - held_metrics["component_macro_ap_bidir"]) <= 1e-14
            and abs(replay["component_macro_bce"]
                    - held_metrics["component_macro_bce"]) <= 1e-14
        )
        reload_rows.append({
            "tag": tag, "update": update,
            "checkpoint_state_sha256": model_sha,
            "checkpoint_reload_identical": checkpoint_match,
            "prediction_reload_metrics_identical": metric_match,
            "replayed": replay,
        })
        trajectory.append({
            "update": update, "last_batch_loss": last_loss,
            "checkpoint": str(checkpoint.relative_to(ROOT)).replace("\\", "/"),
            "checkpoint_sha256": sha_file(checkpoint),
            "checkpoint_state_sha256": model_sha,
            "prediction_artifacts": pred_artifacts,
            "factor_balance_residual": factor_balance(*head_arrays(head)),
            "train": train_metrics,
            "heldout": held_metrics,
        })

    capture(0, None)
    last_loss = None
    for update, (_epoch, _batch_index, batch) in enumerate(stream, start=1):
        last_loss = train_one_update(ctx, head, optimizer, batch, teacher, hcache)
        if update in CHECKPOINTS:
            capture(update, last_loss)
    if len(stream) != N_UPDATES or trajectory[-1]["update"] != N_UPDATES:
        raise ContractError("trajectory did not reach registered update count")
    return {
        **stamp(), "trajectory": tag,
        "initial_product_sha256": hashlib.sha256(
            (init_U.T @ init_V).astype(np.float64).tobytes()).hexdigest(),
        "checkpoints": trajectory,
    }, reload_rows, artifact_paths


def bootstrap_delta(comp0, comp100):
    keys = sorted(set(comp0) & set(comp100))
    if len(keys) != EXPECTED["heldout_components"]:
        raise ContractError("bootstrap component count mismatch")
    delta = np.asarray([comp100[k] - comp0[k] for k in keys], dtype=np.float64)
    rng = np.random.default_rng(BOOT_SEED)
    draws = rng.integers(0, len(keys), size=(N_BOOT, len(keys)))
    boot = delta[draws].mean(axis=1)
    return {
        "units": len(keys), "replicates": N_BOOT, "seed": BOOT_SEED,
        "delta_ap_point": float(delta.mean()),
        "lcb95_one_sided": float(np.percentile(boot, 5.0)),
        "ucb95_one_sided": float(np.percentile(boot, 95.0)),
        "ci95": [float(np.percentile(boot, 2.5)),
                 float(np.percentile(boot, 97.5))],
    }


def trajectory_decision(traj):
    by_update = {x["update"]: x for x in traj["checkpoints"]}
    x0, x100 = by_update[0], by_update[100]
    inference = bootstrap_delta(
        x0["heldout"]["component_ap"], x100["heldout"]["component_ap"]
    )
    train_bce_down = (
        x100["train"]["component_macro_bce"]
        < x0["train"]["component_macro_bce"]
    )
    held_drop = (
        x0["heldout"]["component_macro_ap_bidir"]
        - x100["heldout"]["component_macro_ap_bidir"]
    )
    misaligned = bool(train_bce_down and held_drop >= AP_DROP
                      and inference["ucb95_one_sided"] < 0.0)
    return {
        "train_bce_decreased": train_bce_down,
        "heldout_ap_drop": held_drop,
        "registered_minimum_drop": AP_DROP,
        "component_inference": inference,
        "MISALIGNED": misaligned,
    }


def write_report(verdict, input_manifest, candidate, ray, decisions):
    path = OUT / "PHASE2B_S0R_COMPLETE_PANEL_REPORT.md"
    lines = [
        "# Phase 2B S0R complete-panel replay",
        "",
        f"Terminal verdict: `{verdict}`",
        "",
        "This was a synthetic-only replay over a pair universe constructed from",
        "metadata fields only. No MONN residue-edge or affinity value was read.",
        "",
        "## Contract",
        "",
        f"- Train hash panel: {input_manifest['panels']['train']['pairs']:,} pairs / "
        f"{input_manifest['panels']['train']['components']} components",
        f"- Complete held-out A: {input_manifest['panels']['heldout']['pairs']:,} pairs / "
        f"{input_manifest['panels']['heldout']['components']} components",
        f"- Candidate-path teacher AP: {candidate['teacher_ap_bidir']:.6f}",
        f"- Balanced train-only ray scale: {ray['a_star']:.8f}",
        "",
        "## Decision",
        "",
        "| trajectory | AP(0) | AP(100) | drop | UCB95(delta) | misaligned |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for tag in ("original_gauge", "balanced_gauge"):
        d = decisions[tag]
        tr = d["trajectory"]
        by_u = {x["update"]: x for x in tr["checkpoints"]}
        lines.append(
            f"| {tag} | {by_u[0]['heldout']['component_macro_ap_bidir']:.6f} | "
            f"{by_u[100]['heldout']['component_macro_ap_bidir']:.6f} | "
            f"{d['decision']['heldout_ap_drop']:.6f} | "
            f"{d['decision']['component_inference']['ucb95_one_sided']:.6f} | "
            f"{str(d['decision']['MISALIGNED']).lower()} |"
        )
    lines += [
        "",
        "The result localizes only the synthetic control. It does not identify",
        "biology, affinity, few-shot adaptation or a biological z coordinate.",
        "Real Phase 2B and the frozen law operator remain untouched.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    started = time.time()
    EXEC_DIR.mkdir(parents=True, exist_ok=True)
    generated_reports = []
    execution_artifacts = []
    try:
        ctx = prepare_label_blind_context()
        hcache = HCache(ctx)
        teacher = Teacher(ctx, hcache, DEV_SEED)
        panels = build_panels(ctx, teacher, hcache)
        panel_counts = {
            "train": {
                "pairs": len(panels["train"]),
                "constructs": len({r["sk"] for r in panels["train"]}),
                "components": len({ctx["construct_component"][r["sk"]]
                                   for r in panels["train"]}),
                "synthetic_panel_sha256": panel_hash(panels["train"]),
            },
            "heldout": {
                "pairs": len(panels["heldout"]),
                "constructs": len({r["sk"] for r in panels["heldout"]}),
                "components": len({ctx["construct_component"][r["sk"]]
                                   for r in panels["heldout"]}),
                "synthetic_panel_sha256": panel_hash(panels["heldout"]),
            },
        }
        if panel_counts["train"]["pairs"] != EXPECTED["train_panel"]:
            raise ContractError("train hash panel count mismatch")
        if panel_counts["train"]["components"] != EXPECTED["train_panel_components"]:
            raise ContractError("train hash panel component count mismatch")
        if panel_counts["heldout"]["pairs"] != EXPECTED["heldout"]:
            raise ContractError("complete heldout pair count mismatch")
        if panel_counts["heldout"]["components"] != EXPECTED["heldout_components"]:
            raise ContractError("complete heldout component count mismatch")

        input_manifest = {
            **stamp(),
            "whitelisted_runner_fields": [
                "source_key", "pdb_id", "uniprot_id", "uniprot_sequence",
                "ligand_ccd", "seq_key", "graph_key", "scaffold", "n_atoms",
                "n_res", "component",
            ],
            "forbidden_fields_present": 0,
            "input_hashes": INPUT_HASHES,
            "panels": panel_counts,
            "sources_opened": [
                "metadata_only_records.jsonl", "train_pairs.jsonl",
                "heldoutA_pairs.jsonl", "frozen ESM2 residue cache",
                "frozen B5 protein-only prior", "MONN ligand mol_dict files",
            ],
            "explicitly_not_opened": [
                "MONN interaction edge corpus", "Phase2A residue masks",
                "ChEMBL", "BindingDB", "DAVIS", "KIBA", "recipient labels",
            ],
        }
        p = OUT / "S0R_INPUT_AND_FIREWALL_MANIFEST.json"
        write_json(p, input_manifest)
        generated_reports.append(p)

        stream, stream_info = freeze_stream(ctx)
        execution_artifacts.append(EXEC_DIR / "frozen_stream.jsonl")
        stream_manifest = {**stamp(), "stream": stream_info,
                           "training_reads_serialized_stream": True}
        p = OUT / "S0R_FROZEN_STREAM_MANIFEST.json"
        write_json(p, stream_manifest)
        generated_reports.append(p)

        candidate = candidate_path_witness(ctx, teacher, panels, hcache)
        p = OUT / "S0R_CANDIDATE_PATH_WITNESS.json"
        write_json(p, candidate)
        generated_reports.append(p)
        if candidate["failures"]:
            raise ContractError("candidate path witness failed")

        a_star, ray = estimate_ray_scale(ctx, teacher, panels["train"])
        sqrt_a = math.sqrt(a_star)
        original_U, original_V = teacher.U.copy(), teacher.V.copy()
        balanced_U, balanced_V = sqrt_a * teacher.U, sqrt_a * teacher.V
        product_rel = float(np.linalg.norm(
            balanced_U.T @ balanced_V - a_star * teacher.W
        ) / np.linalg.norm(a_star * teacher.W))
        ray.update({
            "balanced_parameterisation": "U=sqrt(a*)U*, V=sqrt(a*)V*",
            "relative_product_equality_error": product_rel,
            "original_factor_balance_residual": factor_balance(original_U, original_V),
            "balanced_factor_balance_residual": factor_balance(balanced_U, balanced_V),
            "factor_norm_ratio_original": float(np.linalg.norm(original_U)
                                                / np.linalg.norm(original_V)),
            "factor_norm_ratio_balanced": float(np.linalg.norm(balanced_U)
                                                / np.linalg.norm(balanced_V)),
        })
        if product_rel > 1e-12:
            raise ContractError("balanced product equality failed")

        print("S0R original-gauge trajectory", flush=True)
        original, reload_o, artifacts_o = run_trajectory(
            "original_gauge", ctx, panels, stream, teacher, hcache,
            original_U, original_V,
        )
        execution_artifacts.extend(artifacts_o)
        print("S0R balanced-gauge trajectory", flush=True)
        balanced, reload_b, artifacts_b = run_trajectory(
            "balanced_gauge", ctx, panels, stream, teacher, hcache,
            balanced_U, balanced_V,
        )
        execution_artifacts.extend(artifacts_b)

        p_original = OUT / "S0R_ORIGINAL_GAUGE_TRAJECTORY.json"
        p_balanced = OUT / "S0R_BALANCED_GAUGE_TRAJECTORY.json"
        write_json(p_original, original)
        write_json(p_balanced, {**balanced, "train_only_ray_audit": ray})
        generated_reports.extend([p_original, p_balanced])

        reload_rows = reload_o + reload_b
        reload_ok = all(x["checkpoint_reload_identical"]
                        and x["prediction_reload_metrics_identical"]
                        for x in reload_rows)
        reload_audit = {**stamp(), "checks": reload_rows,
                        "all_reload_checks_pass": reload_ok}
        p = OUT / "S0R_ARTIFACT_RELOAD_AUDIT.json"
        write_json(p, reload_audit)
        generated_reports.append(p)
        if not reload_ok:
            raise ContractError("checkpoint or prediction reload mismatch")

        decisions = {
            "original_gauge": {"trajectory": original,
                               "decision": trajectory_decision(original)},
            "balanced_gauge": {"trajectory": balanced,
                               "decision": trajectory_decision(balanced)},
        }
        inference = {
            **stamp(),
            "registered_rule": {
                "train_bce_100_less_than_bce_0": True,
                "heldout_ap_drop_at_least": AP_DROP,
                "one_sided_ucb95_delta_ap_less_than": 0.0,
            },
            "original_gauge": decisions["original_gauge"]["decision"],
            "balanced_gauge": decisions["balanced_gauge"]["decision"],
        }
        p = OUT / "S0R_COMPLETE_PANEL_COMPONENT_INFERENCE.json"
        write_json(p, inference)
        generated_reports.append(p)

        original_bad = decisions["original_gauge"]["decision"]["MISALIGNED"]
        balanced_bad = decisions["balanced_gauge"]["decision"]["MISALIGNED"]
        if balanced_bad:
            verdict = "SURROGATE_AP_MISALIGNMENT_FULL_PANEL"
        elif original_bad:
            verdict = "SCALE_PARAMETERIZATION_MISMATCH"
        else:
            verdict = "SUBSET_SELECTION_ARTIFACT"
        verdict_artifact = {
            **stamp(), "TERMINAL_VERDICT": verdict,
            "original_gauge_misaligned": original_bad,
            "balanced_gauge_misaligned": balanced_bad,
            "complete_heldout_pairs": len(panels["heldout"]),
            "complete_heldout_components": panel_counts["heldout"]["components"],
            "elapsed_seconds": round(time.time() - started, 3),
            "remaining_frozen": [
                "S1 repaired control", "budget scaling", "real Phase2B labels",
                "affinity", "independent structural confirmation",
                "few-shot section", "biological z", "A(F,z)=K(B(z)F(z))",
            ],
        }
        p_verdict = OUT / "PHASE2B_S0R_VERDICT.json"
        write_json(p_verdict, verdict_artifact)
        generated_reports.append(p_verdict)
        p_report = write_report(verdict, input_manifest, candidate, ray, decisions)
        generated_reports.append(p_report)

        not_run = {
            **stamp(), "status": "NOT_RUN",
            "reason": f"S0R terminal verdict {verdict} selects a separate next contract",
            "not_run": ["S1", "budget scaling", "real Phase2B", "affinity",
                        "few-shot section", "z admission"],
        }
        p = OUT / "PHASE2B_POST_S0R_NOT_RUN.json"
        write_json(p, not_run)
        generated_reports.append(p)

        manifest = {
            **stamp(),
            "reports": file_manifest(generated_reports),
            "execution_artifacts": file_manifest(
                [EXEC_DIR / "frozen_stream.jsonl", *execution_artifacts]
            ),
        }
        write_json(OUT / "S0R_DETACHED_ARTIFACT_MANIFEST.json", manifest)
        print(json.dumps({"TERMINAL_VERDICT": verdict,
                          "elapsed_seconds": round(time.time() - started, 1)},
                         indent=2), flush=True)
        return 0
    except Exception as exc:
        verdict = {
            **stamp(), "TERMINAL_VERDICT": "SYNTHETIC_CONTRACT_INVALID",
            "error_type": type(exc).__name__, "error": str(exc),
            "elapsed_seconds": round(time.time() - started, 3),
            "all_downstream": "NOT_RUN",
        }
        write_json(OUT / "PHASE2B_S0R_VERDICT.json", verdict)
        write_json(OUT / "PHASE2B_POST_S0R_NOT_RUN.json", {
            **stamp(), "status": "NOT_RUN", "reason": str(exc),
            "not_run": ["S1", "budget scaling", "real Phase2B", "affinity",
                        "few-shot section", "z admission"],
        })
        print(json.dumps(verdict, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
