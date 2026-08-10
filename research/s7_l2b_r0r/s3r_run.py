"""S3R real structural gauge-free residue differential.

This runner preserves the historical Phase 2B artifacts.  It creates separate
train and held-out label views, trains the S2R direct-W ordinal estimator in
child processes, and scores R1-R5 once on the development holdout.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "s7_l2b_r0r"
sys.path.insert(0, str(HERE))

import p2b_run as P  # noqa: E402
from p2b_residue_residual import (  # noqa: E402
    BATCH_COMPONENTS, CLIP, D_ATOM, D_ESM, EPOCHS, LR, MIN_PARAM_MOVEMENT,
    ORTHO_TOL, SEED_BOOT, SEED_PARAM, SEED_SAMPLER, S4, S5, aggregate,
    ap_exact, build_pairs, component_bootstrap, context_shuffle, g_of,
    nuisance_basis, pair_metrics, project_np, sha_file,
)
from s1r_run import rank_loss_torch  # noqa: E402
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import FeatureStore  # noqa: E402
from s7_run import load_mols  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S3R_REAL_STRUCTURAL_DIRECT_W.md"
AMENDMENT = HERE / "PREREG_PHASE2B_S3R_AMENDMENT_01.md"
AMENDMENT_02 = HERE / "PREREG_PHASE2B_S3R_AMENDMENT_02.md"
PREREG_SHA = "16cef1a477dfe74cc4bc8194065accdccf5936983f918730e6a312196a2fc7ad"
AMENDMENT_SHA = "8f93cdca8f4ea4b0818d911c40bca771e3bd77851b4056ab15d80a2f59107c3e"
AMENDMENT_02_SHA = "27d1cb7ff86f6ef2f9e92f7a352dc773a8ba832bea868aefefac7af19a96ca08"
S2R_VERDICT_SHA = "003261b92b193fab80ea04bdb08da40d768fe37587726de36710d0a3360a3ef7"
S2R_CAL_SHA = "f84ada8e9e887de59c6263e191164f67ad850ba4a91eb6078a33be5af13aad2f"
S2R_SEALED_SHA = "610694caf136bc4a348c6a466a6ae8a82a5c5be64d053f95b7d55e5284ad1340"
CONTROL_SHA = "e187a5f00f0b66328877bacd93b22471fe607e382e811f2674ecfc4a9dec9c33"

OUT = ROOT / "report" / "s7_l2b_r0r"
EXEC = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "phase2b_s3r"
ESM = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "esm2_650M"
OLD_P2B = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "phase2b"
CONTROL_PATH = OLD_P2B / "control_maps.json"
NORM_EPS = 1e-12
N_UPDATES = EPOCHS * 35
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
EXPECTED = {
    "train_pairs": 226765,
    "train_components": 554,
    "heldoutA_pairs": 46818,
    "heldoutA_components": 112,
    "heldoutB_pairs": 30661,
}
LABEL_READ_LOG = []


class S3RContractError(RuntimeError):
    pass


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str),
                    encoding="utf-8")


def write_jsonl(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def read_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            yield json.loads(line)


def write_gzip_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(payload)


def read_gzip_json(path: Path):
    if path.name.endswith("_residue_masks.json.gz"):
        LABEL_READ_LOG.append(str(path.resolve()))
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


def git_execution_surface_clean() -> bool:
    paths = [str(PREREG.relative_to(ROOT)), str(AMENDMENT.relative_to(ROOT)),
             str(AMENDMENT_02.relative_to(ROOT)),
             "research/s7_l2b_r0r/s3r_run.py",
             "tests/test_s7_l2b_phase2b_s3r.py"]
    output = subprocess.check_output(["git", "status", "--porcelain", "--", *paths],
                                     cwd=ROOT, text=True)
    return not output.strip()


def require_absent(paths) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise S3RContractError(f"no-clobber contract: outputs already exist: {existing[:5]}")


def state_hash(head) -> str:
    return hashlib.sha256(
        head.W.detach().cpu().contiguous().numpy().tobytes()).hexdigest()


def normalized(score):
    return score / torch.sqrt(torch.mean(score * score) + NORM_EPS)


class DirectW(nn.Module):
    def __init__(self, seed=SEED_PARAM):
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.W = nn.Parameter(torch.randn(D_ESM, D_ATOM, generator=generator) * 1e-3)
        self.project_norm()

    @torch.no_grad()
    def project_norm(self):
        norm = self.W.norm()
        if not torch.isfinite(norm) or float(norm) <= 0:
            raise S3RContractError("invalid direct-W norm")
        self.W.div_(norm)


def configure_determinism() -> None:
    torch.manual_seed(SEED_PARAM)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED_PARAM)
    torch.use_deterministic_algorithms(True)


def pair_rows(pairs):
    return [{"sk": sk, "a": a, "b": b} for sk, a, b in pairs]


def build_stream(pairs, construct_component):
    stream, update = [], 0
    for epoch in range(EPOCHS):
        chosen = P.hierarchical_sample(pairs, construct_component, epoch)
        by_component = defaultdict(list)
        for component, sk, pair_list in chosen:
            by_component[component].append({
                "sk": sk,
                "pairs": [[a, b] for _same_sk, a, b in pair_list],
            })
        components = sorted(by_component)
        for start in range(0, len(components), BATCH_COMPONENTS):
            update += 1
            stream.append({
                "update": update,
                "epoch": epoch,
                "components": [
                    {"component": component,
                     "constructs": by_component[component]}
                    for component in components[start:start + BATCH_COMPONENTS]
                ],
            })
    if update != N_UPDATES:
        raise S3RContractError(f"real stream has {update} updates, expected {N_UPDATES}")
    return stream


def semantic_stream_hash(stream) -> str:
    payload = json.dumps(stream, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def sanitize_record(record, split, component):
    fields = (
        "source_key", "seq_key", "graph_key", "scaffold", "n_atoms", "n_res",
        "uniprot_sequence", "ligand_ccd", "cohort", "pdb_id",
    )
    out = {key: record[key] for key in fields if key in record}
    out["split"] = split
    out["component"] = component
    return out


def prepare_views() -> dict:
    require_absent([
        OUT / "PHASE2B_S3R_INPUT_AND_FIREWALL_MANIFEST.json",
        OUT / "PHASE2B_S3R_REAL_STREAM_MANIFEST.json",
        OUT / "PHASE2B_S3R_GATE.json",
        OUT / "PHASE2B_S3R_FAIL_CLOSED.json",
    ])
    if EXEC.exists() and any(EXEC.iterdir()):
        raise S3RContractError("no-clobber contract: S3R execution directory is not empty")
    if not git_execution_surface_clean():
        raise S3RContractError("S3R execution surface is dirty or untracked")
    if (sha_file(PREREG) != PREREG_SHA or sha_file(AMENDMENT) != AMENDMENT_SHA or
            sha_file(AMENDMENT_02) != AMENDMENT_02_SHA):
        raise S3RContractError("S3R preregistration hash mismatch")
    if sha_file(CONTROL_PATH) != CONTROL_SHA:
        raise S3RContractError("frozen control-map hash mismatch")
    checks = {
        "s2r_verdict_sha": sha_file(OUT / "PHASE2B_S2R_VERDICT.json"),
        "s2r_calibration_sha": sha_file(OUT / "S2R_CALIBRATION_GATE.json"),
        "s2r_sealed_sha": sha_file(OUT / "S2R_SEALED_VERIFICATION_GATE.json"),
    }
    if checks != {"s2r_verdict_sha": S2R_VERDICT_SHA,
                  "s2r_calibration_sha": S2R_CAL_SHA,
                  "s2r_sealed_sha": S2R_SEALED_SHA}:
        raise S3RContractError("S2R authorizing artifacts changed")
    s2r = json.loads((OUT / "PHASE2B_S2R_VERDICT.json").read_text())
    if s2r.get("TERMINAL_VERDICT") != "BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED":
        raise S3RContractError("S2R did not authorize real structural transfer")

    kept, quarantine, contract, _failures = build()
    comp_of = protein_components(kept)
    train, _held_all, held_a, held_b = make_split(kept, comp_of)
    masks_all = {record["source_key"]: frozenset(i for i, _j in record["edges"])
                 for record in kept}
    construct_component = {}
    for record in kept:
        construct_component.setdefault(record["seq_key"], comp_of[record["source_key"]])
    store, mols, gvec = FeatureStore(), load_mols(), {}
    for record in kept:
        graph_key = record["graph_key"]
        if graph_key not in gvec:
            gvec[graph_key] = g_of(store.atoms(record, mols))
    train_pairs, train_excluded = build_pairs(train, masks_all)
    helda_pairs, helda_excluded = build_pairs(held_a, masks_all)
    heldb_pairs, heldb_excluded = build_pairs(held_b, masks_all)
    ctx = {
        "kept": kept, "quarantine": quarantine, "contract": contract,
        "comp_of": comp_of, "train": train, "held_A": held_a,
        "held_B": held_b, "masks": masks_all,
        "construct_component": construct_component, "gvec": gvec,
        "train_pairs": train_pairs, "heldA_pairs": helda_pairs,
        "heldB_pairs": heldb_pairs,
        "excluded": {"train": train_excluded, "heldoutA": helda_excluded,
                     "heldoutB": heldb_excluded},
    }
    control = json.loads(CONTROL_PATH.read_text(encoding="utf-8"))
    splits = {"train": ctx["train"], "heldoutA": ctx["held_A"],
              "heldoutB": ctx["held_B"]}
    pairs = {"train": ctx["train_pairs"], "heldoutA": ctx["heldA_pairs"],
             "heldoutB": ctx["heldB_pairs"]}
    masks = {}
    metadata = []
    for split, records in splits.items():
        masks[split] = {r["source_key"]: sorted(ctx["masks"][r["source_key"]])
                        for r in records}
        metadata.extend(sanitize_record(r, split, ctx["comp_of"][r["source_key"]])
                        for r in records)

    census = {
        "train_pairs": len(pairs["train"]),
        "train_components": len({ctx["construct_component"][p[0]]
                                  for p in pairs["train"]}),
        "heldoutA_pairs": len(pairs["heldoutA"]),
        "heldoutA_components": len({ctx["construct_component"][p[0]]
                                     for p in pairs["heldoutA"]}),
        "heldoutB_pairs": len(pairs["heldoutB"]),
    }
    if census != EXPECTED:
        raise S3RContractError(f"real pair census changed: {census}")

    train_components = {ctx["construct_component"][p[0]] for p in pairs["train"]}
    held_components = {ctx["construct_component"][p[0]] for p in pairs["heldoutA"]}
    train_graphs = {r["graph_key"] for r in splits["train"]}
    held_graphs = {r["graph_key"] for r in splits["heldoutA"]}
    if train_components & held_components or train_graphs & held_graphs:
        raise S3RContractError("closure or ligand graph firewall failed")

    EXEC.mkdir(parents=True, exist_ok=True)
    metadata_path = EXEC / "metadata_view.jsonl"
    write_jsonl(metadata_path, sorted(metadata, key=lambda row: row["source_key"]))
    for split in splits:
        write_gzip_json(EXEC / f"{split}_residue_masks.json.gz", masks[split])
        write_jsonl(EXEC / f"{split}_pairs.jsonl", pair_rows(pairs[split]))

    graph_keys = sorted(ctx["gvec"])
    np.save(EXEC / "gvec_f64.npy", np.stack([ctx["gvec"][key] for key in graph_keys]))
    write_json(EXEC / "gvec_index.json", graph_keys)
    write_json(EXEC / "construct_component.json", ctx["construct_component"])

    stream = build_stream(pairs["train"], ctx["construct_component"])
    stream_path = EXEC / "real_stream.jsonl"
    write_jsonl(stream_path, stream)
    reloaded = list(read_jsonl(stream_path))
    if semantic_stream_hash(reloaded) != semantic_stream_hash(stream):
        raise S3RContractError("real stream reload mismatch")

    runtime_inputs = [
        ESM / "esm2_650M_index.json", ESM / "esm2_650M_residues.fp16.dat",
        OLD_P2B / "b_prior_f64.npy", OLD_P2B / "b_prior_index.json",
        S5 / "B5_checkpoint.pt", S5 / "heldoutA_B5.f16.dat",
        S4 / "heldoutA_index.json", CONTROL_PATH, PREREG, AMENDMENT, AMENDMENT_02,
        HERE / "s3r_run.py", HERE / "p2b_residue_residual.py",
        HERE / "s1r_run.py", HERE / "s7_dataset.py", HERE / "s7_localizer.py",
        HERE / "s7_run.py", ROOT / "tests" / "test_s7_l2b_phase2b_s3r.py",
        OUT / "PHASE2B_S2R_VERDICT.json", OUT / "S2R_CALIBRATION_GATE.json",
        OUT / "S2R_SEALED_VERIFICATION_GATE.json",
    ]
    if any(not path.is_file() for path in runtime_inputs):
        missing = [str(path) for path in runtime_inputs if not path.is_file()]
        raise S3RContractError(f"missing runtime input: {missing}")
    artifacts = [metadata_path, EXEC / "train_residue_masks.json.gz",
                 EXEC / "heldoutA_residue_masks.json.gz",
                 EXEC / "heldoutB_residue_masks.json.gz",
                 EXEC / "train_pairs.jsonl", EXEC / "heldoutA_pairs.jsonl",
                 EXEC / "heldoutB_pairs.jsonl", EXEC / "gvec_f64.npy",
                 EXEC / "gvec_index.json", EXEC / "construct_component.json",
                 stream_path, *runtime_inputs]
    manifest = {
        "schema": "MetaSieve.S7L2B.P2B.S3R.InputAndFirewall.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "amendment_sha256": AMENDMENT_SHA,
        "amendment_02_sha256": AMENDMENT_02_SHA,
        "preparation_commit": git_head(),
        "execution_surface_clean": True,
        "census": census,
        "component_overlap": 0,
        "ligand_graph_overlap": 0,
        "real_structural_edge_label_reads": sum(len(v) for view in masks.values()
                                                  for v in view.values()),
        "affinity_value_reads": 0,
        "training_process_allowed_label_view": "train_residue_masks.json.gz only",
        "scoring_process_allowed_label_views": ["heldoutA_residue_masks.json.gz",
                                                  "heldoutB_residue_masks.json.gz"],
        "stream": {"updates": len(stream), "semantic_sha256": semantic_stream_hash(stream),
                   "file_sha256": sha_file(stream_path)},
        "controls_sha256": sha_file(CONTROL_PATH),
        "s2r_authorization": checks,
        "artifacts": {str(path.relative_to(ROOT)).replace("\\", "/"): sha_file(path)
                      for path in artifacts},
    }
    write_json(OUT / "PHASE2B_S3R_INPUT_AND_FIREWALL_MANIFEST.json", manifest)
    write_json(OUT / "PHASE2B_S3R_REAL_STREAM_MANIFEST.json", {
        "schema": "MetaSieve.S7L2B.P2B.S3R.Stream.v1",
        "updates": len(stream), "epochs": EPOCHS,
        "semantic_sha256": semantic_stream_hash(stream),
        "file_sha256": sha_file(stream_path),
        "candidate_repeat_and_R5_share_this_exact_stream": True,
    })
    return manifest


class RuntimeContext:
    def __init__(self, label_split: str):
        records = list(read_jsonl(EXEC / "metadata_view.jsonl"))
        self.records = {row["source_key"]: row for row in records}
        self.split_records = defaultdict(list)
        for row in records:
            self.split_records[row["split"]].append(row)
        self.pairs = list(read_jsonl(EXEC / f"{label_split}_pairs.jsonl"))
        self.masks = {key: frozenset(value) for key, value in
                      read_gzip_json(EXEC / f"{label_split}_residue_masks.json.gz").items()}
        self.label_split = label_split
        keys = json.loads((EXEC / "gvec_index.json").read_text())
        matrix = np.load(EXEC / "gvec_f64.npy", mmap_mode="r")
        self.gvec = {key: np.asarray(matrix[i], dtype=np.float64) for i, key in enumerate(keys)}
        self.construct_component = json.loads((EXEC / "construct_component.json").read_text())
        self.esm_index = json.loads((ESM / "esm2_650M_index.json").read_text())
        total = max(offset + length for offset, length in self.esm_index.values())
        self.esm = np.memmap(ESM / "esm2_650M_residues.fp16.dat", dtype=np.float16,
                             mode="r", shape=(total, D_ESM))
        self.prior_index = json.loads((OLD_P2B / "b_prior_index.json").read_text())
        self.prior = np.load(OLD_P2B / "b_prior_f64.npy", mmap_mode="r")
        needed = {row["sk"] for row in self.pairs}
        self.Qs = {sk: nuisance_basis(self.b_prior(sk)) for sk in needed}
        self.control = json.loads(CONTROL_PATH.read_text(encoding="utf-8"))
        self._hcache = {}

    def h(self, sk):
        if sk not in self._hcache:
            offset, length = self.esm_index[sk]
            self._hcache[sk] = np.asarray(self.esm[offset:offset + length],
                                          dtype=np.float32)
            if len(self._hcache) > 96:
                self._hcache.pop(next(iter(self._hcache)))
        return self._hcache[sk]

    def b_prior(self, sk):
        offset, length = self.prior_index[sk]
        return np.asarray(self.prior[offset:offset + length], dtype=np.float64)


def validate_runtime_inputs(ctx: RuntimeContext, expected_split: str) -> None:
    if ctx.label_split != expected_split:
        raise S3RContractError("wrong label view mounted")
    manifest = json.loads((OUT / "PHASE2B_S3R_INPUT_AND_FIREWALL_MANIFEST.json").read_text())
    for relative, expected in manifest["artifacts"].items():
        if expected_split == "train" and ("heldoutA_" in relative or
                                           "heldoutB_" in relative):
            continue
        if sha_file(ROOT / relative) != expected:
            raise S3RContractError(f"input artifact changed: {relative}")


def labels_for(ctx: RuntimeContext, a, b, permuted=False):
    if permuted:
        mapping = ctx.control["within_construct_derangement"]
        a, b = mapping.get(a, a), mapping.get(b, b)
    left, right = ctx.masks[a], ctx.masks[b]
    return left - right, right - left


def train(kind: str) -> dict:
    directory = EXEC / kind
    checkpoint = directory / "final_checkpoint.pt"
    result_path = OUT / f"PHASE2B_S3R_TRAIN_{kind.upper()}.json"
    require_absent([checkpoint, result_path])
    configure_determinism()
    ctx = RuntimeContext("train")
    validate_runtime_inputs(ctx, "train")
    if len(LABEL_READ_LOG) != 1 or not LABEL_READ_LOG[0].endswith(
            "train_residue_masks.json.gz"):
        raise S3RContractError(f"training label-view firewall failed: {LABEL_READ_LOG}")
    stream = list(read_jsonl(EXEC / "real_stream.jsonl"))
    stream_manifest = json.loads((OUT / "PHASE2B_S3R_REAL_STREAM_MANIFEST.json").read_text())
    if semantic_stream_hash(stream) != stream_manifest["semantic_sha256"]:
        raise S3RContractError("training stream changed")

    head = DirectW().to(DEVICE)
    initial = head.W.detach().cpu().clone()
    optimizer = torch.optim.Adam(head.parameters(), lr=LR)
    trace, sampled, near_zero_presentations = [], [], 0
    permuted = kind == "permuted"
    for item in stream:
        optimizer.zero_grad(set_to_none=True)
        component_losses = []
        for component in item["components"]:
            construct_losses = []
            for construct in component["constructs"]:
                sk = construct["sk"]
                h = torch.from_numpy(ctx.h(sk).copy()).float().to(DEVICE)
                z = h @ head.W
                q = torch.from_numpy(ctx.Qs[sk]).double().to(DEVICE)
                pair_losses = []
                for a, b in construct["pairs"]:
                    ga = ctx.gvec[ctx.records[a]["graph_key"]]
                    gb = ctx.gvec[ctx.records[b]["graph_key"]]
                    gd = torch.from_numpy((ga - gb).copy()).float().to(DEVICE)
                    raw = (z @ gd).double()
                    score = raw - q @ (q.T @ raw)
                    if float(torch.sqrt(torch.mean(score * score))) <= 1e-12:
                        near_zero_presentations += 1
                    gain, loss = labels_for(ctx, a, b, permuted=permuted)
                    pair_losses.append(rank_loss_torch(normalized(score), gain, loss))
                    sampled.append(f"{item['update']}|{component['component']}|{sk}|{a}|{b}")
                if pair_losses:
                    construct_losses.append(torch.stack(pair_losses).mean())
            if construct_losses:
                component_losses.append(torch.stack(construct_losses).mean())
        if not component_losses:
            raise S3RContractError("empty training update")
        loss = torch.stack(component_losses).mean()
        loss.backward()
        grad = float(head.W.grad.norm())
        if not np.isfinite(grad) or grad <= 0:
            raise S3RContractError("nonfinite or zero W gradient")
        torch.nn.utils.clip_grad_norm_(head.parameters(), CLIP)
        optimizer.step()
        head.project_norm()
        norm = float(head.W.detach().norm())
        if abs(norm - 1.0) > 1e-5:
            raise S3RContractError("unit-norm contract failed")
        trace.append({"update": item["update"], "epoch": item["epoch"],
                      "loss": float(loss.detach()), "grad_W": grad, "W_norm": norm})
        if item["update"] % 35 == 0:
            print(f"[{kind}] update {item['update']}/{N_UPDATES} loss={trace[-1]['loss']:.6f}",
                  flush=True)

    if len(trace) != N_UPDATES:
        raise S3RContractError("wrong number of optimizer updates")
    movement = float((head.W.detach().cpu() - initial).norm() / initial.norm())
    directory.mkdir(parents=True, exist_ok=True)
    torch.save({"model": head.state_dict(), "optimizer": optimizer.state_dict(),
                "trace": trace, "movement": movement}, checkpoint)
    sampled_path = directory / "sampled_pairs.txt"
    sampled_path.write_text("\n".join(sampled) + "\n", encoding="utf-8")
    result = {
        "kind": kind, "device": DEVICE, "execution_commit": git_head(),
        "checkpoint": str(checkpoint.relative_to(ROOT)).replace("\\", "/"),
        "checkpoint_sha256": sha_file(checkpoint), "W_state_sha256": state_hash(head),
        "sampled_pairs_sha256": sha_file(sampled_path),
        "stream_semantic_sha256": semantic_stream_hash(stream),
        "updates": len(trace), "relative_W_movement": movement,
        "min_grad_W": min(row["grad_W"] for row in trace),
        "final_loss": trace[-1]["loss"], "final_W_norm": trace[-1]["W_norm"],
        "near_zero_score_presentations_retained": near_zero_presentations,
        "label_view_paths_opened": list(LABEL_READ_LOG),
        "label_view_opened": "train", "heldout_label_view_reads": 0,
        "affinity_value_reads": 0,
    }
    write_json(result_path, result)
    write_json(directory / "trace.json", trace)
    return result


def load_head(kind: str) -> DirectW:
    checkpoint = EXEC / kind / "final_checkpoint.pt"
    head = DirectW().to("cpu")
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    head.load_state_dict(payload["model"])
    return head


def delta(ctx: RuntimeContext, W, sk, gk_a, gk_b, h_override=None):
    h = ctx.h(sk) if h_override is None else h_override
    raw = (h.astype(np.float64) @ W) @ (ctx.gvec[gk_a] - ctx.gvec[gk_b])
    return project_np(ctx.Qs[sk], raw)


def evaluate_arm(ctx, W, mode, b5mean=None, score_prefix=None):
    values, meta, scores_flat, score_index = {}, {}, [], []
    offset, near_zero = 0, 0
    by_sk = defaultdict(list)
    for row in ctx.pairs:
        by_sk[row["sk"]].append(row)
    for sk in sorted(by_sk):
        context_h = None
        if mode == "context":
            record = ctx.records[by_sk[sk][0]["a"]]
            context_h = context_shuffle(ctx.h(sk), record["uniprot_sequence"],
                                        int(hashlib.sha256(sk.encode()).hexdigest()[:8], 16))
        z = None
        if mode != "b5diff":
            h = ctx.h(sk) if context_h is None else context_h
            z = h.astype(np.float64) @ W
        for row in by_sk[sk]:
            a, b = row["a"], row["b"]
            pid = f"{a}|{b}"
            ra, rb = ctx.records[a], ctx.records[b]
            if mode == "b5diff":
                score = b5mean[a] - b5mean[b]
            else:
                if mode == "foreign":
                    gk_a, gk_b = ctx.control["foreign_pair_map"][pid]
                elif mode == "chem":
                    mapping = ctx.control["within_construct_derangement"]
                    pa, pb = mapping.get(a, a), mapping.get(b, b)
                    gk_a = ctx.records[pa]["graph_key"]
                    gk_b = ctx.records[pb]["graph_key"]
                else:
                    gk_a, gk_b = ra["graph_key"], rb["graph_key"]
                score = project_np(ctx.Qs[sk], z @ (ctx.gvec[gk_a] - ctx.gvec[gk_b]))
            if not np.isfinite(score).all():
                raise S3RContractError(f"nonfinite prediction in {mode}")
            if float(np.sqrt(np.mean(score * score))) <= 1e-12:
                near_zero += 1
            gain, loss = labels_for(ctx, a, b)
            metrics = pair_metrics(score, gain, loss, int(ra["n_res"]))
            if metrics is None:
                raise S3RContractError("eligible pair lost its real label difference")
            values[pid], meta[pid] = metrics, sk
            if score_prefix is not None:
                scores_flat.append(np.asarray(score, dtype=np.float64))
                score_index.append({"pair_id": pid, "sk": sk, "offset": offset,
                                    "length": int(score.size)})
                offset += int(score.size)
    artifacts = {}
    if score_prefix is not None:
        path = EXEC / f"{score_prefix}_scores_f64.npy"
        index = EXEC / f"{score_prefix}_score_index.jsonl"
        np.save(path, np.concatenate(scores_flat))
        write_jsonl(index, score_index)
        artifacts = {"scores": sha_file(path), "index": sha_file(index)}
    return values, meta, artifacts, {"near_zero_pairs_retained": near_zero}


def component_field(values, metadata, construct_component, field):
    return aggregate({key: value[field] for key, value in values.items()}, metadata,
                     construct_component)


def assert_common_masks(arms, construct_component):
    reference_values, reference_meta = arms["candidate"]
    reference = set(reference_values)
    for name, (values, _metadata) in arms.items():
        if set(values) != reference:
            raise S3RContractError(f"paired Gate mask differs for arm {name}")
        if _metadata != reference_meta:
            raise S3RContractError(f"pair-to-construct mapping differs for arm {name}")
    rows = [f"{pair_id}|{reference_meta[pair_id]}|{construct_component[reference_meta[pair_id]]}"
            for pair_id in sorted(reference)]
    return hashlib.sha256("\n".join(rows).encode()).hexdigest()


def score() -> dict:
    require_absent([OUT / "PHASE2B_S3R_GATE.json",
                    OUT / "PHASE2B_S3R_COMPONENT_TABLES.json",
                    EXEC / "candidate_heldoutA_scores_f64.npy",
                    EXEC / "repeat_heldoutA_scores_f64.npy"])
    configure_determinism()
    ctx = RuntimeContext("heldoutA")
    validate_runtime_inputs(ctx, "heldoutA")
    if len(LABEL_READ_LOG) != 1 or not LABEL_READ_LOG[0].endswith(
            "heldoutA_residue_masks.json.gz"):
        raise S3RContractError(f"scoring label-view firewall failed: {LABEL_READ_LOG}")
    candidate = load_head("candidate")
    repeat = load_head("repeat")
    permuted = load_head("permuted")
    W = candidate.W.detach().double().numpy()
    W_repeat = repeat.W.detach().double().numpy()
    W_perm = permuted.W.detach().double().numpy()

    sidx = json.loads((S4 / "heldoutA_index.json").read_text())
    total = max(offset + length * atoms for offset, length, atoms in sidx.values())
    b5 = np.memmap(S5 / "heldoutA_B5.f16.dat", dtype=np.float16, mode="r", shape=(total,))
    b5mean = {}
    for key, (offset, length, atoms) in sidx.items():
        b5mean[key] = np.asarray(b5[offset:offset + length * atoms], dtype=np.float64).reshape(
            length, atoms).mean(1)

    arms = {}
    candidate_eval = evaluate_arm(ctx, W, "candidate", score_prefix="candidate_heldoutA")
    arms["candidate"] = candidate_eval[:2]
    arms["b5diff"] = evaluate_arm(ctx, W, "b5diff", b5mean=b5mean)[:2]
    arms["foreign"] = evaluate_arm(ctx, W, "foreign")[:2]
    arms["context"] = evaluate_arm(ctx, W, "context")[:2]
    arms["permuted"] = evaluate_arm(ctx, W_perm, "candidate")[:2]
    arms["chemistry_shuffle"] = evaluate_arm(ctx, W, "chem")[:2]
    zero_eval = evaluate_arm(ctx, np.zeros_like(W), "candidate")
    arms["zero_W"] = zero_eval[:2]
    common_mask_sha = assert_common_masks({key: arms[key]
                                           for key in ("candidate", "b5diff", "foreign",
                                                       "context", "permuted")},
                                          ctx.construct_component)

    repeat_eval = evaluate_arm(ctx, W_repeat, "candidate",
                               score_prefix="repeat_heldoutA")
    repeat_values, repeat_meta = repeat_eval[:2]
    if set(repeat_values) != set(arms["candidate"][0]):
        raise S3RContractError("repeat prediction mask differs")
    max_repeat_metric_diff = max(
        abs(repeat_values[key]["ap_bidir"] - arms["candidate"][0][key]["ap_bidir"])
        for key in repeat_values)
    max_repeat_W_diff = float(np.max(np.abs(W - W_repeat)))
    candidate_scores = np.load(EXEC / "candidate_heldoutA_scores_f64.npy", mmap_mode="r")
    repeat_scores = np.load(EXEC / "repeat_heldoutA_scores_f64.npy", mmap_mode="r")
    if candidate_scores.shape != repeat_scores.shape:
        raise S3RContractError("repeat raw prediction shape differs")
    max_repeat_prediction_diff = float(np.max(np.abs(candidate_scores - repeat_scores)))

    component, macro = {}, {}
    for name, (values, metadata) in arms.items():
        for field in ("ap_bidir", "chance_bidir", "ap_gain", "ap_loss"):
            component[(name, field)], macro[(name, field)] = component_field(
                values, metadata, ctx.construct_component, field)

    def gate(left, right, margin, left_field="ap_bidir", right_field=None):
        result = component_bootstrap(component[(left, left_field)],
                                     component[(right, right_field or left_field)],
                                     n_boot=10000, seed=SEED_BOOT)
        result["margin"] = margin
        result["pass"] = bool(result["delta"] >= margin and
                              result["lcb95_one_sided"] > 0)
        return result

    gates = {
        "R1_vs_chance": gate("candidate", "candidate", 0.05,
                              "ap_bidir", "chance_bidir"),
        "R2_vs_frozen_B5_differential": gate("candidate", "b5diff", 0.03),
        "R3_vs_two_ligand_foreign_pair": gate("candidate", "foreign", 0.03),
        "R4_vs_residue_context_corruption": gate("candidate", "context", 0.03),
        "R5_vs_trained_permuted_label_learner": gate("candidate", "permuted", 0.05),
    }

    candidate_train = json.loads((OUT / "PHASE2B_S3R_TRAIN_CANDIDATE.json").read_text())
    repeat_train = json.loads((OUT / "PHASE2B_S3R_TRAIN_REPEAT.json").read_text())
    perm_train = json.loads((OUT / "PHASE2B_S3R_TRAIN_PERMUTED.json").read_text())
    stream_same = len({candidate_train["stream_semantic_sha256"],
                       repeat_train["stream_semantic_sha256"],
                       perm_train["stream_semantic_sha256"]}) == 1
    score_variance = float(np.var(candidate_scores))
    chance_macro = macro[("candidate", "chance_bidir")]
    module = {
        "min_grad_W": candidate_train["min_grad_W"],
        "grad_nonzero": candidate_train["min_grad_W"] > 0,
        "relative_W_movement": candidate_train["relative_W_movement"],
        "movement_pass": candidate_train["relative_W_movement"] >= MIN_PARAM_MOVEMENT,
        "W_norm": float(np.linalg.norm(W)),
        "unit_norm_pass": abs(np.linalg.norm(W) - 1.0) <= 1e-5,
        "raw_score_variance": score_variance,
        "variance_pass": score_variance >= 1e-8,
        "zero_W_ap": macro[("zero_W", "ap_bidir")],
        "zero_W_chance": chance_macro,
        "zero_W_chance_pass": abs(macro[("zero_W", "ap_bidir")] - chance_macro) <= 0.005,
        "context_degrades": macro[("context", "ap_bidir")] < macro[("candidate", "ap_bidir")],
        "foreign_ligands_degrade": macro[("foreign", "ap_bidir")] < macro[("candidate", "ap_bidir")],
        "same_stream_hash": stream_same,
        "repeat_W_max_abs_diff": max_repeat_W_diff,
        "repeat_prediction_max_abs_diff": max_repeat_prediction_diff,
        "repeat_metric_max_abs_diff": max_repeat_metric_diff,
        "repeat_pass": (max_repeat_W_diff <= 1e-7 and
                         max_repeat_prediction_diff <= 1e-7 and
                         max_repeat_metric_diff <= 1e-7),
        "candidate_prediction_sha256": sha_file(EXEC / "candidate_heldoutA_scores_f64.npy"),
        "repeat_prediction_sha256": sha_file(EXEC / "repeat_heldoutA_scores_f64.npy"),
    }
    module["pass"] = bool(module["grad_nonzero"] and module["movement_pass"] and
                           module["unit_norm_pass"] and module["variance_pass"] and
                           module["zero_W_chance_pass"] and module["context_degrades"] and
                           module["foreign_ligands_degrade"] and module["same_stream_hash"] and
                           module["repeat_pass"])

    if not module["pass"]:
        verdict = "S3R_DIRECT_W_TRAINING_OR_PARTICIPATION_FAILED"
    elif not gates["R1_vs_chance"]["pass"]:
        verdict = "REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED"
    elif not gates["R2_vs_frozen_B5_differential"]["pass"]:
        verdict = "NO_INCREMENT_BEYOND_FROZEN_B5_DIFFERENTIAL"
    elif not all(gates[name]["pass"] for name in (
            "R3_vs_two_ligand_foreign_pair", "R4_vs_residue_context_corruption",
            "R5_vs_trained_permuted_label_learner")):
        verdict = "REAL_RESIDUE_STATISTIC_SHORTCUT_DEPENDENT"
    else:
        verdict = "STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT"

    u, singular, vt = np.linalg.svd(W, full_matrices=False)
    W8 = (u[:, :8] * singular[:8]) @ vt[:8]
    rank8_values, rank8_meta = evaluate_arm(ctx, W8, "candidate")[:2]
    _rank8_component, rank8_macro = component_field(rank8_values, rank8_meta,
                                                    ctx.construct_component, "ap_bidir")

    heldout_b = {"status": "NOT_OPENED_PRIMARY_DID_NOT_PASS", "non_gating": True}
    if verdict == "STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT":
        ctx_b = RuntimeContext("heldoutB")
        validate_runtime_inputs(ctx_b, "heldoutB")
        if len(LABEL_READ_LOG) != 2 or not LABEL_READ_LOG[-1].endswith(
                "heldoutB_residue_masks.json.gz"):
            raise S3RContractError(f"heldout-B label-view firewall failed: {LABEL_READ_LOG}")
        heldb_values, heldb_meta = evaluate_arm(ctx_b, W, "candidate")[:2]
        _heldb_component, heldb_macro = component_field(
            heldb_values, heldb_meta, ctx_b.construct_component, "ap_bidir")
        _heldb_chance_component, heldb_chance = component_field(
            heldb_values, heldb_meta, ctx_b.construct_component, "chance_bidir")
        heldout_b = {"status": "OPENED_AFTER_PRIMARY_PASS", "ap_bidir": heldb_macro,
                     "chance_bidir": heldb_chance,
                     "candidate_minus_chance": heldb_macro - heldb_chance,
                     "sign_reversal": bool(heldb_macro - heldb_chance < 0),
                     "pairs": len(heldb_values), "non_gating": True}

    component_tables = {
        name: {field: component[(name, field)]
               for field in ("ap_bidir", "chance_bidir", "ap_gain", "ap_loss")}
        for name in arms
    }
    write_json(OUT / "PHASE2B_S3R_COMPONENT_TABLES.json", component_tables)
    write_json(OUT / "PHASE2B_S3R_MODULE_PARTICIPATION.json", module)
    result = {
        "schema": "MetaSieve.S7L2B.P2B.S3R.Gate.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA, "amendment_sha256": AMENDMENT_SHA,
        "amendment_02_sha256": AMENDMENT_02_SHA,
        "primary_panel": {"pairs": len(arms["candidate"][0]),
                          "components": len(component[("candidate", "ap_bidir")]),
                          "common_mask_sha256": common_mask_sha},
        "macro_ap_bidir": {name: macro[(name, "ap_bidir")] for name in arms},
        "gates": gates, "module_participation": module,
        "rank8_secondary": {"ap_bidir": rank8_macro,
                            "energy_fraction": float((singular[:8] ** 2).sum() /
                                                     (singular ** 2).sum()),
                            "non_gating": True},
        "heldoutB_secondary": heldout_b,
        "near_zero_pairs_retained": candidate_eval[3]["near_zero_pairs_retained"],
        "label_view_paths_opened": list(LABEL_READ_LOG),
        "R6_status": "NOT_RUN_ABSOLUTE_SCALE_ORIGIN_AND_DIFFERENCE_NULLSPACE_NOT_IDENTIFIED",
        "TERMINAL_VERDICT": verdict,
        "affinity_value_reads": 0,
        "claims_not_made": ["exact residue-atom coupling", "physical energy",
                            "affinity", "selectivity", "few-shot section",
                            "biological z", "validated end-to-end DTA"],
        "authorized_next_action": (
            "preregister S3R-I0 positive-scalar B5 integration and an independent structural confirmation design"
            if (verdict == "STRUCTURAL_LIGAND_CONDITIONED_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT"
                and not heldout_b.get("sign_reversal", False))
            else "none; stop at the earliest failed boundary or heldout-B escalation block"),
    }
    write_json(OUT / "PHASE2B_S3R_GATE.json", result)
    report = ["# Phase 2B S3R real structural ordinal result", "",
              f"Terminal verdict: `{verdict}`", "",
              "| arm | component-macro AP_bidir |", "|---|---:|"]
    report.extend(f"| {name} | {value:.6f} |" for name, value in result["macro_ap_bidir"].items())
    report += ["", "| Gate | delta | LCB95 | PASS |", "|---|---:|---:|:---:|"]
    report.extend(f"| {name} | {gate['delta']:.6f} | {gate['lcb95_one_sided']:.6f} | {gate['pass']} |"
                  for name, gate in gates.items())
    report += ["", "R6 was not opened because pairwise ordinal training does not identify",
               "absolute amplitude, ligand-feature origin or the difference-design nullspace.",
               "No affinity value was read and the frozen law operator was not modified."]
    (OUT / "PHASE2B_S3R_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return result


def run_all() -> int:
    prepare_views()
    script = str(Path(__file__).resolve())
    for mode in ("train-candidate", "train-repeat", "train-permuted", "score"):
        subprocess.run([sys.executable, script, mode], cwd=ROOT, check=True)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "train-candidate", "train-repeat",
                                         "train-permuted", "score", "all"))
    args = parser.parse_args(argv)
    started = time.time()
    try:
        if args.mode == "prepare":
            output = prepare_views()
        elif args.mode == "train-candidate":
            output = train("candidate")
        elif args.mode == "train-repeat":
            output = train("repeat")
        elif args.mode == "train-permuted":
            output = train("permuted")
        elif args.mode == "score":
            output = score()
        else:
            return run_all()
        print(json.dumps({"mode": args.mode, "result": output,
                          "elapsed_seconds": round(time.time() - started, 3)},
                         indent=2, default=str), flush=True)
        return 0
    except Exception as exc:
        failure = {
            "schema": "MetaSieve.S7L2B.P2B.S3R.FailClosed.v1",
            "created_utc": "2026-08-10", "mode": args.mode,
            "error_type": type(exc).__name__, "error": str(exc),
            "TERMINAL_VERDICT": "S3R_CONTRACT_OR_LABEL_FIREWALL_FAIL_CLOSED",
            "affinity_value_reads": 0,
        }
        write_json(OUT / "PHASE2B_S3R_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
