"""S4R single-axis graph-aware ligand representation repair of the S3R stage.

Registered by
  research/s7_l2b_r0r/PREREG_PHASE2B_S4R_GRAPH_AWARE_LIGAND_DIRECT_W.md
  (sha256 b0630cc3..., commit 123ed22), authorized by
  research/s7_l2b_r0r/PREREG_PHASE2B_S4R_LIGAND_REPRESENTATION_AUDIT.md
  and its amendment 01, both committed BEFORE this file existed.

Exactly one axis differs from `s3r_run.py`: the ligand statistic. The protein
branch, gauge, estimator, loss, sampler, split, seeds, optimizer, 210-update
budget, control maps and R1-R5 margins are byte-identical. Heldout-B is never
created and never read. No affinity value is opened.
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
    SEED_BOOT, SEED_PARAM, S4, S5, aggregate, build_pairs, component_bootstrap,
    context_shuffle, g_of, nuisance_basis, pair_metrics, project_np, sha_file,
)
from s1r_run import rank_loss_torch  # noqa: E402
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import FeatureStore  # noqa: E402
from s7_run import load_mols  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S4R_GRAPH_AWARE_LIGAND_DIRECT_W.md"
AUDIT_PREREG = HERE / "PREREG_PHASE2B_S4R_LIGAND_REPRESENTATION_AUDIT.md"
AUDIT_AMENDMENT = HERE / "PREREG_PHASE2B_S4R_AUDIT_AMENDMENT_01.md"
PREREG_SHA = "b0630cc34beeba874de78f25fc269b9c2a45b62bb5286e9c7b294c346b670ac3"
AUDIT_PREREG_SHA = "8c3be16973957c6d1e7e735a7c3214d1d7e4f3b5d59f791e3f72271894130138"
AUDIT_AMENDMENT_SHA = "5210197fee00d7b15288f89925d8f282791305bd35b082f17bb5bd38b96745f4"
CONTROL_SHA = "e187a5f00f0b66328877bacd93b22471fe607e382e811f2674ecfc4a9dec9c33"

# ---- frozen S4R contract (preregistration sections 3, 4.1, 4.2, 5, 6) ------
D_LIG = 128
MORGAN_RADIUS = 1
VOCAB_SHA = "a200a4b986af1850fdb1d244f2e002c9b5ae707a114d8a3635053edb215ed877"
S3R_STREAM_SEMANTIC_SHA = "4bc68d54884437ded999fbb5f8fc8997b47b456f2bdde4e0fbafd9df3dcdc3ef"
S3R_CANDIDATE_MACRO_AP = 0.03588006089257408
BASELINE_REPLICATION_TOL = 1e-12
DEGENERACY_TOL = 5e-3

OUT = ROOT / "report" / "s7_l2b_r0r"
PROC = ROOT / "dataset" / "processed" / "s7_l2b_r0r"
EXEC = PROC / "phase2b_s4r"
AUDIT_EXEC = PROC / "phase2b_s4r_audit"
ESM = PROC / "esm2_650M"
OLD_P2B = PROC / "phase2b"
CONTROL_PATH = OLD_P2B / "control_maps.json"
NORM_EPS = 1e-12
N_UPDATES = EPOCHS * 35
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
EXPECTED = {
    "train_pairs": 226765,
    "train_components": 554,
    "heldoutA_pairs": 46818,
    "heldoutA_components": 112,
}
TRAINED_ARMS = ("candidate", "repeat", "permuted", "baseline41")
ARM_DIM = {"candidate": D_LIG, "repeat": D_LIG, "permuted": D_LIG,
           "baseline41": D_ATOM}
LABEL_READ_LOG = []


class S4RContractError(RuntimeError):
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
    paths = [str(PREREG.relative_to(ROOT)), str(AUDIT_PREREG.relative_to(ROOT)),
             str(AUDIT_AMENDMENT.relative_to(ROOT)),
             "research/s7_l2b_r0r/s4r_run.py",
             "research/s7_l2b_r0r/s4r_audit.py",
             "tests/test_s7_l2b_phase2b_s4r.py"]
    output = subprocess.check_output(["git", "status", "--porcelain", "--", *paths],
                                     cwd=ROOT, text=True)
    return not output.strip()


def require_absent(paths) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise S4RContractError(f"no-clobber contract: outputs already exist: {existing[:5]}")


def state_hash(head) -> str:
    return hashlib.sha256(
        head.W.detach().cpu().contiguous().numpy().tobytes()).hexdigest()


def normalized(score):
    return score / torch.sqrt(torch.mean(score * score) + NORM_EPS)


class DirectW(nn.Module):
    """The only trainable object: one matrix, no bias, unit Frobenius norm."""

    def __init__(self, d_lig: int, seed: int = SEED_PARAM):
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.W = nn.Parameter(torch.randn(D_ESM, d_lig, generator=generator) * 1e-3)
        self.project_norm()

    @torch.no_grad()
    def project_norm(self):
        norm = self.W.norm()
        if not torch.isfinite(norm) or float(norm) <= 0:
            raise S4RContractError("invalid direct-W norm")
        self.W.div_(norm)


def configure_determinism() -> None:
    torch.manual_seed(SEED_PARAM)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED_PARAM)
    torch.use_deterministic_algorithms(True)


# ---------------------------------------------- frozen graph-aware statistic
def load_vocabulary() -> list:
    path = AUDIT_EXEC / "selected_ligand_vocabulary.json"
    if not path.is_file():
        raise S4RContractError("frozen S4R-A ligand vocabulary is missing")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (payload["radius"] != MORGAN_RADIUS or
            payload["vocabulary_size"] != D_LIG or
            len(payload["vocabulary"]) != D_LIG):
        raise S4RContractError("ligand vocabulary does not match the registration")
    digest = hashlib.sha256(
        json.dumps(payload["vocabulary"], separators=(",", ":")).encode()).hexdigest()
    if digest != VOCAB_SHA:
        raise S4RContractError("frozen ligand vocabulary hash mismatch")
    return [int(x) for x in payload["vocabulary"]]


def g_graph_of(mol, vocabulary_index: dict, n_heavy: int) -> np.ndarray:
    """Mean over heavy atoms of the radius-1 Morgan neighbourhood one-hot
    descriptor, truncated to the frozen train-derived vocabulary. Same pooling
    and same size normalization as `g_of`; only the descriptor is graph-aware."""
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS)
    working = Chem.Mol(mol)
    Chem.SanitizeMol(working)
    out = np.zeros(D_LIG, dtype=np.float64)
    counts = generator.GetSparseCountFingerprint(working).GetNonzeroElements()
    scale = float(max(n_heavy, 1))
    for environment, count in counts.items():
        column = vocabulary_index.get(int(environment))
        if column is not None:
            out[column] = count / scale
    return out


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
        raise S4RContractError(f"real stream has {update} updates, expected {N_UPDATES}")
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


# --------------------------------------------------------------- preparation
def prepare_views() -> dict:
    require_absent([
        OUT / "PHASE2B_S4R_INPUT_AND_FIREWALL_MANIFEST.json",
        OUT / "PHASE2B_S4R_REAL_STREAM_MANIFEST.json",
        OUT / "PHASE2B_S4R_GATE.json",
        OUT / "PHASE2B_S4R_FAIL_CLOSED.json",
    ])
    if EXEC.exists() and any(EXEC.iterdir()):
        raise S4RContractError("no-clobber contract: S4R execution directory is not empty")
    if not git_execution_surface_clean():
        raise S4RContractError("S4R execution surface is dirty or untracked")
    if (sha_file(PREREG) != PREREG_SHA or
            sha_file(AUDIT_PREREG) != AUDIT_PREREG_SHA or
            sha_file(AUDIT_AMENDMENT) != AUDIT_AMENDMENT_SHA):
        raise S4RContractError("S4R preregistration hash mismatch")
    if sha_file(CONTROL_PATH) != CONTROL_SHA:
        raise S4RContractError("frozen control-map hash mismatch")
    audit = json.loads((OUT / "PHASE2B_S4R_REPRESENTATION_AUDIT.json").read_text())
    if audit["TERMINAL_VERDICT"] != "GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE":
        raise S4RContractError("S4R-A did not authorize a training stage")
    if (audit["selected"]["radius"] != MORGAN_RADIUS or
            audit["selected"]["vocabulary_size"] != D_LIG or
            audit["selected"]["vocabulary_sha256"] != VOCAB_SHA):
        raise S4RContractError("S4R-A selection does not match the registration")
    vocabulary = load_vocabulary()
    vocabulary_index = {environment: i for i, environment in enumerate(vocabulary)}

    kept, quarantine, contract, _failures = build()
    comp_of = protein_components(kept)
    train, _held_all, held_a, _held_b = make_split(kept, comp_of)
    masks_all = {record["source_key"]: frozenset(i for i, _j in record["edges"])
                 for record in kept}
    construct_component = {}
    for record in kept:
        construct_component.setdefault(record["seq_key"], comp_of[record["source_key"]])

    store, mols = FeatureStore(), load_mols()
    gvec, gvec_graph = {}, {}
    for record in kept:
        graph_key = record["graph_key"]
        if graph_key in gvec:
            continue
        atom_matrix = store.atoms(record, mols)
        gvec[graph_key] = g_of(atom_matrix)
        gvec_graph[graph_key] = g_graph_of(mols[record["ligand_ccd"]],
                                           vocabulary_index,
                                           int(atom_matrix.shape[0]))
    for name, table, width in (("gvec", gvec, D_ATOM), ("gvec_graph", gvec_graph, D_LIG)):
        stacked = np.stack([table[key] for key in sorted(table)])
        if stacked.shape[1] != width or not np.isfinite(stacked).all():
            raise S4RContractError(f"invalid {name} matrix")

    train_pairs, train_excluded = build_pairs(train, masks_all)
    helda_pairs, helda_excluded = build_pairs(held_a, masks_all)

    splits = {"train": train, "heldoutA": held_a}
    pairs = {"train": train_pairs, "heldoutA": helda_pairs}
    masks, metadata = {}, []
    for split, records in splits.items():
        masks[split] = {r["source_key"]: sorted(masks_all[r["source_key"]])
                        for r in records}
        metadata.extend(sanitize_record(r, split, comp_of[r["source_key"]])
                        for r in records)

    census = {
        "train_pairs": len(pairs["train"]),
        "train_components": len({construct_component[p[0]] for p in pairs["train"]}),
        "heldoutA_pairs": len(pairs["heldoutA"]),
        "heldoutA_components": len({construct_component[p[0]]
                                    for p in pairs["heldoutA"]}),
    }
    if census != EXPECTED:
        raise S4RContractError(f"real pair census changed: {census}")

    train_components = {construct_component[p[0]] for p in pairs["train"]}
    held_components = {construct_component[p[0]] for p in pairs["heldoutA"]}
    train_graphs = {r["graph_key"] for r in splits["train"]}
    held_graphs = {r["graph_key"] for r in splits["heldoutA"]}
    if train_components & held_components or train_graphs & held_graphs:
        raise S4RContractError("closure or ligand graph firewall failed")

    EXEC.mkdir(parents=True, exist_ok=True)
    metadata_path = EXEC / "metadata_view.jsonl"
    write_jsonl(metadata_path, sorted(metadata, key=lambda row: row["source_key"]))
    for split in splits:
        write_gzip_json(EXEC / f"{split}_residue_masks.json.gz", masks[split])
        write_jsonl(EXEC / f"{split}_pairs.jsonl", pair_rows(pairs[split]))

    graph_keys = sorted(gvec)
    np.save(EXEC / "gvec_f64.npy", np.stack([gvec[key] for key in graph_keys]))
    np.save(EXEC / "gvec_graph_f64.npy",
            np.stack([gvec_graph[key] for key in graph_keys]))
    write_json(EXEC / "gvec_index.json", graph_keys)
    write_json(EXEC / "construct_component.json", construct_component)

    stream = build_stream(pairs["train"], construct_component)
    stream_path = EXEC / "real_stream.jsonl"
    write_jsonl(stream_path, stream)
    reloaded = list(read_jsonl(stream_path))
    semantic = semantic_stream_hash(stream)
    if semantic_stream_hash(reloaded) != semantic:
        raise S4RContractError("real stream reload mismatch")
    if semantic != S3R_STREAM_SEMANTIC_SHA:
        raise S4RContractError(
            "stream identity contract failed: a surface other than the ligand "
            f"representation moved ({semantic})")

    runtime_inputs = [
        ESM / "esm2_650M_index.json", ESM / "esm2_650M_residues.fp16.dat",
        OLD_P2B / "b_prior_f64.npy", OLD_P2B / "b_prior_index.json",
        S5 / "B5_checkpoint.pt", S5 / "heldoutA_B5.f16.dat",
        S4 / "heldoutA_index.json", CONTROL_PATH,
        AUDIT_EXEC / "selected_ligand_vocabulary.json",
        PREREG, AUDIT_PREREG, AUDIT_AMENDMENT,
        HERE / "s4r_run.py", HERE / "s4r_audit.py",
        HERE / "p2b_residue_residual.py", HERE / "s1r_run.py",
        HERE / "s7_dataset.py", HERE / "s7_localizer.py", HERE / "s7_run.py",
        ROOT / "tests" / "test_s7_l2b_phase2b_s4r.py",
        OUT / "PHASE2B_S4R_REPRESENTATION_AUDIT.json",
        OUT / "PHASE2B_S3R_GATE.json",
    ]
    missing = [str(path) for path in runtime_inputs if not path.is_file()]
    if missing:
        raise S4RContractError(f"missing runtime input: {missing}")
    artifacts = [metadata_path, EXEC / "train_residue_masks.json.gz",
                 EXEC / "heldoutA_residue_masks.json.gz",
                 EXEC / "train_pairs.jsonl", EXEC / "heldoutA_pairs.jsonl",
                 EXEC / "gvec_f64.npy", EXEC / "gvec_graph_f64.npy",
                 EXEC / "gvec_index.json", EXEC / "construct_component.json",
                 stream_path, *runtime_inputs]
    manifest = {
        "schema": "MetaSieve.S7L2B.P2B.S4R.InputAndFirewall.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "audit_preregistration_sha256": AUDIT_PREREG_SHA,
        "audit_amendment_01_sha256": AUDIT_AMENDMENT_SHA,
        "preparation_commit": git_head(),
        "execution_surface_clean": True,
        "census": census,
        "excluded": {"train": train_excluded, "heldoutA": helda_excluded},
        "component_overlap": 0,
        "ligand_graph_overlap": 0,
        "changed_axis": "ligand statistic only",
        "ligand_baseline": {"name": "mean-pooled 41-D atom-local features",
                            "dimension": D_ATOM,
                            "W_parameters": D_ESM * D_ATOM},
        "ligand_candidate": {"name": "per-heavy-atom radius-1 Morgan environment "
                                     "counts over a frozen train vocabulary",
                             "dimension": D_LIG,
                             "vocabulary_sha256": VOCAB_SHA,
                             "W_parameters": D_ESM * D_LIG},
        "real_structural_edge_label_reads": sum(len(v) for view in masks.values()
                                                for v in view.values()),
        "residue_label_views_created": sorted(splits),
        "heldoutB_created": False,
        "heldoutB_reads": 0,
        "affinity_value_reads": 0,
        "training_process_allowed_label_view": "train_residue_masks.json.gz only",
        "scoring_process_allowed_label_views": ["heldoutA_residue_masks.json.gz"],
        "stream": {"updates": len(stream), "semantic_sha256": semantic,
                   "file_sha256": sha_file(stream_path),
                   "matches_s3r_stream": True},
        "controls_sha256": sha_file(CONTROL_PATH),
        "artifacts": {str(path.relative_to(ROOT)).replace("\\", "/"): sha_file(path)
                      for path in artifacts},
    }
    write_json(OUT / "PHASE2B_S4R_INPUT_AND_FIREWALL_MANIFEST.json", manifest)
    write_json(OUT / "PHASE2B_S4R_REAL_STREAM_MANIFEST.json", {
        "schema": "MetaSieve.S7L2B.P2B.S4R.Stream.v1",
        "updates": len(stream), "epochs": EPOCHS,
        "semantic_sha256": semantic,
        "file_sha256": sha_file(stream_path),
        "all_trained_arms_share_this_exact_stream": True,
        "identical_to_s3r_stream": semantic == S3R_STREAM_SEMANTIC_SHA,
    })
    return manifest


# --------------------------------------------------------------- runtime
class RuntimeContext:
    def __init__(self, label_split: str):
        records = list(read_jsonl(EXEC / "metadata_view.jsonl"))
        self.records = {row["source_key"]: row for row in records}
        self.pairs = list(read_jsonl(EXEC / f"{label_split}_pairs.jsonl"))
        self.masks = {key: frozenset(value) for key, value in
                      read_gzip_json(EXEC / f"{label_split}_residue_masks.json.gz").items()}
        self.label_split = label_split
        keys = json.loads((EXEC / "gvec_index.json").read_text())
        base = np.load(EXEC / "gvec_f64.npy", mmap_mode="r")
        graph = np.load(EXEC / "gvec_graph_f64.npy", mmap_mode="r")
        self.gvec = {key: np.asarray(base[i], dtype=np.float64)
                     for i, key in enumerate(keys)}
        self.gvec_graph = {key: np.asarray(graph[i], dtype=np.float64)
                           for i, key in enumerate(keys)}
        self.construct_component = json.loads(
            (EXEC / "construct_component.json").read_text())
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

    def ligands(self, kind: str):
        return self.gvec if kind == "baseline41" else self.gvec_graph

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
        raise S4RContractError("wrong label view mounted")
    manifest = json.loads((OUT / "PHASE2B_S4R_INPUT_AND_FIREWALL_MANIFEST.json").read_text())
    for relative, expected in manifest["artifacts"].items():
        if expected_split == "train" and "heldoutA_" in relative:
            continue
        if sha_file(ROOT / relative) != expected:
            raise S4RContractError(f"input artifact changed: {relative}")


def labels_for(ctx: RuntimeContext, a, b, permuted=False):
    if permuted:
        mapping = ctx.control["within_construct_derangement"]
        a, b = mapping.get(a, a), mapping.get(b, b)
    left, right = ctx.masks[a], ctx.masks[b]
    return left - right, right - left


# --------------------------------------------------------------- training
def train(kind: str) -> dict:
    if kind not in TRAINED_ARMS:
        raise S4RContractError(f"unregistered trained arm {kind}")
    directory = EXEC / kind
    checkpoint = directory / "final_checkpoint.pt"
    result_path = OUT / f"PHASE2B_S4R_TRAIN_{kind.upper()}.json"
    require_absent([checkpoint, result_path])
    configure_determinism()
    ctx = RuntimeContext("train")
    validate_runtime_inputs(ctx, "train")
    if len(LABEL_READ_LOG) != 1 or not LABEL_READ_LOG[0].endswith(
            "train_residue_masks.json.gz"):
        raise S4RContractError(f"training label-view firewall failed: {LABEL_READ_LOG}")
    stream = list(read_jsonl(EXEC / "real_stream.jsonl"))
    stream_manifest = json.loads((OUT / "PHASE2B_S4R_REAL_STREAM_MANIFEST.json").read_text())
    if semantic_stream_hash(stream) != stream_manifest["semantic_sha256"]:
        raise S4RContractError("training stream changed")

    ligands = ctx.ligands(kind)
    head = DirectW(ARM_DIM[kind]).to(DEVICE)
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
                    ga = ligands[ctx.records[a]["graph_key"]]
                    gb = ligands[ctx.records[b]["graph_key"]]
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
            raise S4RContractError("empty training update")
        loss = torch.stack(component_losses).mean()
        loss.backward()
        grad = float(head.W.grad.norm())
        if not np.isfinite(grad) or grad <= 0:
            raise S4RContractError("nonfinite or zero W gradient")
        torch.nn.utils.clip_grad_norm_(head.parameters(), CLIP)
        optimizer.step()
        head.project_norm()
        norm = float(head.W.detach().norm())
        if abs(norm - 1.0) > 1e-5:
            raise S4RContractError("unit-norm contract failed")
        trace.append({"update": item["update"], "epoch": item["epoch"],
                      "loss": float(loss.detach()), "grad_W": grad, "W_norm": norm})
        if item["update"] % 35 == 0:
            print(f"[{kind}] update {item['update']}/{N_UPDATES} "
                  f"loss={trace[-1]['loss']:.6f}", flush=True)

    if len(trace) != N_UPDATES:
        raise S4RContractError("wrong number of optimizer updates")
    movement = float((head.W.detach().cpu() - initial).norm() / initial.norm())
    directory.mkdir(parents=True, exist_ok=True)
    torch.save({"model": head.state_dict(), "optimizer": optimizer.state_dict(),
                "trace": trace, "movement": movement}, checkpoint)
    sampled_path = directory / "sampled_pairs.txt"
    sampled_path.write_text("\n".join(sampled) + "\n", encoding="utf-8")
    result = {
        "kind": kind, "device": DEVICE, "execution_commit": git_head(),
        "ligand_representation": ("mean-pooled 41-D" if kind == "baseline41"
                                  else "radius-1 Morgan d=128 per heavy atom"),
        "ligand_dimension": ARM_DIM[kind],
        "trainable_parameters": D_ESM * ARM_DIM[kind],
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
    head = DirectW(ARM_DIM[kind]).to("cpu")
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    head.load_state_dict(payload["model"])
    return head


# --------------------------------------------------------------- evaluation
def evaluate_arm(ctx, W, mode, ligand_kind="candidate", b5mean=None,
                 score_prefix=None):
    """`mode` selects the corruption; `ligand_kind` selects which frozen ligand
    table the arm reads. Every arm walks the same pair list in the same order."""
    ligands = ctx.ligands(ligand_kind)
    values, meta, scores_flat, score_index = {}, {}, [], []
    offset, near_zero = 0, 0
    max_ligand_only_field = 0.0
    by_sk = defaultdict(list)
    for row in ctx.pairs:
        by_sk[row["sk"]].append(row)
    for sk in sorted(by_sk):
        z = None
        if mode != "b5diff":
            h = ctx.h(sk)
            if mode == "context":
                record = ctx.records[by_sk[sk][0]["a"]]
                h = context_shuffle(h, record["uniprot_sequence"],
                                    int(hashlib.sha256(sk.encode()).hexdigest()[:8], 16))
            elif mode == "ligand_only":
                # One residue-mean state, so the residue field is constant by
                # construction rather than by a float64 coincidence.
                h = h.mean(axis=0, keepdims=True)
            z = np.asarray(h, dtype=np.float64) @ W
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
                raw = z @ (ligands[gk_a] - ligands[gk_b])
                if mode == "ligand_only":
                    # raw is a single scalar spread over every residue, so the
                    # field lies in span{1} subset span{Q_P} and the registered
                    # gauge projection annihilates it exactly. Evaluating the
                    # analytic value avoids ranking 1-ulp cancellation noise,
                    # which would otherwise leak the b^P direction.
                    if raw.shape != (1,):
                        raise S4RContractError("ligand-only field is not residue-constant")
                    max_ligand_only_field = max(max_ligand_only_field, abs(float(raw[0])))
                    score = np.zeros(int(ra["n_res"]), dtype=np.float64)
                else:
                    score = project_np(ctx.Qs[sk], raw)
            if not np.isfinite(score).all():
                raise S4RContractError(f"nonfinite prediction in {mode}")
            if float(np.sqrt(np.mean(score * score))) <= 1e-12:
                near_zero += 1
            gain, loss = labels_for(ctx, a, b)
            metrics = pair_metrics(score, gain, loss, int(ra["n_res"]))
            if metrics is None:
                raise S4RContractError("eligible pair lost its real label difference")
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
    return values, meta, artifacts, {
        "near_zero_pairs_retained": near_zero,
        "max_abs_ligand_only_field": max_ligand_only_field}


def component_field(values, metadata, construct_component, field):
    return aggregate({key: value[field] for key, value in values.items()}, metadata,
                     construct_component)


def assert_common_masks(arms, construct_component):
    reference_values, reference_meta = arms["candidate"]
    reference = set(reference_values)
    for name, (values, metadata) in arms.items():
        if set(values) != reference:
            raise S4RContractError(f"paired Gate mask differs for arm {name}")
        if metadata != reference_meta:
            raise S4RContractError(f"pair-to-construct mapping differs for arm {name}")
    rows = [f"{pair_id}|{reference_meta[pair_id]}|{construct_component[reference_meta[pair_id]]}"
            for pair_id in sorted(reference)]
    return hashlib.sha256("\n".join(rows).encode()).hexdigest()


def score() -> dict:
    require_absent([OUT / "PHASE2B_S4R_GATE.json",
                    OUT / "PHASE2B_S4R_COMPONENT_TABLES.json",
                    EXEC / "candidate_heldoutA_scores_f64.npy",
                    EXEC / "repeat_heldoutA_scores_f64.npy"])
    configure_determinism()
    ctx = RuntimeContext("heldoutA")
    validate_runtime_inputs(ctx, "heldoutA")
    if len(LABEL_READ_LOG) != 1 or not LABEL_READ_LOG[0].endswith(
            "heldoutA_residue_masks.json.gz"):
        raise S4RContractError(f"scoring label-view firewall failed: {LABEL_READ_LOG}")
    W = load_head("candidate").W.detach().double().numpy()
    W_repeat = load_head("repeat").W.detach().double().numpy()
    W_perm = load_head("permuted").W.detach().double().numpy()
    W_base = load_head("baseline41").W.detach().double().numpy()

    sidx = json.loads((S4 / "heldoutA_index.json").read_text())
    total = max(offset + length * atoms for offset, length, atoms in sidx.values())
    b5 = np.memmap(S5 / "heldoutA_B5.f16.dat", dtype=np.float16, mode="r", shape=(total,))
    b5mean = {}
    for key, (offset, length, atoms) in sidx.items():
        b5mean[key] = np.asarray(b5[offset:offset + length * atoms],
                                 dtype=np.float64).reshape(length, atoms).mean(1)

    arms = {}
    candidate_eval = evaluate_arm(ctx, W, "candidate", score_prefix="candidate_heldoutA")
    arms["candidate"] = candidate_eval[:2]
    arms["baseline41"] = evaluate_arm(ctx, W_base, "candidate",
                                      ligand_kind="baseline41")[:2]
    arms["b5diff"] = evaluate_arm(ctx, W, "b5diff", b5mean=b5mean)[:2]
    arms["foreign"] = evaluate_arm(ctx, W, "foreign")[:2]
    arms["context"] = evaluate_arm(ctx, W, "context")[:2]
    ligand_only_eval = evaluate_arm(ctx, W, "ligand_only")
    arms["ligand_only"] = ligand_only_eval[:2]
    arms["permuted"] = evaluate_arm(ctx, W_perm, "candidate")[:2]
    arms["chem_shuffle"] = evaluate_arm(ctx, W, "chem")[:2]
    arms["zero_W"] = evaluate_arm(ctx, np.zeros_like(W), "candidate")[:2]
    common_mask_sha = assert_common_masks(
        {key: arms[key] for key in ("candidate", "baseline41", "b5diff", "foreign",
                                    "context", "ligand_only", "permuted")},
        ctx.construct_component)

    repeat_values, repeat_meta = evaluate_arm(ctx, W_repeat, "candidate",
                                              score_prefix="repeat_heldoutA")[:2]
    if set(repeat_values) != set(arms["candidate"][0]):
        raise S4RContractError("repeat prediction mask differs")
    max_repeat_metric_diff = max(
        abs(repeat_values[key]["ap_bidir"] - arms["candidate"][0][key]["ap_bidir"])
        for key in repeat_values)
    max_repeat_W_diff = float(np.max(np.abs(W - W_repeat)))
    candidate_scores = np.load(EXEC / "candidate_heldoutA_scores_f64.npy", mmap_mode="r")
    repeat_scores = np.load(EXEC / "repeat_heldoutA_scores_f64.npy", mmap_mode="r")
    if candidate_scores.shape != repeat_scores.shape:
        raise S4RContractError("repeat raw prediction shape differs")
    max_repeat_prediction_diff = float(np.max(np.abs(candidate_scores - repeat_scores)))

    component, macro = {}, {}
    for name, (values, metadata) in arms.items():
        for field in ("ap_bidir", "chance_bidir", "ap_gain", "ap_loss"):
            component[(name, field)], macro[(name, field)] = component_field(
                values, metadata, ctx.construct_component, field)

    baseline_delta = abs(macro[("baseline41", "ap_bidir")] - S3R_CANDIDATE_MACRO_AP)
    if baseline_delta > BASELINE_REPLICATION_TOL:
        raise S4RContractError(
            "baseline replication contract failed: baseline41 macro AP "
            f"{macro[('baseline41', 'ap_bidir')]!r} differs from the S3R candidate "
            f"by {baseline_delta}")

    def contrast(left, right, margin=None, left_field="ap_bidir", right_field=None):
        result = component_bootstrap(component[(left, left_field)],
                                     component[(right, right_field or left_field)],
                                     n_boot=10000, seed=SEED_BOOT)
        if margin is None:
            result["gating"] = False
            return result
        result["margin"] = margin
        result["gating"] = True
        result["pass"] = bool(result["delta"] >= margin and
                              result["lcb95_one_sided"] > 0)
        return result

    gates = {
        "R1_vs_chance": contrast("candidate", "candidate", 0.05,
                                 "ap_bidir", "chance_bidir"),
        "R2_vs_frozen_B5_differential": contrast("candidate", "b5diff", 0.03),
        "R3_vs_two_ligand_foreign_pair": contrast("candidate", "foreign", 0.03),
        "R3b_vs_ligand_only": contrast("candidate", "ligand_only", 0.03),
        "R4_vs_residue_context_corruption": contrast("candidate", "context", 0.03),
        "R5_vs_trained_permuted_label_learner": contrast("candidate", "permuted", 0.05),
    }
    non_gating = {
        "C1_candidate_minus_baseline41": contrast("candidate", "baseline41"),
        "C2_baseline41_minus_chance": contrast("baseline41", "baseline41", None,
                                               "ap_bidir", "chance_bidir"),
    }

    trained = {kind: json.loads(
        (OUT / f"PHASE2B_S4R_TRAIN_{kind.upper()}.json").read_text())
        for kind in TRAINED_ARMS}
    stream_same = len({row["stream_semantic_sha256"] for row in trained.values()}) == 1
    score_variance = float(np.var(candidate_scores))
    chance_macro = macro[("candidate", "chance_bidir")]
    module = {
        "min_grad_W": trained["candidate"]["min_grad_W"],
        "grad_nonzero": trained["candidate"]["min_grad_W"] > 0,
        "relative_W_movement": trained["candidate"]["relative_W_movement"],
        "movement_pass": trained["candidate"]["relative_W_movement"] >= MIN_PARAM_MOVEMENT,
        "W_norm": float(np.linalg.norm(W)),
        "unit_norm_pass": bool(abs(np.linalg.norm(W) - 1.0) <= 1e-5),
        "raw_score_variance": score_variance,
        "variance_pass": bool(score_variance >= 1e-8),
        "zero_W_ap": macro[("zero_W", "ap_bidir")],
        "zero_W_chance": chance_macro,
        "zero_W_chance_pass": bool(abs(macro[("zero_W", "ap_bidir")] - chance_macro)
                                   <= DEGENERACY_TOL),
        "ligand_only_ap": macro[("ligand_only", "ap_bidir")],
        "ligand_only_is_chance_pass": bool(
            abs(macro[("ligand_only", "ap_bidir")] - chance_macro) <= DEGENERACY_TOL),
        "ligand_only_max_abs_preprojection_field":
            ligand_only_eval[3]["max_abs_ligand_only_field"],
        "ligand_only_branch_alive": bool(
            ligand_only_eval[3]["max_abs_ligand_only_field"] > 0.0),
        "context_degrades": bool(macro[("context", "ap_bidir")] <
                                 macro[("candidate", "ap_bidir")]),
        "foreign_ligands_degrade": bool(macro[("foreign", "ap_bidir")] <
                                        macro[("candidate", "ap_bidir")]),
        "same_stream_hash": stream_same,
        "baseline_replication_abs_delta": baseline_delta,
        "baseline_replication_pass": bool(baseline_delta <= BASELINE_REPLICATION_TOL),
        "repeat_W_max_abs_diff": max_repeat_W_diff,
        "repeat_prediction_max_abs_diff": max_repeat_prediction_diff,
        "repeat_metric_max_abs_diff": max_repeat_metric_diff,
        "repeat_pass": bool(max_repeat_W_diff <= 1e-7 and
                            max_repeat_prediction_diff <= 1e-7 and
                            max_repeat_metric_diff <= 1e-7),
        "candidate_prediction_sha256": sha_file(EXEC / "candidate_heldoutA_scores_f64.npy"),
        "repeat_prediction_sha256": sha_file(EXEC / "repeat_heldoutA_scores_f64.npy"),
        "common_mask_sha256": common_mask_sha,
    }
    module["pass"] = bool(
        module["grad_nonzero"] and module["movement_pass"] and
        module["unit_norm_pass"] and module["variance_pass"] and
        module["zero_W_chance_pass"] and module["ligand_only_is_chance_pass"] and
        module["context_degrades"] and module["foreign_ligands_degrade"] and
        module["same_stream_hash"] and module["baseline_replication_pass"] and
        module["repeat_pass"])

    shortcut = ("R3_vs_two_ligand_foreign_pair", "R3b_vs_ligand_only",
                "R4_vs_residue_context_corruption",
                "R5_vs_trained_permuted_label_learner")
    if not module["pass"]:
        verdict = "GRAPH_LIGAND_REPRESENTATION_TRAINING_FAILED"
    elif not (gates["R1_vs_chance"]["pass"] and
              gates["R2_vs_frozen_B5_differential"]["pass"]):
        verdict = "REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED"
    elif not all(gates[name]["pass"] for name in shortcut):
        verdict = "GRAPH_LIGAND_STATISTIC_SHORTCUT_DEPENDENT"
    else:
        verdict = "GRAPH_AWARE_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT"

    singular = np.linalg.svd(W, compute_uv=False)
    component_tables = {
        name: {field: component[(name, field)]
               for field in ("ap_bidir", "chance_bidir", "ap_gain", "ap_loss")}
        for name in arms
    }
    write_json(OUT / "PHASE2B_S4R_COMPONENT_TABLES.json", component_tables)
    write_json(OUT / "PHASE2B_S4R_MODULE_PARTICIPATION.json", module)
    result = {
        "schema": "MetaSieve.S7L2B.P2B.S4R.Gate.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "audit_preregistration_sha256": AUDIT_PREREG_SHA,
        "audit_amendment_01_sha256": AUDIT_AMENDMENT_SHA,
        "changed_axis": "ligand statistic only",
        "ligand_candidate": {"encoder": "frozen RDKit Morgan radius 1, unfolded",
                             "dimension": D_LIG, "vocabulary_sha256": VOCAB_SHA,
                             "W_parameters": D_ESM * D_LIG},
        "ligand_baseline": {"encoder": "mean-pooled 41-D atom-local features",
                            "dimension": D_ATOM, "W_parameters": D_ESM * D_ATOM},
        "primary_panel": {"pairs": len(arms["candidate"][0]),
                          "components": len(component[("candidate", "ap_bidir")]),
                          "common_mask_sha256": common_mask_sha},
        "macro_ap_bidir": {name: macro[(name, "ap_bidir")] for name in arms},
        "gates": gates,
        "non_gating_contrasts": non_gating,
        "module_participation": module,
        "singular_spectrum_head": [float(x) for x in singular[:8]],
        "near_zero_pairs_retained": candidate_eval[3]["near_zero_pairs_retained"],
        "label_view_paths_opened": list(LABEL_READ_LOG),
        "heldoutB_status": "NOT_CREATED_AND_NOT_READ_BY_REGISTRATION",
        "R6_status": "NOT_RUN_ABSOLUTE_SCALE_ORIGIN_AND_DIFFERENCE_NULLSPACE_NOT_IDENTIFIED",
        "TERMINAL_VERDICT": verdict,
        "affinity_value_reads": 0,
        "claims_not_made": ["exact residue-atom coupling", "physical energy",
                            "affinity", "selectivity", "few-shot section",
                            "biological z", "validated end-to-end DTA"],
        "authorized_next_action": (
            "preregister an independent structural confirmation on a panel that is "
            "not heldout-A"
            if verdict == "GRAPH_AWARE_RESIDUE_DIRECTION_IDENTIFIED_IN_DEVELOPMENT"
            else "none; stop at the earliest failed boundary"),
    }
    write_json(OUT / "PHASE2B_S4R_GATE.json", result)
    write_report(result)
    return result


def write_report(result) -> None:
    lines = ["# Phase 2B S4R graph-aware ligand representation result", "",
             f"Terminal verdict: `{result['TERMINAL_VERDICT']}`", "",
             "| arm | component-macro AP_bidir |", "|---|---:|"]
    lines += [f"| {name} | {value:.6f} |"
              for name, value in result["macro_ap_bidir"].items()]
    lines += ["", "| Gate | delta | LCB95 | margin | PASS |",
              "|---|---:|---:|---:|:---:|"]
    lines += [f"| {name} | {gate['delta']:.6f} | {gate['lcb95_one_sided']:.6f} | "
              f"{gate['margin']:.2f} | {gate['pass']} |"
              for name, gate in result["gates"].items()]
    lines += ["", "| non-gating contrast | delta | LCB95 |", "|---|---:|---:|"]
    lines += [f"| {name} | {row['delta']:.6f} | {row['lcb95_one_sided']:.6f} |"
              for name, row in result["non_gating_contrasts"].items()]
    lines += ["",
              "Heldout-B was neither created nor read. R6 was not opened. No affinity",
              "value was read and the frozen law operator was not modified.", ""]
    (OUT / "PHASE2B_S4R_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def run_all() -> int:
    prepare_views()
    script = str(Path(__file__).resolve())
    for mode in ("train-candidate", "train-repeat", "train-permuted",
                 "train-baseline41", "score"):
        subprocess.run([sys.executable, script, mode], cwd=ROOT, check=True)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "train-candidate", "train-repeat",
                                         "train-permuted", "train-baseline41",
                                         "score", "all"))
    args = parser.parse_args(argv)
    started = time.time()
    try:
        if args.mode == "prepare":
            output = prepare_views()
        elif args.mode.startswith("train-"):
            output = train(args.mode.split("-", 1)[1])
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
            "schema": "MetaSieve.S7L2B.P2B.S4R.FailClosed.v1",
            "created_utc": "2026-08-10", "mode": args.mode,
            "error_type": type(exc).__name__, "error": str(exc),
            "TERMINAL_VERDICT": "S4R_CONTRACT_OR_LABEL_FIREWALL_FAIL_CLOSED",
            "affinity_value_reads": 0,
        }
        write_json(OUT / "PHASE2B_S4R_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
