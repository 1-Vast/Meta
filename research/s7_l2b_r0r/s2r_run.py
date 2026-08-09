"""S2R gauge-free direct-W binary ordinal witness."""
from __future__ import annotations

import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "s7_l2b_r0r"
sys.path.insert(0, str(HERE))

import s0r_run as S  # noqa: E402
import s1r_run as R  # noqa: E402
from p2b_residue_residual import CLIP, D_ATOM, D_ESM, LR, project_np, sha_file  # noqa: E402
from s0_synth import HCache, Teacher, build_panels, hier_macro, pair_ap  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S2R_GAUGE_FREE_DIRECT_W.md"
PREREG_SHA = "b10bb815e47c8b33de49653168bcd06bb0cb1a793711888f7c5c704732b83577"
CALIBRATION_SEEDS = (20260931, 20260932, 20260933)
SEALED_SEED = 20260997
PARAM_SEED = 20260901
AP_THRESHOLD = 0.50
NORM_EPS = 1e-12
EXEC_DIR = S.META_DIR / "s2r_execution"
OUT = S.OUT
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


class S2RContractError(RuntimeError):
    pass


class DirectW(nn.Module):
    def __init__(self, seed=PARAM_SEED):
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.W = nn.Parameter(torch.randn(D_ESM, D_ATOM, generator=generator) * 1e-3)
        self.project_norm()

    @torch.no_grad()
    def project_norm(self):
        norm = self.W.norm()
        if not torch.isfinite(norm) or norm <= 0:
            raise S2RContractError("invalid direct-W norm")
        self.W.div_(norm)


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def stamp():
    x = S.stamp()
    x.update({
        "stage": "P1R2B-PHASE2B-S2R_GAUGE_FREE_DIRECT_W",
        "preregistration_sha256": PREREG_SHA,
        "training_device": DEVICE,
        "trainable_object": "one direct W matrix (1280 x 41)",
        "trainable_parameters": D_ESM * D_ATOM,
    })
    return x


def normalized(d):
    return d / torch.sqrt(torch.mean(d * d) + NORM_EPS)


def state_hash(head):
    return S.hashlib.sha256(head.W.detach().cpu().contiguous().numpy().tobytes()).hexdigest()


def train_update(ctx, head, optimizer, batch, teacher, hcache):
    optimizer.zero_grad(set_to_none=True)
    batch_loss = torch.zeros((), device=DEVICE)
    for _component, constructs in batch:
        component_loss = torch.zeros((), device=DEVICE)
        for sk, pair_list in constructs:
            h = torch.from_numpy(hcache.get(sk).copy()).float().to(DEVICE)
            z = h @ head.W
            q = torch.from_numpy(ctx["Qs"][sk]).double().to(DEVICE)
            pair_losses = []
            for _same_sk, a, b in pair_list:
                ga = ctx["gvec"][ctx["rec_of"][a]["graph_key"]]
                gb = ctx["gvec"][ctx["rec_of"][b]["graph_key"]]
                gd = torch.from_numpy((ga - gb).copy()).float().to(DEVICE)
                raw = (z @ gd).double()
                d = raw - q @ (q.T @ raw)
                gain, loss = teacher.labels(sk, a, b)
                pair_losses.append(R.rank_loss_torch(normalized(d), gain, loss))
            if pair_losses:
                component_loss = component_loss + torch.stack(pair_losses).mean() / len(constructs)
        batch_loss = batch_loss + component_loss / len(batch)
    batch_loss.backward()
    torch.nn.utils.clip_grad_norm_(head.parameters(), CLIP)
    optimizer.step()
    head.project_norm()
    return float(batch_loss.detach())


def score_panel(ctx, W, rows, hcache, prefix: Path | None = None):
    W = np.asarray(W, dtype=np.float64)
    by_sk = defaultdict(list)
    for row in rows:
        by_sk[row["sk"]].append(row)
    values, ordered = [], []
    flat, index = [], []
    offset = 0
    for sk in sorted(by_sk):
        z = hcache.get(sk).astype(np.float64) @ W
        q = ctx["Qs"][sk]
        for row in by_sk[sk]:
            gd = ctx["gvec"][row["gk_a"]] - ctx["gvec"][row["gk_b"]]
            d = project_np(q, z @ gd)
            ap = pair_ap(d, row)[0]
            values.append(ap)
            ordered.append(row)
            if prefix is not None:
                flat.append(d)
                index.append({
                    "sk": sk, "a": row["a"], "b": row["b"],
                    "component": ctx["construct_component"][sk],
                    "offset": offset, "length": int(d.size),
                    "gain_indices": row["gi"].tolist(),
                    "loss_indices": row["li"].tolist(),
                    "ap_bidir": ap,
                })
                offset += int(d.size)
    comp, macro = hier_macro(values, ordered, ctx)
    by_construct = defaultdict(list)
    for value, row in zip(values, ordered):
        if value is not None:
            by_construct[row["sk"]].append(value)
    result = {
        "pairs": len(rows), "constructs": len(by_construct),
        "components": len(comp), "component_macro_ap_bidir": macro,
        "component_ap": comp,
        "construct_ap": {k: float(np.mean(v)) for k, v in by_construct.items()},
    }
    artifacts = {}
    if prefix is not None:
        scores = prefix.with_name(prefix.name + "_scores_f64.npy")
        indices = prefix.with_name(prefix.name + "_index.jsonl")
        np.save(scores, np.concatenate(flat).astype(np.float64, copy=False))
        S.write_jsonl(indices, index)
        artifacts = S.file_manifest([scores, indices])
    return result, artifacts


def replay_ap(scores_path, index_path, construct_component):
    scores = np.load(scores_path, mmap_mode="r")
    values, rows = [], []
    for item in S.load_jsonl(index_path):
        offset, length = int(item["offset"]), int(item["length"])
        yg = np.zeros(length, dtype=np.int8)
        yl = np.zeros(length, dtype=np.int8)
        yg[np.asarray(item["gain_indices"], dtype=np.int64)] = 1
        yl[np.asarray(item["loss_indices"], dtype=np.int64)] = 1
        row = {"sk": item["sk"], "L": length,
               "gi": np.flatnonzero(yg), "li": np.flatnonzero(yl),
               "yg": yg, "yl": yl}
        values.append(pair_ap(np.asarray(scores[offset:offset + length]), row)[0])
        rows.append(row)
    comp, macro = hier_macro(values, rows, {"construct_component": construct_component})
    return {"pairs": len(rows), "components": len(comp),
            "component_macro_ap_bidir": macro}


def rank8_diagnostic(ctx, W, panels, hcache):
    u, singular, vt = np.linalg.svd(W, full_matrices=False)
    W8 = (u[:, :8] * singular[:8]) @ vt[:8]
    train, _ = score_panel(ctx, W8, panels["train"], hcache)
    held, _ = score_panel(ctx, W8, panels["heldout"], hcache)
    return {"singular_values": singular.tolist(),
            "rank8_energy_fraction": float((singular[:8] ** 2).sum()
                                            / (singular ** 2).sum()),
            "rank8_train_ap": train["component_macro_ap_bidir"],
            "rank8_heldout_ap": held["component_macro_ap_bidir"]}


def run_seed(seed, ctx, stream, hcache):
    teacher = Teacher(ctx, hcache, seed)
    panels = build_panels(ctx, teacher, hcache)
    if len(panels["train"]) != S.EXPECTED["train_panel"] or len(panels["heldout"]) != S.EXPECTED["heldout"]:
        raise S2RContractError("panel count mismatch")
    head = DirectW(PARAM_SEED).to(DEVICE)
    optimizer = torch.optim.Adam(head.parameters(), lr=LR)
    losses = []
    for update, (_epoch, _batch_index, batch) in enumerate(stream, start=1):
        losses.append(train_update(ctx, head, optimizer, batch, teacher, hcache))
    if update != S.N_UPDATES:
        raise S2RContractError("update count mismatch")
    norm = float(head.W.detach().norm())
    if abs(norm - 1.0) > 1e-5:
        raise S2RContractError("direct-W norm contract failed")
    directory = EXEC_DIR / str(seed)
    directory.mkdir(parents=True, exist_ok=True)
    checkpoint = directory / "u0210_checkpoint.pt"
    torch.save({"model": head.state_dict(), "optimizer": optimizer.state_dict(),
                "update": update, "last_batch_loss": losses[-1]}, checkpoint)
    expected = state_hash(head)
    loaded = DirectW().to("cpu")
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    loaded.load_state_dict(payload["model"])
    if state_hash(loaded) != expected:
        raise S2RContractError("checkpoint reload mismatch")
    W = head.W.detach().cpu().double().numpy()
    prefix = directory / "heldout"
    train, _ = score_panel(ctx, W, panels["train"], hcache)
    held, artifacts = score_panel(ctx, W, panels["heldout"], hcache, prefix)
    replay = replay_ap(prefix.with_name(prefix.name + "_scores_f64.npy"),
                       prefix.with_name(prefix.name + "_index.jsonl"),
                       ctx["construct_component"])
    if abs(replay["component_macro_ap_bidir"] - held["component_macro_ap_bidir"]) > 1e-14:
        raise S2RContractError("prediction artifact reload mismatch")
    compression = rank8_diagnostic(ctx, W, panels, hcache)
    return {
        "seed": seed, "updates": update, "last_batch_loss": losses[-1],
        "W_norm": norm, "W_state_sha256": expected,
        "checkpoint": str(checkpoint.relative_to(ROOT)).replace("\\", "/"),
        "checkpoint_sha256": sha_file(checkpoint),
        "prediction_artifacts": artifacts,
        "train": train, "heldout": held, "rank8_diagnostic": compression,
        "train_pass": bool(train["component_macro_ap_bidir"] >= AP_THRESHOLD),
        "heldout_pass": bool(held["component_macro_ap_bidir"] >= AP_THRESHOLD),
    }


def preflight(ctx, hcache):
    head = DirectW().to("cpu")
    W = head.W.detach().double().numpy()
    sk = sorted(ctx["Qs"])[0]
    h = hcache.get(sk).astype(np.float32)
    q = ctx["Qs"][sk]
    gd = next(iter(ctx["gvec"].values())) - list(ctx["gvec"].values())[1]
    raw_parent = (torch.from_numpy(h) @ head.W.detach() @ torch.from_numpy(gd).float()).double()
    q_t = torch.from_numpy(q).double()
    parent = (raw_parent - q_t @ (q_t.T @ raw_parent)).numpy()
    direct = project_np(q, (h.astype(np.float64) @ W) @ gd)
    # Evaluation uses float64 accumulation, while the parent training path uses
    # float32 raw followed by float64 projection. The invariant concerns the
    # projection operation for an identical raw vector.
    same_raw = project_np(q, raw_parent.numpy())
    projection_error = float(np.max(np.abs(parent - same_raw)))
    d = torch.from_numpy(parent.copy()).double()
    normalized_error = float(torch.max(torch.abs(normalized(7.0 * d) - normalized(d))))
    failures = []
    if projection_error > 1e-12:
        failures.append("float64 projection parity failed")
    if normalized_error > 1e-8:
        failures.append("positive-scale normalization invariance failed")
    if not np.isfinite(direct).all():
        failures.append("non-finite direct score")
    return {"float64_projection_max_error": projection_error,
            "positive_scale_normalization_max_error": normalized_error,
            "initial_W_norm": float(head.W.norm()), "failures": failures}


def write_terminal(verdict, details, started):
    artifact = {**stamp(), "TERMINAL_VERDICT": verdict, "details": details,
                "elapsed_seconds": round(time.time() - started, 3),
                "remaining_frozen": ["real Phase2B", "affinity",
                                     "independent structural confirmation",
                                     "few-shot section", "biological z",
                                     "A(F,z)=K(B(z)F(z))"]}
    write_json(OUT / "PHASE2B_S2R_VERDICT.json", artifact)
    write_json(OUT / "PHASE2B_POST_S2R_NOT_RUN.json", {
        **stamp(), "status": "NOT_RUN", "terminal_verdict": verdict,
        "not_run": artifact["remaining_frozen"]})
    return artifact


def main():
    started = time.time()
    try:
        if sha_file(PREREG) != PREREG_SHA:
            raise S2RContractError("preregistration hash mismatch")
        ctx = S.prepare_label_blind_context()
        hcache = HCache(ctx)
        stream_path = S.EXEC_DIR / "frozen_stream.jsonl"
        s0_manifest = json.loads((OUT / "S0R_FROZEN_STREAM_MANIFEST.json").read_text())
        stream = S.deserialize_stream(stream_path)
        if sha_file(stream_path) != s0_manifest["stream"]["file_sha256"] or S.stream_hash(stream) != s0_manifest["stream"]["semantic_sha256"]:
            raise S2RContractError("frozen stream mismatch")
        checks = preflight(ctx, hcache)
        write_json(OUT / "S2R_PREFLIGHT_AND_INPUT_MANIFEST.json", {
            **stamp(), "checks": checks, "calibration_seeds": list(CALIBRATION_SEEDS),
            "sealed_seed_not_instantiated": SEALED_SEED,
            "stream_file_sha256": sha_file(stream_path)})
        if checks["failures"]:
            raise S2RContractError("; ".join(checks["failures"]))

        calibration = []
        for seed in CALIBRATION_SEEDS:
            print(f"S2R calibration seed {seed}", flush=True)
            calibration.append(run_seed(seed, ctx, stream, hcache))
        train_pass = all(x["train_pass"] for x in calibration)
        heldout_pass = all(x["heldout_pass"] for x in calibration)
        gate = {**stamp(), "runs": calibration, "train_PASS": train_pass,
                "heldout_PASS": heldout_pass, "PASS": train_pass and heldout_pass,
                "sealed_seed_instantiated": False}
        write_json(OUT / "S2R_CALIBRATION_GATE.json", gate)
        if not train_pass:
            write_terminal("BINARY_ORDINAL_TRAIN_FIT_NOT_IDENTIFIED", gate, started)
            return 0
        if not heldout_pass:
            write_terminal("FINITE_DESIGN_GENERALIZATION_NOT_IDENTIFIED", gate, started)
            return 0

        print(f"S2R sealed verification seed {SEALED_SEED}", flush=True)
        sealed = run_seed(SEALED_SEED, ctx, stream, hcache)
        sealed_pass = sealed["train_pass"] and sealed["heldout_pass"]
        write_json(OUT / "S2R_SEALED_VERIFICATION_GATE.json", {
            **stamp(), "run": sealed, "PASS": sealed_pass})
        verdict = ("BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED" if sealed_pass
                   else "GAUGE_FREE_BINARY_ORDINAL_VERIFICATION_FAILED")
        terminal = write_terminal(verdict, sealed, started)

        report = OUT / "PHASE2B_S2R_GAUGE_FREE_DIRECT_W_REPORT.md"
        lines = ["# Phase 2B S2R gauge-free direct-W witness", "",
                 f"Terminal verdict: `{verdict}`", "",
                 "| seed | train AP | held-out AP | rank-8 held-out AP |",
                 "|---:|---:|---:|---:|"]
        for run in [*calibration, sealed]:
            lines.append(f"| {run['seed']} | {run['train']['component_macro_ap_bidir']:.6f} | "
                         f"{run['heldout']['component_macro_ap_bidir']:.6f} | "
                         f"{run['rank8_diagnostic']['rank8_heldout_ap']:.6f} |")
        lines += ["", "No real structural edge or affinity value was read.",
                  "The result concerns a bounded binary residue-ranking statistic only;"]
        lines += ["biology admission, affinity, few-shot sectioning and the frozen law remain closed."]
        report.write_text("\n".join(lines) + "\n", encoding="utf-8")

        report_paths = [OUT / x for x in (
            "S2R_PREFLIGHT_AND_INPUT_MANIFEST.json", "S2R_CALIBRATION_GATE.json",
            "S2R_SEALED_VERIFICATION_GATE.json", "PHASE2B_S2R_VERDICT.json",
            "PHASE2B_POST_S2R_NOT_RUN.json", "PHASE2B_S2R_GAUGE_FREE_DIRECT_W_REPORT.md")]
        write_json(OUT / "S2R_DETACHED_ARTIFACT_MANIFEST.json", {
            **stamp(), "reports": S.file_manifest(report_paths),
            "execution_artifacts": S.file_manifest(list(EXEC_DIR.rglob("*"))),
            "terminal_verdict": terminal["TERMINAL_VERDICT"]})
        print(json.dumps({"TERMINAL_VERDICT": verdict,
                          "elapsed_seconds": round(time.time() - started, 1)}, indent=2), flush=True)
        return 0
    except Exception as exc:
        terminal = write_terminal("S2R_CONTRACT_INVALID",
                                  {"error_type": type(exc).__name__, "error": str(exc)}, started)
        print(json.dumps(terminal, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
