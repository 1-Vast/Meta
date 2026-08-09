"""S1R: pairwise-ranking objective repair for the Phase 2B synthetic control."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "s7_l2b_r0r"
sys.path.insert(0, str(HERE))

import s0r_run as S  # noqa: E402
from p2b_residue_residual import CLIP, LR, WD, Head, project_np, sha_file  # noqa: E402
from s0_synth import HCache, Teacher, build_panels  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S1R_PAIRWISE_RANK_OBJECTIVE.md"
PREREG_SHA = "b130f5316c6b962e3a9508d2cb8628c0f3eab587c348f0dbf7d07b3fcd987557"
CALIBRATION_SEEDS = (20260921, 20260922, 20260923)
SEALED_SEED = 20260998
PARAM_SEED = 20260901
AP_THRESHOLD = 0.50
ALIGN_DROP = 0.05
EXEC_DIR = S.META_DIR / "s1r_execution"
OUT = S.OUT


class S1RContractError(RuntimeError):
    pass


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def stamp():
    x = S.stamp()
    x.update({
        "stage": "P1R2B-PHASE2B-S1R_PAIRWISE_RANK_OBJECTIVE_REPAIR",
        "preregistration_sha256": PREREG_SHA,
        "only_changed_object": "group-balanced BCE -> all-residue bidirectional pairwise logistic ranking loss",
    })
    return x


def rank_loss_torch(d: torch.Tensor, gain_set, loss_set) -> torch.Tensor:
    """Exact registered all-residue bidirectional pairwise logistic loss."""
    length = int(d.numel())
    terms = []
    for positives, score in ((gain_set, d), (loss_set, -d)):
        pos_idx = torch.tensor(sorted(positives), dtype=torch.long, device=d.device)
        if pos_idx.numel() == 0:
            continue
        keep = torch.ones(length, dtype=torch.bool, device=d.device)
        keep[pos_idx] = False
        negative = score[keep]
        positive = score[pos_idx]
        if negative.numel():
            terms.append(F.softplus(-(positive[:, None] - negative[None, :])).mean())
    if not terms:
        return d.sum() * 0.0
    return torch.stack(terms).mean()


def rank_loss_np(d: np.ndarray, row) -> float:
    length = int(d.size)
    terms = []
    for indices, score in ((row["gi"], d), (row["li"], -d)):
        idx = np.asarray(indices, dtype=np.int64)
        if idx.size == 0:
            continue
        keep = np.ones(length, dtype=bool)
        keep[idx] = False
        diff = score[idx, None] - score[keep][None, :]
        terms.append(float(np.logaddexp(0.0, -diff).mean()))
    return float(np.mean(terms)) if terms else 0.0


def train_update(ctx, head, optimizer, batch, teacher, hcache):
    optimizer.zero_grad(set_to_none=True)
    batch_loss = torch.zeros((), dtype=torch.float32)
    for _component, constructs in batch:
        component_loss = torch.zeros((), dtype=torch.float32)
        for sk, pair_list in constructs:
            h = torch.from_numpy(hcache.get(sk).copy()).float()
            uh = head.uh(h)
            q = torch.from_numpy(ctx["Qs"][sk]).float()
            pair_losses = []
            for _same_sk, a, b in pair_list:
                ga = ctx["gvec"][ctx["rec_of"][a]["graph_key"]]
                gb = ctx["gvec"][ctx["rec_of"][b]["graph_key"]]
                raw = uh @ head.vg(torch.from_numpy((ga - gb).copy()).float())
                d = raw - q @ (q.T @ raw) if q.shape[1] else raw
                gain, loss = teacher.labels(sk, a, b)
                pair_losses.append(rank_loss_torch(d, gain, loss))
            if pair_losses:
                component_loss = component_loss + torch.stack(pair_losses).mean() / len(constructs)
        batch_loss = batch_loss + component_loss / len(batch)
    batch_loss.backward()
    torch.nn.utils.clip_grad_norm_(head.parameters(), CLIP)
    optimizer.step()
    if not all(torch.isfinite(p).all() for p in head.parameters()):
        raise S1RContractError("non-finite model parameter")
    return float(batch_loss.detach())


def panel_rank_loss(ctx, head, rows, hcache):
    U, V = S.head_arrays(head)
    by_sk = defaultdict(list)
    for row in rows:
        by_sk[row["sk"]].append(row)
    values, ordered = [], []
    for sk in sorted(by_sk):
        uh = hcache.get(sk).astype(np.float64) @ U.T
        q = ctx["Qs"][sk]
        for row in by_sk[sk]:
            gd = ctx["gvec"][row["gk_a"]] - ctx["gvec"][row["gk_b"]]
            d = project_np(q, uh @ (V @ gd))
            values.append(rank_loss_np(d, row))
            ordered.append(row)
    component, macro = S.hier_macro(values, ordered, ctx)
    return {"component_macro_rank_loss": macro, "component_rank_loss": component}


def replay_rank_loss(scores_path, index_path, construct_component):
    scores = np.load(scores_path, mmap_mode="r")
    values, rows = [], []
    for item in S.load_jsonl(index_path):
        offset, length = int(item["offset"]), int(item["length"])
        row = {
            "sk": item["sk"], "L": length,
            "gi": np.asarray(item["gain_indices"], dtype=np.int64),
            "li": np.asarray(item["loss_indices"], dtype=np.int64),
        }
        values.append(rank_loss_np(np.asarray(scores[offset:offset + length]), row))
        rows.append(row)
    _component, macro = S.hier_macro(
        values, rows, {"construct_component": construct_component}
    )
    return macro


def save_capture(tag, update, ctx, head, optimizer, panels, hcache, last_loss):
    directory = EXEC_DIR / tag
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / f"u{update:04d}_heldout"
    held, pred_artifacts = S.score_panel(ctx, head, panels["heldout"], hcache, prefix)
    train_rank = panel_rank_loss(ctx, head, panels["train"], hcache)
    held_rank = replay_rank_loss(
        prefix.with_name(prefix.name + "_scores_f64.npy"),
        prefix.with_name(prefix.name + "_index.jsonl"),
        ctx["construct_component"],
    )
    checkpoint = directory / f"u{update:04d}_checkpoint.pt"
    torch.save({"update": update, "model": head.state_dict(),
                "optimizer": optimizer.state_dict(), "last_batch_loss": last_loss}, checkpoint)
    expected_state = S.state_hash(head)
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    loaded = Head()
    loaded.load_state_dict(payload["model"])
    state_ok = S.state_hash(loaded) == expected_state
    replay = S.replay_prediction_artifacts(
        prefix.with_name(prefix.name + "_scores_f64.npy"),
        prefix.with_name(prefix.name + "_index.jsonl"),
        ctx["construct_component"],
    )
    replay_ok = (
        replay["pairs"] == held["pairs"]
        and replay["components"] == held["components"]
        and abs(replay["component_macro_ap_bidir"]
                - held["component_macro_ap_bidir"]) <= 1e-14
    )
    return {
        "update": update, "last_batch_loss": last_loss,
        "checkpoint": str(checkpoint.relative_to(ROOT)).replace("\\", "/"),
        "checkpoint_sha256": sha_file(checkpoint),
        "checkpoint_state_sha256": expected_state,
        "checkpoint_reload_identical": state_ok,
        "prediction_reload_identical": replay_ok,
        "prediction_artifacts": pred_artifacts,
        "train": train_rank,
        "heldout": {**held, "component_macro_rank_loss": held_rank},
    }


def run_trajectory(tag, ctx, panels, stream, teacher, hcache, start, checkpoints):
    head = Head(PARAM_SEED)
    if start == "teacher":
        with torch.no_grad():
            head.U.copy_(torch.from_numpy(teacher.U).float())
            head.V.copy_(torch.from_numpy(teacher.V).float())
    elif start != "random":
        raise ValueError(start)
    optimizer = torch.optim.AdamW(head.parameters(), lr=LR, weight_decay=WD)
    captures = []
    if 0 in checkpoints:
        captures.append(save_capture(tag, 0, ctx, head, optimizer, panels, hcache, None))
    last_loss = None
    target = max(checkpoints)
    for update, (_epoch, _batch_index, batch) in enumerate(stream[:target], start=1):
        last_loss = train_update(ctx, head, optimizer, batch, teacher, hcache)
        if update in checkpoints:
            captures.append(save_capture(
                tag, update, ctx, head, optimizer, panels, hcache, last_loss
            ))
    if captures[-1]["update"] != target:
        raise S1RContractError("trajectory stopped before registered checkpoint")
    if not all(x["checkpoint_reload_identical"] and x["prediction_reload_identical"]
               for x in captures):
        raise S1RContractError("checkpoint/prediction reload failed")
    return {**stamp(), "tag": tag, "teacher_seed": teacher.seed,
            "initialization": start, "checkpoints": captures}


def alignment_decision(trajectory):
    by_update = {x["update"]: x for x in trajectory["checkpoints"]}
    x0, x100 = by_update[0], by_update[100]
    inference = S.bootstrap_delta(
        x0["heldout"]["component_ap"], x100["heldout"]["component_ap"]
    )
    ap_drop = (x0["heldout"]["component_macro_ap_bidir"]
               - x100["heldout"]["component_macro_ap_bidir"])
    loss_down = (x100["train"]["component_macro_rank_loss"]
                 < x0["train"]["component_macro_rank_loss"])
    passed = bool(loss_down and ap_drop <= ALIGN_DROP
                  and inference["ucb95_one_sided"] >= -ALIGN_DROP)
    return {"train_rank_loss_decreased": loss_down, "heldout_ap_drop": ap_drop,
            "component_inference": inference, "PASS": passed}


def student_decision(trajectory):
    final = trajectory["checkpoints"][-1]
    ap = final["heldout"]["component_macro_ap_bidir"]
    return {"updates": final["update"], "heldout_ap_bidir": ap,
            "threshold": AP_THRESHOLD, "PASS": bool(ap >= AP_THRESHOLD)}


def write_terminal(verdict, details, started):
    artifact = {**stamp(), "TERMINAL_VERDICT": verdict,
                "details": details, "elapsed_seconds": round(time.time() - started, 3),
                "remaining_frozen": ["real Phase2B", "affinity",
                                     "independent structural confirmation",
                                     "few-shot section", "biological z",
                                     "A(F,z)=K(B(z)F(z))"]}
    write_json(OUT / "PHASE2B_S1R_VERDICT.json", artifact)
    write_json(OUT / "PHASE2B_POST_S1R_NOT_RUN.json", {
        **stamp(), "status": "NOT_RUN", "terminal_verdict": verdict,
        "not_run": artifact["remaining_frozen"],
    })
    return artifact


def main():
    started = time.time()
    try:
        if sha_file(PREREG) != PREREG_SHA:
            raise S1RContractError("preregistration hash mismatch")
        ctx = S.prepare_label_blind_context()
        hcache = HCache(ctx)
        stream_path = S.EXEC_DIR / "frozen_stream.jsonl"
        s0_manifest = json.loads(
            (OUT / "S0R_FROZEN_STREAM_MANIFEST.json").read_text(encoding="utf-8")
        )
        if sha_file(stream_path) != s0_manifest["stream"]["file_sha256"]:
            raise S1RContractError("reused stream file hash mismatch")
        stream = S.deserialize_stream(stream_path)
        if len(stream) != S.N_UPDATES or S.stream_hash(stream) != s0_manifest["stream"]["semantic_sha256"]:
            raise S1RContractError("reused stream semantic mismatch")

        input_artifact = {**stamp(), "calibration_seeds": list(CALIBRATION_SEEDS),
                          "sealed_seed_committed_but_not_instantiated": SEALED_SEED,
                          "stream_file_sha256": sha_file(stream_path),
                          "stream_semantic_sha256": S.stream_hash(stream),
                          "complete_train_pairs": S.EXPECTED["train_panel"],
                          "complete_heldout_pairs": S.EXPECTED["heldout"],
                          "complete_heldout_components": S.EXPECTED["heldout_components"]}
        write_json(OUT / "S1R_INPUT_AND_STREAM_MANIFEST.json", input_artifact)

        alignments = []
        calibration_cache = []
        for seed in CALIBRATION_SEEDS:
            print(f"S1R calibration alignment seed {seed}", flush=True)
            teacher = Teacher(ctx, hcache, seed)
            panels = build_panels(ctx, teacher, hcache)
            if len(panels["train"]) != S.EXPECTED["train_panel"] or len(panels["heldout"]) != S.EXPECTED["heldout"]:
                raise S1RContractError("synthetic panel count changed")
            trajectory = run_trajectory(
                f"cal_{seed}_teacher", ctx, panels, stream, teacher, hcache,
                "teacher", (0, 100),
            )
            decision = alignment_decision(trajectory)
            alignments.append({"seed": seed, "decision": decision,
                               "trajectory": trajectory})
            calibration_cache.append((seed, teacher, panels))

        alignment_pass = all(x["decision"]["PASS"] for x in alignments)
        if not alignment_pass:
            gate = {**stamp(), "alignment": alignments,
                    "student_trajectories": "NOT_RUN",
                    "PASS": False,
                    "verdict": "PAIRWISE_SURROGATE_STILL_MISALIGNED"}
            write_json(OUT / "S1R_CALIBRATION_GATE.json", gate)
            write_terminal("PAIRWISE_SURROGATE_STILL_MISALIGNED", gate, started)
            return 0

        students = []
        for seed, teacher, panels in calibration_cache:
            print(f"S1R calibration student seed {seed}", flush=True)
            trajectory = run_trajectory(
                f"cal_{seed}_student", ctx, panels, stream, teacher, hcache,
                "random", (210,),
            )
            students.append({"seed": seed, "decision": student_decision(trajectory),
                             "trajectory": trajectory})
        student_pass = all(x["decision"]["PASS"] for x in students)
        calibration_gate = {**stamp(), "alignment": alignments,
                            "students": students, "PASS": student_pass,
                            "sealed_seed_instantiated": False,
                            "verdict": ("CALIBRATION_PASS" if student_pass else
                                        "PAIRWISE_OBJECTIVE_ALIGNED_BUDGET_NOT_IDENTIFIED")}
        write_json(OUT / "S1R_CALIBRATION_GATE.json", calibration_gate)
        if not student_pass:
            write_terminal("PAIRWISE_OBJECTIVE_ALIGNED_BUDGET_NOT_IDENTIFIED",
                           calibration_gate, started)
            return 0

        # The sealed teacher is first instantiated below, after the persisted gate.
        print(f"S1R sealed verification seed {SEALED_SEED}", flush=True)
        teacher = Teacher(ctx, hcache, SEALED_SEED)
        panels = build_panels(ctx, teacher, hcache)
        sealed_alignment_traj = run_trajectory(
            f"sealed_{SEALED_SEED}_teacher", ctx, panels, stream, teacher, hcache,
            "teacher", (0, 100),
        )
        sealed_student_traj = run_trajectory(
            f"sealed_{SEALED_SEED}_student", ctx, panels, stream, teacher, hcache,
            "random", (210,),
        )
        sealed_alignment = alignment_decision(sealed_alignment_traj)
        sealed_student = student_decision(sealed_student_traj)
        sealed_pass = bool(sealed_alignment["PASS"] and sealed_student["PASS"])
        sealed = {**stamp(), "seed": SEALED_SEED,
                  "alignment": sealed_alignment,
                  "student": sealed_student,
                  "alignment_trajectory": sealed_alignment_traj,
                  "student_trajectory": sealed_student_traj,
                  "PASS": sealed_pass}
        write_json(OUT / "S1R_SEALED_VERIFICATION_GATE.json", sealed)
        verdict = ("SYNTHETIC_IDENTIFIABILITY_REPAIRED" if sealed_pass
                   else "PAIRWISE_OBJECTIVE_VERIFICATION_FAILED")
        terminal = write_terminal(verdict, sealed, started)

        report = OUT / "PHASE2B_S1R_PAIRWISE_OBJECTIVE_REPORT.md"
        lines = ["# Phase 2B S1R pairwise-objective repair", "",
                 f"Terminal verdict: `{verdict}`", "",
                 "The only changed object was the synthetic training loss.", "",
                 "## Calibration", "",
                 "| seed | alignment | student AP | student pass |",
                 "|---:|---|---:|---|"]
        for a, student in zip(alignments, students):
            lines.append(f"| {a['seed']} | {a['decision']['PASS']} | "
                         f"{student['decision']['heldout_ap_bidir']:.6f} | "
                         f"{student['decision']['PASS']} |")
        lines += ["", "## Sealed verification", "",
                  f"- Alignment PASS: {sealed_alignment['PASS']}",
                  f"- Student AP: {sealed_student['heldout_ap_bidir']:.6f}",
                  f"- Student PASS: {sealed_student['PASS']}", "",
                  "No real structural edge or affinity label was read. The frozen",
                  "law operator and biological z remain untouched."]
        report.write_text("\n".join(lines) + "\n", encoding="utf-8")

        reports = [OUT / x for x in (
            "S1R_INPUT_AND_STREAM_MANIFEST.json", "S1R_CALIBRATION_GATE.json",
            "S1R_SEALED_VERIFICATION_GATE.json", "PHASE2B_S1R_VERDICT.json",
            "PHASE2B_POST_S1R_NOT_RUN.json", "PHASE2B_S1R_PAIRWISE_OBJECTIVE_REPORT.md",
        )]
        execution = list(EXEC_DIR.rglob("*"))
        write_json(OUT / "S1R_DETACHED_ARTIFACT_MANIFEST.json", {
            **stamp(), "reports": S.file_manifest(reports),
            "execution_artifacts": S.file_manifest(execution),
            "terminal_verdict": terminal["TERMINAL_VERDICT"],
        })
        print(json.dumps({"TERMINAL_VERDICT": verdict,
                          "elapsed_seconds": round(time.time() - started, 1)},
                         indent=2), flush=True)
        return 0
    except Exception as exc:
        terminal = write_terminal("S1R_CONTRACT_INVALID",
                                  {"error_type": type(exc).__name__,
                                   "error": str(exc)}, started)
        print(json.dumps(terminal, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
