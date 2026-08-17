"""Stage Q trainer: decoupled frozen-feature level head (single stage)."""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    centered_task_error, compact_episode, learning_rate_factor,
    normalized_episode, ranking_term, training_label_scale,
)
from tools.research.stageB_complementary.train_stageb import (
    draw_fit_episode, internal_validation_bank, partition_components,
)
from tools.research.stageJ_assay.train_stagej import journal_vocabulary
from tools.research.stageQ_frozenhead.model import FrozenHeadModel

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
PANEL_DIM = 35


def panel_stats(data, spec):
    cells = data.cells
    ligand_mw = {}
    for row in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
        ligand_mw[str(row["drug_key"])] = float(row["molecular_weight"])
    mw, atom_means = [], []
    for i in spec.query:
        ligand = cells[int(i)]["ligand_id"]
        mw.append(ligand_mw[ligand])
        atoms, _, mask = data.ligand_bank.get(ligand)
        n = int(mask.sum())
        if n:
            atom_means.append(atoms[:n].mean(0))
    mw = np.asarray(mw, dtype=np.float32)
    atom_vec = np.asarray(atom_means, dtype=np.float32).mean(0) if atom_means \
        else np.zeros(32, dtype=np.float32)
    return np.concatenate([
        np.asarray([len(spec.query), mw.mean(), mw.std()], dtype=np.float32),
        atom_vec]).astype(np.float32)


def episode_journal_ids(data, spec, index, max_codes=8):
    ids = set()
    for i in spec.query:
        for pid in data.cells[int(i)]["panel_ids"]:
            body = str(pid).split("|")[0]
            parts = body.split("/")
            if len(parts) >= 2 and parts[0].startswith("doi:"):
                publisher = parts[0].split(".")[-1]
                journal = "".join(c for c in parts[1] if c.isalpha())[:4]
                for code in (("pub_" + publisher,), ("jnl_" + journal,)):
                    if code in index:
                        ids.add(index[code])
    ordered = sorted(ids)[:max_codes]
    padded = ordered + [-1] * (max_codes - len(ordered))
    return torch.as_tensor(padded, dtype=torch.long)


def forward_q(model, episode, ids, stats, gate):
    import torch.nn.functional as F_pad
    if ids is not None:
        ids = ids.to(next(model.parameters()).device)
    stats = torch.as_tensor(stats, dtype=torch.float32,
                            device=next(model.parameters()).device)
    if episode.support_atoms.ndim == 3:
        episode = compact_episode(episode)
    support_active = (int(episode.support_mask.sum(-1).max())
                      if episode.support_mask.numel() else 0)
    active = max(support_active, int(episode.query_mask.sum(-1).max()))
    if episode.support_atoms.ndim == 4:
        support_atoms = episode.support_atoms[:, :, :active]
        support_bonds = episode.support_bonds[:, :, :active, :active]
        support_mask = episode.support_mask[:, :, :active]
        query_atoms = episode.query_atoms[:, :, :active]
        query_bonds = episode.query_bonds[:, :, :active, :active]
        query_mask = episode.query_mask[:, :, :active]
    else:
        support_atoms = episode.support_atoms[:, :active]
        support_bonds = episode.support_bonds[:, :active, :active]
        support_mask = episode.support_mask[:, :active]
        query_atoms = episode.query_atoms[:, :active]
        query_bonds = episode.query_bonds[:, :active, :active]
        query_mask = episode.query_mask[:, :active]
    support_atoms = F_pad.pad(support_atoms, (0, 0, 0, active - support_atoms.shape[-2]))
    support_bonds = F_pad.pad(support_bonds, (0, 0, 0, active - support_bonds.shape[-2],
                                              0, active - support_bonds.shape[-2]))
    support_mask = F_pad.pad(support_mask, (0, active - support_mask.shape[-1]))
    query_atoms = F_pad.pad(query_atoms, (0, 0, 0, active - query_atoms.shape[-2]))
    query_bonds = F_pad.pad(query_bonds, (0, 0, 0, active - query_bonds.shape[-2],
                                          0, active - query_bonds.shape[-2]))
    query_mask = F_pad.pad(query_mask, (0, active - query_mask.shape[-1]))
    return model(
        episode.protein_pooled, episode.protein_tokens, episode.protein_mask,
        support_atoms, support_bonds, support_mask, episode.support_y,
        query_atoms, query_bonds, query_mask, adapt=True,
        protein_chemistry=episode.protein_chemistry,
        support_fingerprint=episode.support_fingerprint,
        query_fingerprint=episode.query_fingerprint,
        journal_ids=ids, panel_stats=stats, level_gate=gate,
        task_state_override=None, geometry_coordinates=None,
        geometry_edge_index=None, geometry_available=None,
        geometry_common_frame=None)


def build_model(config, data, vocab):
    base = config.base
    return FrozenHeadModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=base.hidden_dim, task_dim=base.task_dim,
        ligand_layers=base.ligand_layers, pair_dim=base.pair_dim,
        pair_blocks=base.pair_blocks, pair_latents=base.pair_latents,
        pair_heads=base.pair_heads, pair_chunk_size=base.pair_chunk_size,
        support_hidden_dim=base.support_hidden_dim,
        support_blocks=base.support_blocks, adapter_rank=base.adapter_rank,
        adaptive_blocks=base.adaptive_blocks, adapter_scale=base.adapter_scale,
        journal_vocab=len(vocab), panel_dim=PANEL_DIM, use_learned_key=False,
        dtype=torch.float32)


def gpu_probe(data, config, model, vocab, output_dir):
    import subprocess
    device = config.base.device
    model = model.to(device).train()
    torch.manual_seed(0)
    rng = np.random.default_rng(0)
    fit, _ = partition_components(data)
    spec = draw_fit_episode(data, fit, 0, 8, rng)
    label_scale = training_label_scale(data)
    episode = compact_episode(normalized_episode(data.materialize(spec),
                                                 label_scale))
    if device.startswith("cuda"):
        from dataclasses import replace
        episode = replace(episode, **{
            field: getattr(episode, field).to(device)
            for field in episode.__dataclass_fields__
            if isinstance(getattr(episode, field), torch.Tensor)})
    index = {code: i for i, code in enumerate(vocab)}
    ids = episode_journal_ids(data, spec, index).to(device)
    stats = panel_stats(data, spec)
    probe = {"torch_cuda_is_available": bool(torch.cuda.is_available()),
             "torch_version": torch.__version__,
             "cuda_runtime": torch.version.cuda,
             "device_count": torch.cuda.device_count(),
             "configured_device": device,
             "model_parameter_devices": sorted({str(p.device)
                                                for p in model.parameters()}),
             "nvidia_smi_samples": []}
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-5)
    parsed = []
    for step in range(4):
        output = forward_q(model, episode, ids, stats, 1.0)
        loss = F.smooth_l1_loss(output.prediction, episode.query_y.to(
            device=output.prediction.device, dtype=output.prediction.dtype))
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        if step >= 1:
            try:
                text = subprocess.run(
                    ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",
                     "--format=csv,noheader,nounits"],
                    capture_output=True, text=True, timeout=10).stdout.strip()
                parts = text.split(",")
                parsed.append((int(parts[0].strip()), int(parts[1].strip())))
                probe["nvidia_smi_samples"].append(text)
            except Exception:  # noqa: BLE001
                pass
    probe["max_gpu_utilization_percent"] = max((v for v, _ in parsed),
                                                default=None)
    probe["utilization_nonzero"] = any(v > 0 for v, _ in parsed)
    probe["batch_device_check"] = bool(device.startswith("cuda"))
    (output_dir / "GPU_PROBE.json").write_text(
        json.dumps(probe, indent=1, sort_keys=True), encoding="utf-8")
    print(json.dumps(probe, indent=1), flush=True)
    return probe


def internal_score(model, data, banks, label_scale, index):
    model.eval()
    total, count = 0.0, 0
    with torch.no_grad():
        for k, bank in banks.items():
            gate = 1.0 if k == 0 else 0.0
            for episode in bank:
                ids = episode_journal_ids(data, episode.spec, index)
                stats = panel_stats(data, episode.spec)
                output = forward_q(model, episode, ids, stats, gate)
                error = float(F.mse_loss(output.prediction, episode.query_y.to(
                    device=output.prediction.device,
                    dtype=output.prediction.dtype)))
                total += error * label_scale.scale ** 2
                count += 1
    model.train()
    return total / max(count, 1)


def train(data, config, output_dir, progress_path):
    base = config.base
    torch.manual_seed(base.seed)
    rng = np.random.default_rng(base.seed)
    label_scale = training_label_scale(data)
    vocab = journal_vocabulary(data)
    model = build_model(config, data, vocab)
    probe = gpu_probe(data, config, model, vocab, output_dir)
    if not probe["torch_cuda_is_available"] or not probe["batch_device_check"]:
        raise RuntimeError("GPU verification failed")
    model = model.to(base.device)
    index = {code: i for i, code in enumerate(vocab)}
    trainable = [p for p in model.parameters() if p.requires_grad]
    fast_prefixes = ("meta.term.", "transport.")
    fast_parameters = [p for n, p in model.named_parameters()
                       if p.requires_grad and n.startswith(fast_prefixes)]
    fast_ids = {id(p) for p in fast_parameters}
    slow_parameters = [p for p in trainable if id(p) not in fast_ids]
    optimizer = torch.optim.AdamW([
        {"params": slow_parameters,
         "lr": base.learning_rate * base.backbone_lr_scale},
        {"params": fast_parameters, "lr": base.learning_rate},
    ], weight_decay=base.weight_decay)
    base_lrs = [group["lr"] for group in optimizer.param_groups]
    fit_components, internal_components = partition_components(data)
    banks = internal_validation_bank(data, internal_components, label_scale)
    print("arm=" + config.arm + " vocab " + str(len(vocab)) + " fit "
          + str(len(fit_components)) + " internal " + str(len(internal_components)),
          flush=True)
    best_state, best_value, best_step = None, float("inf"), 0
    trace = []
    started = time.monotonic()
    for step in range(1, base.steps + 1):
        model.train()
        factor = learning_rate_factor(step, base)
        for group, start in zip(optimizer.param_groups, base_lrs):
            group["lr"] = start * factor
        optimizer.zero_grad(set_to_none=True)
        support_size = SUPPORT_SIZES[(step - 1) % len(SUPPORT_SIZES)]
        gate = 1.0 if support_size == 0 else (1.0 if config.arm == "Q-UNGATED"
                                              else 0.0)
        selected = []
        for _ in range(base.episodes_per_step):
            requested = int(rng.integers(base.min_query_size,
                                         base.query_size + 1))
            spec = draw_fit_episode(data, fit_components, support_size,
                                    requested, rng)
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            ids = episode_journal_ids(data, spec, index)
            stats = panel_stats(data, spec)
            selected.append((episode, ids, stats))
        loss_value = 0.0
        for episode, ids, stats in selected:
            output = forward_q(model, episode, ids, stats, gate)
            query_y = episode.query_y.to(device=output.prediction.device,
                                         dtype=output.prediction.dtype)
            loss_post = F.smooth_l1_loss(output.prediction, query_y)
            loss_pre = F.smooth_l1_loss(output.zero_shot, query_y)
            loss_rank = ranking_term(output.prediction, query_y, base,
                                     label_scale)
            loss_centered = centered_task_error(output.prediction, query_y)
            loss = (loss_post + base.zero_shot_loss_weight * loss_pre
                    + base.ranking_loss_weight * loss_rank
                    + base.shape_loss_weight * loss_centered
                    + base.support_match_loss_weight
                    * output.support_match_loss)
            loss_value += float(loss.detach())
            (loss / len(selected)).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), base.grad_clip)
        optimizer.step()
        trace.append(loss_value / len(selected))
        if step % base.val_interval == 0 or step == base.steps:
            value = internal_score(model, data, banks, label_scale, index)
            progress = {"step": step, "arm": config.arm,
                        "internal_val_mse_pk": value,
                        "best": min(best_value, value),
                        "elapsed_seconds": time.monotonic() - started}
            print(json.dumps(progress), flush=True)
            if progress_path is not None:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(progress) + "\n")
            if value < best_value:
                best_state = copy.deepcopy(model.state_dict())
                best_value, best_step = value, step
    if best_state is None:
        raise RuntimeError("training produced no checkpoint")
    model.load_state_dict(best_state)
    report = {"arm": config.arm, "best_internal_val_mse_pk": best_value,
              "best_step": best_step, "loss_trace": trace,
              "label_scale": asdict(label_scale),
              "checkpoint_selection": "meta_train internal components only",
              "journal_vocab": len(vocab),
              "fit_components": len(fit_components),
              "internal_val_components": len(internal_components),
              "optimization_steps": base.steps,
              "wall_time_seconds": time.monotonic() - started,
              "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                     if base.device.startswith("cuda") else 0.0),
              "trainable_parameters": int(sum(p.numel() for p in trainable)),
              "gpu_probe": probe}
    payload = {"model_state": model.state_dict(),
               "journal_vocab": vocab, "config": asdict(base),
               "arm": config.arm, "label_scale": asdict(label_scale),
               "stageQ_version": 1}
    torch.save(payload, output_dir / "checkpoint.pt")
    (output_dir / "RESULT.json").write_text(json.dumps(
        {"arm": config.arm, "report": report,
         "meta_test": data.seal_record()}, indent=1), encoding="utf-8")
    print("wrote " + str(output_dir))
    return model, report, label_scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", required=True, choices=("Q", "Q-UNGATED"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()
    progress_path = args.output / "progress.jsonl"
    if progress_path.exists() and not args.force:
        raise SystemExit(f"{progress_path} exists; pass --force to overwrite")
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    config = TrainConfig(arch="similarity_only", steps=args.steps,
                         seed=args.seed, split_directory=str(SPLIT),
                         device=args.device, amp=False)
    from tools.research.stageD_level_panel.train_staged import StageEConfig
    args.output.mkdir(parents=True, exist_ok=True)
    train(data, StageEConfig(base=config, arm=args.arm), args.output,
          progress_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
