"""Fit the preregistered R0 exact-distance arms and score fresh R0-C once."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM, MAX_ATOMS
from research.correspondence_router.r0_exact_distance import (
    component_macro,
    deterministic_derangement,
    negative_log_likelihood,
    paired_component_bootstrap,
    ranked_probability_score,
)
from research.correspondence_router.run_r0b_prefit import (
    FROZEN_CHECKPOINT_SHA256,
    _load_exact_geometry,
    _load_exact_proteins,
)
from research.e0_identifiability.run_tdir_pilot import (
    RESIDUE_CLASSES,
    _load_frozen_model,
)
from research.meta_fewshot.affinity_pair_field import (
    AffinityDirectedPairField,
    coarse_interaction_compatibility,
    exact_distance_loss_per_system,
)
from scripts.build_holo_complex_index import _ccd_molecule
from scripts.cache_r0b_exact_geometry import (
    FROZEN_LIGAND_BANK_SHA256,
    FROZEN_PANEL_SHA256,
    R0C_LIGAND_BANK_SHA256,
    R0C_PANEL_SHA256,
    verify_exact_geometry_bank,
)
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


SEEDS = (1701, 1702, 1703)
EPOCHS = 20
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-4
ORDINAL_WEIGHT = 0.25
RANK = 8
PAIR_BUDGET = 1_500_000
MAX_BATCH = 32
GRADIENT_CLIP = 5.0
BOOTSTRAP_DRAWS = 10_000
BOOTSTRAP_SEED = 1700
@dataclass
class System:
    entry: str
    component: str
    atom_state: np.ndarray
    residue_state: np.ndarray
    distance_prior_slot: np.ndarray
    atom_class: np.ndarray
    residue_class: np.ndarray
    slots: np.ndarray
    labels: np.ndarray
    n1_residue_state: np.ndarray | None = None
    n1_compatibility: np.ndarray | None = None

    @property
    def atoms(self) -> int:
        return len(self.atom_state)

    @property
    def residues(self) -> int:
        return len(self.residue_state)


def _validate_inputs(
    r0b_panel: Path,
    r0b_geometry: Path,
    r0b_proteins: Path,
    r0b_ligands: Path,
    r0c_panel: Path,
    r0c_geometry: Path,
    r0c_proteins: Path,
    r0c_ligands: Path,
    checkpoint: Path,
    prefit: Path,
) -> None:
    expected = {
        r0b_panel: FROZEN_PANEL_SHA256,
        r0b_ligands: FROZEN_LIGAND_BANK_SHA256,
        r0c_panel: R0C_PANEL_SHA256,
        r0c_ligands: R0C_LIGAND_BANK_SHA256,
        checkpoint: FROZEN_CHECKPOINT_SHA256,
    }
    for path, digest in expected.items():
        if sha256_file(path) != digest:
            raise ValueError(f"frozen training input hash differs: {path}")
    verify_exact_geometry_bank(r0b_geometry, contract="r0b")
    verify_exact_geometry_bank(r0c_geometry, contract="r0c")
    for root, panel_sha in ((r0b_proteins, FROZEN_PANEL_SHA256),
                            (r0c_proteins, R0C_PANEL_SHA256)):
        manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        if (manifest.get("records_sha256") != panel_sha or
                manifest.get("checkpoint_sha256") != FROZEN_CHECKPOINT_SHA256):
            raise ValueError(f"exact protein bank is not frozen: {root}")
    prefit_report = json.loads(prefit.read_text(encoding="utf-8"))
    if (prefit_report.get("verdict") != "R0C_PREFIT_ADMISSION_PASS" or
            prefit_report.get("panel_sha256") != R0C_PANEL_SHA256 or
            prefit_report.get("ligand_bank_sha256") != R0C_LIGAND_BANK_SHA256 or
            not all(prefit_report.get("gates", {}).values())):
        raise ValueError("R0-C prefit did not authorize training")


def _atom_classes(record: dict, factory) -> np.ndarray:
    chemistry = _ccd_molecule(Path(record["ccd_path"]))
    molecule = chemistry["molecule"]
    name_to_heavy = {
        name: index for index, name in enumerate(chemistry["heavy_atom_names"])
    }
    atom_to_heavy = {
        atom.GetIdx(): name_to_heavy[atom.GetProp("_CCDAtomName")]
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetProp("_CCDAtomName") in name_to_heavy
    }
    channels = np.zeros((len(name_to_heavy), 8), dtype=np.uint8)
    family_to_index = {
        "Hydrophobe": 0, "LumpedHydrophobe": 0, "Donor": 2,
        "Acceptor": 3, "PosIonizable": 4, "NegIonizable": 5,
    }
    for feature in factory.GetFeaturesForMol(molecule):
        channel = family_to_index.get(feature.GetFamily())
        if channel is None:
            continue
        for atom_index in feature.GetAtomIds():
            if atom_index in atom_to_heavy:
                channels[atom_to_heavy[atom_index], channel] = 1
    for atom_index, heavy_index in atom_to_heavy.items():
        atom = molecule.GetAtomWithIdx(atom_index)
        channels[heavy_index, 1] = int(atom.GetIsAromatic())
        channels[heavy_index, 6] = int(atom.GetAtomicNum() in {9, 17, 35, 53})
    channels[:, 7] = channels[:, :7].sum(axis=1) == 0
    return channels.argmax(axis=1).astype(np.int8)


def _frozen_batches(model, rows: list[tuple[dict, dict, dict]], device: str,
                    batch_size: int = 16) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    result = {}
    for start in range(0, len(rows), batch_size):
        batch = rows[start:start + batch_size]
        size = len(batch)
        X = torch.zeros(size, MAX_ATOMS, ATOM_FEAT_DIM, device=device)
        A = torch.zeros(size, MAX_ATOMS, MAX_ATOMS, BOND_FEAT_DIM, device=device)
        mask = torch.zeros(size, MAX_ATOMS, device=device)
        residue = torch.zeros(size, 128, 128, device=device)
        residue_mask = torch.zeros(size, 128, device=device)
        for index, (_, graph, protein) in enumerate(batch):
            X[index] = graph["X"].float().to(device)
            mask[index] = graph["mask"].float().to(device)
            edge = graph["edge_index"].long().to(device)
            A[index, edge[0], edge[1]] = graph["edge_attr"].float().to(device)
            residue[index] = torch.from_numpy(
                protein["slot_projected"]).float().to(device)
            residue_mask[index] = torch.from_numpy(
                protein["slot_mask"]).float().to(device)
        with torch.inference_mode():
            _, atom_state = model.ligand(X, A, mask)
            bridge = model.bridge(atom_state, mask, residue, residue_mask)
            prior = torch.softmax(bridge.distance_logits, dim=-1)
        for index, (record, _, _) in enumerate(batch):
            atoms = int(mask[index].sum().item())
            result[record["source_entry_id"]] = (
                atom_state[index, :atoms].half().cpu().numpy(),
                prior[index, :atoms].half().cpu().numpy(),
            )
    return result


def _load_systems(panel_path: Path, geometry_dir: Path, protein_dir: Path,
                  ligand_bank_path: Path, frozen_model, device: str,
                  allowed_splits: set[str]) -> list[System]:
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures

    panel = [row for row in read_jsonl(panel_path)
             if row["r0_split"] in allowed_splits]
    by_entry = {row["source_entry_id"]: row for row in panel}
    index = [row for row in read_jsonl(geometry_dir / "index.jsonl")
             if row["source_entry_id"] in by_entry]
    if len(index) != len(panel):
        raise ValueError("panel and exact geometry counts differ")
    geometry = _load_exact_geometry(geometry_dir, index)
    proteins = _load_exact_proteins(
        protein_dir, {row["sequence_sha256"] for row in panel})
    ligands = torch.load(ligand_bank_path, map_location="cpu", weights_only=False)
    frozen_rows = [(row, ligands[row["ccd_sha256"]],
                    proteins[row["sequence_sha256"]]) for row in panel]
    frozen = _frozen_batches(frozen_model, frozen_rows, device)
    factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    systems = []
    for row in panel:
        exact = geometry[row["source_entry_id"]]
        protein = proteins[row["sequence_sha256"]]
        atom_state, prior = frozen[row["source_entry_id"]]
        if len(protein["exact_projected"]) != len(row["sequence"]):
            raise ValueError("exact protein length differs from panel")
        component = row.get(
            "homology_group_id", row.get("r0c_final_component_id"))
        systems.append(System(
            entry=row["source_entry_id"],
            component=str(component),
            atom_state=atom_state,
            residue_state=protein["exact_projected"],
            distance_prior_slot=prior,
            atom_class=_atom_classes(row, factory),
            residue_class=np.asarray([
                RESIDUE_CLASSES.get(code, 5) for code in row["sequence"]
            ], dtype=np.int8),
            slots=exact["slot_of_residue"].astype(np.int64),
            labels=exact["distance_bin"].astype(np.int64),
        ))
    return sorted(systems, key=lambda value: value.entry)


def _batches(systems: list[System], order: np.ndarray) -> list[list[System]]:
    result, current = [], []
    max_atoms = max_residues = 0
    for index in order.tolist():
        system = systems[index]
        next_atoms = max(max_atoms, system.atoms)
        next_residues = max(max_residues, system.residues)
        if current and (len(current) == MAX_BATCH or
                        (len(current) + 1) * next_atoms * next_residues > PAIR_BUDGET):
            result.append(current)
            current, max_atoms, max_residues = [], 0, 0
        current.append(system)
        max_atoms = max(max_atoms, system.atoms)
        max_residues = max(max_residues, system.residues)
    if current:
        result.append(current)
    return result


def _prepare_n1(system: System) -> None:
    if system.n1_residue_state is not None:
        return
    residue = system.residue_state.astype(np.float32).copy()
    compatibility = coarse_interaction_compatibility(
        torch.from_numpy(system.atom_class.astype(np.int64))[None],
        torch.from_numpy(system.residue_class.astype(np.int64))[None],
    )[0].numpy()
    for slot in np.unique(system.slots):
        members = np.flatnonzero(system.slots == slot)
        residue[members] = residue[members].mean(axis=0, keepdims=True)
        compatibility[:, members] = compatibility[:, members].mean(
            axis=1, keepdims=True)
    system.n1_residue_state = residue
    system.n1_compatibility = compatibility


def _pack(systems: list[System], arm: str, device: str,
          intervention: str | None = None) -> dict[str, torch.Tensor]:
    batch = len(systems)
    atoms, residues = max(s.atoms for s in systems), max(s.residues for s in systems)
    atom_state = torch.zeros(batch, atoms, 128, device=device)
    residue_state = torch.zeros(batch, residues, 128, device=device)
    atom_mask = torch.zeros(batch, atoms, dtype=torch.bool, device=device)
    residue_mask = torch.zeros(batch, residues, dtype=torch.bool, device=device)
    atom_class = torch.full((batch, atoms), -1, dtype=torch.long, device=device)
    residue_class = torch.full((batch, residues), -1, dtype=torch.long, device=device)
    distance_prior = torch.zeros(batch, atoms, residues, 5, device=device)
    labels = torch.zeros(batch, atoms, residues, dtype=torch.long, device=device)
    atom_indices = torch.full((batch, atoms), -1, dtype=torch.long, device=device)
    residue_indices = torch.full((batch, residues), -1, dtype=torch.long, device=device)
    for row, system in enumerate(systems):
        a, r = system.atoms, system.residues
        atom_values = system.atom_state.astype(np.float32)
        if arm == "n1":
            _prepare_n1(system)
            residue_values = system.n1_residue_state
        else:
            residue_values = system.residue_state.astype(np.float32)
        atom_classes = system.atom_class.copy()
        residue_classes = system.residue_class.copy()
        if intervention == "residue":
            mapping = deterministic_derangement(
                system.slots.tolist(), namespace=f"R0C-N3|{system.entry}")
            residue_values = residue_values[mapping]
            residue_classes = residue_classes[mapping]
        elif intervention == "atom":
            mapping = deterministic_derangement(
                [0] * a, namespace=f"R0C-N4|{system.entry}")
            atom_values = atom_values[mapping]
            atom_classes = atom_classes[mapping]
        atom_state[row, :a] = torch.from_numpy(atom_values).to(device)
        residue_state[row, :r] = torch.from_numpy(residue_values).to(device)
        atom_mask[row, :a] = True
        residue_mask[row, :r] = True
        atom_class[row, :a] = torch.from_numpy(atom_classes.astype(np.int64)).to(device)
        residue_class[row, :r] = torch.from_numpy(residue_classes.astype(np.int64)).to(device)
        distance_prior[row, :a, :r] = torch.from_numpy(
            system.distance_prior_slot[:, system.slots].astype(np.float32)).to(device)
        labels[row, :a, :r] = torch.from_numpy(system.labels).to(device)
        atom_indices[row, :a] = torch.arange(a, device=device)
        residue_indices[row, :r] = torch.arange(r, device=device)
    compatibility = coarse_interaction_compatibility(atom_class, residue_class)
    if arm == "n1":
        for row, system in enumerate(systems):
            compatibility[row, :system.atoms, :system.residues] = torch.from_numpy(
                system.n1_compatibility).to(device)
    prior_total = distance_prior.sum(dim=-1, keepdim=True)
    normalized_prior = torch.where(
        prior_total > 0,
        distance_prior / prior_total.clamp_min(torch.finfo(distance_prior.dtype).tiny),
        distance_prior,
    )
    return {
        "atom_states": atom_state,
        "atom_mask": atom_mask,
        "residue_states": residue_state,
        "residue_mask": residue_mask,
        "contact_prior": normalized_prior[..., :2].sum(dim=-1).clamp_(0.0, 1.0),
        "distance_prior": distance_prior,
        "compatibility": compatibility,
        "atom_indices": atom_indices,
        "residue_indices": residue_indices,
        "labels": labels,
    }


def _new_model(seed: int, arm: str, device: str) -> AffinityDirectedPairField:
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = AffinityDirectedPairField(
        128, 128, rank=RANK,
        interaction_mode="additive" if arm == "n2" else "bilinear",
    ).to(device)
    with torch.no_grad():
        model.distance_residual.weight.zero_()
    return model


def _predict(model: AffinityDirectedPairField, systems: list[System], arm: str,
             device: str, intervention: str | None = None) -> dict[str, np.ndarray]:
    model.eval()
    output = {}
    order = np.arange(len(systems))
    with torch.inference_mode():
        for batch in _batches(systems, order):
            packed = _pack(batch, arm, device, intervention)
            labels = packed.pop("labels")
            prediction = model(**packed)
            for index, system in enumerate(batch):
                output[system.entry] = prediction.distance_prob[
                    index, :system.atoms, :system.residues].cpu().numpy()
    return output


def _metrics(probability: dict[str, np.ndarray], systems: list[System]) -> dict:
    system_rps, system_nll, component_of = {}, {}, {}
    for system in systems:
        value = probability[system.entry]
        system_rps[system.entry] = float(
            ranked_probability_score(value, system.labels).mean())
        system_nll[system.entry] = float(
            negative_log_likelihood(value, system.labels).mean())
        component_of[system.entry] = system.component
    component_rps = component_macro(system_rps, component_of)
    component_nll = component_macro(system_nll, component_of)
    return {
        "system_rps": system_rps,
        "system_nll": system_nll,
        "component_rps": component_rps,
        "component_nll": component_nll,
        "rps": float(np.mean(list(component_rps.values()))),
        "nll": float(np.mean(list(component_nll.values()))),
    }


def _train(seed: int, arm: str, train: list[System], validation: list[System],
           output: Path, device: str, selected_epoch: int | None = None) -> dict:
    model = _new_model(seed, arm, device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    history, selected_state = [], None
    for epoch in range(1, EPOCHS + 1):
        model.train()
        rng = np.random.default_rng(seed * 10_000 + epoch)
        losses = []
        for batch in _batches(train, rng.permutation(len(train))):
            packed = _pack(batch, arm, device)
            labels = packed.pop("labels")
            optimizer.zero_grad(set_to_none=True)
            prediction = model(**packed)
            loss = exact_distance_loss_per_system(
                prediction, labels, ordinal_weight=ORDINAL_WEIGHT).mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), GRADIENT_CLIP)
            optimizer.step()
            losses.append(float(loss.detach().cpu()))
        validation_metrics = _metrics(
            _predict(model, validation, arm, device), validation)
        history.append({
            "epoch": epoch,
            "train_loss": float(np.mean(losses)),
            "validation_rps": validation_metrics["rps"],
            "validation_nll": validation_metrics["nll"],
        })
        checkpoint = {
            "schema": "MetaSieve.R0C.DistanceResidualCheckpoint.v1",
            "seed": seed,
            "arm": arm,
            "epoch": epoch,
            "model_state": {key: value.detach().cpu() for key, value in model.state_dict().items()},
            "config": {"rank": RANK, "interaction_mode": model.interaction_mode},
        }
        if arm == "full":
            path = output / f"full_seed{seed}_epoch{epoch:02d}.pt"
            torch.save(checkpoint, path)
        if selected_epoch == epoch:
            selected_state = checkpoint
    if arm != "full":
        if selected_state is None:
            raise RuntimeError("selected epoch was not reached")
        path = output / f"{arm}_seed{seed}_selected.pt"
        torch.save(selected_state, path)
    return {"history": history}


def _resume_full(output: Path, validation: list[System], device: str) -> dict:
    expected = {
        f"full_seed{seed}_epoch{epoch:02d}.pt"
        for seed in SEEDS for epoch in range(1, EPOCHS + 1)
    }
    present = {path.name for path in output.iterdir()}
    if present != expected:
        raise ValueError("resume-full requires exactly the frozen Full checkpoints")
    training = {"full": {}}
    for seed in SEEDS:
        history = []
        for epoch in range(1, EPOCHS + 1):
            path = output / f"full_seed{seed}_epoch{epoch:02d}.pt"
            checkpoint = torch.load(path, map_location="cpu", weights_only=False)
            if (checkpoint.get("schema") !=
                    "MetaSieve.R0C.DistanceResidualCheckpoint.v1" or
                    checkpoint.get("seed") != seed or
                    checkpoint.get("arm") != "full" or
                    checkpoint.get("epoch") != epoch):
                raise ValueError(f"invalid Full checkpoint for resume: {path}")
            model = _load_checkpoint(path, device)
            value = _metrics(_predict(model, validation, "full", device), validation)
            history.append({
                "epoch": epoch,
                "train_loss": None,
                "validation_rps": value["rps"],
                "validation_nll": value["nll"],
            })
        training["full"][str(seed)] = {"history": history}
    return training


def _load_checkpoint(path: Path, device: str) -> AffinityDirectedPairField:
    checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    model = AffinityDirectedPairField(
        128, 128, rank=int(checkpoint["config"]["rank"]),
        interaction_mode=checkpoint["config"]["interaction_mode"],
    ).to(device)
    model.load_state_dict(checkpoint["model_state"])
    return model


def _ensemble(predictions: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    return {
        entry: np.mean([prediction[entry] for prediction in predictions], axis=0)
        for entry in predictions[0]
    }


def run(args) -> dict:
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("R0-C residual training is registered for CUDA")
    output = Path(args.output)
    if output.exists() and not args.resume_full:
        raise FileExistsError(f"R0-C training output already exists: {output}")
    paths = {name: Path(getattr(args, name)) for name in (
        "r0b_panel", "r0b_geometry", "r0b_proteins", "r0b_ligands",
        "r0c_panel", "r0c_geometry", "r0c_proteins", "r0c_ligands",
        "checkpoint", "prefit",
    )}
    _validate_inputs(**paths)
    output.mkdir(parents=True, exist_ok=args.resume_full)
    start = time.perf_counter()
    checkpoint_data = torch.load(paths["checkpoint"], map_location="cpu", weights_only=False)
    frozen, _ = _load_frozen_model(
        paths["checkpoint"], int(checkpoint_data["protein_dim"]), args.device)
    train = _load_systems(
        paths["r0b_panel"], paths["r0b_geometry"], paths["r0b_proteins"],
        paths["r0b_ligands"], frozen, args.device, {"train"})
    validation = _load_systems(
        paths["r0b_panel"], paths["r0b_geometry"], paths["r0b_proteins"],
        paths["r0b_ligands"], frozen, args.device, {"val"})
    confirmation = _load_systems(
        paths["r0c_panel"], paths["r0c_geometry"], paths["r0c_proteins"],
        paths["r0c_ligands"], frozen, args.device, {"heldout_b"})
    del frozen

    if args.resume_full:
        training = _resume_full(output, validation, args.device)
    else:
        training = {"full": {}}
        for seed in SEEDS:
            training["full"][str(seed)] = _train(
                seed, "full", train, validation, output, args.device)
    mean_validation = {
        epoch: float(np.mean([
            training["full"][str(seed)]["history"][epoch - 1]["validation_rps"]
            for seed in SEEDS
        ])) for epoch in range(1, EPOCHS + 1)
    }
    selected_epoch = min(mean_validation, key=lambda epoch: (mean_validation[epoch], epoch))
    for arm in ("n1", "n2"):
        training[arm] = {}
        for seed in SEEDS:
            training[arm][str(seed)] = _train(
                seed, arm, train, validation, output, args.device,
                selected_epoch=selected_epoch)

    predictions: dict[str, list[dict[str, np.ndarray]]] = defaultdict(list)
    per_seed = defaultdict(dict)
    for arm in ("full", "n1", "n2"):
        for seed in SEEDS:
            path = (output / f"full_seed{seed}_epoch{selected_epoch:02d}.pt"
                    if arm == "full" else output / f"{arm}_seed{seed}_selected.pt")
            model = _load_checkpoint(path, args.device)
            value = _predict(model, confirmation, arm, args.device)
            predictions[arm].append(value)
            per_seed[arm][str(seed)] = _metrics(value, confirmation)
            if arm == "full":
                predictions["n3"].append(_predict(
                    model, confirmation, arm, args.device, intervention="residue"))
                predictions["n4"].append(_predict(
                    model, confirmation, arm, args.device, intervention="atom"))

    prior_probability = {}
    for system in confirmation:
        probability = system.distance_prior_slot[:, system.slots].astype(np.float64)
        probability /= probability.sum(axis=-1, keepdims=True)
        prior_probability[system.entry] = probability
    metrics = {"prior": _metrics(prior_probability, confirmation)}
    for arm in ("full", "n1", "n2", "n3", "n4"):
        metrics[arm] = _metrics(_ensemble(predictions[arm]), confirmation)
    free_arm = min(("prior", "n1", "n2"), key=lambda arm: metrics[arm]["rps"])
    delta_star = float(json.loads(paths["prefit"].read_text(encoding="utf-8"))["delta_star"])
    g = metrics[free_arm]["rps"] - metrics["full"]["rps"]
    contrasts = {
        "g1_prior": paired_component_bootstrap(
            metrics["prior"]["component_rps"], metrics["full"]["component_rps"],
            seed=BOOTSTRAP_SEED, draws=BOOTSTRAP_DRAWS),
        "g2_free": paired_component_bootstrap(
            metrics[free_arm]["component_rps"], metrics["full"]["component_rps"],
            seed=BOOTSTRAP_SEED + 1, draws=BOOTSTRAP_DRAWS),
        "g3_residue": paired_component_bootstrap(
            metrics["n3"]["component_rps"], metrics["full"]["component_rps"],
            seed=BOOTSTRAP_SEED + 2, draws=BOOTSTRAP_DRAWS),
        "g3_atom": paired_component_bootstrap(
            metrics["n4"]["component_rps"], metrics["full"]["component_rps"],
            seed=BOOTSTRAP_SEED + 3, draws=BOOTSTRAP_DRAWS),
        "nll_full_minus_free": paired_component_bootstrap(
            metrics["full"]["component_nll"], metrics[free_arm]["component_nll"],
            seed=BOOTSTRAP_SEED + 4, draws=BOOTSTRAP_DRAWS),
    }
    g1_seed_direction = all(
        metrics["prior"]["rps"] > per_seed["full"][str(seed)]["rps"]
        for seed in SEEDS)
    if free_arm == "prior":
        g2_seed_direction = g1_seed_direction
    else:
        g2_seed_direction = all(
            per_seed[free_arm][str(seed)]["rps"] > per_seed["full"][str(seed)]["rps"]
            for seed in SEEDS)
    gates = {
        "g1_prior_incremental": (
            contrasts["g1_prior"]["delta"] >= delta_star and
            contrasts["g1_prior"]["lcb95_one_sided"] > 0 and g1_seed_direction),
        "g2_exact_pair_incremental": (
            contrasts["g2_free"]["delta"] >= delta_star and
            contrasts["g2_free"]["lcb95_one_sided"] > 0 and g2_seed_direction),
        "g3_residue_identity": (
            contrasts["g3_residue"]["delta"] >= 0.5 * g and
            contrasts["g3_residue"]["lcb95_one_sided"] > 0),
        "g3_atom_identity": (
            contrasts["g3_atom"]["delta"] >= 0.5 * g and
            contrasts["g3_atom"]["lcb95_one_sided"] > 0),
        "nll_guard": (
            contrasts["nll_full_minus_free"]["ucb95_one_sided"] <=
            0.01 * metrics[free_arm]["nll"]),
    }
    if not gates["g1_prior_incremental"]:
        verdict = "DISTANCE_RESIDUAL_NOT_INCREMENTAL"
    elif not gates["g2_exact_pair_incremental"]:
        verdict = "MARGINAL_OR_SLOT_RECALIBRATION_ONLY"
    elif not (gates["g3_residue_identity"] and gates["g3_atom_identity"]):
        verdict = "ONE_SIDED_IDENTITY_SHORTCUT"
    elif not gates["nll_guard"]:
        verdict = "ORDERED_SCORE_GAIN_MISCALIBRATED"
    else:
        verdict = "R0C_EXACT_DISTANCE_RESIDUAL_CONFIRMED"

    system_rows = []
    for system in confirmation:
        row = {"source_entry_id": system.entry, "component": system.component}
        for arm in metrics:
            row[f"{arm}_rps"] = metrics[arm]["system_rps"][system.entry]
            row[f"{arm}_nll"] = metrics[arm]["system_nll"][system.entry]
        system_rows.append(row)
    write_jsonl(output / "confirmation_systems.jsonl", system_rows)
    checkpoint_hashes = {
        path.name: sha256_file(path) for path in sorted(output.glob("*.pt"))
    }
    report = {
        "schema": "MetaSieve.R0C.ExactDistanceTraining.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "device": args.device,
        "elapsed_seconds": time.perf_counter() - start,
        "execution": {
            "resumed_full_checkpoints": bool(args.resume_full),
            "elapsed_excludes_original_full_training": bool(args.resume_full),
        },
        "config": {
            "seeds": list(SEEDS), "epochs": EPOCHS,
            "learning_rate": LEARNING_RATE, "weight_decay": WEIGHT_DECAY,
            "ordinal_weight": ORDINAL_WEIGHT, "rank": RANK,
            "pair_budget": PAIR_BUDGET, "max_batch": MAX_BATCH,
            "gradient_clip": GRADIENT_CLIP,
        },
        "counts": {
            "train": len(train), "validation": len(validation),
            "confirmation": len(confirmation),
        },
        "selected_epoch": selected_epoch,
        "mean_full_validation_rps_by_epoch": mean_validation,
        "training_history": training,
        "confirmation_metrics": {
            arm: {"rps": value["rps"], "nll": value["nll"]}
            for arm, value in metrics.items()
        },
        "per_seed_confirmation": {
            arm: {seed: {"rps": value["rps"], "nll": value["nll"]}
                  for seed, value in seeds.items()}
            for arm, seeds in per_seed.items()
        },
        "free_arm": free_arm,
        "delta_star": delta_star,
        "g": g,
        "contrasts": contrasts,
        "gates": gates,
        "verdict": verdict,
        "affinity_labels_used": False,
        "production_migration_authorized": False,
        "checkpoint_sha256": checkpoint_hashes,
        "confirmation_systems_sha256": sha256_file(output / "confirmation_systems.jsonl"),
        "input_sha256": {name: (sha256_file(path / "manifest.json")
                                 if path.is_dir() else sha256_file(path))
                         for name, path in paths.items()},
        "runner_sha256": sha256_file(Path(__file__)),
        "model_sha256": sha256_file(Path(__file__).parents[1] / "meta_fewshot" / "affinity_pair_field.py"),
    }
    (output / "RESULT.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r0b-panel", required=True)
    parser.add_argument("--r0b-geometry", required=True)
    parser.add_argument("--r0b-proteins", required=True)
    parser.add_argument("--r0b-ligands", required=True)
    parser.add_argument("--r0c-panel", required=True)
    parser.add_argument("--r0c-geometry", required=True)
    parser.add_argument("--r0c-proteins", required=True)
    parser.add_argument("--r0c-ligands", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--prefit", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--resume-full", action="store_true")
    args = parser.parse_args()
    result = run(args)
    print(json.dumps({
        "selected_epoch": result["selected_epoch"],
        "confirmation_metrics": result["confirmation_metrics"],
        "gates": result["gates"], "verdict": result["verdict"],
        "elapsed_seconds": result["elapsed_seconds"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
