# S0–S4 — Structural Self-Supervision Programme: Terminal Report

Date: 2026-08-08. Branch `research/ssl-b2-structural`.
Environment: `drug`, Python 3.11.15, numpy 1.26.4, gemmi 0.7.5, rdkit 2023.09.6,
torch 2.6.0+cu124, transformers 4.46.3, MMseqs2 (repo-pinned binary).
Label reads: DAVIS `0`, recipient `0`, ChEMBL37 affinity `0`, **any affinity
value `0`** — S0–S4 are entirely label-free.
GPU: frozen ESM-2 t30 **inference only**. **No GPU training was performed.**

---

## 0. Terminal verdict

```text
POSE_FREE_DEPLOYMENT_INPUTS_INSUFFICIENT
```

Two of the six named 3D mechanism channels are genuinely observable from
deployment-available inputs and clearly beat capacity-matched random features.
**None of them is protein-specific.** The attribution control shows the ligand
alone explains everything the full model explains; the protein embedding
contributes nothing, and for hydrophobic burial it measurably hurts.

`POSE_FREE_STRUCTURAL_MECHANISM_IDENTIFIED` is refused.
`STRUCTURAL_TEST_UNDERPOWERED` is ruled out by the S3 power analysis.
`MODEL_REALIZATION_FAILURE` is not indicated: the probe reaches `R2 = 0.30`
where the information exists, so it is not underfitting.

---

## 1. What each stage established

### S0 — freeze
Commits `3281780`, `b659f23`, `12a2765` verified. Regression `73 passed`.
`model/`, production `scripts/`, `theory/`, CSMO, Band and the fixed mesh
untouched throughout.

### S1 — data roles, licences, exposure
`DATASET_ROLE_REGISTRY.json`, `LICENSE_AND_PROVENANCE_AUDIT.md`,
`STRUCTURAL_EXPOSURE_AUDIT.json`.

All 10,468 `pilot20k` PDB ids are treated as **exposed**, including P1B's own
val and test partitions, because those were consumed as P1B's evaluation. The
RCSB query (X-ray, `<= 2.5 A`, bound non-polymer, one protein entity, released
`>= 2024-01-01`) returned 15,003 entries with **exposed overlap `0`**, so S1 did
not fail closed. 1,476 acquired under CC0-1.0; 1,162 usable; 1,118 with a
parsable CCD ligand.

**Final independent block: 1,118 complexes, 621 MMseqs40/80cov protein
clusters, 586 Bemis-Murcko scaffolds.** For comparison, every protein-side
interval in XP1/XP2 rested on **8** kinase groups.

PLINDER **not used**: its annotation licence is separate from its GPL-2.0
software and policy forbids the full 400k pull, so annotations are computed
locally from raw CC0 coordinates instead. PDBbind **not used**: redistribution
terms unverified. Tier-C DTI/KG sources **not downloaded**.

### S2 — frozen 3D teacher
`TEACHER_CONTRACT.md`, `s2_teacher.py`, `TEACHER_REPRODUCIBILITY_AUDIT.json`.

Six named channels from raw holo coordinates: directional H-bond, signed
electrostatics, hydrophobic burial, aromatic orientation, steric overlap, pocket
burial. Audited on 78 complexes:

| property | measured | tolerance |
|---|---|---|
| rotation + translation invariance | `7.6e-15` | `1e-9` |
| atom permutation invariance | `1.6e-14` | `1e-9` |
| determinism | `0.0` exactly | `1e-9` |
| channel degeneracy | none of six | — |

**Defect found and fixed after the contract was frozen.** The RCSB filter used
`polymer_entity_count_protein == 1`, which admits homo-oligomers (one *entity*,
many chains); the `O(nP^2)` bonded-neighbour step then reached 10.5 GB. The
teacher now restricts to whole protein residues within `10 A` of the ligand.
Since every channel is defined at `<= 8 A` and bonded neighbours are `<= 1.8 A`,
this is **exactly equivalent**, not an approximation. The full invariance audit
was re-run afterwards rather than assumed, and channel statistics were unchanged
to three decimals.

### S3 — split and power
`STRUCTURAL_SPLIT_MANIFEST.json`, `STRUCTURAL_POWER_ANALYSIS.json`.

Closures: protein homology (MMseqs 40%/80cov), exact sequence, ligand CCD
identity, Bemis-Murcko scaffold, PDB entry, P1B exposure (zero by construction).
Inference unit is the homology cluster or scaffold, never an atom-residue pair.

**124 effective independence units per fold; minimum detectable `R2 = 0.02` at
100% detection** under the registered rule (95% cluster-bootstrap lower bound
`> 0`). The S6 floor was frozen at `0.02` before the test set was opened. The
nulls below are therefore real nulls.

### S4 — observability

**Registered deviation.** Rather than running the P1B checkpoint (which would
require rebuilding its full preprocessing contract for 1,476 new structures),
S4 asks whether the six channels are reachable from ESM-2 + ECFP **at all**.
P1B's predicted contact/distance geometry is itself a function of exactly those
inputs, so a negative bounds the entire sequence+2D class rather than one model.
A reader should therefore interpret the negative as *"no sequence+2D route
reaches these channels protein-specifically"*, not *"P1B specifically fails"*.

| channel | `R2` vs mean | vs random | **vs deranged protein** |
|---|---|---|---|
| hbond_directional | **+0.268 [+0.166, +0.378]** | **+0.366 [+0.222, +0.505]** | +0.037 [−0.015, +0.084] |
| hydrophobic_burial | **+0.299 [+0.162, +0.454]** | **+0.307 [+0.167, +0.431]** | −0.006 [−0.035, +0.024] |
| steric_overlap | +0.079 [−0.012, +0.137] | +0.058 [−0.005, +0.132] | +0.070 [−0.009, +0.184] |
| pocket_burial | +0.055 [−0.029, +0.138] | +0.052 [−0.029, +0.137] | +0.029 [−0.045, +0.098] |
| aromatic_orientation | +0.033 [−0.013, +0.069] | +0.071 [+0.014, +0.133] | +0.010 [−0.026, +0.041] |
| electrostatic_signed | +0.027 [−0.015, +0.048] | +0.044 [−0.004, +0.103] | −0.006 [−0.031, +0.018] |

Two channels are strongly observable. **No channel beats the deranged-protein
control** — every interval spans zero.

### S4b — attribution

The deranged null could mean "the protein matters but ESM cannot express it" or
"the protein does not matter here". The attribution control separates them on
the same cells and split, with identical hyperparameter selection.

| channel | LIG-ONLY | PROT-ONLY | BOTH | **BOTH − LIG** |
|---|---|---|---|---|
| hbond_directional | **+0.266 [+0.160, +0.391]** | +0.009 [−0.026, +0.047] | +0.268 | +0.0015 [−0.0379, +0.0363] |
| hydrophobic_burial | **+0.331 [+0.204, +0.485]** | −0.017 [−0.044, +0.004] | +0.299 | **−0.0321 [−0.0621, −0.0052]** |
| pocket_burial | +0.077 [+0.027, +0.125] | −0.027 [−0.080, +0.023] | +0.055 | −0.0200 [−0.0732, +0.0327] |
| aromatic_orientation | +0.044 [+0.006, +0.078] | +0.000 [−0.010, +0.010] | +0.033 | −0.0137 [−0.0401, +0.0107] |
| steric_overlap | +0.043 [+0.009, +0.114] | +0.012 [−0.047, +0.052] | +0.079 | +0.0385 [−0.0265, +0.1425] |
| electrostatic_signed | +0.036 [+0.006, +0.054] | +0.001 [−0.013, +0.012] | +0.027 | −0.0101 [−0.0286, +0.0065] |

**Protein-only is ~0 for all six channels. The ligand alone explains everything.**
Adding the protein buys `+0.0015` for the best channel (CI spans zero) and
**costs `−0.032` for hydrophobic burial with a CI excluding zero**.

The mechanism is not mysterious: aggregate hydrophobic burial and aggregate
directional H-bonding are largely set by the ligand's own composition — a large
apolar ligand buries apolar surface in whatever pocket it occupies, and a ligand
with many polar atoms makes many polar contacts. These aggregates are ligand
descriptors wearing structural clothing.

---

## 2. Why GPU training was not authorised

S4's registered authorisation required four conditions. Three are met:

| condition | status |
|---|---|
| teacher reproducible | **met** — machine precision |
| at least one channel non-degenerate | **met** — all six |
| frozen inputs retain information above random controls | **met** — two channels |
| deep head has an identifiable target the lightweight probe underfits | **not met** |

The fourth fails twice over. The probe already reaches `R2 = 0.30` where
information exists, so it is not underfitting; and the information it reaches is
ligand-side, so a Mechanistic Distillation Network trained on this teacher from
these inputs would learn **ligand chemistry**. That is precisely the population
shortcut the programme forbids, and it would produce a statistic that fails the
S6 correct-vs-deranged criterion by construction.

Training it anyway would have produced an impressive-looking reconstruction
number and no biology. The honest action is to not train it.

---

## 3. The ten questions

| # | question | answer |
|---|---|---|
| 1 | Which 3D channels are reproducible from experimental structure? | **All six.** Invariances at machine precision, none degenerate |
| 2 | Which are observable from sequence+2D? | **Two** — directional H-bond (`+0.268`) and hydrophobic burial (`+0.299`), both beating random features |
| 3 | Which require protein structure and a predicted pose? | Unresolved, but the S4b attribution says the *protein-specific* part of all six is unreachable from sequence+2D. If it is reachable at all, it needs a pose |
| 4 | Which contain affinity direction beyond ligand-only? | **Not tested** — no channel passed the structural Gate, so S8 was never entered and no affinity value was read |
| 5 | Which are correct-protein specific? | **None.** Every deranged interval spans zero; protein-only `R2 ~ 0` |
| 6 | Which survive protein-group and scaffold closure? | The ligand-side signal does; the protein-specific signal does not exist to survive |
| 7 | Which replicate across sources? | Not applicable — nothing was admitted to replicate |
| 8 | Can `k <= 5` identify a target-specific section? | Not reached; requires an admitted channel |
| 9 | Deep biological meta-learning, or a shortcut? | The observable part is a **ligand descriptor**, explicitly a shortcut. No deep biological model is claimed |
| 10 | Admitted into the probability-law operator? | **No.** `K(B(z)F(z))` was never scored, deliberately |

---

## 4. Honest limitations

1. **S4 is an upper bound on the input class, not a P1B test.** A reader wanting
   "does P1B specifically carry this" does not have that answer.
2. **Aggregate channels only.** The teacher stores pair-local contributions, but
   S4 scored the six complex-level aggregates. A pair-local target could in
   principle be more protein-discriminative; that is a genuine untested variant,
   though the protein-only `R2 ~ 0` makes it unpromising.
3. **Mean-pooled ESM.** A residue-local protein representation (rung `B1`) was
   not tested. This is the most defensible remaining pose-free variant.
4. **One ligand per complex** (the largest non-solvent component), so multi-ligand
   and cofactor-mediated mechanisms are outside scope.
5. **Waters, metals and protonation are excluded**, so channels 1 and 2 describe
   direct polar contacts only and cannot represent water-mediated bridges.
6. **1,118 complexes** is a moderate corpus; the power analysis says it is ample
   for `R2 >= 0.02`, but rarer mechanisms may be under-represented.

---

## 5. What the decision tree indicates next

The pose-free branch is closed for protein specificity. The tree's conditional
`S7` pose-aware pilot is the indicated continuation, and its precondition is
satisfied in the qualified sense that the *protein-specific* part of the teacher
is not observable from sequence+2D.

`S7` would need, and this programme did not attempt: a frozen pose predictor
(DiffDock / Uni-Mol class) with pinned version and licence, PoseBusters
validity checks, 1,000–3,000 protein-ligand pairs chosen without affinity
information, and the native-pose / predicted-pose / pose-free / randomised-pose /
deranged-protein comparison. Its central risk is already quantified by this
programme: XP2 showed that noise in a *predicted* input degrades the
target-specific component faster than the shared one, so pose error is the
principal threat to `S7`, not compute.

Before `S7`, the cheaper and more informative variant is rung `B1`:
**residue-local** protein states against **atom-local** ligand states with
pair-local channel targets, which is the one pose-free representation this
programme has not tested and the one whose failure would be genuinely decisive
for the whole sequence+2D class.

---

## 6. Artifacts

| file | content |
|---|---|
| `DATASET_ROLE_REGISTRY.json` | every source, tier, role, licence, exposure |
| `LICENSE_AND_PROVENANCE_AUDIT.md` | licence, redistributability, attribution, non-use reasons |
| `STRUCTURAL_EXPOSURE_AUDIT.json` | P1B exposure policy and counts |
| `TEACHER_CONTRACT.md` | frozen channel definitions and invariance claims |
| `TEACHER_REPRODUCIBILITY_AUDIT.json` | measured invariances and degeneracy |
| `STRUCTURAL_SPLIT_MANIFEST.json` | closures and independence unit |
| `STRUCTURAL_POWER_ANALYSIS.json` | MDE and the frozen S6 floor |
| `S4_OBSERVABILITY_AUDIT.json` | the six-channel observability table |
| `S4B_ATTRIBUTION.json` | ligand-side vs protein-side attribution |
