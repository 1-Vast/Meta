# S7_L2B — development report: exact-residue localization

Date: 2026-08-09.
Preregistration: `research/s7_l2b_r0r/PREREG_S7_L2B_UNIFIED.md`,
SHA-256 `2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92`,
**committed as `ce186f4` before any model existed and before any AP was computed.**

```text
Data / provenance ........... PASS   (R0R-1, hash-verified, bit-identical)
Closure construction ........ RESOLVED (protein partition + ligand filter)
Evaluator contract .......... PASS
Trainability control ........ PASS
Matched baseline B4 ......... TRAINED AND FROZEN
Registered Gate outcome ..... FAIL CLOSED (1 of 5 gates met the effect size)
Escalation rule ............. FIRED — B5 authorised
B5 (frozen ESM2-650M) ....... blocked on weight acquisition, see preflight audit
Confirmation cohort ......... SEALED, never opened
```

## 1. The closure correction that made this measurable

R0R-2 measured that union-**merging** ligand identity into the inference
partition produces a **93.19 %** giant component over 14,589 complexes, making
component inference impossible. Ablation localised the cause to the ligand side:
protein closure alone gives 1,994 components with a 3.11 % largest, while adding
exact ligand graph alone drives it to 76.22 %.

The correction is not a relaxation. Ligand closure is enforced as a
**disjointness filter between train and held-out**, which is *stricter* between
the two sets than a merge would have been, while leaving the partition usable:

| Set | Complexes | Components | Largest | Positives |
|---|---:|---:|---:|---:|
| train | 9,758 | — | — | 151,065 |
| **held-out A** (protein- and exact-ligand-graph-disjoint) | **2,415** | **196** | 15.2 % | 36,046 |
| held-out B (additionally scaffold-disjoint) | 1,881 | 160 | 16.5 % | 29,073 |

The claim scope is stated up front: generalisation is over **proteins**; the
ligand-disjoint strata additionally test unseen ligands.

## 2. Contract checks that had to pass first

**Atom mapping.** MONN's `atom_name` lists deposited ligand atoms *including
hydrogens*, while `mol_dict` is a heavy-atom molecule, and `atom_idx` is a plain
identity range — **not** a mapping into the molecule. The correct mapping is the
rank of a slot among the heavy positions. Validated per record: **14,586 of
14,589 pass**, 3 quarantined on heavy-count mismatch, **zero** positive edges on
hydrogens, **zero** residue indices out of range.

**Evaluator.** AP is computed in float64 with ties broken by a fixed
`(residue_index, atom_slot)` order. Self-test: AP is bit-identical under a random
permutation of the flattened matrix, returns exactly 1.0 when all rows are
positive, and returns undefined (excluded) when a complex has no positives.

**Trainability.** The identical pipeline recovers a *known* function of the
frozen inputs at macro-AP **0.7588** against a prevalence of 0.0081. Optimisation
and the objective are therefore **not** the defect for anything below.

## 3. Result — held-out A, complete-matrix AP, 196 components

| Arm | What it is | macro-AP |
|---|---|---:|
| `B0` | prevalence | 0.00250 |
| `BL` | **ligand-only** | 0.00450 |
| `BP` | wrong protein | 0.00485 |
| `BM` | motif shuffle | 0.00530 |
| `BX` | wrong ligand | 0.00768 |
| **`B4`** | **non-PLM residue baseline** | **0.02295** |

| Gate | Contrast | Δ | LCB95 | Threshold | |
|---|---|---:|---:|---:|---|
| G1 | B4 − B0 | +0.02045 | +0.01719 | 0.02 | **PASS** |
| G2 | B4 − BL | +0.01845 | +0.01523 | 0.02 | FAIL |
| G3 | B4 − BP | +0.01810 | +0.01487 | 0.02 | FAIL |
| G4 | B4 − BM | +0.01765 | +0.01408 | 0.02 | FAIL |
| G5 | B4 − BX | +0.01527 | +0.01238 | 0.02 | FAIL |

**Registered outcome: FAIL CLOSED.** Four of five contrasts fall below the
preregistered 0.02 absolute-AP effect size. The threshold was frozen before any
model existed and is **not** relaxed. Statistical significance does not override
a preregistered practical-effect requirement.

Held-out B (scaffold-strict, 160 components) reproduces the ordering:
B4 0.02153 against BL 0.00452, Δ = +0.01701 [LCB +0.01399]. **No sign reversal**,
so the result is not an artifact of ligand-scaffold overlap.

## 4. What the controls say — read carefully

**Every one of the five contrasts is directionally positive with a lower bound
well above zero.** The identifiability requirement is met in *direction* on every
control; it is the *effect size* that is not met.

* Substituting a **wrong protein** collapses B4 from 0.02295 to **0.00485** —
  essentially back to the ligand-only level. Almost all of B4's advantage over
  ligand-only requires the correct protein.
* **Motif shuffle** collapses it to 0.00530, so residue order and position are
  load-bearing.
* **Wrong ligand** collapses it to 0.00768 — well above ligand-only but far below
  B4, so the correct ligand matters substantially.
* **Ligand-only** sits at 0.00450 against a prevalence of 0.00250.

### The auto-assigned label was wrong and is withdrawn

The runner emitted `S7L2B_LIGAND_ONLY_SHORTCUT` from a decision-logic defect: it
mapped "G2 did not pass" straight to the ligand-shortcut verdict, but G2 can fail
on effect size while being strongly directional, which is exactly what happened.
That label is contradicted by its own data — B4 is **5.1×** ligand-only, and
ligand-only barely clears prevalence.

Corrected: `S7L2B_PROTEIN_SIGNAL_IDENTIFIED_BELOW_PREREGISTERED_EFFECT_SIZE`.
Adding a descriptive label does **not** convert the outcome into a pass; the Gate
outcome remains FAIL CLOSED. Full adjudication in
`S7L2B_DEVELOPMENT_GATE_ADJUDICATED.json`; the raw run output is retained
unedited.

## 5. Failure localization

| Candidate | Verdict | Evidence |
|---|---|---|
| data | **not the cause** | six source hashes verified; corpus bit-identical across runs |
| closure | **not the cause** | 196 components, largest 15.2 %, ligand-disjoint from train |
| optimization | **not the cause** | trainability control recovers a known function at 0.759 |
| identifiability support | **not the cause** | all five contrasts have LCB > 0 |
| **representation** | **THE CAUSE** | explicit non-PLM residue features carry real but weak protein information (0.0229) |
| transfer | **not tested** | confirmation cohort sealed |

## 6. Escalation — preregistered, and now evidence-backed

PREREG §7 authorises `B5` only if `B4 − BL` fails G2, **or** `B4` macro-AP stays
below 0.10. **Both conditions hold** (G2 fails at +0.01845; B4 = 0.0229).

This is exactly what constraint 9 requires before adding a larger encoder: prior
experiments must *explicitly indicate missing information*. That evidence now
exists and is measured — the residue representation is the binding constraint,
and it is the only component B5 changes. `B5` must still clear G1–G5 **plus**
G6 (`B5 − B4 ≥ 0.02`, LCB > 0).

## 7. B5 status

Not run. Blocked solely on acquiring the frozen `esm2_t33_650M_UR50D` weights;
see `GPU_PREFLIGHT_AUDIT.json` for the measured cause. `PLM_RUNTIME_NOT_FEASIBLE`
is **not** claimed — the GPU contract was never exercised. Nothing here is
evidence about the PLM hypothesis.

## 8. Boundaries held

No affinity, DAVIS, recipient or few-shot support label was read. `theory/`,
`model/`, `contracts/`, production `scripts/`, CSMO, Band, simplex, positive
ridge, mesh and production `z` are unmodified. All new code is under `research/`.
The confirmation cohort was never scored, and R0R-3 publication/time closure
remains unbuilt, so it could not be opened even on a development pass. Only the
preregistration was committed; nothing was pushed.
