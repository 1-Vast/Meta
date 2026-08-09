# S7/L2B Phase 2A — audit-only attribution of B5: synthesis

Date: 2026-08-10.
Preregistration `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`,
SHA-256 `4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e`,
frozen before any Phase 2A metric was computed, with amendments 01–03 each
frozen before the phase it governs. Repository commit `623602e`.

```text
TERMINAL VERDICT

    LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING

AUTHORIZED NEXT ACTION

    preregister one ligand-conditioned residue residual head
```

Nothing was trained. No affinity value, DAVIS, KIBA or recipient label was read.
No frozen surface was modified. Nothing was committed or pushed.

## The one-sentence result

**The MONN labels are strongly ligand-conditioned at the residue level, B5 is
not, and neither the labels nor B5 show identifiable coupling at the
residue–atom edge level** — so the correct repair is a residue-differential
head, not a pair-coupling head.

---

## 1. What changed relative to the Phase 1 reading

Phase 1 reported that a wrong ligand retained 92.5% of B5's residue AP and
concluded "mostly generic pocket". That conclusion was about **B5**. It was
silently being carried as if it were also a statement about the **corpus**. It
is not.

The wrong-ligand control used an *arbitrary foreign* ligand. Phase 2A replaced
it with the scientifically correct comparison — a **real alternative ligand of
the same exact construct**, against a **noise floor measured from the data
itself**: two crystal structures of the same construct with the *same* ligand.

| comparison (component-macro Jaccard of residue masks) | value |
|---|---:|
| same construct, **same** ligand, different crystal (noise floor) | 0.636 |
| same construct, **different** scaffold-distinct ligand | 0.416 |
| **T1 ΔJ, paired over 292 closure components** | **+0.258 [LCB95 +0.234]** |

The registered minimum meaningful effect was 0.05. The observed effect is **five
times** that, and it survives on held-out A alone (+0.191 [+0.117], 27
components) and under the unpaired fallback (+0.220 [+0.198]).

The variation is also **chemistry-associated**, not noise: within a construct,
mask dissimilarity rises with ligand chemical distance at Spearman
**ρ = +0.322 [LCB +0.299]**, and the ligand-permutation control gives a median
per-construct p of 0.03 with 84.8% of constructs on the positive side.

Both registered criteria (`T1` and `T6`) pass, so `TEACHER_GENERIC_POCKET_ONLY`
is refuted on this corpus.

A calibration worth keeping: the replicate floor is **0.636, not 1.0**. Roughly
44% of the observed alternative-ligand mask difference is ligand-attributable;
the remaining ~56% sits at the same level as pure crystallographic/annotation
variation between replicates of the identical system. Any differential objective
must be built expecting that noise.

## 2. Where B5's ligand information actually lives

Decomposing the sealed logits on the complete, uniformly weighted mask
(`G = Proj_W(1, α_r, β_a) + C`, orthogonality achieved at `1.2e-9` against a
registered tolerance of `1e-8`):

| arm | full | residue marginal | atom marginal | additive | coupling `C` |
|---|---:|---:|---:|---:|---:|
| **B5** | **0.0698** | 0.0404 | 0.0051 | 0.0398 | **0.0113** |
| B4 | 0.0232 | 0.0162 | 0.0055 | 0.0151 | 0.0064 |
| BX5 wrong ligand | 0.0197 | 0.0360 | 0.0033 | 0.0277 | 0.0035 |
| BP5 wrong protein | 0.0046 | 0.0046 | 0.0049 | 0.0056 | 0.0035 |
| BL ligand-only | 0.0057 | 0.0033 | 0.0057 | 0.0057 | 0.0030 |

This localises the ligand dependence precisely, and the answer is not the one
Phase 1's headline implied:

* B5's **residue marginal** barely notices the ligand — a wrong ligand retains
  `0.0360/0.0404 = 89%` of it. That is the generic-pocket term.
* B5's **coupling residual** is where the ligand matters — a wrong ligand
  retains only `0.0035/0.0113 = 31%` of it, a **3.3×** drop.

So B5 does use ligand identity, but only through the pair term, and that term is
small in absolute size.

## 3. Why the coupling still fails its registered bar

| contrast | Δ | LCB95 | margin | |
|---|---:|---:|---:|---|
| B5 coupling − degree-preserving rewiring null | +0.0060 | +0.0046 | 0.01 | **FAIL** |
| B5 coupling − BX5 coupling | +0.0079 | +0.0062 | 0.01 | **FAIL** |

Both are **statistically clearly above zero** and both are **below the
preregistered practical margin of 0.01**. That is the honest reading: the effect
is real but small, and the bar was set before the number was seen. It is not
rounded up into a pass.

The teacher's own edge coupling fails too. Against the registered rewiring
specification (checkerboard swaps, 100×E burn-in, 30×E between samples, 20
independent rewires, **zero** degree-preservation violations), the teacher's
marginal-orthogonal coupling sits at **median z = +0.413** with 63.4% of
complexes above their own null, against a registered threshold of z ≥ 2.0. This
reproduces I-2's earlier median of +0.41 under a different, stricter rewiring
implementation — an independent confirmation, not a repetition.

Mixing was checked rather than assumed: edge overlap with the original decays
1.000 → 0.334 → 0.298 → 0.292 → 0.292 at 0, 1, 5, 10, 30 swaps per edge, and
successive independent samples overlap at 0.292 — i.e. the chain reaches its
degree-constrained plateau well before the first sample is taken.

`TC = false` and `BC = false` therefore give rule 7:
**`LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`**.

## 4. The headroom number that reframes the programme

The well-posed label-fitted additive ceiling — the AP obtainable by recovering
the *true* residue and atom margins exactly — is **0.389**.

```text
B5 full            0.0698   = 17.9% of the additive ceiling
B5 additive part   0.0398   = 10.2% of the additive ceiling
B5 residue margin  0.0404   = 19.8% of the true residue-margin ceiling (0.204)
```

**The bottleneck is not coupling. It is the residue marginal itself**, which is
recovered at roughly one fifth of what the labels permit. Spending the next
stage on a pair-coupling head would be optimising a term worth 0.011 while
leaving 0.32 of additive AP unclaimed. This is registered here as the reason the
authorized repair is residue-side.

*(The logistic Rasch additive null registered in section 7 is reported but
flagged `rasch_converged: false`. The design is completely separated — the
matrix is 0.07% positive, so almost every residue row and atom column has no
positive at all and its coefficient diverges. Its AP is **not** a valid ceiling
and is not used as one; the least-squares additive projection above is the
well-posed object.)*

## 5. Label semantics: audited, and not ambiguous

| check | result |
|---|---|
| interaction-type census | 212,556 typed edges: hydrophobic 68,153; H-bond 57,647; π-stacking 34,917; salt bridge 25,447; water bridge 17,427; π-cation 7,754; halogen 1,211 |
| indirect (water-mediated) fraction | **8.2%**, below the 20% registered threshold |
| metal-mediated | 0.0% — absent from this taxonomy |
| registered sensitivity: T1 with water bridges removed | ΔJ **+0.278 [+0.253]** — the conclusion **strengthens**, it does not reverse |
| dense-distance comparator | built locally on **1,909** complexes from already-governed mmCIF, median mapped sequence identity **1.000** |
| PLIP positives within 5.0 Å of a ligand heavy atom | **88.1%**; only 9.0% lie beyond, consistent with the water-bridge fraction |
| geometry ⊄ PLIP | 46% at 5.0 Å — expected: PLIP applies chemical *and* geometric criteria, so the labels are a strict subset of proximity. This is **not** evidence of missing positives |
| second frozen interaction tool | **ABSENT** — that specific comparison remains **UNRESOLVED**, recorded rather than substituted |

`LABEL_SEMANTICS_AMBIGUOUS` is not demonstrated. PU learning and a soft teacher
remain **unauthorized**, exactly as the registration requires when ambiguity is
not positively shown.

The comparator deserves one note: the registration allowed recording this
question as unresolved *if no comparator existed*. One did — 2,068 MONN entries
already had local coordinates. Declaring it unresolved would have been a
fabrication in the convenient direction, so amendment 03 required it to be
built.

## 6. Contract and data integrity

Phase 0 passed fail-closed on all seven checks over 26 hashed artifacts. The
load-bearing one was `C3`: Phase 1's marginal decomposition indexed the
B5-family memmaps with the B4-family offset table. Both tables were rebuilt
independently and proved **identical key-for-key and offset-for-offset** (the
ESM availability filter dropped zero records), so the Phase 1 B5/BX5/BP5
marginal numbers were correctly aligned. Every sealed prediction hash matched
its recorded manifest. The evaluation mask is the complete `n_res × n_atoms`
matrix with uniform weights — 52,062,975 cells, 36,237 positives, density
`6.96e-4` — which is what makes classical double centering admissible here; the
general weighted-ALS solver was used regardless and agreed with it to `1.1e-14`.

The corpus is data-identifiable with large margins: 14,585 records, 2,846 exact
constructs, 1,994 closure components (largest 3.1% of records), **1,093**
constructs carrying scaffold-distinct ligand pairs across **779** components,
**292** components carrying both replicate and alternative pairs, **323,410**
scaffold-distinct within-construct pairs. Label-blind power for ΔJ = 0.05 is
≥ 0.81 at the paired-component count even under the most pessimistic registered
σ of 0.30.

## 7. Boundary classification

| class | status |
|---|---|
| DATA/LABEL INSUFFICIENCY | **excluded** for residue-level ligand conditionality; **confirmed** for edge-level coupling — the teacher itself does not carry identifiable pairing beyond its margins |
| BIOLOGICAL REPRESENTATION FAILURE | **confirmed and localised**: the labels are ligand-conditioned, B5's residue marginal is not (89% wrong-ligand retention) |
| OBJECTIVE/OPTIMIZATION FAILURE | not tested here; Phase 2A trains nothing |
| SUPPORT-SECTION NON-IDENTIFIABILITY | not reached |

## 8. What this does not license

No affinity, ranking, transfer, few-shot or `z`-admission claim is made or
implied. A structural PASS does not establish affinity semantics. The frozen
operator `A(F,z) = K(B(z)F(z))` is unchanged and no bridge theorem is claimed.
The confirmation cohort was not opened. Real ChEMBL/BindingDB training, DAVIS,
KIBA, recipient data, few-shot adaptation, production `z`, CSMO, Band, mesh and
P2–P4 all remain frozen.

## 9. Machine-readable artifacts

`PHASE2A_PREREGISTRATION_HASH.json`, `PHASE2A_INPUT_MANIFEST.json`,
`PHASE2A_DATA_IDENTIFIABILITY_CENSUS.json`, `PHASE2A_CONSTRUCT_GROUPS.json`,
`PHASE2A_TEACHER_CONDITIONALITY.json`,
`PHASE2A_MARGINAL_COUPLING_AUDIT.json`, `PHASE2A_COMPONENT_TABLES.json`,
`PHASE2A_LABEL_SEMANTICS.json`, `PHASE2A_VERDICT.json`, and the consoles
`pa3_console.txt`, `pa4_console.txt`. Code under `research/s7_l2b_r0r/pa*.py`.

## 10. Chronology limitation, stated plainly

Commit authorization was not granted for this run. The preregistration and its
three amendments are anchored by SHA-256 and are embedded in every output
artifact, but they carry **no git commit timestamp**. That is a strictly weaker
chronological guarantee than the one behind `ce186f4` / `139effd`, and the
Phase 2A verdict should be read with that attached until the files are
committed.
