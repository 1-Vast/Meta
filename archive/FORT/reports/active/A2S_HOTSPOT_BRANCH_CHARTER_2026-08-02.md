# A2S-HOTSPOT — branch charter

**Branch:** `research/a2s-hotspot-sparse-20260802` (created 2026-08-02, forked from
`research/openmut-medip-audit-20260728`; nothing committed).
**Opened because:** Gate R0 established that the k ≤ 5 selection ceiling is real (+0.074 CI) but
unreachable from any observable derived from the current representation, and that the gap
(≈ 5 AUC points against ≈ 25 needed) cannot be closed within the k ≤ 5 setting by tuning.
**Predecessors:** `A2S_RIP_GATE_R0_DECISION_2026-08-02.md`,
`A2S_MODE_GENERALIZATION_DECISION_2026-08-02.md`.

---

# 1. The general biological principle

Four independent measurements in this programme have failed in the same way, and one principle
retrodicts all four.

## 1.1 The principle

> **Binding hot spots.** The free energy of a protein–ligand interaction is not distributed
> evenly over the contact surface. A small subset of contacts contributes most of the binding
> energy, and the rest contributes comparatively little.

Established by alanine-scanning of the human growth hormone–receptor interface
(Clackson & Wells, *Science* 267:383, 1995), generalised across interfaces by
Bogan & Thorn (*JMB* 280:1, 1998), and standard in fragment-based and hot-spot-mapping practice
since. The medicinal-chemistry statement of the same structure is older:

> **Free–Wilson additivity.** Within a congeneric series, substituent contributions to potency
> are approximately additive and position-specific (Free & Wilson, *J. Med. Chem.* 7:395, 1964).

**INFERENCE — the transferable consequence.** If a target discriminates among ligands through a
few interaction determinants, then in a basis whose coordinates *are* interaction determinants a
target's response head is **sparse**, and different targets are sparse on **different**
coordinates. It says nothing about the head being low-rank, low-dimensional, or shared.

## 1.2 What the principle retrodicts

| Measurement | Previously read as | Predicted by the hot-spot principle |
|---|---|---|
| **G2** — source-head spectrum nearly flat (top 3 = 34.7 %), rank-2 projection retains −6 % | "no structure exists" | **A set of sparse vectors with different supports has a flat covariance spectrum.** A flat spectrum is the *signature* of heterogeneous sparsity, not its absence |
| **G3** — pooled ESM-2 predicts the head at −0.019 [−0.073, +0.019] | "protein is uninformative" | A sequence-averaged embedding cannot express *which few residues* form the hot spot. The principle predicts pooled features fail and residue/pocket-level features are required |
| **G4** — dense empirical-Bayes learning-curve knee at **k ≈ 10** | an unexplained empirical constant | Sparse recovery needs `O(s·log(d/s))` measurements. Measured effective sparsity in the current basis is `s ≈ 8`, `d = 26` ⇒ `8·log(26/8) ≈ 9–10`. **The theory produces the measured number** |
| **R0b** — per-compound margin AUC 0.555, posterior covariance adds 0.004 | "no usable uncertainty" | A dense estimate spreads its mass over all coordinates, so every per-compound quantity derived from it is blunt. Sharpness requires concentration, which requires sparsity |

**FACT — the direct test, already run** (`A2S_RIP_GATE_R0_DECISION` §5, 52 targets / 50 components).
Truncating each target's own head to its top-`s` **coordinates** versus its top-`s` **source
principal directions**:

| `s` | coordinate-truncated | rank-truncated |
|---:|---:|---:|
| 8 | **+0.0343 [+0.0152, +0.0523]** | +0.0085 [−0.0043, +0.0234] |
| 26 | +0.0542 | +0.0542 |

Eight coordinates retain 63 % of the gain; eight principal directions retain 16 % with an interval
crossing zero. Across 52 targets the top-weighted coordinate takes **20 distinct values of 26**.

## 1.3 The forward prediction that defines this branch

> **HYPOTHESIS H0.** There exists a representation — interaction-determinant space rather than
> generic ligand-descriptor space — in which the per-target response head is **2–3 sparse**. In
> that representation the learning-curve knee moves from `k ≈ 10` to `k ≈ 5`, and the k ≤ 5
> deployment budget becomes sufficient.

This is quantitative and falsifiable: `s·log(d/s)` must fall to ≈ 5, and the G4 curve must shift
left. **The branch does not need a new mechanism; it needs a new basis, and it already has the
estimator, the certification layer and the controls.**

---

# 2. What carries over

| Retained | Why |
|---|---|
| The frozen cross-fitted base, the v2 balanced lock, the episode machinery | unchanged substrate; results stay comparable |
| The empirical-Bayes head estimator with a source-measured prior | no free parameters; only its basis changes |
| **Cross-task threshold transfer (R0c PASS)** | the certification layer is sound and reusable the moment a sharp statistic exists |
| The **selection ceiling** (+0.074 CI at 40 % coverage, k=5) | the prize this branch is trying to reach |
| **Magnitude-matched wholesale** and **random-selection-at-matched-coverage** controls | both fired in R0; mandatory in every successor |
| The G4 learning curve | the baseline to shift left; the primary success criterion |
| Stratification by similarity to the *estimator's own* training rows | G1 showed "scaffold-cold" ≠ "chemically distant" |

## What is retired

A2S-RIP as a mechanism (P4/P5 fired); A2S-MODE (G2); the rank-`m` code family incl. A2S-IDA (G2);
per-pair transport reliability (Q2). All remain as named baselines.

---

# 3. Branch gates, in order

**H1 — sparsity in the current basis. DONE (§1.2).** `s ≈ 8`, supports heterogeneous, coordinate
truncation beats rank truncation. The principle is supported in-basis; the sparsity is not yet
tight enough.

**H2 — do the supports themselves transfer?** Are the coordinates a target uses predictable from
anything label-free (protein family, pocket residues, assay context)? If yes, the sparsity prior
becomes target-conditioned and `s log(d/s)` drops further. Cheap: one pass on the existing heads.

**H3 — build the interaction-determinant basis.** Candidates, in order of local availability:
pharmacophore-feature counts × pocket-property descriptors; Free–Wilson substituent–position
indicators derived from the target's own series; residue-level (not pooled) protein features from
the existing `target_esm2_segments32.npz`, KLIFS pocket residues where a target has them. Repeat H1
in each. **Success criterion: effective `s` falls to 2–3 with top-`s` truncation retaining ≥ 60 %.**

**H4 — sparse recovery at k ≤ 5.** LASSO / OMP / spike-and-slab with a meta-learned sparsity prior
versus the dense empirical-Bayes head, on the identical G4 protocol. **Success criterion: the
learning-curve knee moves left — a k=5 lower bound above the 0.005 MDE where G4 measured −0.003.**

**H5 — re-run R0 in the new basis.** The selection ceiling and the margin AUC are recomputed. The
mechanism is admitted only if the margin AUC rises materially above 0.555 *and* the implementable
rule beats the magnitude-matched control.

**Stop rule.** H3 fails to reduce `s` below ≈ 6, or H4 fails to move the knee ⇒ report the measured
sparsity bound and the representation requirement, and stop. The honest deliverable would then be:
*the k ≤ 5 budget is insufficient for target-specific ranking adaptation on open ChEMBL pKi, and the
measured requirement is `s·log(d/s)` ≈ 10 labels in the best available representation* — a
quantitative statement no prior work in this area has produced.

---

# 4. Risks

1. **The principle may be right and unreachable.** Hot spots are a property of the *complex*;
   this programme has no complex structures for most targets, and G3 showed pooled sequence carries
   nothing. If the determinant basis cannot be built from available open data, H3 fails on data, not
   on theory. That is the most likely failure and it is diagnosable in one pass.
2. **Free–Wilson requires congeneric series.** Its additivity holds *within* a series. The
   admissible-stratum work (Q1) showed the corpus does contain local series; the scaffold-disjoint
   protocol deliberately evaluates across them. A basis that only works within series would reproduce
   the transport class's distance limit and must be tested against it explicitly.
3. **Sparsity is not identifiability.** `O(s log(d/s))` assumes incoherent measurements. Support
   compounds drawn from one series are highly coherent. The support policy therefore becomes part of
   the mechanism: a diversity-aware draw is predicted to matter here, and that prediction is
   registered now.
4. **Power.** 50–54 probe components, MDE80 ≈ 0.005–0.010. Unchanged.

---

# 5. Standing constraints (unchanged)

`locked` role and the A2S recipient roster stay sealed. `probe` is a development role. Every result
names its support policy and its relation stratum. Every ranking result is reported beside the
frozen base, the magnitude-matched control and a random-selection control. No `model/` promotion
without a frozen protocol and a one-time `locked` evaluation.
