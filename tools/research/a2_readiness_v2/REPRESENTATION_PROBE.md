# Phase 2: the A2 premise, tested on the representations A2 actually consumes

Authority: `REPRESENTATION_PROBE_meta_val.json`. Frozen A0 (seed 20260815),
double-cold `meta_val` scored **once**, after every choice was fixed on
`meta_train` component folds. `meta_test` unreachable.

## Why this stage exists

The v1 audit measured the final scalar and concluded A2's premise was
falsified. That was not a valid inference: **A2-min never consumes the scalar.**
It builds `z = A_φ(e0)` from an internal representation. A trunk whose scalar
readout is protein-inert could still carry protein-conditioned information
upstream, in which case A2 should be *amended* — pointed at the readout —
rather than rejected.

This stage tests the internal tensors directly.

## Part A — the ligand-differential of every internal representation

For each representation, remove the within-target mean (the level, which
cannot order anything) and ask how much the remainder — the part that *can*
order ligands — changes when the protein is replaced by the nearest legal
cross-component donor.

| representation | width | level shift | **differential cosine** | differential shift | 95% CI on cosine |
|---|---:|---:|---:|---:|---|
| `occupancy` | 24 | 0.0027 | **1.0000** | 0.0052 | [1.0000, 1.0000] |
| `mean_state` | 192 | 0.0021 | **1.0000** | 0.0031 | [1.0000, 1.0000] |
| `max_state` | 192 | 0.0016 | 0.9984 | 0.0503 | [0.9976, 0.9991] |
| `embed` (**the A2 plan's `e0`**) | 96 | 0.0011 | 0.9998 | 0.0146 | [0.9997, 0.9999] |
| `section` | 48 | 0.0006 | 0.9998 | 0.0156 | [0.9997, 0.9999] |
| `interaction` (scalar) | 1 | 0.0219 | 1.0000 | 0.0044 | [1.0000, 1.0000] |
| `ligand` *(protein-blind reference)* | 192 | 0.0000 | 1.0000 | 0.0000 | [1.0000, 1.0000] |

The `ligand` row is the calibration of the instrument: a representation that is
protein-blind *by construction* reads exactly 1.0000 / 0.0000. Every
protein-conditioned representation reads 0.998–1.0000 — closer to the
protein-blind reference than to anything that could be called
protein-conditioned ordering.

`max_state` is the exception worth naming: at cosine 0.9984 and differential
shift 0.0503 it retains **16× more** protein-differential than `mean_state`.
It is the best surviving candidate, and Part B tests it.

**Part A does not by itself reject A2.** A 1.5% change in a 96-dimensional
representation could in principle be a low-dimensional, highly informative
direction. That is what Part B is for.

## Part B — is any of it a *usable, transferable* SAR coordinate?

The smallest thing that could work, trained by ordinary gradient descent
(Adam) — no closed form, no pseudoinverse, no inner loop:

```text
delta_hat(P, L_i, L_j) = s · ⟨ g(P), U (e_i − e_j) ⟩
```

`U: R^D → R^R` is A2-min's `A_φ`; `g` maps the pooled protein into the same
space. Trained on `meta_train` within-target ligand pairs, rank selected on
`meta_train` component folds **separately for each configuration**, scored once
on `meta_val`. The metric is the within-target correlation between predicted
and true pair differences, averaged per target then per component.

The decisive comparison is against the identical model whose `g(P)` is a
learned **constant**: if conditioning on the protein buys nothing, there is no
protein-conditioned SAR coordinate for k≤5 labels to identify, whatever
operator is wrapped around it.

### `embed` — the A2 plan's `e0`

meta_train fold curves (held-out components):

| rank | protein-conditioned | shared direction |
|---:|---:|---:|
| 4 | 0.041 | 0.179 |
| 8 | **0.131** | 0.178 |
| 16 | 0.063 | **0.182** |

meta_val, scored once:

| arm | rank | trainable params | Δ-affinity `r` | 95% CI | sign acc |
|---|---:|---:|---:|---|---:|
| `protein_conditioned` | 8 | 42,313 | +0.0862 | [−0.0805, +0.2560] | 0.519 |
| **`shared_direction`** (no protein) | 16 | **1,553** | **+0.2623** | **[+0.1295, +0.4055]** | 0.564 |
| `wrong_protein` | 8 | 42,313 | +0.1235 | [−0.0540, +0.3030] | 0.536 |
| `reference_protein` | 8 | 42,313 | −0.1771 | [−0.3224, −0.0288] | 0.441 |
| `label_shuffled` | 8 | 42,313 | +0.0392 | [−0.0990, +0.1666] | 0.490 |
| `random_feature` | 8 | 42,313 | −0.0331 | [−0.1309, +0.0559] | 0.480 |
| `random_projection` | 8 | 41,545 | +0.0944 | [−0.0657, +0.2654] | 0.516 |
| **`protein_permuted`** (capacity-matched null) | 8 | 42,313 | +0.0695 | [−0.0696, +0.2182] | 0.536 |

**Protein conditioning gain: +0.0167 over a capacity-matched permuted
protein; −0.1761 against a protein-free direction.**

`protein_permuted` is the control that closes the obvious objection. The
protein head carries ~42k parameters against the constant's 16, so a
protein-conditioned arm that loses could be losing to overfitting rather than
to the absence of signal. `protein_permuted` keeps the identical architecture,
parameter count, optimiser and schedule and destroys only the *correspondence*
between protein and target. It scores +0.0695 against the real protein's
+0.0862. The real protein is worth 0.017 in correlation — indistinguishable
from a protein assigned at random, and indistinguishable from shuffled labels
(+0.0392).

### Every other representation

| representation | protein-conditioned | shared direction | protein-permuted | **gain vs permuted** |
|---|---:|---:|---:|---:|
| `embed` | +0.0862 | **+0.2623** | +0.0695 | +0.0167 |
| `mean_state` | −0.1320 | −0.0476 | −0.0039 | **−0.1281** |
| `max_state` *(highest retention)* | +0.0459 | +0.0783 | +0.0860 | **−0.0401** |
| `ligand` *(protein-blind control)* | +0.1562 | +0.1188 | +0.1464 | +0.0098 |

`max_state` and `ligand` were added after the first run and are labelled as
post-hoc. Both were added in the conservative direction: `max_state` is the
representation Phase 3 identified as retaining the *most* protein-differential,
so testing it can only help A2; `ligand` is protein-blind, so testing it can
only take credit away from my own positive finding. Neither is treated as a
selected representation.

**No representation at any stage of the trunk shows a protein-conditioning gain
larger than 0.017, and two show a negative one.**

## What A2 required, and what was measured

| A2 admission requirement (`NEXT_RESEARCH_PLAN_A2_MOMENT`, §3 S1) | measured |
|---|---|
| an internal representation predicts ligand-pair affinity differences on held-out protein components | ✅ `embed`, +0.2623 [+0.130, +0.406] |
| that prediction **weakens under wrong-protein substitution** | ❌ `wrong_protein` +0.1235 vs `protein_conditioned` +0.0862 — it *improves* |
| it **exceeds** ligand-only, shuffled-label and random-feature controls | ❌ +0.0862 vs shuffled-label +0.0392, permuted-protein +0.0695 — inside the noise |
| it is not confined to repeated or highly similar ligands | zero scaffold overlap by split construction; the ordering that exists is concentrated in the two least-novel terciles (Phase 1 §5) |

**A2's premise is falsified on the representations A2 consumes, not merely on
the endpoint scalar. It is not amendable: there is no protein-conditioned
information anywhere in the trunk to point a better readout at.**

## The positive finding, stated at its true scope

`embed` carries a **protein-independent** transferable SAR direction:

* Δ-affinity `r` = **+0.2623 [+0.1295, +0.4055]** on held-out protein
  components — resolved;
* from **1,553 trainable parameters**, a rank-16 linear projection and one
  constant vector;
* stable across every rank tested on `meta_train` folds (0.178–0.182), so it is
  not a hyperparameter artifact;
* **better than the raw ligand encoder**: `ligand`'s shared direction is
  +0.1188 [−0.0511, +0.2915] (unresolved), fold curve 0.133–0.138 against
  `embed`'s 0.178–0.182.

That last point is the interesting one. `embed` is built from
`cat(ligand, mean_state, max_state, wide_summary, occupancy)`. The contact
channels *do* improve within-target SAR ordering over the ligand encoder alone
— but the direction that reads them is protein-independent. The interaction
machinery is functioning as a **better ligand descriptor**, not as a
protein-conditioned one.

Four things this is **not**, and must never be reported as:

1. not meta-learning — no support label enters it;
2. not protein-conditioned — that is precisely what was ruled out;
3. not a performance result — it is a probe on frozen features on `meta_val`,
   not a trained DTA model, and it is not comparable to any k=0 MSE on record;
4. not necessarily new capability — fixed Morgan/Tanimoto residual transport
   already exploits ligand-side SAR continuity at k≥2 (Stage R0/R6). Whether
   this direction adds anything to that comparator is **untested**.

## Effective resolution of this probe

On the protein-blind `ligand` representation, `protein_permuted` (+0.1464,
interval excluding zero) outscored `shared_direction` (+0.1188, interval
including zero) — on a representation where the protein provably carries
nothing. Differences of ~0.03–0.05 between arms here are noise. Only the
`embed` shared-direction result (+0.2623, more than 3× the noise band, stable
across ranks and folds) is treated as a finding.

## Commands

```bash
conda run -n drug python -m tools.research.a2_readiness_v2.extract_features --split meta_train --output tools/research/a2_readiness_v2/features
conda run -n drug python -m tools.research.a2_readiness_v2.extract_features --split meta_val   --output tools/research/a2_readiness_v2/features
conda run -n drug python -m tools.research.a2_readiness_v2.representation_probe --features tools/research/a2_readiness_v2/features --output tools/research/a2_readiness_v2/REPRESENTATION_PROBE_meta_val.json
```
