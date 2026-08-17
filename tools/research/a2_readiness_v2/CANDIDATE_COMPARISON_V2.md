# Candidate comparison v2 (2026-08-16)

Every candidate is scored against what Phases 1-3 actually measured, not
against what the previous cycle inferred. Two of the four premises the v1
package relied on did not survive.

## What changed under the candidates

| premise the v1 package used | status after Phases 1-3 |
|---|---|
| ordering is interaction-borne, not ligand-borne | ✅ **confirmed** (+0.1855 [+0.0566, +0.3236]) |
| the interaction branch's ordering is protein-inert | ✅ **strengthened** (all 5 donor strata, floor exactly 0) |
| the architecture can express protein-conditioned ordering; training removed it | ❌ **withdrawn** — random-init sensitivity is undirected (pairwise cosine −0.003) |
| the collapse is in the readout, so a better-aimed gradient reaches it | ❌ **relocated** to the fusion/pooling inside `ContactGrammar` |
| *(new)* some internal representation carries protein-conditioned SAR | ❌ **falsified** — gain over a capacity-matched permuted protein is +0.017 |
| *(new)* a protein-**independent** transferable SAR direction exists | ✅ **measured**, +0.2623 [+0.1295, +0.4055] |

## The candidates

### A2-min — protein-conditioned SAR moment update

`z = A_φ(e0)`, `c_S = mean_i r_i z_i`, `δ = η(k)⟨c_S, z_q⟩`.

**Rejected, by direct measurement of the representation it consumes.**

Phase 2 trained exactly the readout A2-min needs — a low-rank projection of
`e0` scored against protein conditioning — on `meta_train` component folds and
evaluated once on `meta_val`. The protein contributes +0.0167 in Δ-affinity
correlation over an identical model given a *permuted* protein, and the
protein-free variant beats the protein-conditioned one by 0.176. The same holds
on `mean_state` (−0.128) and on `max_state` (−0.040), the stage that retains the
most protein-differential of any.

A2-min's S1 gate required the benefit to *weaken* under wrong-protein
substitution. It strengthens (+0.1235 wrong vs +0.0862 correct). This is not a
tuning failure: the operator would be identifying a coordinate in a space where
the protein carries nothing identifiable.

**Not amendable.** The v1 package left open that A2 might be repaired by
pointing it at a better readout. Phase 2 tested four representations spanning
the whole trunk, including the two that Phase 3 identified as retaining the
most protein-differential. There is nothing to point at.

### CPC — centered protein counterfactual

Compute the existing protein contrast on the within-target **centered**
prediction, so `∂L/∂(protein_head) ≡ 0` and the level cannot satisfy it.

**Structurally valid; scientifically unsupported as drafted. Not authorized.**

The mechanism holds — the centering identity is verified on the real model
class (11 probes, `a2_readiness/tests/`), and F6/F8 confirm the incumbent's
uncentered form is satisfiable by a level shift. What has gone is the reason to
expect it to help:

* its motivating evidence was E3, that the architecture demonstrably expresses
  protein-conditioned ordering and training removes it. E3 is withdrawn;
* its target was E4, a readout invariance. The loss is upstream of every
  readout (F10), so a gradient aimed at the readout aims at the wrong place;
* the objective is satisfiable by **degrading the donor prediction** (F8), and
  by inventing arbitrary protein-dependent ordering uncorrelated with truth —
  which is exactly what an untrained network already does (Phase 1 §7);
* the effect it would need (Δ`r` ≈ 0.087) sits at 1.7× the measured
  same-configuration retraining noise on aggregate `r` (0.051, R14 screening).

One argument survives in its favour and should be stated plainly: **Phase 2
measured a *frozen* trunk.** A representation that was never trained to carry
protein-conditioned ordering not carrying it is weak evidence that it could not
be trained to. No objective in R0-R14 ever asked for it — every protein control
in the project is uncentered and therefore level-only (F7). So CPC addresses a
genuinely open question. It just no longer has a measured premise, a located
target, or a shortcut-proof gate, and `PREREGISTRATION_V2.md` supplies the last
of those but cannot supply the first two.

### Attention variants (support self-attention, query→support cross-attention, protein-conditioned support attention, low-rank/kernelised, uncertainty-gated, SAR-coordinate attention)

**Still rejected, and now for a second, stronger reason.**

The v1 argument was that all seven act at k≥1 while the bottleneck is at k=0
where the support set is empty. That stands. Phases 1-3 add: the protein path
is exactly invariant to residue-slot order (F9), so no attention over this
input can be biologically justified as reading a pocket, a contact, or any
ordered structural feature — and the attention that exists already changes by
146% under a protein swap while contributing nothing usable downstream. More
attention over the same unordered bag is not a mechanism, it is more capacity
pointed at a channel measured to be empty.

### Set aggregation without attention (the simpler control)

Unchanged as the correct baseline for any k≥1 mechanism. Fixed
Morgan/Tanimoto residual transport remains the strongest reproducible k≥2
comparator and is parameter-free apart from two scalars.

### Capacity increases on the trunk

**Not indicated.** Phase 3 shows the trunk *does* deliver a large protein
signal into `context` (Jacobian 1.13e+02) and loses it at the pooling. More
capacity upstream of a 3,400× attenuation adds parameters to a channel that is
already saturated with information nothing downstream uses.

### M0 / MSA protein-side calibration probe

Independent lane, untouched by this cycle. It addresses target-level
calibration, which is roughly half of the k=0 error, and it makes no
within-target ordering claim. Nothing in Phases 1-3 supports or weakens it. It
still requires the recorded MMseqs2 executable and a governed UniRef snapshot
before it may run.

### The one thing Phases 1-3 newly support

A **protein-independent** transferable SAR readout on `embed`: Δ-affinity `r`
+0.2623 [+0.1295, +0.4055] on held-out protein components, from 1,553
parameters, stable across every rank tested, and better than the raw ligand
encoder (+0.1188, unresolved).

This is not a candidate for the current objective as it stands, and must not be
promoted as one without answering the question that decides it: **does it add
anything to fixed Morgan/Tanimoto transport?** That comparator already exploits
ligand-side SAR continuity at k≥2 and is the incumbent's strongest few-shot
mechanism. The measurement needed is a matched head-to-head on the same pairs,
and it has not been made. It is also a probe on frozen features, not a trained
DTA model, so it is not comparable to any k=0 MSE on record.

It would also violate the standing objective in its current form: the mandate
requires protein-conditioned few-shot learning, and this is explicitly and
measurably not that.

## Scoreboard

| candidate | premise measured? | target located? | shortcut-proof? | verdict |
|---|---|---|---|---|
| A2-min | ❌ falsified on 4 representations | — | — | **rejected** |
| CPC (drafted) | ❌ premise withdrawn | ❌ aimed downstream of the loss | ❌ donor destruction unguarded | **not authorized** |
| CPC (revised, `PREREGISTRATION_V2`) | ❌ still absent | ⚠️ would act on the trunk end-to-end, not the readout | ✅ gated | **blocked on a prerequisite** |
| attention family (7 variants) | ❌ wrong k, empty channel | — | — | **rejected** |
| non-attention set aggregation | ✅ (Tanimoto, R0/R6) | n/a | ✅ | retained as comparator |
| trunk capacity | ❌ | ❌ | — | not indicated |
| M0 / MSA | independent lane | n/a | n/a | unchanged |
| protein-independent SAR direction | ✅ measured | ✅ `embed` | untested vs Tanimoto | **new evidence, not a candidate** |

## The strongest competing explanation for everything above

That A0's trunk was never asked for protein-conditioned ordering, so its
absence is a property of the training objective rather than of the data or the
architecture — and that the whole double-cold protocol (346 training targets,
9-21 ligands per target, no complex geometry, sequence + 2D inputs) does not
contain enough within-target signal for *any* model to learn protein-specific
SAR at k=0.

Phases 1-3 cannot separate these. Both predict exactly what was measured. The
prerequisite in `PREREGISTRATION_V2.md` §1 is designed to separate them, and it
is cheap: it asks whether protein-conditioned within-target ordering is
learnable *at all* on this data, with the trunk free and the objective
explicitly demanding it, before any admission-grade arm is trained.
