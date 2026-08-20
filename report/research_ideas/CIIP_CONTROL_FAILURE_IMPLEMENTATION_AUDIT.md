# CIIP-1A Control Stage: Implementation and Failure-Analysis Record

Date: 2026-08-20  
Scope: `tools/research/stageCIIP_potential_bridge/`  
Status: research analysis only; production `model/` and `scripts/` are unchanged.

## 1. Purpose

This document records exactly what the current CIIP-1A control experiment computes,
what its result establishes, and which conclusions remain invalid. It is not a new
preregistration and must not override `PREREGISTRATION_STAGE1_CONTROLS.md`.

The experiment asks a narrow question:

> Does a mutation-centered local ESM representation produce a transferable,
> ligand-dependent WT-to-variant contrast beyond matched protein and annotation
> controls?

It does **not** directly test general cold-target DTA, deployable mutation-free
protein representations, or binding affinity in the Ki/Kd sense.

## 2. Frozen data contract

- Source: Duong-Ly kinase inhibitor panel.
- Endpoint: raw percent inhibition; values are not relabeled as Ki, Kd, pK, or
  DeltaDeltaG.
- Admitted single-point WT/variant pairs: 65.
- Matched ESM-covered subset: 49 pairs.
- Original split on the 49 pairs: 32 train, 8 validation, 9 test.
- Test coverage: 6 parent kinases and 9 mutation pairs.
- Common ligands are retained independently for every WT/variant pair.
- Target for pair `i`:

  ```text
  d_i,l = y_variant_i,l - y_WT_i,l
  c_i,l = d_i,l - mean_l(d_i,l)
  ```

  The centering removes a mutation-wide scalar shift and retains only the
  ligand-dependent profile. Consequently, a mutation that changes every ligand
  by the same amount is intentionally invisible to this target.

## 3. Production of protein and ligand inputs

The control trainer loads `DATA1A.json`, `DATA1A.npz`, `DATA2X2.json`,
`DATA2X2.npz`, and the frozen `q1_esm_cache.npz`.

Ligand input is a 2048-bit ECFP4 vector. The positive protein input is a
640-dimensional radius-6 mean of residue-level ESM-2 embeddings extracted at
the verified mutation coordinate for WT and mutant sequences.

The implemented potential is:

```text
eP = ReLU(p_enc(P))
eL = ReLU(l_enc(L))
s(P,L) = sum_r alpha(eP)_r * psi(eL)_r
f(P,L) = mu + bP(eP) + bL(eL) + s(P,L)
g(Pwt, Pvar, L) = s(Pvar,L) - s(Pwt,L)
chat = g - mean_L(g)
```

The centered arm trains only the squared contrast loss on `chat`. The additive
`mu`, `bP`, and `bL` branches are therefore not expected to receive gradients;
this is mathematically correct for the frozen centered-only objective.

## 4. Control arms

The same potential, optimizer, budget, split, checkpoint rule, and keyed random
streams are used for every arm.

1. `oracle_local_esm_correct`: WT and variant windows at the verified mutation
   coordinate.
2. `family_preserving_shuffle`: variant windows permuted within the same parent;
   WT windows remain fixed.
3. `random_local_window`: a keyed non-mutation position with `abs(q-true_pos)>6`
   is used for both WT and variant windows.
4. `ligand_only`: both protein inputs are zero, so the antisymmetric contrast is
   identically zero.
5. `ligand_invariant_shift`: an explicit zero predictor; this isolates a
   protein-main-effect shortcut.
6. `random_protein`: the variant window is replaced by a window from another
   parent while the WT window and ligand rows remain recipient-side.
7. `free_pairwise`: an antisymmetric pair MLP, retained only as an expressivity
   diagnostic; it is not an integrable potential and is not deployable.

## 5. Observed single-seed result

The result artifact is `CONTROL_RESULT.json`. The test aggregate is:

| arm | R2 | Spearman | sign accuracy | nonconstant test pairs |
|---|---:|---:|---:|---:|
| correct local ESM | 0.007 | 0.331 | 0.702 | 9/9 |
| family shuffle | -0.004 | -0.007 | 0.498 | 9/9 |
| random local window | **0.129** | 0.320 | **0.705** | 9/9 |
| random protein | 0.000 | 0.001 | 0.493 | 9/9 |
| ligand-only | 0.000 | undefined | 0.494 | 0/9 |
| free pairwise | 0.027 | 0.125 | 0.562 | 9/9 |

Correct minus random-window observed pair-mean R2 is approximately `-0.122`,
with parent-cluster bootstrap interval approximately `[-0.457, 0.033]`.
Correct minus family and correct minus random-protein intervals also cross zero.
Correct and random-window nonconstant coverage are both 9/9.

The feature-level annotation audit is positive: mean correct WT/variant delta
norm is approximately `0.531`, versus `0.0267` for the matched random window,
with all 49 pairs larger for the correct site. This proves that the verified
mutation position changes the ESM representation more strongly; it does not
prove that the change predicts ligand-specific activity.

## 6. What the result means

Under the frozen verdict rules, the current evidence supports a provisional
`ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED` outcome, pending formal adjudicator output.
The defensible scientific statement is narrower:

> On the oracle-covered subset, the correct mutation-centered local ESM window
> did not outperform the matched random-window, family-shuffle, or random-protein
> controls on transferable centered ligand-response prediction.

This is not a biological falsification. It is not evidence that protein-conditioned
interaction is absent in BindingDB or in other assay systems.

## 7. Root-cause hypotheses, ranked

### H1: the estimand is too aggressive for the available signal

Centering removes mutation-wide effects and leaves only the ligand-specific
reordering. The Duong-Ly endpoint is a single-concentration functional assay,
not a direct binding-affinity measurement. ATP competition, catalytic activity,
construct effects, and noise can remain after centering.

### H2: random-window is not a pure null

ESM is contextual. A residue token far from the annotated mutation can still
reflect a sequence change through the transformer context. Therefore the random
window tests explicit coordinate dependence, not complete mutation independence.

### H3: the split is pair-cold, not parent/protein-cold

The same parent kinases occur in train, validation, and test. The experiment
therefore measures mutation-pair generalization within partially observed parent
families, not strict unseen-protein cold-start transfer.

### H4: statistical power is small

The test contains only 9 pairs and 6 parents, and only one seed was run. A
parent-cluster bootstrap cannot compensate for six independent clusters. Several
correct-arm gains are concentrated in ABL1/KIT, while other parents are near zero
or negative.

### H5: the current model lacks an explicit mutation/effect decomposition

The model supplies a local window difference only through `alpha(Pvar)-alpha(Pwt)`.
It has no explicit mutation-site delta channel, mutation-agnostic main-effect
head, or auxiliary task that stabilizes the mutation representation. The centered
loss alone can encourage small interaction outputs.

## 8. Required follow-up diagnostics before any successor training

1. Add a contextual-leakage control: construct WT/variant embeddings with the
   mutation token masked or locally replaced outside the target site, and compare
   full-context, site-only, and context-masked deltas.
2. Add a ligand-agnostic mutation head and retain centered ligand-specific output
   as a separate head. Report both components; do not add them implicitly.
3. Re-run the controls with parent-disjoint and scaffold-disjoint splits.
4. Report per-parent and per-mutation-class results, not only the nine-pair mean.
5. Use multiple keyed seeds only after the above diagnostic is frozen.
6. Separate functional inhibition from direct binding/DeltaDeltaG endpoints; never
   merge their labels into one regression target.
7. Add a true same-parent matched negative: preserve WT/variant sequence context
   while independently permuting ligand labels, rather than relying only on a
   synthetic random-protein pair.

No successor should be promoted to `CIIP-1B`, BindingDB Bridge, or production
integration until a preregistered control matrix passes on a valid split and the
representation is available without mutation-coordinate oracle metadata.

## 9. Literature anchors

- eSIG-Net uses residue-level WT/MT ESM embeddings at the mutation site and
  combines discrepancy learning with the original interaction prediction to avoid
  a trivial discrepancy-only solution:
  https://www.nature.com/articles/s41592-026-03086-x
- ESM documentation distinguishes residue-level (`per_tok`) from mean-pooled
  embeddings and provides ESM-1v mutation-effect scoring:
  https://github.com/facebookresearch/esm
- PremPLI predicts mutation-induced protein-ligand affinity changes using
  structurally defined complexes and DeltaDeltaG-style targets:
  https://www.nature.com/articles/s42003-021-02826-3
- CS-DTA documents entity-disjoint and similarity-controlled cold-start splits:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC13161074/
- Duong-Ly reports a functional HotSpot kinase assay and percent remaining
  activity, not a direct Ki/Kd endpoint:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4740242/

## 10. Governance

This file is an analysis record. It does not change:

- `PREREGISTRATION_STAGE1_CONTROLS.md`;
- `CONTROL_RESULT.json`;
- the CIIP-1A verdict rules;
- `model/` or production `scripts/`;
- authorization for CIIP-1B, BindingDB Bridge, or production integration.

The next formal action is to run the frozen adjudicator, write the stage report,
and synchronize `history.md`, `task.md`, and `report/EVIDENCE_LEDGER.md`.
