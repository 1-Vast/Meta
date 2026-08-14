# L-CIPF + ELMT: Failure Diagnosis and Replacement Record

## Scope

This record supersedes the former TERM active-path proposal. It is a code and
data-flow audit, not a claim of SOTA performance. The BindingDB Ki corpus is
still the only main-task dataset. Cartesian geometry is disabled because the
production episode materializer has no coordinates or packed geometry edges.

## Why the former model failed

1. **CIPF was globally dense.** Every ligand atom was paired with every protein
   residue and normalized by global softmax. With no pocket/contact mask, true
   interaction evidence was diluted by non-contact residues; this is consistent
   with the near-zero wrong-protein gap.
2. **Primitive slots were not identifiable.** Slots were freely permutable and
   sign/scale could be absorbed by downstream weights. The additive atom and
   residue bias did not encode complementary cross terms. Label permutation was
   therefore nearly as good as the real assignment.
3. **TERM had a label shortcut.** Evidence was one input among task identity,
   ligand-change and primitive features. The router could generate correction
   without using support labels, explaining `real ~= permuted`.
4. **Residual accounting was inconsistent.** `y-f0` was used as TERM evidence,
   while the correction was added after level calibration, so level and SAR
   could explain the same residual.
5. **Confidence was a dead path.** Detached entropy removed gradients from the
   confidence/prior branch; reliability was effectively global rather than
   query-specific.
6. **The training budget and selection criterion were misaligned.** The old
   short run trained several towers from scratch, while checkpoint selection
   minimized full MSE and could select a strong level-only model with a dead
   meta branch.
7. **Geometry was not part of the result.** BindingDB materialization supplies
   sequence/2D inputs only, so Cartesian experiments could not explain or repair
   the observed failure.

## Replacement

**Localized CIPF (L-CIPF)** uses ligand-conditioned top-k residue localization,
sequence-slot chemistry, and bounded primitive responses. Four channels are
weakly anchored by reliable sequence/ligand complements (charge, aromatic and
hydrophobic); remaining channels are learned residual primitives. No H-bond or
3D claim is made from unavailable features.

**Evidence-Locked Meta-Transport (ELMT)** makes support residual values the only
label-dependent path. Keys, attention and signed transport directions are
label-blind; the result is linear in support residual values. Reliability uses
evidence coherence, query-support coverage and support count. There is no
Transformer confidence branch, ridge, matrix inverse, inner loop or deployment
optimization.

Residuals are conserved: a scalar target-level shrinkage is removed first, and
ELMT sees only the leftover residual. Binding supervision is a contrastive
likelihood over correct versus counterfactual assignments; k=1 uses an equal-
magnitude residual flip.

Training now has a short representation warm-up followed by ELMT meta-training.
Validation uses an admission score that penalizes full MSE worse than the
level-only cut and penalizes absent permutation/counterfactual binding gap.

## Verification

- Synthetic held-out shared/level/private ELMT gates: **3 passed**.
- Focused regression tests: **38 passed**.
- `compileall`: passed.
- Real BindingDB CUDA smoke (`2` steps, Phase A -> Phase B, k=0/1): passed.
  The smoke is implementation validation only, not a performance result.
- Cartesian main path: intentionally disabled; no geometry performance claim.

The next admissible experiment is a one-seed BindingDB A/B/C/D comparison. A
three-seed or sealed result is not authorized until correct-label binding,
level-cut improvement, and wrong-protein controls pass on that development run.

## BindingDB development run

Run: `lcipf_elmt_bindingdb_dev_small_seed20260814` (CUDA, 20 steps, reduced
development width, k sampled from 0/1/2/3/5). The run completed, but failed the
admission criteria:

| metric | value |
|---|---:|
| full MSE | 1.8660 |
| level-only MSE | 1.4296 |
| zero-shot MSE | 3.0800 |
| binding gap (permuted - full) | -0.0127 |
| wrong-protein gap | +0.0160 |
| complete foreign gap | +1.8856 |

The negative binding gap means the correct support labels did not beat the
counterfactual assignment. Full prediction was also worse than the scalar
level-only baseline. The foreign gap is not treated as proof of label binding,
because that intervention replaces the complete donor episode. Consequently,
this model is not authorized for a three-seed or sealed evaluation. The full
machine-readable result is in
`report/meta_fewshot/lcipf_elmt_bindingdb_dev_small_seed20260814/RESULT.json`.

## Medium-model power test

The default model was increased to 3,191,427 trainable parameters:
`hidden_dim=192`, `task_dim=48`, `ligand_layers=4`, `pair_dim=96`, four pair
blocks, 24 latent slots, three episodes per optimization step, and 20 query
examples. Cartesian remains disabled. A 10-step CUDA run completed with a
peak allocation of approximately 3.99 GB.

The resulting test metrics were:

| metric | value |
|---|---:|
| full MSE | 1.2431 |
| level-only MSE | 1.9248 |
| binding gap (permuted - full) | -0.0319 |
| wrong-protein gap | -0.0043 |

The larger model uses substantially more compute and improves aggregate MSE
in this short run, but both binding and wrong-protein controls fail. This is a
capacity/optimization diagnostic, not an admission result. The raw output is
`report/meta_fewshot/lcipf_elmt_medium_power_smoke_20260814/RESULT.json`.
