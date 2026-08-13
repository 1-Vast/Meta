# HyperSAR / PBCNet2-inspired optimization record

## Outcome

The active cold-target few-shot path no longer uses ridge regression, a linear solve, or test-time
gradient adaptation. The retained candidate is a single-stage episodic neural model with:

- atom--residue bipartite interaction fields and task-gated low-rank pair updates;
- permutation-invariant residual-bound Set2Code inference;
- query-specific Siamese reference comparison over ligand and protein--ligand endpoints;
- support-anchored absolute level prediction;
- neural query--support matching SAR with a support-only leave-one-out reliability gate.

The PBCNet2.0 connection is an architectural inspiration: reference--query relative recognition and
fine-grained interaction comparison. It is not a reproduction of its coordinate-based Cartesian
rank-2 TensorNet because this deployment corpus has no complex coordinates for every target--ligand
pair.

## Verification

- Environment: `drug`, CUDA.
- Final repository verification: 264 tests passed.
- Active-path static scan: no `ridge`, `torch.linalg.solve`, Cholesky, or pseudoinverse calls.
- Independent review: the two P1 findings (override shape ambiguity and external nested-manifest
  metadata validation) were fixed. A final review also renamed the partial foreign-code control so
  it cannot be mistaken for complete foreign-support replacement; no P0/P1 remains.

`foreign_code_state` replaces only the transient target code. Recipient matching SAR, reliability,
level, protein, and query remain fixed. It diagnoses the code branch only; `permuted_state` is the
complete label-binding counterfactual for the recipient support set.

## Development evidence

All numbers below are development experiments, not confirmatory gate authorization.

| Candidate | Best validation MSE (pK^2) | Test full MSE | SAR-cut MSE | SAR gain | Binding permutation gap |
|---|---:|---:|---:|---:|---:|
| global HyperSAR, 40 steps | 1.030 | 2.080 | 2.080 | +0.0001 | -0.0023 |
| Siamese endpoint, 40 steps | 1.165 | 1.910 | 1.917 | +0.0068 | +0.0001 |
| residual-bound stop-gradient, 120 steps | 0.920 | 1.703 | 1.699 | -0.0037 | +0.0008 |
| ungated matching, 80 steps | 1.129 | 1.780 | 1.766 | -0.0143 | +0.0907 |
| support-anchored matching, 80 steps | 1.046 | 1.396 | 1.367 | -0.0287 | -0.0166 |
| LOO reliability + matching loss, 60 steps | 1.081 | 1.375 | 1.367 | -0.0079 | -0.0009 |

Positive SAR gain means `MSE(SAR-cut) - MSE(full) > 0`. The best positive development gain was
small and did not persist after the support-label binding and stronger evaluation changes. The
support-anchored level equation substantially improved total few-shot MSE and approached the
support-mean baseline, but the final structural SAR branch still had a negative marginal effect.

## Scientific status

Engineering and protocol checks pass. The requested performance claim does not: under the executed
small/medium single-stage budgets, the structural modulation branch has not been shown to be an
important positive performance source, and the model has not established outstanding cold-target
few-shot performance. The result artifacts retain `G2`, `G3a`, and `G3b` authorization as false.

The next scientifically valid step is a predeclared multi-seed nested-k comparison against strong
baselines and `SAR-cut`, using the existing evaluator and unchanged meta-test manifest. It should not
be replaced by more post-hoc architecture search on the same test observations.

## Detailed defect and risk register

### D1 — the core SAR branch has not produced transferable positive utility

**Evidence.** In the final LOO-matching run, full MSE is 1.3749 while SAR-cut MSE is 1.3670, so
`SAR-cut - full = -0.0079`. The reliability-only run is also negative (-0.0104), as are the
support-anchored matching (-0.0287) and 120-step residual-bound (-0.0037) runs.

**Interpretation.** The model can fit useful task level information, but the ligand-specific
support-conditioned correction is not reliably transported to held-out target components. This is
the principal unmet requirement: the proposed innovation is not yet an important positive source
of test performance.

**Consequence.** HyperSAR must remain a candidate module. It must not be described as SOTA,
outstanding, or causally validated.

### D2 — support-mean anchoring, not structural modulation, explains most few-shot improvement

Zero-shot MSE is 3.5119 in the final run, whereas support-mean level-only MSE is 1.3659 and full MSE
is 1.3749. Most of the apparent few-shot gain therefore comes from observing the target's affinity
level. The structural SAR magnitude is only 0.0152 pK on average and slightly worsens MSE.

This is not leakage: the support labels are legal deployment inputs. It is nevertheless a
mechanistic attribution problem. Reporting only full-versus-zero-shot would drastically overstate
the contribution of the novel interaction branch.

### D3 — support sets are too small and noisy to identify local SAR consistently

The model is trained with k in {1,2,3,5}. For k=1, centered residual evidence is exactly zero, so
only absolute level adaptation is possible. At k=2 or 3, a single noisy Ki value can reverse the
centered ordering. BindingDB measurements also mix laboratories and experimental protocols even
after restricting to standardized Ki. The LOO matching loss has only k-1 references per held-out
support point and becomes a high-variance learning signal.

The reliability gate attenuates this failure but cannot create missing information. Its mean final
SAR scale is about 0.066, showing that the model largely learns to suppress the novel branch.

### D4 — the structure branch is not genuinely atomic 3D recognition

The ligand tower uses a 2D molecular graph. The protein tower uses cached sequence-model residue
tokens. The rectangular pair field learns atom--residue relations without pair-specific complex
coordinates, distances, bond angles, or orientations. Consequently it cannot express the
Cartesian scalar/vector/tensor decomposition, rotation equivariance, or subtle geometry used by
PBCNet2.0.

This fallback is deployment-compatible for cold targets but limits atomic claims. The implementation
is PBCNet2.0-inspired reference--query recognition, not PBCNet2.0 replication. A true 3D extension
requires legally available predicted/docked complex coordinates plus explicit missing-geometry
controls; silently inserting predicted poses would change the input contract and evaluation claim.

### D5 — reference compression and matching use weak supervision

Set2Code compresses all support observations into one target reference before producing
query-specific low-rank adapter codes. A residual-bound moment prevents complete label ignorance,
but information can still be lost. The direct matching path avoids that bottleneck, yet it is
supervised only by affinity values: there is no pocket-contact or interaction-type label in the
active single-stage corpus.

Geometry pretraining is optional and disabled in the reported runs. Enabling a geometry checkpoint
would require a separate, explicitly declared comparison because it changes the available prior and
could be confused with multi-stage training.

### D6 — counterfactual evidence is incomplete and currently weak

`permuted_state` is the complete recipient-support label-binding counterfactual: it recomputes the
task code, matching correction, reliability, and level from permuted support labels. Its final gap
is -0.0009 MSE, meaning the correct binding is not measurably superior in this development run.

`foreign_code_state` is narrower. It replaces only the transient target code while keeping recipient
matching SAR, reliability, level, protein, and query fixed. It diagnoses the code/adapter path, not
complete foreign-support adaptation. A true foreign-support control would need to override the
entire support-conditioned state, including label-bound residuals and matching references.

### D7 — validation and test behavior are not stable enough for a performance claim

The 120-step experiment reached validation MSE 0.9196 at step 40 but test MSE 1.7026. Later variants
with better anchoring reached test MSE near 1.37 but worse validation minima. This rank reversal
indicates component heterogeneity and selection variance. The reported test population has only 24
episodes under the compact development setting, so changes of a few thousandths are not convincing
effects.

Required remedy: freeze one model, one manifest, one metric hierarchy, and run multiple model seeds
with component-level paired bootstrap confidence intervals. Architecture selection must end before
that confirmatory evaluation.

### D8 — current metrics emphasize squared error and do not establish ranking excellence

The smoke result contains MSE-based counterfactual summaries, but a medicinal-chemistry SAR claim
also needs per-target CI/Spearman, pairwise sign accuracy, calibration, and nested-k curves. These
are supported by the dedicated evaluator but have not been run at confirmatory scale for the final
candidate. Low global MSE can coexist with poor within-series ranking.

### D9 — memory headroom is narrow on the target GPU

The final medium configuration peaks near 7378 MB, and earlier 128/64 configurations exceeded 8 GB
or reached about 10.6 GB. On an 8 GB RTX 4060, allocator fragmentation or longer ligands can cause
OOM. The second query pass through the interaction trunk is the dominant cost.

Mitigations already present are compact ligand shards, active-atom cropping, pair chunking, AMP, and
96/48 hidden/pair widths. A future implementation should cache the non-adaptive prefix of the pair
trunk rather than recompute it, provided numerical equivalence is tested.

### D10 — the evaluation evidence is development-only and affected by repeated search

Several architectures were selected after observing the same development test diagnostics. Even
though the code labels all gates false, the final test numbers are no longer an untouched estimate
of generalization. They are useful for falsification and debugging, not publication claims.

A new sealed confirmation manifest or an external cold-target benchmark is required. The current
meta-test labels must not be used for further architecture or loss decisions.

### D11 — assay scope is narrow

The active implementation specializes to standardized Ki in one governed BindingDB corpus. It has
no assay-context encoder and makes no validated claim for Kd, IC50, EC50, or cross-assay transport.
Combining these endpoints without an explicit observation model would add label noise and invalidate
the current scalar interpretation.

### D12 — model complexity has increased faster than confirmed utility

The final candidate contains Set2Code attention, residual-bound moments, Siamese query conditioning,
low-rank pair adapters, direct support matching, support-anchored level transfer, and a LOO
reliability gate. Each component is testable, but their combined code and optimization surface are
larger than the measured positive effect.

The next iteration should be an ablation decision, not another additive module. Retain a component
only if multi-seed component-level evidence shows a positive marginal contribution over the
support-anchored baseline.

## Required acceptance gates before claiming success

1. Full must beat support-mean and support-anchored SAR-cut baselines with a positive component-level
   confidence interval on at least three fixed model seeds.
2. Correct support binding must beat permuted labels, and complete recipient support must beat a
   fully replaced foreign-support state.
3. Positive utility must appear across nested k={1,2,3,5}, with k=1 reporting zero structural SAR by
   construction rather than being pooled with k>1.
4. Within-target CI/Spearman and pairwise sign accuracy must improve alongside MSE.
5. A declared 3D variant must separately report coordinate availability, pose provenance, missing-
   geometry fallback, and rotation/translation invariance tests.
6. Peak memory must remain below a predeclared safe ceiling (recommended 7.2 GB on an 8 GB device),
   or the model width/chunking must be reduced before formal runs.
