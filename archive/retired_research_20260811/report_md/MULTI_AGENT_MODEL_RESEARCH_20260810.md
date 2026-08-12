# Multi-agent Meta-Section research and validation

Date: 2026-08-10

## Verdict

The minimal support-identifiable Meta-Section implementation passes its
synthetic positive control. Real BindingDB M0/M1 training must not start yet:
the O1 PASS is structural, while the executable numeric Ki corpus falls below
the frozen evaluation-depth Gate and lacks matching full-corpus T-BASIS rows.

```text
SYNTHETIC_META_SECTION_POSITIVE_CONTROL_PASS
NUMERIC_FEWSHOT_CORPUS_GATE_REQUIRED
REAL_META_SECTION_TRAINING_NOT_AUTHORIZED
```

Three independent analyses covered primary literature, repository/data
interfaces, and experimental/statistical design. The implementation and final
decision were integrated against the current execution authority.

## Architecture amendment review

The later Theory-Projected Q-PMA prompt did not pass architecture review. Its
penalized query-specific support weights can retain a component in
`null(M_S^T)` and use residual variation outside the declared function family.
With zero attention weight it collapses to the existing ridge leverage; with a
nonzero weight it can violate the off-coverage zero-correction contract and no
longer gives one shared task coefficient across queries.

A repaired row-space attention coordinate is mathematically safer but changes
the family from `m(P,L)` to `m_section(q; X_S)`. It is therefore a future,
separately preregistered hypothesis rather than a refactor of the frozen core.
The complete derivation is in
`research/meta_fewshot/ARCHITECTURE_GATE_20260810.md`. No `model/` or `scripts/`
file was changed.

## Model decision

The first falsification remains the frozen minimal family:

```text
phi(P,L) -> U^T phi, d<=5 -> positive-ridge support section -> query correction
```

Q-PMA or another learned attention module is not introduced now. The closed-form
section is already query-dependent through the leverage of `m_q` against the
support row space. Attention would add nonlinear learned support weighting and
would confound the untested coefficient-family hypothesis with a second model
hypothesis. It remains a replacement candidate only after the linear section's
support and biology effects are identified and query misspecification is the
localized blocker.

## Numeric corpus blocker

The frozen O1 artifact has 21,473 metadata-defined cells. A read-only,
hash-matched coverage audit reconstructed the exact 66 dependency roots and
found:

| availability | source k=5 targets | evaluation k=5 targets |
|---|---:|---:|
| structural O1 cells | 336 | 33 |
| exact numeric Ki available | 313 | 29 |
| exact numeric Ki and old 288D bank | 224 | 12 |

Only 18,509 of 21,473 structural pairs have an exact numeric Ki row; 2,964 do
not. The existing 288D bank belongs to the older 12,457-cell quotient corpus
and overlaps only 11,817 O1 pairs. A partial join would recreate the
estimand-mismatch failure.

The value `29` is diagnostic, not a new formal Gate result: numeric admission
and the aggregation of 1,389 replicated pairs have not been preregistered.
Nevertheless, value availability alone is enough to show that the current
claim of 33 executable k=5 tasks is unsupported. The next corpus stage must
freeze aggregation, persist cell/target-to-component membership, seal label
views, and re-run the unchanged `>=30` and MDE thresholds before feature
generation or training.

The diagnostic is now reproducible through
`research/meta_fewshot/numeric_availability_audit.py`. It reconstructs and
verifies all 66 frozen dependency roots, uses only identity fields from the
exact-label artifact, performs no aggregation and writes
`report/meta_fewshot/FS_NUMERIC_AVAILABILITY_AUDIT.json`.

## Paper-locked data and experiment status

AdaMBind and MetaSieve experiments are frozen as separate lanes in
`research/meta_fewshot/PAPER_PROTOCOL_MATRIX_20260810.md`. The official public
AdaMBind CSVs were downloaded from repository commit `01a169a6...` and hashed:
BindingDB 42,203 rows, Davis 30,056 and KIBA 118,254. The manifest is
`dataset/raw/adambind_public_01a169a6/acquisition_manifest.json`.

The public snapshot lacks the paper's CD-HIT 40% split manifests, while its
defaults conflict with the paper on headline shots and outer iterations. The
versioned Zenodo archive `10.5281/zenodo.18595084` is restricted. Consequently
the data are admitted for audit only: neither exact AdaMBind reproduction nor
MetaSieve real-data training was started.

## Synthetic validation

The synthetic experiment used 96 source tasks and 48 held-out tasks with a
planted three-dimensional family in 12 observed dimensions. It meta-learned
`U` through query loss while support adaptation used only the positive-ridge
closed form.

| check | result |
|---|---:|
| planted subspace overlap | 0.9789 |
| correct-support MSE | 0.1250 |
| zero-section MSE | 2.6603 |
| foreign-support MSE | 7.0179 |
| permuted-support MSE | 6.4501 |
| maximum section rank at k=5 | 3 |

Unit validation also covers positive ridge, dimension limits, support-order
invariance, rank at most k, off-row zero correction/coverage, outer gradient to
`U`, analytic covariance versus Monte Carlo, and deterministic bounded-noise
radius. Query labels are passed only to the outer loss, never to the section.

## Literature findings

- MetaDTA is the closest attention-based target-as-task precedent: continuous
  pIC50 tasks, ANP query/support conditioning, and repeated random or scaffold
  episodes. It does not use protein features, uses 10--100 support examples,
  and does not validate Ki k<=5 biology. Source:
  https://openreview.net/forum?id=yzlif16IASM
- FS-Mol establishes separate task pools, repeated support sampling,
  task-macro evaluation, and the importance of source-task diversity. It is a
  ChEMBL27 assay-level binary benchmark with 4,938/40/157 train/valid/test
  tasks, not quantitative DTA. Source:
  https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/file/8d3bba7425e7c98c50f52ca1b52d3735-Paper-round2.pdf
- The final FS-CAP paper uses BindingDB with 1,754 training and 41 test targets;
  context sizes are 1/2/4/8. Numbers for 4,807 BindingDB assays and 98,593
  PubChem assays belong to the earlier preprint and must not be cited as final
  paper results. FS-CAP is ligand/context-only and has no protein biology or
  law-valued output. Source: https://pubs.acs.org/doi/10.1021/acs.jcim.4c00485
- AdaMBind validates protein-as-task, 5/40-shot evaluation, and a CD-HIT 40%
  novel-target split. Its predictor uses one-hot protein input plus a 1D CNN;
  ESM-2 is used for later task-similarity analysis. Its MAML/adaptive sampler
  is not evidence that MetaSieve should add unconstrained task parameters.
  Source: https://www.nature.com/articles/s41467-026-70554-5
- R2-D2 supports differentiating through a dual closed-form ridge solve, which
  is computationally natural for k<=5. Its experiments are image
  classification, so its DTA relevance is a transfer hypothesis rather than
  direct evidence. Source: https://arxiv.org/abs/1805.08136

## Frozen next decision

1. Reopen O1 only as a numeric-corpus amendment; do not alter the structural
   split or thresholds to rescue the count.
2. Freeze replicate aggregation and exact-value admission, persist a self-
   contained split mapping, and re-run numeric FS-C0/FS-C1.
3. If it passes, generate complete correct/wrong-protein/ligand-only 288D rows
   for that exact corpus and seal source/validation/evaluation label views.
4. Run source-only M0 with `d=0..5`, target-balanced episodes and nested
   `k=1/2/3/5`, then correct/zero/foreign/permuted support controls.
5. Open the frozen evaluation once only if M0 passes; use dependency-component
   inference. Do not tune Q-PMA or any other replacement on consumed targets.
