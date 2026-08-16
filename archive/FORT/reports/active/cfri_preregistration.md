# CFRI preregistration — Round 1 (registered 2026-07-25, before any CFRI development number was read)

## Hypothesis

A target-conditioned residual `g(t,d)` can predict how an unseen target reorders ligands relative to
the global ligand-only ranking `b(d)`, without learning target identity, scaffold identity or source
shortcuts.

    y(t,d) = b(d) + g(t,d) + eps

Null hypothesis: within-target ligand ordering on strictly dual-cold targets is explained by ligand
marginals alone, and no target-conditioned term adds resolvable ranking information.

## Substrate and contract

* Registry `dataset/public/plinder_2024_06_v2/processed/dualcold/registry.parquet`,
  sha256 `1adde66d3fd8bdb5f1520023b979e3892286ce9015fbf373b9a03eaafcc5e359`, 11,367 rows,
  one row per (cluster, ligand parent connectivity).
* Splits assigned on connected components of the (cluster, UniProt accession) graph, so train /
  development / confirmation are disjoint on cluster AND accession by construction. Drug axis:
  Bemis-Murcko scaffold disjoint, parent connectivity disjoint, Tanimoto < 0.95 to any train ligand.
  All six audited overlap axes are exactly zero.
* Ligand tensor: 1034 = Morgan(r=2, 1024 bit) ++ 10 physicochemical descriptors, descriptors
  standardized on train rows only.
* Target tensor: frozen ESM-2 650M (MIT licence) pooled (1280) + 8 ordered segment tokens
  (8 x 1280), sequences resolved from the local ChEMBL-37 UniProt mirror. No fine-tuning.
  No accession, cluster id, family label, source or split value is ever an input.
* Development: 40 targets with >= 4 ligands and a resolved sequence, in 33 target components.
* Confirmation (PLINDER `test` components) is not touched. `sealed_test_consumed=false`.

## Cross-fitting

5 folds over training target components. `b^(-f)` is fitted without fold `f`, and every training
residual `r = y - b^(-component)(d)` comes from a base that never saw that row's own target
component. Provenance is asserted programmatically (`provenance_violations` must be 0). The base
deployed at evaluation is a separate full-fit base, frozen and shared by every additive/interaction
arm, so `CFRI - B0` differs by exactly the interaction term.

## Objective (weights frozen here, from train-only loss magnitudes; never tuned on development)

    L = L_residual_rank
      + 0.3 * L_affinity        (mean squared error of b + g against y)
      + 1.0 * L_target_centering    (mean_d g(t,d))^2
      + 1.0 * L_base_orthogonality  corr(g, b)^2 inside the episode

## Arms

`B0`, `T0`, `A0`, `I0`, `R0`, `CFRI`, `CFRI-Tshuffle`, `CFRI-Lshuffle`, `CFRI-Tpool`,
`CFRI-Trandom` — identical data, identical optimizer budget, matched interaction capacity.

## Power (frozen before the run, `reports/active/dualcold_power.json`)

Retraining the SAME ligand-only base under different seeds moves per-component development Spearman
by sd 0.23–0.26. On 33–34 components this gives a paired **MDE80 of 0.0954** under the zero-shot
full-query protocol (0.131 under the reserved-pool protocol). The task's nominal 0.03 threshold is
therefore **unresolvable** on this substrate; the task permits stricter thresholds when the power
audit justifies them, so Gate Z is frozen at the empirical MDE80.

## Gate Z (frozen; failure means stop, not retune)

| id | criterion |
|----|-----------|
| Z1 | `CFRI - B0` development target-macro Spearman >= **0.0954** |
| Z2 | grouped target-component bootstrap LCB95 > 0 |
| Z3 | positive direction in every development target component |
| Z4 | RMSE no more than 2% worse than B0 |
| Z5 | significant damage under BOTH target derangement and ligand derangement |
| Z6 | LCB95 > 0 in at least one registered stratum |

If Gate Z fails, Round 2 (three seeds) and Round 3 (BERP-GS) are **not** authorized, and no
threshold, evaluation set or capacity change may be used as a rescue.

## Prior expectation

Recorded honestly in advance: the substrate supplies ~470 training targets with >= 2 ligands and
~20k within-target ligand pairs, of which only ~650 pairs are measured on more than one target. This
is thin supervision for a target-conditioned reordering function, and four previous programs on this
repository (Davis docking, holo receptors, PLINDER native pose, Gate-P2A) found protein-side
information not load-bearing for dual-cold within-target ranking. A negative Gate Z is the expected
outcome; the purpose of the run is to make that verdict falsifiable and quantitative.
