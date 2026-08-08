# E-AFF-H0C Fixed Radial Interaction Residual

Status: registered before H0C feature generation, source-affinity access,
model fitting, or metric computation.

## Question

After removing support-matched ligand SAR and both fixed-tensor marginals, does
the remaining radial chemistry-distance interaction residual add material,
correct-partner-specific affinity ranking information?

This is a research-only confirmation diagnostic. H0A informed its hypothesis,
so H0C excludes every H0A task and ligand row but cannot claim an untouched
closure-family validation panel.

## Selection

- Start from governed E0 label-blind rows only.
- Exclude all 107 tasks used by H0A.
- Require at least 40 distinct ligand states and a deterministic strict
  Bemis-Murcko scaffold partition with 20 support and 20 test ligands.
- Treat all acyclic ligands as one explicit empty-Murcko group.
- Select one eligible task per closure component and all task/scaffold/ligand
  choices by fixed SHA-256 ordering. Affinity values cannot affect selection.
- Support and test scaffold sets must be disjoint. A task that cannot satisfy
  the exact 20/20 contract is ineligible; no random ligand split is allowed.

## Frozen Inputs

P1B, local-state cache, atom/residue chemistry, six radial coordinates,
T-BASIS calibration, D1 closure folds, OOF ligand-prior procedure, Ridge
`alpha=10`, and score-blind `<0.40` derangement are unchanged.

The calibrated tensor is not assumed nonnegative. Define its algebraic
double-centered residual as

```text
phi_null = chemistry_marginal * radial_marginal / total
psi = (phi - phi_null) / total
```

All totals must be positive. `psi` must have chemistry and radial marginals
zero to numerical tolerance. It is not called a probability distribution or a
physical energy.

## Support-Matched Nuisance

For each task, first target the residual from the frozen global closure-OOF
ligand prior. Fit a task-local ligand Ridge on the frozen 128D pooled ligand
state using the same 20 support labels. Its support predictions are five-fold
cross-fitted by fixed ligand-hash folds. Fit the interaction direction only to
pairwise differences of these cross-fitted nuisance residuals. Refit the same
ligand nuisance on all 20 support rows only to score the 20 test rows.

No test label, deranged feature, task ID embedding, alternative fingerprint,
hyperparameter search, or orientation feature enters either fit.

## Arms

- `Global-L`: frozen closure-OOF ligand prior.
- `Local-L`: Global-L plus the support-matched ligand nuisance.
- `Interaction-C`: Local-L plus the direction fitted on `psi(P,L)`.
- `Interaction-D`: Local-L plus the same direction applied to `psi(P',L)`.

The wrong protein is evaluation-only. Report task CI, component-macro CI,
endpoint secondary results, and a 2,000-draw closure-component bootstrap.

## Frozen Criteria

Interaction value requires `CI(Interaction-C) - CI(Local-L) >= 0.03` with a
positive 95% component-bootstrap LCB. Partner specificity requires
`CI(Interaction-C) - CI(Interaction-D) >= 0.03` with a positive LCB.

Terminal verdicts are:

- `FIXED_RADIAL_INTERACTION_RESIDUAL_AND_PARTNER_SIGNAL_OBSERVED`
- `FIXED_RADIAL_INTERACTION_RESIDUAL_WITHOUT_PARTNER_SPECIFICITY`
- `WEAK_PARTNER_CONDITIONED_INCREMENT_BELOW_GATE`
- `FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED`
- `DATA_OR_NUMERICAL_CONTRACT_FAIL_CLOSED`

Even the strongest verdict authorizes only a separately sealed validation of
the same hypothesis. It does not authorize RFSA, DAVIS, production `z`,
orientation, CSMO/Band changes, or P2-P4.
