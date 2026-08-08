# Universal Mechanistic Basis Research Analysis

Updated: 2026-08-07

## Route Change

T-DIR-P0 showed that sparse PLIP event classification is not a stable primary
interface for the current model. PLIP remains useful as an auxiliary semantic
audit, but the main structural object is now a fixed, continuous and
permutation-invariant chemogeometric basis learned from structure as privileged
information.

The deployment contract remains:

```text
protein sequence + ligand 2D graph + support + observable context
```

Holo coordinates may define training targets but cannot become a required
deployment input.

## Evidence-Calibrated Architecture

### 1. UMBD: Universal Mechanistic Basis Distillation

```text
privileged holo coordinates -> fixed Phi*
sequence + 2D graph         -> predicted Phi
```

`Phi` must have a fixed analytic gauge. It cannot be an arbitrary rotated
latent embedding. The basis is built in increasing order:

1. two-body chemistry x radial moments;
2. low-order directional/angular terms;
3. only if required, selected many-body correlations.

T-BASIS-R0 now establishes the first item for the existing P1B frontend. It
does not establish items two or three.

### 2. CDAC: Cross-Dataset Delta-Affinity Calibration

Only after the structural basis is frozen and passes partner-recoverability:

```text
r = y - f_L^OOF(L)
r_hat = theta_gamma^T Phi(P,L)
```

The shared object is `Phi`; endpoint/dataset observation maps remain separate
through observable context `gamma`. Raw Ki, Kd, IC50 and KIBA scales must not be
pooled as one physical label. The E0R2 point-residual plus within-task
residual-difference objective remains the registered numerical template.

### 3. RFSA: Rank-Aware Few-Shot Section Adapter

Only after source affinity and transfer Gates pass. Support may adapt directions
that are identified by the support design. The necessary checks are:

```text
rank(Psi_S), singular values, conditioning, query row-space coverage
```

`d <= k` alone is insufficient. Uncovered query directions must shrink to the
population state or abstain.

## Theory Integration

The frozen operator remains:

```text
z -> F(z) in simplex -> B(z)F(z) -> K(beta)
```

Deep bioinformatics integration requires biology to determine the admitted,
bounded and observable coordinates of `z`. The future candidate state is:

```text
z_bio = [interaction basis groups,
         population delta-affinity contributions,
         identified support-section coordinates,
         rank/coverage diagnostics,
         observable assay context]
```

The dimension of `z` is chosen only after these objects are identified. The old
engineering choice of 28 dimensions is not a biological requirement.

## Gate Order

```text
E0R2 synthetic objective/design/solver                  PASS
T-DIR-P0 sparse PLIP-event learnability                NEGATIVE PILOT
T-BASIS-R0 fixed two-body radial recoverability        PASS
T-BASIS-A angular/many-body privileged distillation    NOT REGISTERED
E-AFF source OOF delta-affinity                        FROZEN
cross-dataset replication                              FROZEN
RFSA few-shot identified section                       FROZEN
biological z admission                                 FROZEN
P2-P4                                                  FROZEN
```

No later stage may compensate for an earlier failure by changing CSMO, Band,
theory, ligand-baseline leakage rules or closure governance.

## Current T-BASIS-R0 Evidence

- fresh 320-complex panel;
- 320 distinct homology groups and exact sequences;
- validation/test held out from P1B training;
- fixed 288D analytic teacher;
- test reconstruction gain `0.5312 [0.4433, 0.5962]`;
- test wrong-partner degradation `0.1561 [0.1070, 0.2007]`;
- all affinity/recipient reads zero.

This supports a structure-as-privileged-information route. It does not yet
support the word `universal`: one corpus, one radial basis and one frozen
frontend have been tested.

## Literature Boundary

- CompBind supports complex-guided pretraining with structure-free affinity
  inference: https://pubmed.ncbi.nlm.nih.gov/41562952/
- 3D Denoisers Are Good 2D Teachers supports distilling 3D representation
  knowledge into a 2D graph encoder:
  https://ojs.aaai.org/index.php/AAAI/article/view/31986
- ACE supports fixed local, systematically refinable invariant bases:
  https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.6.013804
- LUPI provides the general privileged-information precedent:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6238365/

These works motivate the architecture. They do not establish MetaSieve's basis,
affinity or few-shot Gates.
