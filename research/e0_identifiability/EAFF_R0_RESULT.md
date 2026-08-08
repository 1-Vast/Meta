# E-AFF-R0 Readout Diagnosis Result

## Verdict

```text
READOUT_BLIND_TO_TASK_LEVEL_AFFINITY_LOCATION
PERFECT_LEVEL_PREDICTOR_SCORES_CHANCE_AT_EVERY_VARIANCE_SHARE
```

The readout used by every affinity result in this project is exactly invariant
to the affinity channel a protein sets for a chemical series. No affinity label
was read.

## 1. Exact Invariance

Using the repository's own `metrics.concordance` on the real H0C task sizes
(54 tasks, 20 test ligands each), maximum absolute change in within-task
concordance:

| Transform | Max deviation |
|---|---:|
| per-task prediction shift | `0.0` |
| per-task prediction rescale | `0.0` |
| per-task label shift | `0.0` |
| per-task label rescale | `0.0` |

Not approximately zero — exactly zero. The metric compares only the signs of
pairwise differences inside a task, and a task-level affinity location and scale
are exactly such maps.

## 2. Credit Assigned To The Level Channel

Generating `y[t,i] = level[t] + within[t,i]` and scoring oracles under the
registered readout and under a location-sensitive error:

| Level share of variance | CI, level oracle | CI, within oracle | RMSE, level oracle | RMSE, global mean |
|---:|---:|---:|---:|---:|
| 0.059 | **0.5000** | 1.0000 | 1.001 | 1.033 |
| 0.200 | **0.5000** | 1.0000 | 1.006 | 1.166 |
| 0.500 | **0.5000** | 1.0000 | 0.970 | 1.514 |
| 0.800 | **0.5000** | 1.0000 | 0.961 | 2.079 |
| 0.941 | **0.5000** | 1.0000 | 1.019 | 3.412 |
| 0.985 | **0.5000** | 1.0000 | 1.017 | 8.832 |

A predictor that knows a task's affinity level perfectly and nothing else scores
exactly chance under the registered readout, at every variance share. The same
predictor's advantage under a location-sensitive error grows from `1.033/1.001`
to `8.832/1.017` across the sweep. The readout is not merely insensitive to the
level channel; it is algebraically blind to it.

## 3. The Frozen Geometry Was Not Inert

From H0C's published per-task scores:

| Contrast | Mean | Median | Tasks unchanged |
|---|---:|---:|---:|
| correct − local ligand | −0.00391 | −0.00531 | 3 / 54 |
| deranged − local ligand | −0.00543 | −0.00263 | 1 / 54 |
| correct − deranged | +0.00152 | +0.00000 | 7 / 54 |

The geometry term moved 51 of 54 tasks. The null result is therefore not an
inert feature; it is a feature whose *sign* carried no consistent affinity
direction under this readout — which is precisely the failure the project
records as "structural information is not affinity information".

## 4. Compounding Upstream Removal

Independently of the metric, the H0C design removed the level channel twice
before the geometry was consulted. `support_target = y − global_ligand_prior`,
then `interaction_target = support_target − task_local_ligand_nuisance`, where
that nuisance is fit on 20 labelled support examples of the correct protein's
own task. The resulting `local_score` is added to **both** the correct and the
deranged arm, so both arms hold the correct protein's task level before the
contrast is taken. This is a design fact readable from
[run_eaff_h0c.py:241](run_eaff_h0c.py:241) and
[run_eaff_h0c.py:264](run_eaff_h0c.py:264), not an inference from the numbers,
and it explains why replacing the protein cost so little.

## 5. What This Does And Does Not Establish

Establishes: the registered readout, and the H0C residualization feeding it,
cannot detect a protein contribution expressed as a task-level affinity
location, whatever its size.

Does **not** establish: that protein-specific affinity lives in that channel.
That is a hypothesis, and testing it is the job of a separately registered Gate
with held-out proteins, a location-sensitive metric, a sequence-only baseline
and no shared support nuisance. R0 authorizes nothing.

Past results are not overturned. P1C, P1R\*, E-AFF-P0, H0A and H0C remain valid
statements about within-task ranking. What changes is their scope: they are
evidence about ranking, not about affinity level, and they were never evidence
about the functional the frozen theory controls.

## Why This Was A Reasonable Design To Begin With

Within-task concordance was a defensible choice: absolute affinity levels across
tasks are confounded by assay protocol, and a cross-task metric invites exactly
the assay-identifier shortcut the project's requirements prohibit. The metric
defeated that shortcut completely. The cost was that it also discarded the
channel where a protein's affinity contribution primarily sits. The repair is
not to abandon confound control but to control the confound explicitly — assay
stratification through the context map, held-out proteins, and a sequence-only
baseline — while measuring a location-sensitive functional.
