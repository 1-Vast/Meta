# MEDIP S0 synthetic engineering preregistration

**Frozen:** 2026-07-28, before `model/medip.py`,
`research/medip_synthetic.py`, `tests/test_medip.py`, or any MEDIP result
artifact existed.

**Environment:** `D:\anaconda\envs\drug\python.exe`.
**Formal seeds:** `1729, 1730, 1731, 1732, 1733, 1734`.

**Ordering amendment A1 (2026-07-28, before any smoke or formal MEDIP
run):** S0 is allowed before `OMUT-I0` by the same independent
estimator-falsification logic as R-MAON G0 because it uses generated features
and generated outcomes only. This corrects the more conservative ordering in
`open_evidence_pretraining_crossreview` step 6. An S0 pass cannot unlock real
outcomes, grant representation or biological credit, authorize open-data
pretraining, or advance any real-label stage. Those actions remain blocked by
the applicable OMUT information, provenance, topology, firewall, and power
gates.

## 1. Question and claim boundary

S0 asks whether a small evidence-decoupled observation module can recover a
known low-dimensional target-ligand interaction when synthetic measurements
mix exact endpoints, censoring, binary observations, ordinal observations, and
same-ligand cross-target selectivity comparisons.

This is an engineering calibration only. It uses generated Gaussian feature
coordinates and generated outcomes. It reads no registry, dataset, affinity
column, affinity value, development asset, confirmation asset, or sealed-test
asset. A pass is not biological evidence, open-data evidence, predictive
evidence, or authorization to train on real affinity data.

The tested design has two explicit boundaries:

1. The interaction encoder receives only target and ligand coordinates.
   Endpoint, assay, source, document, and provenance metadata may enter only
   observation heads.
2. Observation heads explain measurement processes. They do not alter the
   latent interaction score and cannot be used as inputs to the interaction
   encoder.

## 2. Frozen synthetic design

Each of six replicates contains 36 targets and 48 ligands. Targets have four
observed coordinates and ligands have five observed coordinates. The first 24
targets and first 32 ligands form the training rectangle. The remaining 12
targets and 16 ligands form a strict dual-cold test rectangle; neither a test
target nor a test ligand occurs in training.

The latent score is

```text
s(t,l) = t^T a + l^T b + t^T Theta l,
```

where `Theta` is a seeded rank-two matrix. Its scale is fixed from the training
rectangle before any observation is generated. Target and ligand main effects
are deliberately smaller than the interaction. No protein sequence,
structure, molecule, or biological label is simulated.

Two exact endpoint processes have distinct intercepts and endpoint-by-source
offsets. There are three synthetic sources and endpoint-specific Gaussian
noise scales. Endpoint and source assignment are correlated with observed
target-ligand coordinates so that ignoring measurement metadata can induce a
false interaction rather than merely inflate independent noise.

For the continuous measurement stream, values within endpoint-specific
detection limits are retained exactly. Values below or above a limit are
stored only as left- or right-censored observations. Binary and four-level
ordinal measurements are independently sampled from monotone functions of the
same latent score with source-specific observation offsets. Same-ligand
cross-target pairs are sampled from the training rectangle, and their signed
order is generated from the latent score difference.

## 3. Frozen model

The interaction encoder is the direct low-dimensional score above. The full
model contains:

- a direct target-by-ligand interaction matrix plus target and ligand main
  effects;
- separate exact endpoint intercepts, endpoint-by-source offsets, and
  endpoint-specific noise scales;
- a stable Tobit-style censored likelihood using `torch.special.log_ndtr`;
- a binary logistic head with source-specific observation offsets;
- an ordered-logistic head with strictly ordered learned cut points and
  source-specific observation offsets;
- a same-ligand cross-target logistic selectivity loss.

The exact observation slope is fixed at one. This anchors the latent scale and
prevents a learned endpoint slope from being traded arbitrarily against the
interaction matrix.

All variants use the same optimizer, initialization rule, observations,
training steps, and loss weights. The formal implementation will use
deterministic CPU PyTorch, full-batch Adam, 500 optimization steps, learning
rate `0.03`, and gradient clipping at norm `10`. The loss is the sum of mean
continuous (exact plus censored), binary, ordinal, and selectivity negative log
likelihoods with weights `1.0, 0.35, 0.35, 0.75`, plus `1e-4` times the squared
parameter norm of the interaction encoder.

## 4. Frozen comparison

Five variants are fitted independently in every replicate:

| variant | change from the correct model |
| --- | --- |
| `correct_multifidelity` | no corruption |
| `forced_endpoint_merge` | both exact endpoint IDs are replaced by one shared endpoint ID |
| `endpoint_source_shuffle` | endpoint and source IDs are independently permuted within each observation stream |
| `selectivity_shuffle` | selectivity signs are permuted across fixed same-ligand target pairs |
| `separable_null` | the target-by-ligand interaction matrix is absent; target and ligand main effects remain |

No variant may receive more parameters in the interaction encoder. The
endpoint/source corruption changes only observation metadata. The selectivity
corruption changes only comparison labels.

## 5. Metrics

Let `S` and `S_hat` be the true and predicted score matrices on the strict
dual-cold test rectangle. The reordering estimand is the doubly centered
interaction:

```text
D(S) = S - row_mean(S) - column_mean(S) + grand_mean(S).
```

The primary interaction metrics are Pearson correlation between `D(S)` and
`D(S_hat)` and same-ligand cross-target ordering accuracy. Ties in a separable
prediction count as one half. Observation calibration is evaluated on the
noise-free mean for every test pair, endpoint, and source, using exact-mean
RMSE. Endpoint residual bias is the maximum absolute mean residual over the
two endpoints.

Per-replicate contrasts are computed before taking medians. No replicate may
be removed.

## 6. Frozen gates

The module passes S0 only if all gates pass:

| gate | requirement |
| --- | --- |
| `S0_RECOVERY` | correct median mixed-difference correlation `>= 0.85` and median cross-target ordering accuracy `>= 0.80` |
| `S0_SEPARATE_ENDPOINTS` | forced endpoint merge minus correct median exact-mean RMSE `>= 0.20`, and correct median maximum endpoint residual bias `<= 0.15` |
| `S0_METADATA_DESTRUCTION` | correct minus endpoint/source-shuffle median mixed-difference correlation `>= 0.10`, and shuffle minus correct median exact-mean RMSE `>= 0.15` |
| `S0_SELECTIVITY_DESTRUCTION` | correct minus selectivity-shuffle median ordering accuracy `>= 0.08` |
| `S0_INTERACTION_NULL` | correct minus separable-null median ordering accuracy `>= 0.20`; separable-null doubly centered prediction has maximum absolute magnitude `<= 1e-6` |
| `S0_ARCHITECTURE` | tests prove metadata is absent from the interaction encoder signature, metadata perturbation leaves latent scores bit-identical, and metadata parameters occur only below observation heads |
| `S0_NUMERICS` | all formal losses and metrics are finite; extreme left- and right-censor likelihood tests are finite with finite gradients |
| `S0_REPRODUCIBLE` | a repeated formal run is byte-equivalent after excluding runtime duration |

All gates pass:

```text
MEDIP_S0_ENGINEERING_CALIBRATION_PASS__SYNTHETIC_ONLY
```

Any gate fails:

```text
MEDIP_S0_ENGINEERING_CALIBRATION_STOP
```

## 7. Prohibited rescue and reopening

After the first formal result, no seed, sample size, assignment mechanism,
model dimension, training step count, loss weight, threshold, metric, or gate
may be changed. A failed variant may not be dropped. Increasing rank, adding
an encoder, extending optimization, or weakening a gate is prohibited.

S0 cannot authorize real open-data pretraining. That stage may reopen only
after source and endpoint lineage, censor semantics, inactive retention,
fold-local construction, target/homology/ligand/scaffold/chemical-neighbour
and provenance closure, sufficient independent components, and a
preregistered power target all pass outside this synthetic calibration.

## 8. Planned artifacts

```text
model/medip.py
research/medip_synthetic.py
tests/test_medip.py
reports/active/medip_s0.json
reports/active/medip_s0_decision.md
```
