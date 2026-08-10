# X1A crossed-interaction ICC precondition — evidence consolidation

## Terminal result

```text
X0/X0-B census recovery and hashing ...... PASS (3/3 data files byte-verified)
label firewall ........................... PASS (Ki/Kd read only after prereg commit)
G1 Ki  UCB95(rho) < 0.0915 ............... PASS  (3.88e-07)
G2 Kd  UCB95(rho) < 0.0164 ............... PASS  (1.33e-03)
G3 no cluster dominates (capped) ......... PASS  (Ki 0.0387, Kd 0.2066)
G4 effective n >= 245 .................... PASS  (Ki 827.0, Kd 604.3)
terminal verdict ......................... X1_ICC_PRECONDITION_PASSED
X1B ...................................... HISTORICAL AUTHORIZATION WITHDRAWN
X2 ....................................... NOT AUTHORIZED, NOT TRAINED
```

X1A trains nothing and introduced zero parameters.

> **Authorization correction (2026-08-10):** the values above remain the
> historical amended-audit result. They no longer authorize X1B because the
> residual ICC is structurally attenuated and does not estimate dependence of
> X1B's statistic `q = DD^2 - v_noise`. The current verdict is
> `X1A_ICC_PRECONDITION_NOT_ESTABLISHED`; see
> `X1A_REPAIR_AND_CURRENT_BOUNDARY.md`.

## 1. Recovery and firewall

X0/X0-B artifacts were recovered from commit `24a9ae0^` and byte-verified:

```text
cells.jsonl                  898df882...  VERIFIED   58,668 cells
dependency_components.jsonl  8970d059...  VERIFIED   70 components
panels.jsonl                 f378cdd6...  VERIFIED   631 panels
report.json                  1f89bcdc... on disk vs cbaa4c83... recorded  MISMATCH
```

The `report.json` mismatch is disclosed rather than repaired: the manifest was
written before the report was finalised. Only the three verified data files are
used, and the X0-B design and statistical unit are **not** rebuilt or replaced.

ChEMBL37 values were read only after `PREREG_X1_ICC_AND_DATA_CONTRACT.md` was
committed (`008c82a`). 63,859 activity ids already enumerated by the
label-blind X0 census were joined against the pinned archive; only
`activity_id`, `standard_relation`, `standard_type` and `pchembl_value` were
selected. All 63,859 returned rows carried `standard_relation '='` and a
non-null `pchembl_value`, so no row was dropped for censoring. BindingDB,
DAVIS, KIBA, PDBbind and recipient reads are `0`. No out-of-fold residual was
computed.

## 2. Two estimator defects found and corrected

**The parent ICC estimator was degenerate.** Section 5 of the parent
registration fitted additive target and ligand effects *within panel*. A
per-panel intercept forces every panel's residual mean to exactly zero, so
`var(cluster)` is identically zero for any dataset. Demonstrated on synthetic
panels with injected cluster offsets of 10, 20 and 30 log units, all returning
residual means of order `1e-15`.

The first execution duly returned `rho = 0.0000`, `ucb95 = 0.0000` and would
have reported `X1_ICC_PRECONDITION_PASSED`. **That result is void.** A Gate that
cannot fail is not a Gate. Amendment 01 replaced the within-panel fit with a
global per-endpoint fit and changed no threshold, seed, endpoint, cluster
definition, firewall rule or verdict name. The correction ran against this
stage's convenience: the void result was already a pass.

**G3 and G4 used the wrong unit.** They were computed on measurement counts.
The registered statistical unit is the X0-B cell-disjoint DD unit. Both now use
the frozen X0-B per-cluster unit sizes, which reproduce X0-B's capped totals
exactly — Ki `sum(min(size, 32)) = 827`, matching X0-B's own `units: 827` at
its breakeven cap.

## 3. Results

| | Ki | Kd |
|---|---:|---:|
| cells with usable values | 53,673 | 4,995 |
| measurements | 57,489 | 6,370 |
| replicate-supported cells | 2,968 | 1,363 |
| dependency clusters | 36 | 12 |
| `rho` point estimate | 2.22e-07 | 2.35e-05 |
| `rho` UCB95, cluster bootstrap | 3.88e-07 | 1.33e-03 |
| `rho*` threshold | 0.0915 | 0.0164 |
| largest cluster share, uncapped | 0.4818 | 0.4006 |
| largest cluster share, capped | 0.0387 | 0.2066 |
| design effect | 1.0000 | 1.0012 |
| effective `n` | 827.0 | 604.3 |

Variance components of the additively-adjusted residual:

| component | Ki | Kd |
|---|---:|---:|
| cluster | 8.49e-08 | 8.12e-05 |
| panel | 0.0 (truncated) | 0.0 (truncated) |
| cell | 0.0 (truncated) | 2.19e-03 |
| replicate | 0.38190 | 3.45871 |

Both endpoints clear every registered Gate, so the registered terminal verdict
is `X1_ICC_PRECONDITION_PASSED` for Ki and Kd.

## 4. Three caveats that must travel with this pass

**The pass is biased in the unsafe direction.** The global additive fit is
heavily parameterised: it consumes 42% (Ki) and 32% (Kd) of cells as
parameters, and 12.2% / 14.5% of ligands appear in exactly one cell, so their
residual is forced to exactly zero. Every structural variance component is
therefore shrunk toward zero, and `rho` is a *lower* bound far more credibly
than a point estimate. A gate of the form "`rho` must be small" is easier to
pass under this bias, not harder. This is stated because the direction matters:
the audit's own error mode favours the outcome it produced.

**The nested decomposition is partly non-identified.** `var(panel)` truncated to
zero for both endpoints and `var(cell)` truncated to zero for Ki. After the
additive fit, between-cell and between-panel structure is smaller than
replicate noise, so the moment estimators return negative values. The
truncation is registered and reported, but it means the decomposition resolves
only two of four levels for Ki.

**Replicate noise dominates and sets a hard bar for X1B.** Replicate variance is
99.99998% of the Ki total and 99.93% of the Kd total:

```text
Ki  replicate SD = 0.618 log units  ->  detectable interaction RMS = 0.309
Kd  replicate SD = 1.860 log units  ->  detectable interaction RMS = 0.930
```

at X0's frozen interaction-RMS-to-noise ratio of `0.5`. A crossed interaction
RMS of `0.93` log units for Kd would be a near ten-fold selectivity swing as an
*RMS*, not as an extreme. Kd clears the dependence precondition on paper while
being, on this evidence, close to unusable for detecting anything smaller. Ki's
`0.309` is demanding but not implausible.

This is precisely what X1B is built to adjudicate: its estimand
`I_real^2 = max(0, E[DD^2] - E[v_noise])` subtracts exactly this noise. X1A
tested dependence, not detectability, and only dependence passed.

## 5. What this does and does not authorize

`X1_ICC_PRECONDITION_PASSED` authorizes **X1B only**, for Ki and Kd separately.
It does not authorize X2, any trainable component, any 3D route, support
adaptation, or any claim about affinity energy, selectivity, causal
interaction, few-shot capability or `z` admission.

X1B was not run in this session. It remains the authorized next stage and must
test interaction *variance* against replicate and assay noise — never whether
`mean(DD)` differs from zero, since opposing selectivity effects cancel.

## 6. Governance

- Trains nothing; zero trainable parameters.
- Ki and Kd analysed completely separately, never merged, pooled or averaged.
- X0-B design and statistical unit not rebuilt and not replaced.
- No threshold, margin, seed or Gate was weakened. The one amendment corrected
  a provably degenerate estimator and moved no threshold.
- `model/`, production `scripts/`, `theory/`, CSMO, Band, the mesh, production
  `z` and `A(F,z)=K(B(z)F(z))` are unmodified. All code is under `research/`.
- Historical artifacts preserved; no failed verdict rewritten.
