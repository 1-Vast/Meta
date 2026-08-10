# Preregistration — X1

## Crossed protein-by-ligand affinity interaction: ICC precondition and data contract

Stage identifier: `E-AFF-X1A_ICC_AND_DATA_CONTRACT`

Written 2026-08-10, after the C1 terminal verdict
`EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER` (commit `e68847c`) and **before
any ChEMBL37 Ki or Kd value has been read**. All numbers quoted below come from
the already-published, label-blind X0 and X0-B artifacts.

This document registers **X1A only**. X1B requires X1A to pass; X2 requires X1B
to pass and its own separate preregistration.

## 1. Question

Do the available public source data identify a real target-dependent ligand
preference after removing additive protein and ligand effects?

```text
DD = y(P1, La) - y(P1, Lb) - y(P2, La) + y(P2, Lb)
```

`DD` is the primary biological estimand. Same-target ligand differences,
out-of-fold residuals, wrong-protein penalties and residue contact scores are
**not** substitutes for it and may not be reported as if they were.

X1A does not test `DD`. It tests whether the dependence structure of the source
permits `DD` to be tested at all at the frozen sensitivity.

## 2. Recovered X0-B census — frozen inputs

X0 and X0-B artifacts were recovered from commit `24a9ae0^` and verified
byte-for-byte against the X0 manifest:

```text
cells.jsonl                  sha256 898df88235401a2be2341ae1ab222e6c5903202796c8312d8e9091cf76741562  VERIFIED
dependency_components.jsonl  sha256 8970d059bcbf4dd1ca0f3b7ec9e5ab3f594740ec37bcb5c173fdb66ffb8a7779  VERIFIED
panels.jsonl                 sha256 f378cdd610205fc02850e84e5530830d2255fcf539f4ecca81518eb19bf036bc  VERIFIED
```

`report.json` does not match its manifest entry (`1f89bcdc…` on disk versus
`cbaa4c83…` recorded). The manifest was written before the report was
finalised. The three data files that carry the census all verify, and only they
are used. This discrepancy is disclosed, not repaired.

The X0-B design is **not rebuilt and its statistical unit is not replaced**:

| | Ki | Kd |
|---|---:|---:|
| cell-disjoint units | 11,168 | 1,041 |
| dependency clusters | 36 | 12 |
| distinct target pairs | 205 | 49 |
| distinct targets | 224 | 73 |
| distinct ligands | 19,062 | 1,256 |
| largest uncapped cluster share | 0.4818 | 0.4006 |
| breakeven `rho*` | 0.0915 | 0.0164 |
| cap at `rho*` | 32 | 125 |
| required effective `n` | 245 | 245 |

## 3. Frozen definitions

```text
endpoint        exactly two, Ki and Kd, analysed completely separately and
                never merged, pooled, averaged or substituted for one another
scale           pChEMBL value, -log10(molar). Rows without a pchembl_value are
                excluded rather than converted, so no unit inference occurs
relation        standard_relation '=' only; censored '>' '<' '>=' '<=' '~' are
                excluded and counted
construct       protein_sequence_sha256 as recorded by X0
organism        the X0 panel context assay_organism, unchanged
stereochemistry ligand_connectivity_key, the stereochemistry-stripped InChIKey
                connectivity block, exactly as X0 defined it
document        the D1 homology-document closure component recorded by X0
assay / panel   panel_id and assay_activity_ids as recorded by X0
target pair     an unordered pair of distinct protein_sequence_sha256 in one panel
ligand pair     an unordered pair of distinct ligand_connectivity_key in one panel
cell            (panel_id, protein_sequence_sha256, ligand_connectivity_key)
closure         closure_component_id, the D1 homology-document closure component
cluster         dependency_component_id, the X0 dependency component; this is
                the inference unit and is never replaced
replicate       two or more activity_ids for one cell under one assay_chembl_id
```

Nothing in this list may be redefined after a value is read.

## 4. Label firewall

- ChEMBL37 Ki/Kd values may be read **only after this document is committed**.
- Reads are restricted to `activity_id`, `standard_relation`, `pchembl_value`
  and `standard_type`, joined on the activity ids already enumerated by X0.
  No new selection of rows, assays, targets or documents occurs.
- The pinned archive is
  `dataset/raw/source_affinity/chembl37_sqlite_v1/extracted/chembl_37/chembl_37_sqlite/chembl_37.db`,
  sha256 `4be13df3b68e25dcd0bff44bf094033b5aebe98f415acdc8c1cdf380e0c15142`.
- BindingDB, DAVIS, KIBA, PDBbind, recipient labels and every previously
  consumed confirmation panel remain **unread**.
- No out-of-fold residual is computed and none may be offered as evidence of
  biological interaction.

## 5. The ICC estimand

`DD` cancels every effect additive in target or in ligand. The dependence that
survives into `DD` is therefore the dependence of the **additively-adjusted**
measurement, not of the raw affinity.

Within each endpoint separately, remove additive target and ligand effects by
one two-way least-squares fit on the cell means,

```text
ybar[cell] = mu + a[target] + b[ligand] + r[cell]
```

fitted within panel so that panel-level offsets cannot leak into `r`. Then
decompose the adjusted residual with a nested random-effects model

```text
r = u[cluster] + v[panel in cluster] + w[cell in panel] + e[replicate]
rho = var(u) / (var(u) + var(v) + var(w) + var(e))
```

`rho` is the intra-cluster correlation that enters `DEFF = 1 + (m-1) rho`.
Variance components are estimated by the standard unbalanced one-way
ANOVA/Henderson-III moment estimators applied at each nesting level, with
negative estimates truncated at zero and the truncation reported.
`var(e)` is identified by exact-assay replicates only.

### 5.1 Upper confidence bound

The one-sided 95% upper confidence bound is the 95th percentile of a
**cluster bootstrap**: resample the `G` dependency clusters with replacement,
`10,000` draws, seed `20260903`, recomputing every variance component inside
each draw. Clusters are the inference unit. Rectangles, cells and measurements
are never resampled as IID, and `G` is `36` for Ki and `12` for Kd.

A bootstrap over 12 clusters is coarse by construction. That coarseness is a
property of the source design, not a defect of the estimator, and the interval
width is reported rather than smoothed.

## 6. Frozen Gates

```text
G1  Ki   UCB95(rho) <  0.0915        the X0-B breakeven rho*
G2  Kd   UCB95(rho) <  0.0164        the X0-B breakeven rho*
G3  no cluster dominates: in the capped design actually used for inference at
    the estimated rho, the largest cluster contributes <= 0.25 of total weight
G4  recomputed effective n under the actual nested design >= 245
```

`G1` and `G2` are inherited unchanged from X0-B and may not be weakened,
averaged, or replaced by a point estimate. Endpoints are judged independently:
one may pass while the other fails, and a failing endpoint stops there.

`G3` is evaluated on the **capped** design because X0-B's per-cluster cap is
the mechanism that controls domination; the uncapped share is reported
alongside for transparency. The `0.25` value is the same giant-component cap
this programme already uses for closure components in the C0 correspondence
audit, reused unchanged rather than chosen here.

`G4` uses `DEFF = 1 + (m-1) rho` with `m` the mean capped cluster influence and
`rho` the point estimate, exactly the X0-B formula, evaluated at the measured
`rho` rather than on a grid.

## 7. Terminal verdicts

Exactly one, by earliest failed precondition:

```text
X1_DATA_CONTRACT_INVALID        recovery, hashing, firewall or coverage fails
X1_INTERACTION_UNDERDETERMINED  G1, G2, G3 or G4 fails for both endpoints, or
                                for the only endpoint that reached them
X1_ICC_PRECONDITION_PASSED      at least one endpoint clears G1/G2 and clears
                                G3 and G4; that endpoint alone proceeds
```

Only `X1_ICC_PRECONDITION_PASSED` authorises X1B, and only for the endpoints
that passed.

## 8. Stopping rules

One run. If an endpoint fails its ICC threshold, that endpoint stops. It is
forbidden to weaken the threshold, merge Ki with Kd, average the two, switch to
a point estimate, acquire another dataset, re-pack the design, change the
cluster definition, or train any model in response.

X1A trains nothing. No model, no parameters, no gradients.

## 9. Boundary

`model/`, production `scripts/`, `theory/`, CSMO, Band, the mesh, production
`z` and the operator

```text
A(F, z) = K(B(z) F(z))
```

are unmodified. All experimental code stays under `research/`. Historical
artifacts are preserved and no failed verdict is rewritten.

X1A can identify at most whether the source dependence structure permits a
crossed interaction test. It identifies no interaction, no affinity energy, no
selectivity, no causal mechanism, no few-shot capability and no `z` admission.
