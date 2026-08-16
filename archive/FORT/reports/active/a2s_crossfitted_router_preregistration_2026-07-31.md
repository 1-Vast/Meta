# A2S Cross-Fitted Router Preregistration

Date: 2026-07-31  
Status: protocol executed; target-balanced pKi router NO-GO

## Motivation

The scalar source adapter is not a sufficient universal transfer rule by
itself. After global feature alignment and target-macro source weighting were
fixed, the registered cross-fitted router was executed. It produced pKi gains
of `-1.144`, `-1.213`, and `-1.315` at `k=1,3,5`; this route is NO-GO. The
protocol does not authorize Mamba/Transformer capacity expansion.

## Cross-fitting contract

Let `H` be the abundant source-target set. For each pseudo-recipient `h in H`:

1. fit the pooled source predictor and source adapters on `H minus {h}`;
2. select `h` support/query units covariate-only;
3. compute every candidate-source routing feature on the support set; and
4. use query utility only as the meta-routing label for this pseudo-recipient.

The router is trained on these held-out source-target episodes, never on a
natural recipient target. Natural recipient labels enter the final evaluation
only through the declared support and query. Source and recipient target IDs
must be disjoint in every final artifact.

## Router features and controls

The small router may use only frozen, pre-registered quantities:

* source/recipient protein relation, with continuous similarity bins;
* source/recipient ligand chemical-space compatibility;
* source adapter residual loss on recipient support;
* source data depth, replicate reliability, and provenance quality; and
* endpoint compatibility.

Required controls are all-source averaging, random source, protein-similarity-
only, chemistry-similarity-only, support-compatibility-only, and the current
scalar adapter. A parameter-matched linear/ridge router is the first model;
neural routing is not admitted until this control passes.

## Abstention

The router outputs a source score and a calibrated transfer-risk estimate. It
may set transfer weight to zero and fall back to recipient calibration when the
cross-fitted estimated gain is not positive at the declared confidence level.
The gate threshold is frozen from source pseudo-recipient episodes and cannot
be selected on natural-recipient queries.

## Primary evaluation

Use pKi natural scarce recipients on the target-side single-cold primary track
with `k={1,3,5}`. Report separate chemical similarity, scaffold-cold,
global-drug-cold, and homology-cold tracks. pKd is a secondary replication and
is never pooled with pKi.

The primary estimand is target-macro transfer gain versus recipient calibration
at each `k`, with target-level/component bootstrap intervals. Required outputs
are RMSE, MAE, within-target Spearman, pairwise accuracy, negative-transfer
rate, abstention rate, and source-selection calibration. The router must beat
the current scalar source control and random-source control; a positive result
only on pseudo-recipient episodes is insufficient.

## Stop rules

Stop the router route if any of the following occurs:

* pKi k=3/5 transfer gain is not positive against recipient calibration;
* random source is as good as learned routing;
* the gate reduces coverage without reducing negative transfer;
* gains occur only in one similarity bin or only on pseudo-recipients;
* source/recipient, support/query, document, assay, or chemical closure fails; or
* the result requires unregistered target labels, confirmation/sealed data, or
  a larger neural architecture.
