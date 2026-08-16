# RECRO-DTA v2 Stage L0 — decision

Verdict: **`RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE`**. The R0 cross-document residual signal (+0.334) is
an artifact of provenance duplication (the same measurements recorded under different ChEMBL document
IDs), not leakage-free biological cross-environment ligand reordering. R0-b and R1 are **not
authorized**.

Preregistration `reports/active/recro_v2_preregistration.md`; runner `research/recro_l0.py`; result
`reports/active/recro_l0.json`. Source: licensed raw `chembl37_pKi.jsonl.gz` (ChEMBL CC BY-SA, sha256
`f36dd43d…bbd0`), TRAIN split only; no dev/confirmation/sealed label read. Seed 1729.

## Required outputs

**1. Raw-record lineage (L0-A).** Analysis restricted to train (target,conn); the registry aggregates
target-ligand-endpoint then assigns splits, so every raw replicate of a cell shares one split, and no
development/confirmation raw record enters the train-only replication graph. No split-crossing detected.

**2. Provenance-family audit (L0-B).** Across 8,432 cross-document within-target comparisons:
**79.75% of co-measured (target,ligand) cells are exact-value duplicates** (|ΔpK| < 0.01; median spread
0.000). Clustering the 984 documents by pooled shared-cell identity (union-find, ≥50% identical links a
pair) yields **463 provenance families**, and **91.6% of comparisons fall within a single family** —
i.e. nominal "cross-document" replication is overwhelmingly the same data under different IDs.

**3. Cross-fitted vs non-cross-fitted (L0-C).** Properly potency-removed residual with a well-estimated
global ligand potency reproduces at +0.334 on all doc-pairs but collapses to +0.090 family-disjoint
(below). A leave-target-out variant with noisy per-ligand nuisance *inflates* to +0.873 — it fails to
remove potency and manufactures apparent signal, a positive leakage sentinel that poor cross-fitting is
dangerous here.

**4. Leakage-sensitivity ladder (L0-D), residual component-macro (grouped bootstrap):**

| isolation level | residual ρ | n_comp |
|---|---|---|
| all doc-pairs (= R0) | +0.3345 [+0.219,+0.446] | 131 |
| assay-disjoint | +0.3345 [+0.220,+0.445] | 131 |
| **document-family-disjoint** | **+0.0901 [−0.0557,+0.2358]** | 78 |

Assay-disjointness does not help (same paper uses several assay IDs); only provenance-family
disjointness removes the artifact, and it collapses the signal 73%. Raw (potency-inclusive) ordering
still reproduces family-disjoint (+0.700) — generic ligand potency (B0) is genuinely reproducible
across independent labs; the *target-specific* part is the artifact.

**5. Matched negative controls (L0-E).** Family-disjoint matched wrong-target residual reproducibility
is +0.1645 [−0.064,+0.381] (n=21). Primary contrast **residual_correct − residual_wrong = −0.0744**:
on provenance-disjoint data the correct target does not beat the matched wrong target.

**6/7. R0-b and R1.** Not run / not authorized — L0 failed its frozen pass criterion.

**8. Effect sizes / power.** Family-disjoint residual +0.090, LCB −0.056 (< 0); correct−wrong −0.074;
signal drop vs R0 73.1%; family-disjoint components n=78 (residual), n=21 (wrong-target). The clean set
is smaller and noisier, so a *small* biological effect below detection cannot be strictly excluded — but
there is no positive evidence for one, and the correct-vs-wrong contrast is negative.

**9. Final verdict: `RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE`** (frozen gates: family-disjoint residual
LCB > 0 FALSE; beats matched wrong by 0.03 FALSE).

## What this overturns and what stands

* **Overturned:** the RECRO R0 nuance "target-specific reordering reproduces across documents (+0.334)"
  and the earlier memory refinement claiming within-target cross-document reproducibility is biological.
  It was provenance duplication. AMOB's O0 +0.434 (staged Harmonic data with no document IDs) is now
  strongly suspected to be the same artifact and cannot be trusted as biological.
* **Stands:** generic ligand potency (B0) reproduces across independent provenance families
  (+0.700 raw family-disjoint) — reaffirming B0 as the only genuinely reproducible signal. Target-
  specific protein-conditioned reordering remains unidentified once provenance is controlled.

## Program implication

This is the leakage attribution the whole program needed: the recurring "cross-document / cross-environment
reordering signal" (CROSSDOC, AMOB, RECRO R0) is dominated by ChEMBL document-ID duplication of the same
measurements. Under strict provenance-family isolation the target-specific signal is null. The binding
conclusion tightens: **no admissible open source has yet shown leakage-free, provenance-disjoint,
target-specific ligand reordering.** A credible route now requires provenance metadata (DOI/PMID/patent/
authors) to build true document-family firewalls, or a prospectively measured factorial panel with
independent environments by construction. No model is authorized; the reproducible predictor remains B0.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing; L0 read train only).
