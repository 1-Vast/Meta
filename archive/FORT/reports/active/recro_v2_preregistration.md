# RECRO-DTA v2 — preregistration (Stage L0 leakage attribution; then R0-b/R1)

Date: 2026-07-26. RECRO-DTA v2 is the leakage-attributed refinement of RECRO. The immutable prior
result stands: raw ChEMBL pKi gave a cross-document residual-ranking correlation ≈ +0.334, correct >
wrong-target, and the original R0 power gate failed on single-document-pair variance. That result does
**not** authorize R1 and must not be rescued by relaxing the original gate. L0 runs first and decides
whether the +0.334 is biological or an artifact (contamination / residualization / provenance /
pseudoreplication / med-chem selection / assay reproducibility) before any R0-b or model training.

## Data (admissible)

Licensed raw per-record `chembl37_pKi.jsonl.gz` (target, ligand, pK, assay_id, doc_id; ChEMBL
CC BY-SA), restricted to the registry TRAIN split (zero-overlap with dev/confirmation on all 7 axes).
No dev/confirmation/sealed label read.

## Stage L0 tests and frozen gates

* **L0-A firewall:** analysis restricted to train (target,conn); every raw replicate of a cell shares
  one split; no dev/confirmation raw record enters. Fail -> `RECRO_ANALYTICAL_LEAKAGE_DETECTED`.
* **L0-B provenance:** report exact-duplicate cross-document cell fraction; cluster documents into
  provenance families by pooled shared-cell identity (union-find, identity fraction >= 0.5 links a
  document pair). All decisive tests use **document-family-disjoint** comparisons.
* **L0-C cross-fitting sentinel:** compare properly potency-removed residual (well-estimated global
  ligand potency) against a leave-target-out variant; a signal that appears only under poor/
  non-cross-fitted nuisance is `RECRO_RESIDUALIZATION_LEAKAGE_DETECTED`.
* **L0-D leakage-sensitivity ladder:** identical comparisons at all -> assay-disjoint ->
  document-family-disjoint isolation; report effect + grouped CI at each.
* **L0-E matched wrong-target** on the family-disjoint set (same ligand pairs, matched exposure).
  Primary contrast: residual_correct − residual_matched_wrong.

**L0 passes only if** the component-grouped LCB95 of the family-disjoint residual reproducibility is
> 0 **and** it exceeds the matched wrong-target control by > 0.03. Verdicts:
`RECRO_CROSS_ENVIRONMENT_SIGNAL_CONFIRMED` / `RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE` /
`RECRO_RESIDUALIZATION_LEAKAGE_DETECTED` / `RECRO_ANALYTICAL_LEAKAGE_DETECTED` /
`RECRO_LEAKAGE_AUDIT_UNDERPOWERED`.

## R0-b and R1 (only if L0 passes)

R0-b: environment-BUNDLE consensus ranking per target (not single document-pairs); component-level
inference; freeze the R1 MDE from the actual R1 component-macro endpoint (unseen target + unseen
scaffold), not document-pair variance. R1: cross-fitted conditional ordinal operator with a frozen
out-of-fold ligand-pair offset and a double-centered protein-conditioned interaction; full destructive
control suite. Confirmation set is development-only (already read); a final claim needs the unopened
sealed test or a new disjoint panel. `sealed_test_consumed=false`; `confirmation_labels_read=true`.
