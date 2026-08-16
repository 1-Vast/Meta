# RECRO-DTA — preregistration (Stage R0 replication-graph audit)

Date: 2026-07-26. RECRO-DTA is the assay/document-firewalled refinement of AMOB (not a new module):
distinguish a shared biological target-conditioned ordering rule from document/assay/campaign-specific
correlation, via **replication across independent environments**. This preregisters Stage R0 with
frozen gates before any result is read.

## Admissible data (resolves AMOB's blocker)

AMOB's O0 (+0.434) could not be certified because the staged Harmonic CSVs had no assay/document IDs.
RECRO R0 uses instead the **licensed raw per-record ChEMBL-37 extract**
`dataset/public/chembl_37/processed/chembl37_pKi.jsonl.gz` (fields: target, smiles, pK, assay_id,
doc_id; ChEMBL CC BY-SA 3.0), which the aggregated registry collapsed to medians but which retains
per-document granularity. Restricted to the registry **TRAIN** split (zero-overlap with dev/confirmation
on all 7 firewall axes); no dev/confirmation/sealed label read. Feasibility probe: 248 targets have a
document-pair sharing >=5 ligands (vs CROSSDOC's 11-13 cross-source units).

## Estimand and test

Environment = ChEMBL document. For target t, per document d, per ligand: median pK. Generic ligand
potency `b(conn)` = mean over target-level values; residual `resid = pK - b(conn)` removes generic
potency, leaving target-specific ordering (the double-centered `z^perp` interaction, within-target
form). Primary R0 test: does the **potency-residualized** within-target ligand ordering reproduce across
independent documents (cross-document Spearman over shared ligands), with enough entity-disjoint
independent units to be powered?

Controls (R0 step 7): raw (potency-inclusive) reproducibility; **wrong-target residual** (t's document
vs a different target's document over shared ligands) — protein-specificity; sign agreement before/after
potency removal. Unit = homology component; entity-disjoint packing reuses no target or document.

## Frozen gates

* **Power:** entity-disjoint independent units >= 25 **and** empirical MDE80 <= 0.10. Else
  `RECRO_REPLICATION_GRAPH_INSUFFICIENT`.
* **Cross-environment:** component-level residual reproducibility grouped LCB95 > 0 **and** exceeds the
  wrong-target residual control by > 0.03. Else `RECRO_ORDINAL_SIGNAL_NOT_CROSS_ENVIRONMENT`.
* Both pass -> `RECRO_R0_PASS_AUTHORIZE_R1` (authorizes exactly one short R1; no model trained in R0).

R0 authorizes nothing else (no R1 model, no Bayesian few-shot, no confirmation/sealed access). A pass
does not claim transfer; R1 (simultaneous unseen-target + unseen-scaffold) is a separate gate.
`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing).
