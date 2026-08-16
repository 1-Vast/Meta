# Papyrus F0 preregistration: raw-provenance quantitative-anchor audit

Date: 2026-07-27  
Candidate count: agent-proposed candidate 3 of at most 3.  
Role: final data/supervision gate before any Papyrus pretraining or Mamba architecture comparison.

## Question

Can the CC-BY-SA Papyrus 05.7++ release supply **raw document-resolved quantitative evidence** for a
document-disjoint, protein-conditioned local residual anchor? This is a data-identifiability audit,
not a model benchmark.

Papyrus is an aggregation. Its high-quality label alone is insufficient: a numeric value may count
only if it can be assigned to one exact source/document without ambiguity. The audit must reject a
row rather than infer a per-document value from semicolon-combined metadata.

## Frozen input and rows

Use only `05.7++_combined_set_without_stereochemistry.tsv.xz` and
`05.7_combined_set_protein_targets.tsv.xz` from Zenodo record 13987985, version 05.7, recorded in
`dataset/public/papyrus_05_7/raw/`.

Retain one candidate observation only when all conditions hold:

1. Papyrus quality is `High` or `Very High`;
2. target is human wild type with a mapped UniProt accession and one nonempty ChEMBL classification;
3. relation is exact `=` and `pchembl_value_Mean` is finite;
4. exactly one of `type_Ki` and `type_KD` is one, all other endpoint-type flags are zero;
5. `source`, `doc_id`, `all_doc_ids`, and `pchembl_value` each contain no semicolon;
6. `pchembl_value_N == 1`;
7. standardized SMILES and connectivity are present.

Rows are keyed by `(connectivity, target_id, doc_id)`. Exact duplicate keys collapse only if their
numeric pK agrees to 0.01; otherwise every conflicting key is excluded. Missing/censored records
never become negatives and endpoint types never pool.

## Frozen support units

A document replication is a `(connectivity, target_id)` cell represented by at least two retained,
distinct `doc_id` values. It is the minimal unit required to compare the same local interaction under
independent documents. A target is replication-supported only if it has at least ten such cells.

An anchor-capable ligand is measured on at least five retained targets, and a target-contrast ligand
has at least two retained targets. Target classes are the first nonempty ChEMBL classification level;
classification values containing `;` cannot define a class.

No claim of document-family independence is made at F0. F0 asks the prior question: whether the
release preserves enough raw per-document numeric cells to construct that stricter firewall at F1.

## Frozen gates

All must pass:

1. at least 100 human wild-type exact-pK targets across ten top-level classes;
2. at least 100,000 retained raw observations and at least 20,000 unique parent--target cells;
3. at least 500 document-replicated parent--target cells across at least 30 replication-supported
   targets and ten top-level classes;
4. at least 1,000 anchor-capable ligands and 1,000 target-contrast ligands;
5. source and exact document fields are nonmissing for 100% of retained cells;
6. current-run FORT development/confirmation labels remain unread, historical confirmation remains
   quarantined, and `sealed_test_consumed=false`.

Pass verdict: `PAPYRUS_F0_RAW_PROVENANCE_SUFFICIENT_AUTHORIZE_F1`.

Fail verdict: `PAPYRUS_F0_RAW_PROVENANCE_INSUFFICIENT_STOP`.

On failure, do not lower replication thresholds, use Papyrus aggregated multi-document values as raw
observations, pool Ki/Kd/IC50, or start Mamba/pretraining. First issue the required overall failure
report, then reopen only data/measurement-design exploration rather than an architecture rescue.

