# Papyrus F0 decision: raw provenance is insufficient

Date: 2026-07-27

## Verdict

`PAPYRUS_F0_RAW_PROVENANCE_INSUFFICIENT_STOP`.

Papyrus 05.7++ does not preserve the document-resolved repeated parent--target
observations required for the registered target-conditioned local residual
anchor. No Papyrus pretraining, DTA training, or Mamba architecture comparison
is authorized from this corpus.

## Frozen inputs and firewall

The audit used only the frozen Papyrus activity and target metadata files named
in `papyrus_f0_preregistration.md`. Their SHA-256 values are recorded in
`papyrus_f0.json`:

- activity: `8004e0d1027a760f205b45264386f792e7d49658da39f77f52e660a6f19760dd`;
- target metadata: `832e564fb82daea0e4da79abcb44834d10104229382874e79915a1288d80783c`.

No FORT development or confirmation labels were read. Historical confirmation
remains quarantined and `sealed_test_consumed=false`.

## Provenance interpretation

The release stores both `doc_id` and `all_doc_ids`. The README says that
`doc_id` need not be the document from which an activity value originates.
Consequently, a row counted as a raw document observation had to have a single,
non-semicolon source, `doc_id`, `all_doc_ids`, and numeric value, with
`doc_id == all_doc_ids`. This is an ambiguity exclusion, not a relabeling rule.

Papyrus encodes an absent `type_other` flag as a blank string in this release.
The audit treats that blank only as an absent endpoint flag; mixed endpoint
strings remain excluded. Ki and Kd were retained as separate endpoint labels and
were never pooled.

## Results

| registered requirement | result | status |
| --- | ---: | --- |
| human WT exact-pK targets / top-level classes | 1,468 / 14 | pass |
| resolved observations / parent--target cells | 147,434 / 147,434 | pass |
| document-replicated parent--target cells | 0 / 500 required | fail |
| replication-supported targets / classes | 0 / 30 required; 0 / 10 required | fail |
| anchor-capable ligands (at least 5 targets) | 2,779 | pass |
| target-contrast ligands (at least 2 targets) | 30,636 | pass |
| source and document metadata completeness | 100% | pass |

The single failed registered gate is therefore decisive. There are no duplicate
`(connectivity, target_id, doc_id)` keys to resolve, and the strict replication
graph has zero components.

## Why this is not a threshold artifact

Before the `doc_id == all_doc_ids` ambiguity exclusion, 248,954 rows met every
other frozen activity, target, endpoint, structure, and atomic-metadata rule.
Of these, 101,520 had conflicting document fields. The stricter filter retained
147,434 rows. The result is not caused by that conservative choice: the complete
707,461-row activity table has exactly one aggregated row per
`(connectivity, target_id)` pair. Thus neither document field can produce two
raw, independent observations of the same parent--target cell.

The semicolon-combined values that could contain multi-document information are
explicitly excluded by the preregistration. Splitting them would manufacture
pseudo-observations, so it is not a valid rescue.

## Scientific consequence

This failure concerns raw measurement provenance, not model capacity. Papyrus
has substantial target, chemical, and endpoint coverage, but its aggregate table
cannot identify whether a local target-conditioned effect replicates across
independent documents. F1 through F4 remain closed, including the controlled
Mamba-in-Transformer ladder.

