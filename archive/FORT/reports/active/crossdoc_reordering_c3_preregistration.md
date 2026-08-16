# CROSSDOC candidate-3 preregistration

Date frozen: 2026-07-26, before cross-document affinity values were compared with KIRHub.

## Confirmation-partition contamination correction

Before this candidate was run, a schema-inspection command loaded the full ChEMBL registry and
printed its first five rows. Those rows were from the existing `confirmation` partition and
included affinity values. They were not used for document selection, hypotheses, thresholds, or
any calculation, but strict protocol treats observation itself as contamination. That partition
is therefore permanently disqualified as confirmation. Candidate 3 uses parquet predicate
pushdown to read `train` only and reports `confirmation_labels_read=true`. Any future confirmation
must be a newly isolated independent partition or source.

## Candidate and source

**Candidate 3 -- CROSSDOC (cross-document target-specific reordering replication).** This is the
third and final autonomous candidate in the reopened round (**3/3**). It changes the data condition
and estimand rather than the network. The user's additional literature synthesis motivates
assembling locally validated parts: panel-relative comparisons, endpoint-specific observation
models, and document isolation. No single source is required to have already solved strict
dual-cold prediction.

The label-blind overlap audit used only identifiers and counts. Among ChEMBL-37 `train` rows, 57/92
KIRHub compounds and 118 KIRHub targets overlap exactly. Three single-document panels were frozen
because each has at least 10 target profiles with five exact shared compounds:

1. `CHEMBL1908390`, pKd only: 316 exact-document rows; 17 targets with >=5 shared compounds;
2. `CHEMBL1201862`, pKi only: 243 rows; 15 targets with >=5;
3. `CHEMBL3991601`, pKd only: 160 rows; 10 targets with >=5.

Rows aggregating multiple documents are excluded. ChEMBL is CC BY-SA 3.0; KIRHub supplementary
data are used for internal academic analysis under their source terms and are not redistributed.
The KIRHub workbook and ChEMBL registry hashes must be recorded.

## Hypothesis and estimand

The hypothesis is that KIRHub's target-specific ligand reordering is biological enough to
replicate directionally in independent continuous-affinity documents, rather than being only
single-concentration saturation, global ligand potency, or KLIFS taxonomy.

For each document `p`, endpoint, target `t`, and shared ligand `l`:

1. convert values to within-target ranks separately within KIRHub and that document;
2. for every ligand, subtract the leave-target-out same-KLIFS-group mean rank in the corresponding
   source;
3. correlate the two residual profiles for target `t`, requiring at least five jointly observed
   ligands.

Raw pKi and pKd values are never pooled or converted. Only per-target rank-residual correlations
are combined. The statistical row is `(document, endpoint, target homology component)`;
overlapping target components across documents are additionally collapsed for the pooled summary.

This is external mechanism replication, not dual-cold prediction: exact target and ligand overlap
is required by design. It cannot be reported as affinity generalization performance.

## Controls and confounding firewall

- `direct_rank`: within-target cross-source Spearman before removing global/taxonomy effects;
- `group_residual_rank`: primary estimand;
- `global_residual_rank`: sensitivity using leave-target-out global rather than group mean;
- within-group target-label permutation of the dense KIRHub profiles, 2,000 fixed permutations;
- 32 fixed ligand-label permutations within each target, averaged as a negative control;
- each document and endpoint reported separately; no raw-scale pooling;
- exact single-document rows only; assay identifiers, replicate counts, and document identity
  retained;
- ChEMBL `development`, buffer, and the now-quarantined `confirmation` rows are excluded from this
  experiment; the accidental prior schema display is recorded above;
- no model, feature fitting, target-ID predictor, or selection by observed correlation.

Potential pseudo-replication from many ligands is prevented by target homology-component
aggregation and bootstrap. Three documents do not justify asymptotic document-level claims, so
document-specific intervals remain mandatory.

## Frozen success gate

CROSSDOC establishes independently replicated target-specific reordering only if all hold:

1. at least 30 target-document units, 25 distinct strict homology components, and two documents;
2. pooled homology-component macro group-residual correlation >= +0.10 with LCB95 > 0;
3. at least two of the three documents have positive point estimates and their combined
   target-unit LCB95 > 0;
4. within-group target-label permutation one-sided p <= 0.01;
5. ligand-label permutation mean is within +/-0.03 and observed minus that null has LCB95 > 0;
6. the global-residual sensitivity has positive LCB95.

Failure means current public cross-document overlap is too sparse/noisy to establish external
mechanism replication. No document may be dropped post hoc, no endpoint may be selected after
results, and no multi-document aggregate row may be added. A pass is a data/estimand breakthrough
only; it does not authorize a final predictive model, confirmation access, or multi-seed training.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (quarantined schema-display event).
