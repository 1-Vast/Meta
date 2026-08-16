# CROSSDOC candidate-3 decision

Verdict under the frozen gate: **`CROSSDOC_C3_FAIL_PUBLIC_OVERLAP`**.

The external mechanism signal is large and passes every destruction test, but the independent
coverage requirement fails. Exact single-document, endpoint-separated ChEMBL-37 train rows yield
only 13 usable target-document units and 11 strict homology components, versus the frozen minimum
30 and 25. The severe loss occurs after requiring exact compound overlap, at least five jointly
observed ligands, source-specific rank residualization, and a same-KLIFS-group leave-target-out
baseline.

Across the 11 components that remain:

- direct KIRHub-to-ChEMBL rank correlation: +0.5337 [+0.3101, +0.7336];
- primary group-residual correlation: +0.4946 [+0.3156, +0.6727];
- global-residual sensitivity: +0.4600 [+0.3041, +0.6139];
- ligand-permutation null: +0.0163 [-0.0136, +0.0426];
- observed minus ligand permutation: +0.4783 [+0.3031, +0.6524];
- within-document/group target permutation: null mean +0.1664, one-sided p=0.001499.

Two documents have positive target-unit means: CHEMBL1201862 pKi is +0.5926
[+0.4113,+0.7750] over eight units; CHEMBL1908390 pKd is +0.2500
[-0.0750,+0.5250] over four. CHEMBL3991601 pKd contributes only one usable unit and cannot support
an interval.

Thus target-specific ligand reordering does replicate directionally beyond global ligand potency
and coarse kinase taxonomy in the small exact-overlap subset. It is a credible signal, not adequate
external confirmation. No document is dropped and the coverage gate is not relaxed.

Protocol correction: a prior schema-inspection command displayed five labeled ChEMBL
`confirmation` rows. They did not enter this experiment, which uses parquet predicate pushdown for
`train` only, but that confirmation partition is permanently quarantined.
`confirmation_labels_read=true`; `sealed_test_consumed=false`.

Candidate ledger: 3/3. No predictive model is authorized.
