# Stage P0 report — protein function annotation (GO) probes: REJECTED at identifiability

No training of the DTA model. The local ProteinKG25 corpus (477,262 proteins,
47,229 GO terms, protein sequences + GO triples) was sequence-matched against
the governed targets: 313/387 targets (81%) receive GO annotation bags.
SGD linear probe with component-fold selection on meta_train; meta_val read
once; meta_test never constructed. Authority: P0_GO_PROBES.json.

| quantity | value |
|---|---|
| targets matched | 313 / 387 (81%) |
| GO vocabulary (train) | 3,580 |
| fold-selected weight decay | 1.0 (max shrinkage) |
| GO level MSE (meta_val) | **2.2699** |
| grand-mean baseline on covered targets | **1.4329** |

The GO bag probe is 58% WORSE than the constant: protein function
annotations carry no transferable cross-component level signal, consistent
with the D0 measurement that component identity transfers -1.1% of level
variance (GO terms are a finer encoding of the same family structure).

## Ledger update

The protein-function-annotation family (protein family context) is
falsified for the k=0 level. The external-representation ledger is now
complete for every locally available legal family: sequence LMs
(ESM-150M/650M frozen, LoRA-tuned), ligand LMs (ChemBERTa-77M),
structure/pocket priors, panel composition, assay covariates
(journal/publisher, endpoint, counts, documents) and protein function
annotations (GO). None breaks the k=0 level wall; the bounded conclusion
(report/BOUNDARY_20260817_NIGHT.md) stands.
