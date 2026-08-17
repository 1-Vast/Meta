# Stage M0 report — ChemBERTa ligand-side LM probes: REJECTED at identifiability

No training of the DTA model. Frozen ChemBERTa-77M pooled ligand embeddings
(local snapshot DeepChem/ChemBERTa-77M-MLM, 600-dim, attention-masked mean
over SMILES tokens; manifest: tools/runtime/chemberta_ligand_pooled/).
Probes with component-fold selection on meta_train; meta_val read once;
meta_test never constructed. Authority: M0_CHEMBERTA_PROBES.json.

## Q1: within-target ordering

Frozen linear probe on within-target centered affinity, the Stage C protocol:

| fold decay | train-fold r | meta_val r (component bootstrap) |
|---|---|---|
| 1.0 (selected) | +0.120 | **+0.1472 [-0.0261, +0.3179]** |

The interval crosses zero and the mean sits below the occupancy record
(+0.2182 [+0.0751, +0.3670]): the ligand-side LM carries no NEW
within-target ordering signal beyond the incumbent trunk representations.

## Q2: target level

Panel-mean ChemBERTa embedding -> episode level: the probe collapses to the
grand mean (all weight decays give identical fold MSE; meta_val level MSE
2.1547 = the meta_train-constant baseline exactly). The ligand LM carries
zero measured cross-component level signal.

## Verdict

The ligand-side language-model input family is REJECTED at the
identifiability gate for both level and shape. The external-representation
ledger is now: sequence LMs (ESM-150M/650M frozen, LoRA-tuned), ligand-side
LM (ChemBERTa), structure/pocket priors, panel composition, assay
covariates (journal/publisher, endpoint, counts) — every locally testable
legal family measured, none breaks the k=0 level wall. The bounded
conclusion (report/BOUNDARY_20260817_NIGHT.md) now covers the full
locally-available legal-input space.
