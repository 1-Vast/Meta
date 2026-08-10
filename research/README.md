# Research index

## Active programme

The active research question is one minimal trainable meta-learning model for
unseen-target few-shot DTA. The source-data and linear-witness pipeline lives in
`crossed_interaction/`.

Implemented and retained:

- `bindingdb_cq_r0.py`: governed BindingDB metadata projection;
- `bindingdb_cq_r1.py`: exact endpoint extraction and audit;
- `prepare_bindingdb_cq_corpus.py`: strict conflict closure;
- `quotient_operator.py`: complete-panel quotient;
- `generate_tbasis_features.py`: frozen CUDA T-BASIS feature generation;
- `train_cq_linear.py`: panel-balanced shared-linear witness.

The shared-linear witness is terminal-negative.

`meta_fewshot/` holds the episodic stage. Its Phase 0 label-blind feasibility
census ran and failed closed with `FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE`: the
source split supports episodes (220 targets at `k=5`) but the evaluation split
yields only 16 held-out targets at `k=5`, below the declared 30, with
`MDE_d = 0.622` above the declared `0.600` ceiling. Leakage is zero on all five
audited axes. No model was preregistered and none was trained, so the
`d<=5` target coefficient subspace and the `k=1/2/3/5` support section remain
untested.

## Admission rule

Research code enters `model/` or `scripts/` only after partner, affinity,
independent-transfer and support-identifiability Gates. No failed historical
route is current execution authority. Detailed terminated artifacts are
recoverable from Git and summarized in `history.md` and
`archive/FAILED_RESEARCH_ARCHIVE_20260810.md`.
