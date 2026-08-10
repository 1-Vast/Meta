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

The shared-linear witness is terminal-negative. The next stage is planned but
not implemented: learn one `d<=5` target coefficient subspace through
support/query episodes, then evaluate support-identifiable sections at
`k=1/2/3/5`.

## Admission rule

Research code enters `model/` or `scripts/` only after partner, affinity,
independent-transfer and support-identifiability Gates. No failed historical
route is current execution authority. Detailed terminated artifacts are
recoverable from Git and summarized in `history.md` and
`archive/FAILED_RESEARCH_ARCHIVE_20260810.md`.
