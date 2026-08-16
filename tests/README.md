# Maintained Test Suite

- Active model and meta-learning: `test_bpsf.py`, `test_qpsmp_meta.py`,
  `test_term_synthetic.py`, `test_train_qpsmp_meta.py`, `test_qpsmp_nested.py`,
  `test_interaction_grammar_synthetic.py`.
- `test_interaction_grammar_synthetic.py` holds the Stage 1 falsification
  gates for the active candidate: exact k=0 identity, query-specific and
  trainable k=1 label effect, support permutation invariance, query
  permutation equivariance, level-only abstention, shared-mechanism recovery,
  private-mechanism rejection, label linearity, no query-label input, no dead
  trainable branch, and a protein-conditioned zero-shot trunk gate.
- Data and governance: `test_data_*`, `test_*corpus*`, `test_*seal*`.
- Optional geometry: `test_cartesian.py`, `test_mechanistic_bridge.py`.
- Research lineage: tests named after modules under `research/`.
- Compatibility contracts: `test_*contract_repairs.py` and package-layout tests.

Run all maintained tests with:

```powershell
conda run -n drug python -m pytest tests -q
```

The `test/` package is a manual audit surface, not a second pytest suite.
