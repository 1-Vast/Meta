# Model package

The active model surface is intentionally small.

- `encoders.py`: ligand graph and protein-bank projections.
- `bpsf.py`: bipartite atom/protein interaction field and latent readout.
- `qpsmp_meta.py`: retained QPSMP/BPSF endpoint and episodic adaptation contract.
- `interaction_grammar.py`: interaction-first comparator trunk.
- `similarity_grammar.py`: retained A0 and fixed Morgan/Tanimoto transport.
- `level_shape.py`: R3/R4 level-shape Pareto models.
- `reltransport.py`: R5-R10 formal relative-transport arms.
- `cartesian.py`: O(3) algebra module; tested but not a current performance path
  because the DTA corpus has no common-frame complexes.
- `runtime.py`: shared masked tensor utilities.

Removed families include analytic/frozen operators, relative grammar, locality
grammar and direct shape. Their scientific verdicts remain in the R-series report
and Git history; they are not valid CLI choices.

The next model family is not implemented yet. Its exact contract is in
`report/NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md`.
