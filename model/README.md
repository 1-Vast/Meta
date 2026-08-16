# Model package

## Active 2026-08-15 candidate: interaction grammar

`interaction_grammar.py` implements `InteractionGrammarModel`, selected by
`--arch grammar` in `scripts/train_qpsmp.py` and `scripts/evaluate_qpsmp.py`.
It replaces the bipartite pair trunk with:

- `ResidueEncoder`: pooled ESM slots plus a residue MLP and chemistry bias;
- `ContactGrammar`: atom-to-residue cross attention feeding a globally shared
  contact-type dictionary, giving per-pair type occupancies;
- an interaction embedding that is simultaneously the zero-shot readout and the
  few-shot kernel key;
- `TransferabilityTransport`:
  `f = f0(q) + shrink(n) * sum_k softmax_k(sim) * rho(q,k) * r_k`, where
  `rho in (0,2)` is a per-(query, support) transferability coefficient. A single
  support observation therefore produces a **query-specific** correction, which
  the retained kernel could not do. `rho == 1` with flat weights recovers the
  shrunken support mean, so level-only abstention stays inside the class.

Contracts: `k=0` and `adapt=False` return exactly the zero-shot endpoint;
support labels enter only as residual values; support order does not change the
output; queries are independent; query labels are not an input. Held-out
synthetic gates live in `tests/test_interaction_grammar_synthetic.py`.

`--arch bpsf` retains the previous `QPSMPBioModel` unchanged as the control arm.

## Retained 2026-08-15 baseline architecture

The active cold-target candidate is an interaction-first BPSF endpoint with a
positive label-locked residual kernel. Kernel keys are built from interaction
endpoints and ligand embeddings; support labels enter only through residual
values. It uses no ridge, closed-form solve, inner loop, or test-time gradient.
It jointly trains k={0,1,2,3,5}; k=0 is exact zero-shot, while k=1 is currently
limited to scalar residual calibration rather than ligand-specific SAR.

`cartesian.py` supplies an optional sparse O(3)-equivariant scalar/vector/
symmetric-traceless-rank-2 encoder. It is enabled only for declared coordinate
inputs as a bias on the same primitive field. The current BindingDB main bank
has no coordinates, so its active path is sequence+2D and makes no atomic 3D
claim.

D-MEMT and HyperSAR descriptions below are historical comparator context; they
no longer define the active architecture.

This package retains verified operator primitives and the current research
model implementation:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py` and
  `runtime.py` retain mathematical/operator diagnostics;
- `bpsf.py` implements the localized bipartite pair field, aligned mechanism
  slots, and its source-only contact/distance supervision. Geometry labels are
  never deployment inputs;
- `qpsmp_meta.py` integrates the pair trunk, mechanism-evidence set encoder,
  difference-only query transport, and scalar slot gates. The active package
  exposes no analytic solver;
- `cartesian.py` implements the optional sparse O(3) field encoder. It rejects
  cross-sample edges and must not be called with independently framed protein
  and ligand coordinates as if they were one complex.

Current validation is development-only and does not imply that the
target-specificity or Cold Target performance gates have passed. HyperSAR,
Set2Code, pooled/one-pass atom-residue, and analytic-section variants are
historical ablations only.

The Cartesian branch is PBCNet2.0/TensorNet-inspired, not a reproduction or
performance claim. The governed BindingDB corpus has no complex coordinates
for its deployment pairs.

Removed support encoders, arbitrary biological latents and terminal-negative
frontends remain recoverable from Git. Adding another PLM, attention branch or
adapter without a falsified information Gate is outside the model contract.

## Level-shape factorization (Stage R3, 2026-08-15)

- `level_shape.py` implements the level-shape factorized cold-target predictor,
  `f = ligand_prior(L) + target_level(P) + centered_interaction(P, L)`. The
  centered branch subtracts a per-protein constant computed from **learned
  anchor ligand embeddings**, so it has exactly zero mean in the anchor basis
  and cannot express a target-level offset, while every prediction stays
  inductive — no query-panel mean and no other transductive statistic. Its
  final `mix` layer and readout carry no bias because the centering makes any
  post-nonlinearity constant structurally unidentifiable.
- `similarity_grammar.py`, `relative_grammar.py` and `locality_grammar.py`
  remain as the tested Stage 5-8 mechanisms. `similarity_grammar.py`
  (`--arch similarity_only`) is the retained comparator; `relative_grammar.py`
  and `locality_grammar.py` are retained **rejected** arms.

`tests/test_level_shape.py` holds the 19 structural gates that must pass before
any level-shape training run: exact anchor centering, level-branch constancy
across queries, protein-blind ligand prior, per-query batch independence, exact
k=0 identity, support permutation invariance, query equivariance, label-locked
residuals, and full gradient coverage.
