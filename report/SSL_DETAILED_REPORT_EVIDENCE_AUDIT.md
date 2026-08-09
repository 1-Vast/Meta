# SSL detailed report: evidence audit and failure triage

Updated: 2026-08-09.

## Scope

This audit evaluates `METASIEVE_SSL_DETAILED_ANALYSIS_REPORT_2026-08-08.md`
against the repository, all local Git branches, the GitHub branch tips, the
frozen theory handoffs, and primary dataset/paper sources. It does not import
unverified metrics into the project ledger.

## Evidence verdict

| Report claim | Repository evidence | Verdict |
|---|---|---|
| S0-S4 aggregate structural observability | Recoverable at commit `608decf`; consolidated in `history.md` | VERIFIED HISTORICAL |
| S5 pair-local v1/v2 results | No code, input manifest, output JSON, checkpoint, or Git commit found | EXTERNAL CLAIM, NOT REPRODUCED |
| minimal bridge retraining result | No training trace, confirmation manifest, or checkpoint found | EXTERNAL CLAIM, NOT REPRODUCED |
| S7/S8 pose-aware v3/v4 results | No preregistration, model code, split hash, or result artifact found | EXTERNAL CLAIM, NOT REPRODUCED |
| S9 Kd/Ki risk results | No sealed affinity manifest, fold hash, prediction file, or metric artifact found | EXTERNAL CLAIM, NOT REPRODUCED |
| NCE score depends on the negative proposal | Follows from the Bayes classifier log-density ratio | SUPPORTED THEORY |
| score/force matching leaves a system-dependent additive constant | Follows from gradient invariance under `E(x,P,L)+C(P,L)` | SUPPORTED THEORY |
| the additive constant alone explains the claimed S9 null | Requires the missing S8/S9 artifacts and ablations | PLAUSIBLE, NOT IDENTIFIED |
| MISATO/PLINDER automatically solve affinity identification | Neither supplies complete bound/unbound partition functions and a clean universal affinity gauge | NOT SUPPORTED |

The formal repository state therefore remains `S5 REGISTERED_NOT_RUN`. The
external report is valuable as a hypothesis and code-review input, but it does
not supersede the immutable evidence ledger.

## Reproducibility blocker

Before any S5-S9 conclusion can be promoted, the following exact artifacts are
required:

1. source commit and branch for every stage;
2. immutable input manifests and raw-source hashes;
3. score-blind split, closure, and exposure manifests;
4. model/config/checkpoint hashes and complete training traces;
5. per-unit predictions, controls, bootstrap units, and metric JSON;
6. label-access audit for S9;
7. exact derangement and ligand-relevance maps;
8. environment and dependency lock.

If these artifacts cannot be recovered, the only valid continuation is the
existing pair-local S5 preregistration. Recreating a similar experiment is a new
run, not verification of the reported numbers.

## Scientific assessment

### Supported

- A pose-discrimination score is not automatically a binding free energy.
- Native-decoy NCE identifies a log ratio relative to a declared decoy
  distribution; changing that distribution changes the learned score.
- Coordinate score matching identifies derivatives, not cross-system energy
  offsets.
- Absolute affinity additionally depends on ensemble, solvent, protonation,
  entropy, standard state, and assay context.
- The executable multi-view CSMO preserves simplex/operator legality, but it
  does not inherit the full tensor-sieve approximation and uniform convergence
  theorem.

### Not yet proved

- That pose-free P1B pair-local states fail: current `main` has not run S5.
- That pose-aware S8 succeeds: the reported prospective artifacts are absent.
- That gauge is the dominant S9 failure rather than label noise, endpoint
  mismatch, representation error, ensemble omission, or split shift.
- That five named thermodynamic coordinates are observable, stable, or
  identifiable from public deployment inputs.
- That MISATO trajectories provide an unbound-state or absolute-affinity
  anchor. They provide valuable MD/QM supervision, not a complete free-energy
  identification experiment.

## Dataset feasibility

### MISATO

MISATO is appropriate for a structure/trajectory observability study: it
contains molecular-dynamics and quantum-chemical information for roughly
20,000 experimental protein-ligand complexes. It can test whether ensemble
statistics are learnable. It cannot by itself identify cross-system binding
free energies, because trajectory supervision still lacks a common affinity
zero point and complete bound/unbound thermodynamic cycles.

### PLINDER

PLINDER is appropriate for protein, pocket, ligand, interaction-similarity and
apo-linkage governance. The current public repository warns that its BindingDB
affinity annotation is disabled because of a parsing bug. It must not be used as
an unverified affinity-label source. Version, annotation licence, split version,
and known bugs must be pinned before acquisition.

### Quantitative gauge anchors

A gauge-fixing experiment needs independently governed quantitative edges:
same-target, same-endpoint, assay-compatible matched-series differences, or
validated relative/absolute free-energy calculations. The entire connected
matched-series graph must remain in one split. Kd and Ki remain separate;
IC50 is excluded from the first thermodynamic test.

## Theory compatibility

The general identifiability theory and the executable law-valued operator are
related interfaces, not an already proved reduction. A future biological state
may enter the operator only after it provides:

- a declared compact domain and off-coverage abstention;
- permutation and gauge invariance of every reported quantity;
- a conservative outer uncertainty object, not merely a point embedding;
- no more than `k` support-derived continuous degrees of freedom;
- source closure-OOF affinity increment and a sealed transfer result.

Names such as `DeltaG_solv` or `DeltaS_ensemble` do not create biological
semantics. Each coordinate needs an executable intervention/measurement
contract. Basis-dependent latent coefficients cannot be given physical names.

## Decision tree

1. **R-VERIFY:** recover and reproduce S5-S9 artifacts. No recovery means the
   report remains external evidence.
2. **S5:** if necessary, execute the existing pair-local frozen-P1B
   observability preregistration.
3. **G-DATA:** only after a verified pose-aware structural state and a verified
   affinity null, census ensemble coverage and quantitative gauge-anchor graphs
   without reading evaluation labels.
4. **G0-A:** test single-pose versus ensemble structural observability. This can
   identify an ensemble contribution, not affinity.
5. **G0-B:** add one frozen, low-capacity difference anchor and test fresh
   closure components. Compare against ligand retrieval and capacity-matched
   baselines.
6. **STOP:** if the anchored state fails correct-versus-baseline and
   correct-versus-wrong-protein risk contrasts in independent domains, stop the
   absolute-affinity structural route.

No real-affinity, DAVIS, few-shot, production-`z`, or P2-P4 stage is authorized
by this audit.

## Code hardening completed from the audit

The following defects were independently reproduced and repaired without
changing the frozen mathematical operator:

- `build_band_operator` now uses the local frozen `context_index` instead of a
  deleted `biological.kappa` import;
- statistic vectors outside `[0,1]` fail closed instead of being silently
  clamped by the sieve;
- zero-atom ligands and zero-residue bridge inputs fail closed;
- deployment manifests require every semantic artifact hash and nonempty
  frontend/source hashes.

These are contract repairs. They do not constitute an SSL or affinity result.

## Primary references checked

- NERE, NeurIPS 2023: https://proceedings.neurips.cc/paper_files/paper/2023/file/6a45a1b0697ee086bd8bf494cacc6567-Paper-Conference.pdf
- MISATO, Nature Computational Science 2024: https://www.nature.com/articles/s43588-024-00627-2
- PLINDER repository and known issues: https://github.com/plinder-org/plinder
- PLINDER preprint: https://www.biorxiv.org/content/10.1101/2024.07.17.603955v1
- Leak-Proof PDBBind: https://arxiv.org/abs/2308.09639
