# MetaSieve-DTA

MetaSieve investigates trainable meta-learning for few-shot drug-target affinity prediction on
proteins absent from source training.

## Current Candidate

The active candidate (`--arch grammar`, `model/interaction_grammar.py`) pairs a
protein-conditioned interaction trunk — atom-to-residue cross attention onto a
globally shared contact-type dictionary — with a label-locked residual transport
whose per-support coefficient depends on the query. The formal zero-shot output
is read only from the protein-ligand interaction endpoint. Few-shot support
labels enter through residual values only. Ridge, closed-form adapters, inner
loops, and test-time gradients are excluded. The previous BPSF model is retained
unchanged as `--arch bpsf`.

Three-seed development results reduce MSE by 12-18% at every k in {0,1,2,3,5}
against the retained baseline, but the governed admission gate was **refused**:
the gain is attributable to the zero-shot trunk and target-level calibration
rather than query-specific transfer, and within-target ranking degrades. See
[report/CURRENT_MODEL_EVIDENCE.md](report/CURRENT_MODEL_EVIDENCE.md).

An optional sparse Cartesian scalar/vector/rank-2 encoder feeds the same slot
contract only when legal common-frame coordinates are supplied. The active
BindingDB bank has no coordinates, so current results validate the sequence+2D
fallback rather than atomic 3D recognition.

The interface is trainable and unit/integration tested. Current three-seed
results are development-only; confirmatory Cold Target utility and Cartesian
performance remain unadmitted.

## Canonical Reading Order

1. [Project file organization](docs/PROJECT_FILE_ORGANIZATION.md).
2. [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md): current status and authority map.
3. [Active protocol](task.md).
4. [Current model evidence](report/CURRENT_MODEL_EVIDENCE.md).
5. [Evidence ledger](report/EVIDENCE_LEDGER.md).

`history.md`, dated reports, `archive/theory/`, and the rest of `archive/` are
provenance, not current authority. Failed results are retained and must not be
interpreted as function-class impossibility proofs.

## Inputs

Deployment inputs are protein sequence or a legal cached sequence/structure representation, ligand
molecular graph, declared context, and the unseen target's disjoint support observations. Target IDs
are lookup keys only. Query labels and persistent target-specific parameter memory are prohibited.

## Verification

```powershell
conda run -n drug python main.py verify tests
```

Large data and embedding assets are described in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).
