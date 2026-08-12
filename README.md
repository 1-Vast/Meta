# MetaSieve-DTA

MetaSieve investigates trainable meta-learning for few-shot drug-target affinity prediction on
proteins absent from source training.

## Current Candidate

The candidate innovation is the neural Quotient-Preserving Section Meta-Potential (QPSMP): shared
protein and ligand encoders, ligand-conditioned protein localization, a crossed endpoint scalar
potential, and a centered zero-preserving neural support adapter. The analytic ridge implementation
is a comparator and geometric diagnostic only.

The interface is trainable and unit/integration tested. Protein-specific G2, few-shot G3, biological
interpretation, and end-to-end Cold Target utility remain unadmitted.

## Canonical Reading Order

1. [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md): current status and authority map.
2. [Pure mathematical theory](theory/CURRENT_THEORY/PURE_MATHEMATICAL_THEORY.md).
3. [QPSMP model theory](theory/CURRENT_THEORY/QPSMP_COLD_TARGET_MODEL_THEORY.md).
4. [Implementation contract](docs/MODEL_IMPLEMENTATION_CONTRACT.md).
5. [Active protocol](task.md).
6. [Evidence ledger](report/EVIDENCE_LEDGER.md).

`history.md`, dated reports, `theory/FINAL_FROZEN_THEORY/`, and `archive/` are provenance, not current
authority. Failed results are retained and must not be interpreted as function-class impossibility
proofs.

## Inputs

Deployment inputs are protein sequence or a legal cached sequence/structure representation, ligand
molecular graph, declared context, and the unseen target's disjoint support observations. Target IDs
are lookup keys only. Query labels and persistent target-specific parameter memory are prohibited.

## Verification

```powershell
conda run -n drug python main.py verify tests
```

Large data and embedding assets are described in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).
