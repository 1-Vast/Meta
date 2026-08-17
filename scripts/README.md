# Scripts

## Active workflow

- `qpsmp_data.py`: governed episodic materialization and split seal.
- `train_qpsmp.py`, `evaluate_qpsmp.py`: retained baseline training/evaluation.
- `train_level_shape.py`, `train_reltransport.py`, `train_grammar_shape.py`:
  formal R-series model trainers.
- `stageR0_*`, `stageR2_*`, `stageR3_*`, `stageR6_*`, `stageR9_*`, `r14_*`:
  retained analyses required by the compact evidence ledger.
- `audit_research_record.py`, `verify_project.py`, `project_status.py`:
  evidence and repository integrity.

## Data and structure utilities

Corpus acquisition, homology governance, protein/ligand banks, double-cold split,
structure supervision and geometry audits remain because they produce governed
inputs. Structure utilities do not imply current DTA geometry coverage.

## Removed workflow

Pre-R0 diagnostic/smoke scripts, research-only C/Q observables, analytic solvers,
old relative/locality trainers and direct-shape training were deleted after their
verdicts were consolidated. Restore them only from the commits named in
`archive/README.md` and only on a temporary branch.
