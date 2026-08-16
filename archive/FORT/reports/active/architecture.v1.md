# Active Architecture

Status: strict few-shot target adaptation on train and development roles only.
All numerical paths require CUDA in the `drug` environment.

## Model Modules

| File | Responsibility |
| --- | --- |
| `model/protein.py` | Bidirectional Mamba scans and landmark attention over frozen residue tokens |
| `model/ligand.py` | Morgan/physicochemical fingerprint encoding |
| `model/interaction.py` | Ligand-residue pooling and pair-feature construction |
| `model/posterior.py` | Generic Cholesky residual posterior |
| `model/reorder.py` | Protein-conditioned subspace, calibration, and ranking posterior |
| `model/likelihood.py` | Endpoint/source observation likelihoods |
| `model/ligandbase.py` | Ligand-only comparison baseline |
| `model/gradadapt.py` | Gradient-adaptation comparison baseline |

The primary proposed inference path is the Bayesian reordering posterior.
`main.py train` executes the registered Wave 1 comparison suite: shared B0,
calibration, ligand-only posterior, gradient baseline, protein-conditioned
Bayesian posterior, wrong support, permuted labels, and protein-free basis.
The first Wave 1 result failed the ligand-only and gradient admission gates, so
the flexible-kernel increment remains blocked.

## Data Modules

`scripts/preprocess.py` maps registered rows without reading labels by default.
`scripts/episode.py` performs outcome-blind support selection and closure.
`scripts/audit.py` freezes endpoint-specific k=5 rosters. `scripts/train.py`
reads only train/development rows and records CUDA utilization, power, memory,
and wall time.

No structure, pose, pocket, external model implementation, confirmation row,
or sealed outcome enters the active runtime.
