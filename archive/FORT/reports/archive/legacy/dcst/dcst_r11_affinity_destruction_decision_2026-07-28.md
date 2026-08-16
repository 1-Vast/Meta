# DCST-R11 affinity-destruction identification decision

Date: 2026-07-28  
Decision: `STOP_SOURCE_ENERGY_TRANSFER_ADVANCE_TO_STRUCTURE_PROMPT`

## Result

R11 reproduced the R6 segment mechanism but certified `0/4` privileged PMCE
bands versus `1/4` PMCE-NoPriv. The bounded train-split three-way
identification loss did not generalize to held-source spectral directions.

The privileged candidate had one positive true utility of `0.20219`, but its
target-destroyed utility was larger (`0.27647`). Another direction had
ligand-destroyed utility `0.20746` versus true `0.02091`. All confidences were
zero. The no-privileged active direction had true utility `0.14835` and
confidence `0.67987`.

Wall time was `191.703 s`; peak allocated CUDA memory was `564.8 MiB`. No new
downstream affinity label was loaded.

## Decision

R6 established a privileged pair-specific structural measure, but R7-R11 show
that transporting a source-fitted affinity energy over absolute positions,
learned roles, fixed roles, or continuous content does not yield a
privileged-specific held-source energy certificate.

The next experiment stops transferring source affinity energies. It transfers
only the already accepted R6 structural measure as a frozen Stage-2 prompt,
then fits a downstream observation head. This matches the user's intended
two-stage boundary more directly: Stage 1 learns interaction information from
the high-quality structural source; Stage 2 learns the dataset-specific
affinity mapping.

