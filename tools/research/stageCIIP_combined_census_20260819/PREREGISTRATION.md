# Stage CIIP-0c combined local-panel census preregistration (2026-08-19)

Frozen BEFORE the combined census computation. Completes the CIIP
positive-control data census after KiRHub -> DATA BLOCKER and Davis ->
INSUFFICIENT ALONE. Read-only audit of the remaining LOCAL panels:
Anastassiadis 2011 MOESM23 (% remaining activity matrix) and Duong-Ly
mmc2/mmc3 (mutation metadata + % inhibition matrix). Same frozen
10-item checklist as stageCIIP_davis_census_20260819 applies per panel.

## Admissibility classification (frozen)

Each panel is classified against the CIIP surfaces:

- SAME-PLATFORM / SAME-ENDPOINT: admissible for CIIP-1A (within-parent
  capacity) and, if held-out parents >= 10, for CIIP-1B (transfer).
- CROSS-PLATFORM or CROSS-ENDPOINT: replication-only surface; NEVER
  merged into CIIP training; never relabeled pK/Ki/Kd.
- Davis (Kd binding) and Duong-Ly (% inhibition) are DIFFERENT
  endpoints: their rows are never pooled into one training surface.

## Frozen stop rules

- A panel is CIIP-1A-admissible iff usable pairs >= 20 and median common
  ligands >= 5 (single-point mutations only; multi-mutant tags such as
  V559D/T670I are counted separately and excluded from the single-mutant
  estimand).
- CIIP-1B-admissible iff additionally held-out parents (>=2 single-mutant
  pairs per parent) >= 10.
- If NO panel reaches CIIP-1A: verdict UNRESOLVED/DATA-BLOCKED; no
  CIIP training anywhere; no return to BindingDB exact-MMP retraining.
- Missing/censoring/duplicate/saturation and construct fields are
  quantified per panel, never guessed.

## Evidence marks

[V] primary local file (SHA-pinned); [I] interpretation; [U] unverified.
