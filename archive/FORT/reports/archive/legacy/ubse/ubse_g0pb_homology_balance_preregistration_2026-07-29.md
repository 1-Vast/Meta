# UBSE-G0PB overrepresented-homology removal preregistration

Date: 2026-07-29  
Status: frozen after G0P failure and before balanced recomputation

## Binding correction

G0P passed source scale, resource ceiling, conflict-free packing, residual
training scale, target-domain support, and firewall. It failed only:

- largest conflict component `28.0397% > 20%`;
- largest homology panel share `6.7618% > 5%`.

G0PB performs one removal-only correction. Starting from the exact G0P panel
enumeration and homology assignment, calculate each homology component's panel
share. Remove all panels belonging to a homology component whose **pre-removal
share is greater than 5%**. Do not iteratively identify new blocks and do not
remove a scaffold, PubMed, panel, or component based on any post-removal
result.

Then recompute the exact G0P conflict graph, deterministic packing, first-88
audit closure, residual training substrate, target support, and all seven G0P
gates with unchanged thresholds.

## Decision

Pass:
`REQUEST_UBSE_G1_CENTERED_CONTACT_STUDENT_PREREGISTRATION`.

Failure:
`STOP_UBSE_BALANCED_SAME_SCAFFOLD_PANEL_TOPOLOGY_INADEQUATE`.

This is the only post-G0P balancing correction. Failure closes the current
BioLiP binding-residue-list student route. No further degree threshold,
component pruning, or source-specific rescue is permitted.

All G0P source hashes, allowed fields, firewall rules, clustering constants,
packing order, gates, and compute boundary remain binding.
