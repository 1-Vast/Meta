# B1 preregistration (DRAFT — freezes only after Q0-Q2 pass)

Stage: same-study wild-type -> mutant positive control on the Duong-Ly 2016
panel (98 rows x 183 compounds, single platform, duplicated measurements).
This document is a design draft; the frozen SHA is recorded in
B1_PREREGISTRATION_SHA256.txt when and only when Q0, Q1 and Q2 have all
passed their frozen gates. Until then this file has no legal force.

## Estimand

Selectivity contrast of percent-remaining-activity between a mutation row m
and its wild-type parent p, centred over the 21 wild-type parent rows:

  Delta(m) = y(m) - y(p),   CSC(m | m,p) = Delta(m) - mean_{p' in WT}(Delta(p'))

where y are logit-transformed percent-remaining-activity values, censored
rows excluded pairwise, and the centring set is the train wild-type parents
only (frozen I6 rule).

## Frozen decisions (to be locked at B1 freeze time)

- Parents: the 21 wild-type rows from the typed pair table (CAS row excluded).
- Mutants: the 65 admitted single-point pairs; quarantined rows are reported
  separately and never admitted.
- Endpoint ladder: (1) logit(% remaining) OLS; (2) interval-censored
  (known bounds); (3) % remaining raw. Censoring ladder reported per row.
- Grouping: same-parent mutant blocks; cluster bootstrap over parents
  (frozen 2.5% percentile, 1000 draws, SHA-256 seed).
- Reference leak check: reference statistics computed on train parents only;
  eval parents' labels never enter the centring term.
- Multiplicity: per-mutant CSC with 95% cluster-bootstrap CI; a mutant is
  positive only if its CSC lower bound > 0 AND its uncensored-cell count
  exceeds the per-row minimum (frozen at 30).

## Gate

B1 passes if the platform shows a positive CSC for the designed positive
controls: BRAF V600E vs BRAF (expected selective for BRAF/MEK tool
compounds only — verified positive-control compound set frozen from the
public panel metadata), plus >= 1 additional designed positive pair with
its expected tool compound. Failure of every designed positive control
falsifies the platform for selectivity work and blocks B2/C/D.

## Deliverables

B1_POSITIVE_CONTROL.json: per-pair CSC, CI, censoring ladder, parent
blocks, compound fingerprints (name-level only, no structural fingerprints
committed), bootstrap seeds, and the PASS/FAIL verdict.
