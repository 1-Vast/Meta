# Core Task 1 validation funnel — re-adjudicated short funnel (2026-08-19)

Governing plan for the mechanism line (M-line). This replaces any
open-ended ladder for Core Task 1. Frozen as the funnel; each concrete
run still gets its own prereg + SHA before execution. No synthetic run
may be launched from this document before the Q2d-1e verdict and (if 1e
fails) the span-param diagnostic verdict are recorded.

## The only scientific question

Does a protein-conditioned interaction signal exist that:
(a) transfers across cold protein components,
(b) depends on the ligand,
(c) improves affinity or selectivity prediction, and
(d) cannot be explained by ligand-only, protein identity, family
identity, assay context, random permutation, or plain main effects?

## Funnel (promotion-based; stop early)

1. Structure & data-flow tests (minutes, CPU): module interfaces,
   forward/backward/optimizer, k=0 bypass, permutation invariants,
   label-flow isolation, input-dim documentation vs code.
2. Single-seed, single primary difficulty, FOUR-ARM screen:
   correct-protein / ligand-only / family-preserving shuffled protein /
   oracle-or-positive-control.
   - STOP IMMEDIATELY if correct does not clearly beat ligand-only AND
     shuffled protein (no full ladder, no second seed).
3. Only a clear single-seed pass triggers: 3 seeds, full negative set,
   component/target bootstrap, per-surface reporting.
4. Synthetic results only certify the harness's ability to DETECT the
   signal. Core Task 1 itself must ultimately be answered on REAL data:
   same-platform matched WT/variant pairs with identical ligands, or
   ligand-identical protein-pair panels.
5. Functional-inhibition endpoints are never called pK/Ki/Kd. Censoring,
   ATP/Km context, construct identity, mutation numbering and assay
   context are handled per dataset before any training.
6. Real positive controls require a data census and legal-pair count
   FIRST; below the frozen power floor, no training happens.

## Budget discipline (all new candidates, P-line and M-line)

Promotion ladder only:
- CPU / structure gate;
- single seed, short budget (point estimate);
- point estimate fails -> STOP;
- passes -> full budget, single seed;
- passes -> multi-seed;
- passes -> independent evaluation surface (meta_test), once.
No 5-levels x 3-seeds x 8-arms x 8-restarts x 6000-step first runs.
Every launch records in commands.jsonl: estimated GPU-hours, max memory,
arms, seeds, and the expected stopping point.

## Standing stop / authorization states

- Q2d-1e PASS -> synthetic learner qualification PASS; Q2d-2
  representation matrix allowed only via the same funnel (single seed,
  single difficulty, four-arm screen first).
- Q2d-1e FAIL -> exactly one frozen diagnostic (span-param) runs.
- Diagnostic PASS -> recorded conclusion only: "the original learner
  failed mainly from unidentifiable parameterization / null-space
  drift"; span-parameterized form retained as a synthetic candidate;
  NOT promoted beyond the synthetic gate.
- Diagnostic FAIL -> low-rank bilinear synthetic learner family CLOSED;
  no further synthetic successors; Core Task 1 remains UNRESOLVED (not
  a biological falsification).
- Synthetic qualification never blocks the P-line (practical few-shot
  performance line) and never substitutes for real-data mechanism
  evidence.
