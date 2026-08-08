# Agent handoff

## Start here

Read, in order:

1. `report/CURRENT_RESEARCH_STATUS.md`
2. `task.md`
3. `experiment.md`
4. `research/ssl_b2_structural_observability/PREREG_S5_LOCAL_MECHANISM_OBSERVABILITY.md`
5. `history.md`

## Active work only

Implement S5-A/S5-B in `research/` first.  Do not edit `model/`, production
`scripts/`, frozen theory, CSMO/Band or the probability-law operator.  Do not
read affinity values, DAVIS or recipient labels.

Before training, prove atom/residue mapping, chain compatibility, pseudo-teacher
invariances, exact-residue-to-slot retention and synthetic trainability.  The
1,118 S4 complexes are development-exposed; seal a new RCSB confirmation block
before final scoring.

Use `conda run -n drug` and report progress approximately every 15 minutes for
long-running work.  Keep failed implementations out of production directories;
record terminal results in `history.md` and remove duplicate artifacts after
consolidation.
