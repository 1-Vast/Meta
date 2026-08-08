# Research Archive Boundary

The executable research branches were consolidated on 2026-08-08.

Terminal negative, mixed, synthetic-only and consumed-development code was
removed from the working tree after its scientific conclusions, metrics,
provenance limits and deletion scope were recorded in `history.md` and
`EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`. The exact deleted files remain
recoverable from Git history through commit `8b7789e`.

The only reusable algebra admitted from the F6I branch is the internal,
gauge-separated component statistic in `model/component_statistic.py`, with a
label-safe JSONL interface in `scripts/evaluate_component_statistic.py`.
Neither file is connected to the production operator or exported by `model`.

New experiments must be created here only after a new preregistration defines
the dataset release, endpoint, split closure, controls, Gate and promotion
boundary. Downloading a dataset does not itself authorize affinity-label access
or production integration.
