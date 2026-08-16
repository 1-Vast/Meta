# LEXOR L0 Overall Report: Metadata Certification Stopped

Date: 2026-07-27

## Verdict

`LEXOR_L0_CORPUS_FRAME_INSUFFICIENT_STOP`

The expanded metadata audit cannot authorize L1, a live model call, or model
training. This is a failure of the registered **metadata-only certification
frame**. It is not evidence that no recoverable panel exists anywhere in the
open literature.

## What Ran

| phase | result | scientific meaning |
| --- | --- | --- |
| Frozen local L0 | 3 documents, 3 provisional families, 0 eligible environments | The initial local frame had no explicit post-firewall query-depth declaration. |
| L0B transport diagnostics | Zenodo anonymous `size=100` rejected; `size=25` accepted | A transport defect was isolated and corrected without accepting candidate records. |
| L0B11 discovery | 7 public metadata requests, 766 deduplicated documents | Candidate reachability is much larger than the local frame. |
| L0B11 local audit | 762 provisional metadata families, 64 accepted licenses, 0 explicit query depths | No document can satisfy the registered >=40 post-firewall scaffold-diverse query gate. |

The 762 family count is only a provisional metadata linkage result. The
minimal OpenAlex projection deliberately omitted authors, affiliations, raw
cells, and assay protocol fields; it cannot establish true measurement
independence.

## Why the Gate Failed

All 766 discovery components have
`scaffold_diverse_query_ligands = null`. This is intentional: title, record
size, total compounds, citation count, and a matrix-like repository title are
not substitutes for a count after ligand, scaffold, target, and provenance
firewalls. Of the 766 records, 702 also lack a verified license on the
registered whitelist.

The formal L0 gate therefore has zero eligible environments, no acquisition
list, and no evaluable MDE80 proxy. The correct causal statement is:

> Public bibliographic/repository metadata can discover possible sources, but
> it does not expose the raw table topology needed to certify LEXOR's strict
> dual-cold query-depth condition.

This differs from a corpus-exhaustion claim. It also differs from a positive
feasibility result: neither the required independent provenance graph nor the
required query-depth distribution is observable at this stage.

## Transport Correction

Earlier L0B stop reports initially attributed HTTP 400 too narrowly because
the runner did not persist the failing request. That attribution has been
corrected in `lexor_l0b_transport_attribution_correction.md`. A bounded probe
then established the actual cause: Zenodo rejects anonymous `size=100` and
accepts `size=25`. L0B11 completed after this correction. This operational fix
does not change the scientific verdict.

## Consequences

* L1 remains locked. The `.env` credential was not read and no model API was
  called.
* No raw measurement, full text, supplement, FORT label, or confirmation label
  was read in L0B11.
* No pretraining, Mamba comparison, affinity model, or support-posterior route
  is authorized.
* The existing prospective-panel conclusion remains the only established path
  after a source-level observability audit fails. It is not appropriate to
  infer that failure from metadata absence alone.

## Next Exploration: L0C

The next defensible exploration is a separately preregistered **L0C raw-table
observability audit**. It would be deterministic and no-LLM: select only
license-verified sources, hash the acquired bytes, recover table shape and
query counts in code, and apply the existing firewalls before deciding whether
an L1 fixture/API gate is even meaningful. It must not use an LLM to create
values, train a model, or inspect sealed/confirmation labels.

L0C is a protocol amendment, not an authorized continuation of L0B11. It
requires a new source manifest and an explicit decision about which licensed
raw tables may be opened. Until that is registered, the correct state is
`L1_LOCKED`.
