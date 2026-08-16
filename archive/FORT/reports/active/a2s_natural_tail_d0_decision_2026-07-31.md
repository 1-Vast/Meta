# A2S-DTA Natural-Tail Gate D0 Decision

Date: 2026-07-31  
Final status: `DATA_NOT_READY`  
Training authorization: `false`

## Executed command

```powershell
D:\anaconda\envs\drug\python.exe main.py natural-tail-audit --out dataset\processed\a2s_natural_tail_d0.v1.json
```

The command completed in 37.5 seconds. It read only the canonical global row
index and the registry fields `target`, `conn`, `endpoint`, `scaffold`,
`assays`, `docs`, `hcluster`, and `dual_cold_split`. It did not load affinity,
replicate SD, development, confirmation, or sealed outcomes, and no model fit
or CUDA training ran.

Input hashes:

- registry: `0e754f73f5d75913d61791d6ccd08e05662cd8015fc608ba370d4ee303e6b784`;
- document metadata: `5c920d4b33b88389c5331879ce1b620fe9b70e5006c3b3b845b1ad7102734109`.

## Feasibility audit settings

The executed source rule was `n_eff >= 100` and the recipient-candidate rule was
`n_eff < 30`. These are the existing operational candidate settings, not a
newly validated empirical tail partition. Support budgets are nested
`k={1,3,5}`. The strict diagnostic attempted support-before-query by ChEMBL
release, a held-out document `src_id`, and support/query closure in parent,
scaffold, document, assay, and source family.

`chembl_release` is only an ingestion-order proxy. The local registry and
document metadata do not contain publication or measurement dates, so the
current source cannot establish genuine temporal precedence. This limitation
is a gate failure, not a missing-value imputation opportunity.

## Cutflow

| Metadata-only stage | pKi recipient targets |
| --- | ---: |
| Candidate `n_eff < 30` | 193 |
| Raw `k=5` plus 10-query count envelope | 85 |
| Optimistic parent/document/assay closure, five arbitrary support combinations | 40 |
| Above plus five distinct support scaffolds | 34 |
| ChEMBL-release temporal envelope | 22 |
| Temporal plus parent/scaffold/document/assay closure | 11 |
| Held-out `src_id` source-family envelope | 1 |
| Source-family full closure | 0 |
| Strict time/source/parent/scaffold/document/assay roster | **0** |

The 40-target row is intentionally optimistic: it does not require time,
source-family separation, common query rows, query-scaffold closure, or true
measurement dates. It is already below the frozen minimum of 50, so no
reasonable stricter interpretation can reverse the D0 decision.

## Post-run adversarial claim boundary

The emitted strict roster is not an admission-grade reusable roster. Independent
review found that it adds an unregistered scaffold-cold constraint, does not
guarantee five distinct support sets, treats ChEMBL ingestion release and
document `src_id` as stronger time/lineage variables than they are, and writes
constructed zero overlap rather than recomputing every emitted overlap. The
`n_eff >= 100` and `<30` role thresholds also remain provisional in `task.md`.

These issues cannot turn this run into a PASS. The DATA stop rests only on two
weaker facts: the local source lacks true time/lineage metadata, and under the
executed candidate definition an exhaustive parent/document/assay-closed
upper bound requiring five distinct feasible k=5 support combinations is 40,
below the prespecified target-count floor of 50. The stricter zero-recipient
diagnostic is reported for sensitivity only and is not used as proof that a
valid natural-tail roster was constructed.

## Dependency and power

Candidate recipients share source family `src_id=1` at 92.75% target coverage.
Joining recipients that share a homology component or any source family puts
all 193 candidates in one dependency component. Candidate/source overlap is
also substantial: 12 recipients share source homology, 93 share a parent, 81
share a document, and all 193 share at least one source family with the head
pool. These are reported diagnostics; they are not converted into independent
replicates.

The strict roster contains zero recipients and zero independent components.
Consequently MDE80 is undefined, and the empty-roster SHA-256 is
`4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945`.
Rows, queries, draws, and seeds cannot substitute for the absent component
count.

The executed JSON artifact SHA-256 is
`1c2fb0a8501bc69333931e7ce50b9960b7b2d3cb58c0382d1a23c4434d7d5c27`.

## Decision

Gate D0 is a DATA STOP. Gate S0, the one-seed A2S-MAP kill test, MAML,
AdaMBind, formal baselines, ablations, five-seed training, pKd replication,
and encoder expansion remain blocked. More epochs, relaxed closure, a lower
recipient threshold, pseudo-tail results, or pKd cannot rescue this failure.

The only next action is to acquire and freeze a provenance-rich pKi
natural-tail source with true date/lineage metadata. Before rerunning D0,
preregister the empirical resource thresholds and correct the draw, closure,
overlap, component, and certificate contract. At least 50 recipients must then
survive that common-roster audit before any model score is read.
