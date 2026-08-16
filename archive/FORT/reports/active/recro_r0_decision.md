# RECRO-DTA Stage R0 — decision

Verdict (frozen gate): **`RECRO_REPLICATION_GRAPH_INSUFFICIENT`** — with a load-bearing qualification:
the cross-environment protein-conditioned **signal is strongly confirmed**; the insufficiency is the
strict entity-disjoint power metric, not absence of signal or a document artifact.

Preregistration `reports/active/recro_preregistration.md`; runner `research/recro_r0.py`; result
`reports/active/recro_r0.json`. Source: licensed raw per-record `chembl37_pKi.jsonl.gz` (ChEMBL
CC BY-SA 3.0, sha256 `f36dd43d…bbd0`), restricted to the registry TRAIN split; no
development/confirmation/sealed label read. Seed 1729.

## What RECRO resolved that AMOB could not

AMOB's O0 (+0.434) was uncertifiable because the staged Harmonic CSVs lacked assay/document IDs, so the
signal could not be shown independent of same-campaign protocol correlation. RECRO R0 uses the raw
ChEMBL per-record extract (per-document values, retained before the registry's median aggregation) and
**directly firewalls by document**: the same target's within-target ligand ordering is compared across
**independent documents**.

## Result (train-only, document-firewalled, 131 independent homology components)

Cross-document reproducibility of within-target ligand ordering:

| arm | component-macro ρ | interpretation |
|---|---|---|
| raw (potency-inclusive) | +0.867 [+0.836,+0.895] | ordering reproduces, dominated by generic potency |
| **residual (potency removed)** | **+0.334 [+0.221,+0.447]** | **target-specific ordering reproduces across documents** |
| wrong-target residual control | −0.207 [−0.347,−0.064] | strongly target-specific (real +0.334 vs wrong −0.207) |
| entity-disjoint residual (n=120) | +0.190 [+0.039,+0.345] | survives no-reused-target/document packing |
| sign agreement, residual doc-pairs | 93.7% | 94% of document-pairs agree in sign after potency removal |

Replication graph: 132 targets with a ≥5-shared-ligand document-pair, 131 homology components, 120
entity-disjoint units (from 248,775 train records, 8,326 documents).

## Why the frozen gate returns INSUFFICIENT

The preregistered power gate required entity-disjoint units ≥25 **and** empirical MDE80 ≤ 0.10. Units
= 120 (pass), but each entity-disjoint unit is a **single document-pair Spearman over ~5 shared
ligands**, so per-unit noise is large (SD 0.86) → **MDE80 = 0.195** (component-level MDE ≈ 0.14; both
> 0.10). The graph can resolve the large biological signal (+0.33) but **cannot certify a small R1
model gain** (historically 0.02–0.05 over B0) at the program's 0.10 resolution. The gate is not
relaxed; `RECRO_REPLICATION_GRAPH_INSUFFICIENT` is reported as frozen.

This is explicitly **not** `RECRO_ORDINAL_SIGNAL_NOT_CROSS_ENVIRONMENT`: the ordinal signal *does*
survive document-independent evaluation (residual +0.334, wrong-target −0.207, 94% sign agreement).

## Scientific significance (and correction to a prior conclusion)

This refines the program's "document-overlap binding constraint" ([[document-overlap-binding-constraint]]).
That conclusion was about the **factorial 2×2 (double-difference)** structure (CROSSDOC 11–13 units;
BM2-PIRR P0 160 entity-disjoint blocks). RECRO shows that the **within-target (single-difference),
potency-residualized** reordering reproduces across independent documents **strongly and over 131
components** — so the target-specific reordering is a **real, document-reproducible biological signal**,
not a protocol artifact. AMOB's +0.434 is therefore corroborated as real by an admissible,
document-firewalled source. What remains hard is: (a) certifying a *small model gain* on the noisy
entity-disjoint graph, and (b) the unchanged crux — **transferring the signal to unseen targets via
protein features** (R1), which four mechanisms + TR-0 + PFSC-0 show is where the program's wall lies.

## Decision

R1 is **not authorized** (R0 power gate failed). No model trained, no threshold relaxed, no
confirmation/sealed access. Per the frozen protocol this stops the mechanism at R0.

**Highest-value continuation (requires re-preregistration):** the R0 insufficiency is a power-*metric*
artifact (single-doc-pair noise), so a re-preregistered R0-b should estimate power from the
**component-macro** endpoint R1 would actually use (many query ligands per held target, as on the Metz/
Reinecke panels where MDE80 was 0.016–0.067), not from single-document-pair reproducibility. If R0-b is
powered, R1 (simultaneous unseen-target + unseen-scaffold, full destructive-control suite from the
RECRO spec) becomes the decisive transfer test — and its expected failure mode is the
protein-featurization wall (real, document-reproducible signal that protein features cannot carry to
unseen targets), i.e. the likely eventual verdict is `RECRO_SIGNAL_REAL_BUT_NO_DUAL_COLD_TRANSFER`
unless the aligned active-site coordinate finally beats shuffled/random protein under strict dual-cold.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing; RECRO read train only).
