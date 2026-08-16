# Pretraining-anchor round: three agent-proposed candidates stopped

Date: 2026-07-27

## Overall judgment

All three agent-proposed data/supervision candidates in the reopened round have
reached terminal F0/G0 stops. None supplies an identified, provenance-defensible
target-conditioned supervision signal for strict dual-cold affinity transfer.
This is an overall report before reopening exploration, as required. It does
not authorize a fourth candidate, a model-training rescue, or a Mamba benchmark.

| candidate | terminal verdict | binding failure |
| --- | --- | --- |
| MMP-X (1/3) | `MMPX_F0_ACCESSION_FIREWALL_INSUFFICIENT_STOP` | Exact-accession isolation reduces 39 nominal cross-source transformation-family repetitions to 8 units across 5 families. |
| TCOPA (2/3) | `TCOPA_G0_INSUFFICIENT_TARGET_CONTRAST_SUPPORT_STOP` | ToxCast has target contrasts, but a 18.4% dominant scaffold and insufficient target-fold depth violate the frozen dual-cold pretext design. |
| Papyrus F0-P (3/3) | `PAPYRUS_F0_RAW_PROVENANCE_INSUFFICIENT_STOP` | The aggregate release has zero document-replicated parent--target cells after raw-provenance filtering; the full table has one row per parent--target pair. |

The user-supplied SAFSA direction is not included in the three-candidate quota.
It independently stopped as
`SAFSA_G0_FAMILY_SELECTIVITY_NOT_IDENTIFIABLE_STOP`: its own-family score was
only +0.0047 above ligand promiscuity and had a negative family-bootstrap lower
bound.

## What the failures say together

The failures are different manifestations of one measurement-design problem:

1. A chemically plausible unit is insufficient when its exact proteins recur
   across sources (MMP-X).
2. A dense systematic binary panel is insufficient when its scaffold and
   target-fold allocation cannot support the registered dual-cold contrast
   (TCOPA).
3. A large quantitative aggregate is insufficient when it has discarded the
   raw repeated document records needed to establish independent replication
   (Papyrus).

The result does not say that local interaction primitives, target contrasts, or
Mamba blocks are false in principle. It says that these particular datasets do
not identify a transferable target-conditioned effect under the required
firewalls. More GPU memory or a larger encoder cannot create missing independent
measurements.

## Consequences for the model program

- Do not train SAFSA, MMP-X, TCOPA, or Papyrus pretraining models.
- Do not enter F1, F2, F3, or F4. In particular, do not use the successful
  Mamba E0 component implementation as a reason to bypass the data gate.
- Do not replace document/provenance isolation with random splits, aggregated
  semicolon fields, pooled Ki/Kd/IC50 labels, or a ligand-only score presented
  as protein-conditioned transfer.
- Preserve the ChEMBL confirmation quarantine: current-run labels remain unread
  and `sealed_test_consumed=false`.

## Exploration reopened, but no fourth experiment is started

The candidate limit is exhausted. Exploration is now limited to source and
measurement-design discovery, not a new performance route. A future candidate
must first demonstrate all of the following from raw records before it consumes
another experiment slot:

1. Per-measurement quantitative values with exact target accession, assay and
   document/provenance-family identifiers, rather than an aggregate pair table.
2. At least two genuinely independent provenance families for enough identical
   parent--target cells to construct a powered replication graph.
3. Multi-family target coverage, scaffold-diverse ligand support, and target
   folds that satisfy the dual-cold query-depth requirement after all firewalls.
4. A pre-registered target-correct versus ligand-only, target-shuffle, random
   target, and matched wrong-target mechanism test.

The next permissible work is therefore an audit of raw-data sources or a
prospective factorial measurement design. It is not an architecture change.

