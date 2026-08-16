# OpenMut `OMUT-F0` decision

**Verdict:** `OMUT_F0_DAVIS_ROLE_FROZEN__PRESERVE_CONFIRMATION`

**Date:** 2026-07-28.
**Preregistration:** `reports/active/omut_f0_preregistration.md`,
SHA-256 `dad31cbefc8b0be4be4bf46cf160d01195d694b57ac2cc25153ee7e45dc11eb8`.
**Runner:** `research/omut_f0.py`. **Result:** `reports/active/omut_f0.json`.
**Policy record:** `manifests/omut_f0_davis_policy.v1.json`.
**Predecessor:** `OMUT-D0`, anchor
`1953cff4c1d7301c51d1ef934c0c5c913f7c154022ad515c53e416bbce8f82f9`.

## 1. The decision

The one-way Davis role choice was put to the project owner on 2026-07-28 with both
consequences stated. The recorded answer is:

```text
PRESERVE_CONFIRMATION
```

`panel_davis` remains the sealed, single-use, target-conditioned confirmation panel.
Its manifest digest is unchanged at
`84035ad45e0a47a8708520acb9caec638e53533a246cf7b8af68465f85199e69`, with
`single_use=true`, `role=confirmation`, `consumed=false`, `sealed_test_consumed=false`.

**The stated consequence holds: the DAVIS-Complete WT-to-mutant contrast is foreclosed.**
Every overlapping value is excluded, and the WT arm of a WT-mutant pair overlaps by
construction. DAVIS-Complete is not a usable mutation substrate under this policy.

The trade was made with the arithmetic in view: DAVIS-Complete carries about 10-11
canonical WT base proteins and 37 clean single-substitution constructs, giving an
optimistic MDE80 of `0.089-0.085` at paired SD `0.10` against a frozen material threshold
of `0.03`. Retiring the confirmation gate would have bought a substrate already known to
be underpowered for the mechanism question. The sealed panel was kept instead.

Reversal requires a new explicit human decision recorded as a superseding amendment with
its date and authority. A code change is not a reversal, and the runner refuses to
overwrite a conflicting policy record.

## 2. Gates

| gate | result |
| --- | --- |
| `F0_D0_COMPLETE` | pass |
| `F0_POLICY_RECORDED` | pass |
| `F0_POLICY_IRREVERSIBLE` | pass |
| `F0_PANEL_SEALED_STATE` | pass |
| `F0_OVERLAP_KEYS_FROZEN` | pass |
| `F0_DOI_COLLAPSE_FROZEN` | pass |
| `F0_ACCESS_GUARD` | pass |
| `F0_NO_VALUES_READ` | pass |

Zero Parquet rows read, zero row groups opened, zero affinity fields materialised. F0
read exactly two metadata objects: the `panel_davis` manifest, and the `panel_davis`
Parquet **footer** for its column schema. A footer is not a row. The panel's `affinity`
column was read by *name* only, which the test suite asserts explicitly.

A second execution against the written policy manifest reproduced the verdict and
confirmed the recorded policy is identical, exercising the irreversibility check against
a real prior record rather than an empty slot.

## 3. What is now executable rather than prose

The value of this stage is that two constraints stopped being prose and became importable
objects that later stages cannot satisfy by intending to be careful.

**`davis_access_guard(source_id, reads_affinity=...)`** raises `PolicyViolation` on any
Davis-derived affinity read under the recorded policy. It raises rather than returning a
flag, so a caller cannot proceed by ignoring a return value. It permits Davis *metadata*,
and permits non-Davis affinity sources, so it encodes the policy rather than a blanket
ban. The runner demonstrates all four behaviours in an in-run self-test, and the gate
fails if any of them changes.

**`provenance_unit(record)`** collapses a record to its primary document. Records sharing
any document identifier collapse to one unit; the key is order- and case-independent; and
a record carrying only a database name raises, because a database name is never a
provenance unit. This is the arithmetic the whole programme's independent-`n` claims rest
on, and it is now tested rather than asserted.

**Cell-level Davis overlap key**, frozen against the real panel schema:

```text
(UniProt accession of the canonical WT base protein, ligand parent connectivity, endpoint)
```

with `endpoint = pKd` fixed by the panel manifest. All three fields were checked to exist
in the `panel_davis` Parquet footer. The overlap *computation* is deliberately not done
here: it needs DAVIS-Complete construct canonicalisation, which `task.md` assigns to D1.

**Re-export edges**, asserted from D0 rather than rediscovered downstream: BindingDB
imports ChEMBL; MdrDB aggregates PLATINUM, GDSC, DepMap, AIMMS, KinaseMD, TKI and RET;
DAVIS-Complete contains the original Davis panel, which is ChEMBL document
`CHEMBL1908390`, which is `panel_davis`; the Zenodo archive is the same DAVIS-Complete
release. Document keys are recorded for all nine affinity-bearing D0 sources.

## 4. Claim boundary

F0 froze policy. It is **not** evidence about substitution semantics, about
target-specific ligand reordering, about topology, about power, or about strict dual-cold
transfer. The programme verdict is unchanged:

> **3** - current data do not identify the substitution-geometry or tau-teacher mechanism;
> new prospective measurement conditions or a newly recovered, source-resolved public
> substrate are required.

## 5. What this unlocks, and the honest state of the route

`OMUT-F0` is complete. It unlocks **`OMUT-D1` on non-Davis sources only**.

That restriction is the whole story of this stage, and it should not be softened. After
D0 and F0, the sources admissible for a D1 label-free topology audit are:

| source | usable at D1 | why |
| --- | --- | --- |
| BindingDB 202607 Articles (local, verified) | yes, through a registered blind projector | 640 columns frozen, carries `Article DOI`, `PMID`, `Institution`, per-chain UniProt IDs and per-chain sequences |
| ChEMBL variant layer | yes, bounded API projection | 20,150 variant assays, 119,801 variant activities, 0.49% of ChEMBL_37 |
| PLATINUM | **no** | correct schema, no licence; `blocked_rights` |
| DAVIS-Complete | **no** | foreclosed by this policy |
| MdrDB | index only | re-exports PLATINUM; not an independent lineage |
| Binder2030, ProteinGym, MaveDB | **no** | `blocked_unresolved`; MAVE and fitness scores are not affinity labels |

So D1 must establish whether BindingDB and the ChEMBL variant layer, after DOI-level
collapse, contain enough document-independent WT-to-single-substitution rectangles with
shared ligands to make the mechanism question answerable at all. The prior expectation
from the standing record is that they will not: the binding constraint identified across
three independent audits is document-independent factorial overlap, and nothing in D0
supplied a new independent lineage. D1 is worth running precisely because it settles that
with a measurement rather than an expectation, and because both admissible sources are
now frozen, licensed, and bounded.

If D1 confirms the expectation, the two live routes are unchanged and both sit outside
this runner: resolve PLATINUM's usage and redistribution terms, or commission the
prospective A0 panel.
