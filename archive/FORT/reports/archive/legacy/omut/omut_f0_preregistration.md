# OpenMut `OMUT-F0` preregistration

**Frozen:** 2026-07-28, before `research/omut_f0.py` existed and before any `OMUT-F0`
artefact was written.
**Environment:** `D:\anaconda\envs\drug\python.exe` (Python 3.11.15, torch 2.6.0+cu124,
CUDA 12.4, RTX 4060 Laptop GPU).
**Predecessor:** `OMUT-D0`, verdict `OMUT_D0_SOURCE_FREEZE_COMPLETE`, anchor
`substantive_registry_sha256`
`1953cff4c1d7301c51d1ef934c0c5c913f7c154022ad515c53e416bbce8f82f9`.
**Authority:** `task.md`, Route A; `OMUT-F0` unlock condition "D0 complete and one
irreversible Davis choice recorded".

## 1. The recorded human decision

`OMUT-F0` is not a computation. Its input is a one-way choice that a runner may record
and enforce but must never make. The choice was put to the project owner on 2026-07-28
with both consequences stated, and the answer was:

```text
PRESERVE_CONFIRMATION
```

`panel_davis` remains the sealed single-use target-conditioned confirmation panel. Every
DAVIS-Complete value that overlaps it is excluded. The direct consequence, stated before
the choice and restated here, is that **the DAVIS-Complete WT-to-mutant contrast is
foreclosed**: the WT arm of essentially every candidate pair is an overlapping value.

The supporting arithmetic, from the D0-era record and not re-derived here: DAVIS-Complete
carries about 10-11 canonical WT base proteins and 37 clean single-substitution
constructs, giving an optimistic MDE80 of `0.089-0.085` at paired SD `0.10` against a
frozen material threshold of `0.03`. Consuming the confirmation gate would therefore have
bought a substrate already known to be underpowered for the mechanism question.

## 2. Question and claim boundary

F0 asks one question:

> Can the Davis role policy and the DOI-level source-circularity firewall be frozen as
> *executable* constraints, binding on every later stage, without reading a value?

F0 is a policy freeze. It is not evidence about substitution semantics, ligand
reordering, topology, power, or dual-cold transfer. A pass authorises entry into
`OMUT-D1` and nothing else, and it authorises D1 **only on non-Davis sources**.

F0 does not authorise: reading any affinity, relation, or censor value; reading any row
of `panel_davis`; downloading `davis_complete.tab`; construct canonicalisation or
mutation-topology derivation, which are D1 deliverables; or any coordinate, estimator, or
predictive work.

## 3. The no-value boundary

F0 reads exactly two objects, both metadata:

| object | mode | permitted |
| --- | --- | --- |
| `dataset/public/chembl_37/processed/panel_davis/manifest.json` | whole-file digest and JSON keys | yes; it holds no affinity |
| `dataset/public/chembl_37/processed/panel_davis/registry.parquet` | **footer only**: column names, types, row count | yes; a footer is not a row |

Zero Parquet rows may be read. The runner must not construct a reader that materialises a
row group, and must assert `rows_read == 0`. `davis_complete.tab` remains name-guarded by
the D0 transport layer.

The overlap **key definition** is frozen at F0. The overlap **computation** is not: it
requires DAVIS-Complete construct canonicalisation, which `task.md` assigns to D1.

## 4. Gates

All eight are evaluated; all eight must pass.

| gate | condition |
| --- | --- |
| `F0_D0_COMPLETE` | `reports/active/omut_d0.json` exists, its verdict is `OMUT_D0_SOURCE_FREEZE_COMPLETE`, its preregistration hash matches, and its substantive anchor equals the value recorded in section 1 |
| `F0_POLICY_RECORDED` | exactly one policy from `{preserve_confirmation, retire_confirmation}`, with the decision date and the recording authority |
| `F0_POLICY_IRREVERSIBLE` | if a policy manifest already exists, its recorded policy is identical; a differing prior policy aborts the run |
| `F0_PANEL_SEALED_STATE` | the `panel_davis` manifest digest equals `84035ad45e0a47a8708520acb9caec638e53533a246cf7b8af68465f85199e69` and it reports `single_use=true`, `role=confirmation`, `consumed=false`, `sealed_test_consumed=false` |
| `F0_OVERLAP_KEYS_FROZEN` | every field named in the cell-level overlap key exists in the `panel_davis` Parquet footer schema, and the DAVIS-Complete side names a file frozen in the D0 registry |
| `F0_DOI_COLLAPSE_FROZEN` | every D0 source with `acquisition.affinity_bearing == true` carries a recorded document-key field list, each field drawn from that source's D0-frozen schema where the source has one |
| `F0_ACCESS_GUARD` | the guard refuses DAVIS-Complete affinity access under the recorded policy, and permits a non-Davis source, both demonstrated by an in-run self-test |
| `F0_NO_VALUES_READ` | zero Parquet rows read, zero affinity fields materialised, zero firewall violations |

## 5. The firewall being frozen

**Cell-level Davis overlap key.** A DAVIS-Complete cell overlaps `panel_davis` when it
matches on the tuple

```text
(UniProt accession of the canonical WT base protein, ligand parent connectivity, endpoint)
```

with `endpoint = pKd` fixed by the panel manifest. Under `PRESERVE_CONFIRMATION` every
matching cell is excluded, and because the WT arm of a WT-mutant pair matches by
construction, the contrast is unavailable.

**DOI-level source circularity.** One primary document is one provenance unit, however
many databases re-export it. A database name is never a provenance unit. The collapse key
is the document identifier, resolved per source from the D0-frozen schemas, and it is
applied transitively: if two records share any document identifier they are one unit.

**Known re-export edges**, frozen from D0 and asserted rather than rediscovered:
BindingDB re-exports ChEMBL records; MdrDB re-exports PLATINUM, GDSC, DepMap, AIMMS,
KinaseMD, TKI, and RET; DAVIS-Complete contains the original Davis panel, which is ChEMBL
document `CHEMBL1908390`, which is `panel_davis`.

## 6. Verdicts

```text
OMUT_F0_DAVIS_ROLE_FROZEN__PRESERVE_CONFIRMATION   all eight gates pass
OMUT_F0_INCOMPLETE_STOP                            any gate fails
OMUT_F0_POLICY_CONFLICT_ABORT                      a differing prior policy exists
```

## 7. Stop rules

- The policy is one-way. Reversing it requires a new explicit human decision, recorded as
  a superseding amendment with its date and authority; the runner cannot reverse it, and
  silent reassignment is prohibited.
- A pass does not make DAVIS-Complete usable. Under the recorded policy it stays
  foreclosed for the WT-mutant contrast, and D1 must proceed on non-Davis sources.
- F0 failure is not rescued by weakening the overlap key, by narrowing the collapse rule,
  or by reclassifying a sealed asset.
- No GPU computation is performed at F0.

## 8. Deliverables

- `research/omut_f0.py` — the policy freeze, the enforceable guard, and the gates.
- `tests/test_omut_f0.py` — guard, key, collapse, gate, and irreversibility tests.
- `reports/active/omut_f0.json` — machine-readable result.
- `manifests/omut_f0_davis_policy.v1.json` — the one-way policy record.
- `reports/active/omut_f0_decision.md` — verdict, consequences, and what D1 may now use.
