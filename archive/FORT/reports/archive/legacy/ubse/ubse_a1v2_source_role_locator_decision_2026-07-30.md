# UBSE-A1-v2 source-role, locator, and advice decision

Date: 2026-07-30  
Decision:
`FREEZE_A1V2_CORRECTED_METADATA_AND_HEAD_CERTIFICATES__RETRAIN_P0A__KEEP_EVENTS_LOCKED`

## 1. Outcome

The label-blind A1-v2 source program now has frozen metadata roles, corrected
locators, and a strict HEAD-only availability certificate:

- `A1-R`: 153 targets, 459 instances, 153 complete three-instance units;
- `A1-S`: 1,260 fit, 59 development, and 81 legacy-audit panels after
  whole-panel chemical-neighbour remediation;
- `A1-C`: 512 primary plus 64 ordered reserves;
- locator resolution: 1,035/1,035 strict unique raw joins, zero missing,
  ambiguous, nonscalar, file-instance duplicate, or pair-context duplicate;
- official URL set: 421 A1-R, 512 A1-C primary, and 64 reserve URLs;
- strict H0V: 997/997 HTTP 200 `application/gzip`, no redirect followed, no
  response-body iteration, and zero actual downloaded bytes.

These results freeze identity and addressability only. No coordinate body,
contact/event label, affinity value, confirmation result, or sealed outcome
was read.

## 2. Corrected source certificate

The original three manifests remain unchanged. The deterministic v2
certificate adds checks omitted from the first result:

| Retained A1-S roles | Max target containment | Max ECFP4 Tanimoto | Conflicts |
|---|---:|---:|---:|
| fit / development | 0.380597015 | 0.469387755 | 0 / 0 |
| fit / legacy audit | 0.351102941 | 0.493150685 | 0 / 0 |
| development / legacy audit | 0.078125000 | 0.410714286 | 0 / 0 |

All exact-resource conflicts are also zero. `SR-3` therefore passes under the
full frozen closure, not merely the earlier exact-overlap check.

Accession is corrected to optional annotation. It is absent in 12 A1-R rows,
20 A1-C primary rows, and 2 reserve rows, but all 1,035 targets match their
exact sequence SHA-256 identity and every required topology/locator field is
complete. No role changes.

The P0A result hash is now bound without decoding its values. Current P0A
supervision overlaps 153/153 A1-R, 512/512 primary A1-C, and 64/64 reserves.
`SR-4` remains a hard failure.

## 3. A1-R dependence correction

The current pair-context identity is unique, but target units are not all
independent:

- 32 PDB groups touch 35 units;
- 41 PubMed groups touch 47 units;
- 7 physical-ligand groups cover 14 rows;
- joint PDB/PubMed/physical-ligand closure gives 124 components;
- component sizes are 106 x 1, 9 x 2, 7 x 3, and 2 x 4;
- component assignment SHA-256 is
  `46bbff8c6be38ffce699b03faabd86871261ad660f2d667f987ea6680a14e83c`;
- worst-case component-correlation Kish effective size is `98.7722`.

Thus 153 locator-complete targets pass the literal completeness gate but are
not 153 independent observations. Event inference must use the frozen 124
components, 2,000 component-bootstrap replicates with seed 1729, joint
weights across unit members/extractors/event types, and largest-component
leave-one-out sensitivity. Event-specific MDE and power remain to be frozen
before extraction; no post-result component threshold is accepted.

## 4. New-advice disposition

The UCE, virtual-cell, and AdaMBind suggestions do not add a new observation
channel:

- retain pair/interface-as-object, masked-event reconstruction, and a
  coupling-residual atlas only as equal-budget auxiliary arms after
  `A1-direct`;
- retain only a conditional bound-state/coupling predictor; bound-only
  BioLiP cannot identify causal unbound-to-bound response;
- retain strict episodic sampling only as grouped domain-generalization
  organization, with an equal-budget grouped-ERM control;
- stop AdaMBind/MAML for primary `k=0`; reserve a separately reported
  chemistry-cold `k=5` task after coupling and affinity gates pass.

The defensible conditional contribution remains:

> source-isolated, cross-deposition and cross-extractor certified typed
> events, distilled into deployment-side coupling-only partial transport
> under a pair-conditioned rank-one exact null and true residue-FG
> checkerboards, followed by fresh strict-dual-cold confirmation and a later
> affinity exact-null increment.

Universal embeddings, masked tokens, pair retrieval, OT, MAML, PEFT, or a
two-stage map-to-affinity pipeline are not individually novel claims.

## 5. Gate status

```text
SR0 identity/firewall                    PASS
SR1 A1-R metadata/locator topology       PASS
SR2 A1-C metadata/locator topology       PASS
SR3 A1-S full retained closure           PASS
SR4 inherited supervised membership      FAIL
SR5 locator half                         PASS
SR5 extractor/FG/assembly half           PENDING
SR6 event freshness/power                PENDING
```

Overall pass remains false:

```text
FREEZE_A1V2_CORRECTED_METADATA_AND_HEAD_CERTIFICATES
RETRAIN_P0A
PREREGISTER_COMPONENT_INFERENCE_AND_EXTRACTORS
KEEP_COORDINATE_BODIES_AND_EVENTS_LOCKED
```

## 6. Next authorized work

1. Freeze P0A-v2 membership closure over all 512 primary and 64 reserves,
   retain a no-P0A primary arm, and retrain the three fixed seeds. The
   independent safe-column estimate is 54,868 rows, 32,769 targets, about
   11.13 million residues, and about 33,216 windows; expected CUDA time is
   5.5-6 hours on the RTX 4060 Laptop GPU.
2. Freeze PLIP version/dependencies, an independent second extractor or
   manual chemistry protocol, FG SMARTS/symmetry, altloc/occupancy,
   water/metal/covalent, model/assembly, and missing-residue rules.
3. Freeze predicted-monomer/cofold training, template, PDB, target-homology,
   and pocket-neighbour membership ledgers.
4. Only after those gates may a coordinate-body hash/download contract be
   considered. Affinity, Stage-2, confirmation scoring, and sealed access
   remain locked.

## 7. Authoritative artifacts

- `reports/active/ubse_a1v2_source_roles_preregistration_2026-07-30.md`
- `reports/active/ubse_a1v2_source_roles.json`
- `reports/active/ubse_a1v2_source_role_certificate_v2.json`
- `reports/active/ubse_a1v2_source_role_certificate_correction_2026-07-30.md`
- `reports/active/ubse_a1v2_locator_completion_preregistration_2026-07-30.md`
- `reports/active/ubse_a1v2_locator_completion.json`
- `reports/active/ubse_a1v2_remote_availability.json`
- `reports/active/ubse_a1v2_head_firewall_verification.json`
- `reports/active/ubse_a1v2_head_firewall_verification_preregistration_2026-07-30.md`
