# UBSE-A1-v2 source-role certificate correction

Date: 2026-07-30  
Status: frozen identity and inference correction; event access remains locked

## 1. Supersession and immutable artifacts

The original metadata manifests remain immutable:

| Artifact | SHA-256 |
|---|---|
| A1-R metadata | `e6eca57f15975340540d2c2c0afd4e2775f7026ca8af8375c9b3ef3e4299fee9` |
| A1-C metadata | `6e9b2e24246db1b0853f6f6714bb2f6f2cc9e9bfd22fbbb0d665e586e67add24` |
| A1-S roles | `e56020d2e656f6d63fe29d93757e4d77e3adad2d9b8683eeb0e508171258c9ab` |

The first non-deterministic certificate
`ubse_a1v2_source_role_certificate.json` is preserved for audit history and
superseded by:

```text
reports/active/ubse_a1v2_source_role_certificate_v2.json
```

The v2 file is deterministic. Its file SHA-256 is
`f3c581bb168fdd16ef7ef8314957564f4a657fa55ac780a67472e6fd8e343d96`;
its canonical scientific-payload SHA-256 is
`152d0a3e9b46fd343bd2fcfd0abef77c30e2d3af46c16153148f0962c56d1c3e`.

## 2. A1-S closure correction

The original result recorded exact overlap and pre-remediation chemical
maxima but did not bind retained cross-role homology and ECFP4 checks into
`SR-3`. The read-only v2 certificate adds the missing full comparisons:

| Retained roles | Max target containment | Max ECFP4 Tanimoto | Target/ligand conflicts |
|---|---:|---:|---:|
| fit / development | 0.380597015 | 0.469387755 | 0 / 0 |
| fit / legacy audit | 0.351102941 | 0.493150685 | 0 / 0 |
| development / legacy audit | 0.078125000 | 0.410714286 | 0 / 0 |

All equality-as-conflict rules remain unchanged. The old
`validation.max_similarity=0.6984` and
`legacy_audit.max_similarity=0.7768` values describe the pools before whole
panel removal; they are not retained-role maxima.

## 3. Accession amendment

`accession` is optional descriptive metadata. It does not participate in
target identity, role eligibility, closure, hash ordering, selection, or
inherited-model membership. Canonical target identity is the exact normalized
sequence SHA-256 `target_key`; homology membership uses that sequence and the
frozen unique-4-mer containment.

The phrase "empty identifiers are ineligible" applies to required topology
and coordinate identifiers, including target key, sequence, PDB, PubMed,
connectivity, scaffold, receptor/ligand chains, filename serial, and later
`mmcif_auth_seq_id`. It does not apply to optional accession annotation.

Observed missing accession counts are:

- source: 2,146 rows / 1,714 targets;
- A1-R: 12 rows / 8 targets;
- A1-C primary: 20 rows / 20 targets;
- A1-C reserve: 2 rows / 2 targets.

Every one of the 1,035 frozen A1-R/A1-C rows has a target key matching its
exact sequence hash. All other required source identity fields have zero
empty rows. This amendment changes no selected row or manifest hash.

## 4. A1-R dependency identity

Two identities are now explicit:

```text
pair_context_id =
  pdb | receptor_chain | ligand_ccd | ligand_chain | mmcif_auth_seq_id

physical_ligand_id =
  pdb | ligand_ccd | ligand_chain | mmcif_auth_seq_id
```

The 459 A1-R rows have zero duplicate pair-context identities. Seven physical
ligand groups cover 14 rows because distinct receptor chains can share one
physical ligand residue. A1-R, A1-C, and A0 have zero pair-context and
physical-ligand overlap across roles.

Two A1-R units belong to one inference component if they share an exact PDB,
PubMed, or physical-ligand identity; transitive closure is binding. The fixed
component assignment has SHA-256:

```text
46bbff8c6be38ffce699b03faabd86871261ad660f2d667f987ea6680a14e83c
```

It yields:

| Component size | Count |
|---:|---:|
| 1 | 106 |
| 2 | 9 |
| 3 | 7 |
| 4 | 2 |

There are 124 components over 153 target units. The largest has four units;
47 units occur in 18 non-singleton components. The frozen-order
resource-disjoint subset contains 125 units, but it is not substituted for
the frozen 153-unit role.

The component count is descriptive, not a post-result acceptance threshold.
Under perfect within-component correlation, Kish effective sample size is
`98.7722`; for ICC `0.25/0.50/0.75/1.00`, the corresponding effective sizes
are `134.53/120.05/108.38/98.77`. Therefore 153 locator-complete targets, or
the original 128 target-completeness threshold, must not be called 153 or 128
independent observations.

## 5. Frozen inference correction

Before any event value is read:

1. retain target-unit equal weighting for the point estimate;
2. replace target bootstrap by 2,000 component-bootstrap replicates with
   seed 1729;
3. apply one replicate weight jointly to every unit, correct/wrong member,
   extractor, and event type in the same fixed component;
4. report the realized complete target units, components, component-size
   distribution, Kish effective size, 2.5% lower confidence bound, and
   leave-one-largest-component-out sensitivity;
5. freeze the event-specific material-effect and power rule in the extractor
   preregistration. No component-count threshold is adopted after observing
   event coverage.

This correction preserves the claim boundary already stated in the source
preregistration: one extractor across PDB/PubMed supports only
cross-deposition repeatability. Independent reliability additionally requires
the separately frozen second extractor or manual chemistry audit.

## 6. Current decision

The corrected metadata certificate passes identity, A1-S closure, locator
identity, dependency accounting, and exact-sequence target identity. It does
not pass inherited-model membership or event inference:

```text
FREEZE_A1V2_CORRECTED_METADATA_CERTIFICATE_
RETRAIN_P0A_PREREGISTER_COMPONENT_INFERENCE_KEEP_EVENTS_LOCKED
```

Current P0A overlaps all 153 A1-R targets, all 512 A1-C primary targets, and
all 64 reserves. Coordinate bodies, event extraction, affinity, confirmation
scoring, Stage-2, and sealed access remain locked.
