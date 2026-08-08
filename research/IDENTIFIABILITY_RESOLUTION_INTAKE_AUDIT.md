# Identifiability Resolution Intake Audit

Updated: 2026-08-08

## Decision

The supplied PKIS/KLIFS/KiSSim and section-operator packages are retained under
`research/` as governed development evidence. They are **not** promoted into
`model/` or normal `scripts/`.

The strongest supported update is:

```text
KINASE_PANEL_COMPONENT_IDENTIFIABILITY_OBSERVED_IN_DEVELOPMENT
F6I_REGISTERED_TOTAL_GATE_NOT_ADMISSIBLE
FRESH_ENDPOINT_CONSISTENT_EXTERNAL_ADMISSION_NOT_RUN
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

This is narrower than the title "final resolution". The package establishes a
useful component decomposition on already-consumed kinase activity panels, but
it does not establish transferable Ki/Kd affinity energetics or production
admission.

## Verified intake

- The archive SHA-256 is
  `02d5a3586a90caba6cf6392edc0796d24948b06d49c6a389377d509635d842c1`.
- All 107 archive entries passed a path-traversal check.
- The three standalone files are byte-identical to their archive copies.
- The isolated package suite passes: `39 passed` in the `drug` environment.
- PKIS v1/v2 artifact output hashes match their manifests.
- One test-only defect was repaired: it previously read a non-distributed
  KiSSim CSV by a working-directory-relative path. The repaired test creates a
  temporary file with the same schema. Runtime research code still requires the
  real KiSSim source explicitly.

One historical provenance limitation remains: the F0 manifest binds
`ceiling_probe.py` SHA-256 `2c925e...`, while the supplied current file hashes
to `6ca01a...`. Therefore the included F0 result is retained as historical
evidence but cannot be claimed byte-for-byte rerunnable from the supplied code
alone. Later stage manifests that bind their stage scripts match the supplied
files.

## Evidence accepted

1. The PKIS panels contain strong target-by-ligand interaction structure after
   double centering.
2. A source-atlas biological surface can provide a protein-dependent zero-shot
   interaction component on the consumed kinase panels.
3. Separating that surface from a one-dimensional support-derived location
   statistic removes the earlier multiplicative gauge between protein identity
   and freely adapted section coefficients.
4. The support statistic is bounded and invariant to permutation of the
   support multiset; its dimension is one, hence no greater than `k=5`.
5. The research `law_bridge.py` passes its own mass, barycentre, stochasticity,
   bandedness and permutation-invariance tests.

## Claims not admitted

1. **No affinity-energy identification.** PKIS and Anastassiadis are activity
   panels with endpoint/context differences; these results do not identify a
   physical free-energy component or the governed ChEMBL Ki/Kd estimands.
2. **No untouched external validation.** PKIS2 and Anastassiadis are explicitly
   consumed development panels. The registered F6I total verdict is
   `F6I_COMPONENTS_NOT_ADMISSIBLE` because the Anastassiadis nearest-protein raw
   MSE contrast reverses slightly.
3. **No production-law equivalence.** The research law bridge places a
   categorical law on a seven-point mesh and applies a mean-preserving Markov
   diffusion. The frozen production operator uses a 33-point mesh, valid CDF
   bands, simplex coefficients and `B(z)` band assembly. The notation is
   compatible, but the objects are not type-identical. Passing the research
   tests cannot replace the production operator contract.
4. **No biological `z` admission.** The proposed coordinates `(b, tau, c)` have
   not passed a fresh endpoint-consistent source Gate and sealed transfer Gate.
5. **No full rerun from this repository alone.** The compact artifacts are
   tracked, but upstream Informer/PKIS, KLIFS, KiSSim and Anastassiadis source
   files are not redistributed here.

## Failure triage

| Boundary | Current result | Meaning |
|---|---|---|
| interaction existence on consumed kinase panels | observed | data support a component-level interaction claim |
| zero-shot biological surface | development signal observed | useful research mechanism, not fresh validation |
| one-dimensional task location | numerically and structurally valid | removes the named gauge but is only a location correction |
| total F6I external Gate | not admissible | full component package is not production-qualified |
| research law invariants | tests pass | mathematical compatibility demonstration only |
| frozen production operator equivalence | not established | research bridge must remain isolated |
| fresh endpoint-consistent replication | not run | required next evidence boundary |

## Next admissible step

Freeze the current decomposition and evaluate it once on a new,
endpoint-consistent, target/ligand/document-governed external panel. The test
must compare the correct biological surface against support-free,
protein-null/deranged and wrong-support controls, with the complete Gate frozen
before outcome access. A PASS may authorize a separate production-interface
registration; it must not directly overwrite `model/`, the frozen theory, or
the current CSMO/Band implementation.
