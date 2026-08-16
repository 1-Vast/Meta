# LEXOR L0C-R0 Reinecke Calibration Preregistration

## Scope

`L0C-R0` is a one-run, local, deterministic calibration of raw-table
observability. It replays only the already development-seen Reinecke S1, S2,
and S3 supplements. It is not source discovery, a new-source acquisition,
an L0 pass attempt, an L1 fixture, an entity-resolution pass, or a model
experiment.

## Frozen Input

<!-- LEXOR_L0CR0_SOURCE_MANIFEST_SHA256: c406208350e83c05b6ec372e627e7546f6385b5ae39f925f03f6a14e1acc972d -->

The manifest is calibration-only because the historical acquisition did not
retain exact per-file URLs. The runner must reject any manifest which changes
that role or whose bytes differ from the frozen digest before opening a table.

## Registered Checks

1. Verify the open-license entries, file media types, byte lengths, and SHA-256
   values before parsing a workbook.
2. Verify every registered sheet locator and required header.
3. Rebuild the existing deterministic parent/scaffold and accession-candidate
   linkage, then roundtrip every accepted S2 cell to its original coordinates.
4. Build a per-cell ledger containing parent, scaffold, accession candidate,
   construct state, endpoint, unit, qualifier, provenance family, and source
   coordinates.
5. Rebuild the ledger twice and require identical order-sensitive SHA-256
   digests.

## Decision Rule

The calibration passes mechanically only when all file and table checks pass,
every accepted cell roundtrips, the two ledgers are identical, and no external
network, LLM/API call, model training, confirmation-label read, or sealed-test
read occurs.

The source remains noncountable when its construct is not reported in the
local source table. Its one provenance family is additionally insufficient for
the L0 portfolio threshold. A mechanical pass therefore emits only
`L0CR0_MECHANICAL_OBSERVABILITY_PASS__SOURCE_NONCOUNTABLE` and authorizes no
scientific-stage transition.
