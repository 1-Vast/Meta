# LEXOR L0C Raw-Table Observability Design

Date: 2026-07-27

## Status

`DESIGN_REGISTERED__L0C_R0_CALIBRATION_COMPLETED__NEW_SOURCE_ACQUISITION_NOT_STARTED`

L0C is the exploration opened after the L0B11 metadata-certification stop. It
is not a model experiment, an LLM experiment, or a fourth supervision
candidate. It exists to answer one narrower question that L0B cannot observe:
whether licensed raw source tables can expose the post-firewall query-depth and
provenance topology required by LEXOR L0.

## Why an Amendment Is Needed

L0B11 found 766 candidate documents and 64 records with an accepted repository
license, but no metadata record exposes a query count after parent, scaffold,
target, and provenance firewalls. Computing that number requires source-table
content and deterministic chemical/entity normalization. Metadata-only L0
therefore cannot distinguish a missing corpus from a missing field.

## Proposed Scope

Before any source byte is opened, freeze a source manifest containing for each
candidate:

* persistent record/document identifier and provenance-family hypothesis;
* explicit accepted license and official license URL;
* exact source file URL, media type, SHA-256 after acquisition, and retrieval
  timestamp;
* a declaration that the file is a table/supplement rather than a narrative
  article alone.

Only manifest entries passing the existing open-license whitelist may be read.
The audit then uses deterministic code only to:

1. recover table dimensions, headers, cell locators, and verbatim values;
2. resolve a ligand to a parent/scaffold only when an auditable structure or
   unambiguous deterministic identifier is present;
3. resolve a target to an exact accession and construct state only when the
   source supplies sufficient evidence;
4. construct source/protocol provenance families without treating repository
   record IDs, plates, or technical replicates as independent environments;
5. count query ligands only after all registered firewalls have been applied.

## Decision Rules

L0C may make no LLM/API call and may not train a model. It stops a source when
any required license, source hash, table locator, structure/parent mapping,
target/construct mapping, endpoint, or provenance evidence is missing. It
cannot replace absent fields with title text, LLM inference, a database mirror,
or a pooled activity value.

The resulting source-level inventory would re-run the existing L0 gate without
changing its 30-family, >=40 query-depth, whitelist, or MDE80 thresholds. A
failure after this observability audit would support the prospective factorial
panel conclusion. A pass would authorize only a separately frozen L1 blind
fixture run, not training or downstream affinity evaluation.

## Authorization Boundary

No candidate source has been selected or acquired under L0C. Proceeding from
this design requires an explicit user-approved source manifest because it
expands the work from public metadata to raw licensed table bytes. The `.env`
credential remains out of scope until a valid L0/L0C pass and a frozen L1
fixture chain exist.

## L0C-R0 Calibration Result

On 2026-07-27, the bounded Reinecke `L0C-R0` calibration replayed only the
already development-seen local S1/S2/S3 supplements. Its frozen calibration
manifest hash was
`c406208350e83c05b6ec372e627e7546f6385b5ae39f925f03f6a14e1acc972d`.

The deterministic ledger recovered 9,346 evidence-bound S2 cells, with zero
source-coordinate roundtrip errors and identical hashes across two complete
replays:
`3e0e96ece56dd5f715f1f88d7a151ff51d0121194279129f3d85f6c2250e712a`.
This validates the L0C mechanical observability path: frozen file hashes,
license entries, table locators, headers, parent/scaffold evidence, source
coordinates, endpoint/unit/qualifier, and provenance can be placed in a stable
ledger.

It intentionally did not pass source eligibility. The local source table
contains gene labels but no source-verified accession or construct; the frozen
local cache yields an accession candidate only. All 9,346 ledger records are
therefore noncountable, query depth after all firewalls is `null`, and the
source has only one provenance family. The registered result is
`L0CR0_MECHANICAL_OBSERVABILITY_PASS__SOURCE_NONCOUNTABLE` in
`reports/active/lexor_l0cr0.json`.

This is neither a new discovery source nor a blind fixture. It did not download
or open any new source bytes, call a network/API, train a model, consume a
sealed test, or authorize L0, L1, or downstream training. A prospective source
still needs a complete acquisition manifest with exact file URL and post-
acquisition receipt before its raw bytes may be inspected.
