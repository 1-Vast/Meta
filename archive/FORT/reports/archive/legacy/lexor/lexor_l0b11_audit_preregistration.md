# LEXOR L0B11 Corpus-Frame Audit Preregistration

Date: 2026-07-27

## Frozen input

Candidate inventory: `manifests/lexor_l0b11_candidate_documents.v1.json`

<!-- LEXOR_L0_INVENTORY_SHA256: 6dcd6e10218051e76c7bfd7dd403fe48ffe26e94fadb8da1325d88d927bd8807 -->

Discovery artifact: `manifests/lexor_l0b11_open_metadata_discovery.v1.json`

<!-- LEXOR_L0B_DISCOVERY_SHA256: 96c016e1601cc1ef66565fe4969a4bf151b8652175accb607f3a7f1cc52c13ec -->

The inventory is closed at 766 deduplicated candidate documents. It must not be
amended after this audit begins.

## Registered audit

Run `research/lexor_l0.py` against the frozen inventory using its registered
metadata-only provenance union-find, open-license whitelist, 30-family gate,
explicit post-firewall scaffold-diverse query depth >=40 gate, and high-noise
MDE80 design proxy <=0.03. No count may be inferred from a title, abstract,
record size, total compounds, or citation count.

The discovery artifact explicitly records zero documents with an allowed
post-firewall scaffold-diverse query-ligand count. This is a pre-run
observation, not an audit verdict; the runner determines the formal result.

## Firewall state

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`
* `fort_labels_read=false`
