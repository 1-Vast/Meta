# LEXOR L0 Preregistration: Local Corpus Frame Audit

Date: 2026-07-27

## Scope and firewall

This is the first LEXOR stage. It is a metadata-only audit of the already
registered local open-source inventory. The runner may read only this
preregistration and the frozen candidate-document inventory. It must not read a
raw measurement table, affinity value, FORT development/confirmation/sealed
label, API credential, or external endpoint. It must not acquire a document,
call an LLM, or train a model.

The conclusion is limited to this frozen local inventory. It is not a claim
that all reachable open-access literature has been enumerated. A later
web-discovery expansion would require a new inventory hash and a new L0
preregistration; it cannot be appended after observing this result.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`

## Frozen input

Inventory: `manifests/lexor_l0_candidate_documents.v1.json`

<!-- LEXOR_L0_INVENTORY_SHA256: ae271c1168b57b8d6dae70805da1f3c6e1167c8014cb7a8549bf5d651a656411 -->

The inventory contains source-level bibliographic fields, verified license
status, and explicitly reported query-depth metadata only. A missing
post-firewall scaffold-diverse query count remains missing: it may not be
filled from a matrix, an aggregate summary, or an inferred value.

The source metadata references are frozen as follows:

| source metadata artifact | SHA-256 |
| --- | --- |
| `manifests/open_sources.json` | `8244a05949459d5303f91b239aab4ba3b8c664fa9f3780db560176d2bee1a4da` |
| `dataset/public/spd_2023/manifest.json` | `1d2bfc46aff4e0ecdfa23751b4115d8db538eaf8e48d02ff5fcb63ccc98e45f1` |
| `reports/active/reinecke2024_admissibility.json` | `0ea2b04fc4ccd27a00b524a4f50c10efaf7ea60676dc6b5c4b03d95651a42050` |

## Frozen provenance-family definition

Each candidate source begins as one document. Union-find merges two documents
when any one of the following registered links holds:

1. exact normalized DOI, PMID, or patent identifier;
2. author-set Jaccard at least `0.50`;
3. affiliation-set Jaccard at least `0.50`;
4. a common explicitly declared `shared_cell_group_id` in the metadata
   inventory.

Both-empty author or affiliation sets have Jaccard `0.0` and never link a
pair. The runner does not infer shared cells from numerical values. A
provenance family contributes at most one provisional measurement environment
at L0, regardless of how many targets or rows a single campaign reports.

## Frozen gate

An eligible provisional environment requires all of:

1. a document in the inventory with a recorded, verified open license;
2. an explicit metadata count of at least `40` scaffold-diverse query ligands
   after the intended firewall;
3. a distinct provenance family.

The L0 count gate is at least `30` eligible provenance families. The
acquisition list is the complete list of source documents supporting those
families; every listed document must have a verified open license.

Because L0 reads no outcomes, it cannot run the empirical retraining-noise
functions in `research/panel_power.py` or `research/dualcold_power.py` without
violating its own boundary. It instead reports a transparent design proxy,

```text
MDE80_proxy = 2.49 * noise_MAD / sqrt(n_families * q_min)
```

at both Landrum-Riniker maximal-curation noise floors (`0.27` IC50 and `0.45`
Ki log-MAD). This proxy is only a metadata feasibility screen; it does not
replace the PP-DC power calculation required at L4. The high-noise proxy must
be at most `0.03` to pass L0.

All four conditions are conjunctive: enough eligible families, all selected
families with at least 40 declared scaffold-diverse queries, a verified
acquisition list, and high-noise `MDE80_proxy <= 0.03`.

## Verdict rule

* Pass: `LEXOR_L0_CORPUS_FRAME_SUFFICIENT_CONTINUE`. This authorizes only
  creation of an L1 preregistration and blind fixture construction; it does not
  authorize an API call by itself.
* Fail: `LEXOR_L0_CORPUS_FRAME_INSUFFICIENT_STOP`. L1-L5 remain closed for
  this inventory. The decision must state whether the failure is a local-frame
  limitation or a complete discovery result.

## Runner and outputs

Run `conda run -n drug python research/lexor_l0.py`. It writes:

* `reports/active/lexor_l0.json`
* `reports/active/lexor_l0_decision.md`

Both artifacts record the inventory and preregistration hashes, family graph
diagnostics, admission counts, power proxy, firewall state, and no-call state.
