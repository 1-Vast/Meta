# Specification: a genuine physical label seal for the QPSMP lineage

Status: **specification only. Not implemented. Not authorized.**
Requires explicit user authorization before any artifact is created, because
executing it means writing a new `meta_test`-labelled shard.

## 1. What is wrong today, stated exactly

`QPSMPData` provides **logical exclusion after parsing**. Three properties are
routinely conflated and only two hold:

| property | holds? | evidence |
|---|---|---|
| **fail-closed** — a caller cannot admit the sealed split by omission | ✅ | default `include_meta_test=False`; a bare `True` raises |
| **unreachable after construction** — no index, map or `materialize` call can reach a sealed row | ✅ | `tools/tests/test_meta_test_seal_contract.py`, 18 tests |
| **metric-unconsumed** — no sealed value entered a recorded computation | ✅ | verified bit-identical re-run, 105/105 fields |
| **physically isolated** — the labels were never on this process's read path | ❌ | `cells.jsonl.gz` is one all-label artifact; every sealed label is decompressed and parsed on **every** construction |

So the current guarantee is: *the sealed labels pass through process memory on
every load and are discarded before anything can reach them.* A bug in
`_governed_cells`, an exception traceback carrying a frame reference, a memory
dump, or a future refactor that keeps the pre-filter list would all breach it.
None of those has happened. That is not the same as them being impossible.

**The words to use, and the words never to use again:**

* say "logical exclusion after parsing" — accurate;
* say "population sealed, process unsealed" for the R14 incident artifacts;
* say "metric-unconsumed, verified" for the numerical claim;
* **do not** say "physical seal", "physically sealed", or "physically excluded"
  of this lineage until this specification is implemented.

## 2. The pattern already exists in this repository

`scripts/seal_compiled_dataset.py` implements exactly the required design for
the *other* corpus lineage, and it is stronger than anything proposed here:

```python
PERMITTED_SPLITS = frozenset({"source", "metaval"})
SPLIT_MAP = {"train": "source", "val": "metaval", "test": "recipient"}
...
if mapped in PERMITTED_SPLITS:          # recipient rows are simply never written
    ...
manifest = {..., "recipient_label_artifact_emitted": False}
```

`SealedCompiledDataset.__init__` then opens **exactly one** label file, verifies
its hash against the manifest, refuses a manifest that does not assert
`recipient_label_artifact_emitted is False`, and refuses rows whose `split`
disagrees with the mounted view.

This is a working, tested implementation of the property the QPSMP lineage
lacks. This specification is therefore not a new design — it is **applying the
repository's own proven pattern to the corpus that the entire R0-R14 record
actually uses**.

## 3. Required design

### 3.1 Artifact layout

```
dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_sealed/
  manifest.json                  # bindings and hashes; no labels
  governance.jsonl               # cell_id -> split, target, component. NO pK.
  meta_train/cells.jsonl.gz      # labelled
  meta_val/cells.jsonl.gz        # labelled
  meta_test/                     # ABSENT until an authorized sealing run
```

The `meta_test` directory **does not exist** in the development artifact. Its
labels live outside the development tree entirely, at a path no ordinary loader
constructs. The confirmation run is the only process that ever receives it.

### 3.2 Mount contract

```python
GovernedSplitView(directory, visible=("meta_train",))                  # training
GovernedSplitView(directory, visible=("meta_train", "meta_val"))       # development
GovernedSplitView(directory, visible=("meta_test",), authorization=…)  # confirmation
```

* `visible` is a whitelist. A view opens exactly the label files it names and
  no others — this is the property `include_meta_test=False` cannot provide,
  because the sealed rows are in the file it must open anyway.
* Requesting `meta_test` requires a written `authorization` string **and** a
  path that is not inside the development artifact.
* Mounting `meta_train` and `meta_test` together raises, as
  `SealedCompiledDataset` already does for source+metaval.

### 3.3 Manifest bindings

The manifest must bind, and the loader must verify on every mount:

| field | purpose |
|---|---|
| `schema` | `MetaSieve.GovernedSplitSeal.v1` |
| `source_corpus` | `bindingdb_ki_main_v0` |
| `source_manifest_sha256` | pins the corpus the split was cut from |
| `split_assignment_sha256` | the existing `9d8c7289…` assignment hash |
| `cell_id_index_sha256` | hash of the sorted `cell_id` list per split |
| `artifacts` | `{path: sha256}` for every emitted label file |
| `meta_test_label_artifact_emitted` | **must be `false`** in the development artifact |
| `counts` | rows and targets per split, so a truncated shard is detected |

A mount fails closed on any hash mismatch, on an unexpected `schema`, or if
`meta_test_label_artifact_emitted` is not `false` for a development mount.

### 3.4 Preserving the record

The sealing run **must not change any recorded number**. The verification is
the one already used for the R14 repair: rebuild the `meta_val` nested episode
banks under the sealed loader and require the ligand identities and pK values
to be identical, for every support size, to the current loader's. Episode cell
*indices* will differ — they already do under the seal flag — so the comparison
is on biology, not on integers, exactly as
`test_sealing_changes_indices_but_not_the_meta_val_episodes` does now.

### 3.5 Migration

`QPSMPData` keeps its signature and gains a `view` parameter. When a
`GovernedSplitView` is supplied it never opens `cells.jsonl.gz`; when it is
absent, behaviour is unchanged and `seal_record()["isolation"]["level"]` stays
`logical_exclusion_after_parsing`. Nothing in R0-R14 needs re-running.

## 4. The one-time sealing process — requires authorization

**I have not run this and will not without an explicit instruction.** Executing
it creates a `meta_test`-labelled artifact, which the standing rules forbid me
from creating or inspecting unauthorized.

```
1. Read bindingdb_ki_main_v0/cells.jsonl.gz and the frozen double-cold
   assignment.json (hash 9d8c7289…).
2. Partition rows by assigned split. Never hold two splits in one structure.
3. Write meta_train/ and meta_val/ label shards plus governance.jsonl into the
   development artifact.
4. Write the meta_test label shard to an OUT-OF-TREE path supplied by the user;
   record only its sha256 and row count in the development manifest, never its
   path or contents.
5. Emit manifest.json with meta_test_label_artifact_emitted: false.
6. Verify: rebuild the meta_val nested banks under the sealed loader and assert
   byte-equality of ligand ids and pK against the current loader.
7. Emit a SEALING_RECORD.json: input hashes, output hashes, row counts, the
   authorization string, and the verification result.
```

The process reads `meta_test` labels exactly once, in step 2-4, inside a
dedicated trusted process that produces no metrics and imports no model code.

## 5. What this does and does not buy

**Does:** removes the sealed labels from every development process's read path;
makes a breach require an explicit out-of-tree path rather than a wrong default;
makes the guarantee auditable from the manifest rather than from a code review
of one function.

**Does not:** change a single recorded number; retroactively make the R14
artifacts process-sealed (they stay `process_unsealed` with their correction
block permanently); or reduce the need for the authorization gate, which stays.

## 6. Recommendation

Implement **after** the current scientific question is settled, not before.
The seal defect is a governance debt with a verified-zero numerical impact, and
the R-series record is unaffected. Sealing now would consume a session and
create the one artifact the rules are most careful about, to fix a risk that
the fail-closed default and 18 contract tests already reduce to "a future
refactor could reintroduce it".

The correct trigger is **before the confirmation run** — the process that opens
`meta_test` for real should be the first one that has to obtain it from an
out-of-tree path under written authorization. Until then the honest label is
the one now printed by every artifact and by `audit_research_record`:
*logical exclusion after parsing, population sealed, process sealed since
2026-08-16, physically isolated: no.*
