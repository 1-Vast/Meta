"""One-time governed builder: separately hashed per-split label artifacts.

This is the administrative process specified in
`tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md`. It is the **only**
process permitted to read the all-label corpus together with the sealed split
assignment, and it produces no metric. It exists so that ordinary development
loaders have **no filesystem dependency** on the `meta_test` label artifact.

What it does *not* do, by construction:

* it does not change any target/component assignment. The split is read from
  the frozen `bindingdb_ki_double_cold_v1` assignment whose sha256 is verified
  before anything is written;
* it does not compute, print, or store any statistic of `meta_test` labels
  beyond a row count and a file hash;
* it does not place the `meta_test` label artifact inside the development
  directory. That path is supplied separately and recorded only by hash.

Layout produced:

```
<output>/                       the development surface
  manifest.json                 bindings, counts, hashes. NO labels.
  governance.jsonl              cell_id -> split/target/component. NO pK.
  meta_train/cells.jsonl.gz     labelled
  meta_val/cells.jsonl.gz       labelled
                                (no meta_test directory)
<sealed-output>/               supplied separately, outside the dev tree
  meta_test/cells.jsonl.gz      labelled
  SEALING_RECORD.json
```

`manifest.json` records `meta_test_label_artifact_emitted: false` for the
development surface. `GovernedSplitView` refuses to mount a development
manifest that says otherwise, and refuses to mount `meta_test` at all without a
written authorization plus an explicit out-of-tree path.

Run (requires --authorization; writes nothing without it):

```
python -m scripts.build_governed_split_views \
  --corpus dataset/processed/meta_fewshot/bindingdb_ki_main_v0 \
  --split-directory dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1 \
  --output dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views \
  --sealed-output <path outside the development tree> \
  --authorization "<written reason>"
```
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SCHEMA = "MetaSieve.GovernedSplitSeal.v1"
DEVELOPMENT_SPLITS = ("meta_train", "meta_val")
SEALED_SPLIT = "meta_test"
ALL_SPLITS = (*DEVELOPMENT_SPLITS, SEALED_SPLIT)

# A written reason, not a token. The 2026-08-16 incident was a seal opened by a
# bare flag, so "authorized" has to cost a sentence a reviewer can read back.
AUTHORIZATION_MIN_CHARS = 16

REQUIRED_MANIFEST_FIELDS = (
    "artifacts", "cell_id_index_sha256", "counts", "governance_jsonl_sha256",
    "source_manifest_sha256", "split_assignment_sha256",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_authorization(authorization: object) -> str:
    """Fail closed on a missing or malformed meta_test authorization.

    Three distinct rejections, because they are three distinct mistakes:
    a caller that forgot the argument, a caller that passed something that is
    not written text, and a caller that passed a placeholder.
    """
    if authorization is None:
        raise ValueError(
            "mounting meta_test requires a written authorization; none given")
    if not isinstance(authorization, str):
        raise ValueError(
            "malformed meta_test authorization: expected written text, got "
            f"{type(authorization).__name__}")
    text = authorization.strip()
    if not text:
        raise ValueError(
            "mounting meta_test requires a written authorization; the value "
            "given is blank")
    if len(text) < AUTHORIZATION_MIN_CHARS:
        raise ValueError(
            "malformed meta_test authorization: a written reason of at least "
            f"{AUTHORIZATION_MIN_CHARS} characters is required, got "
            f"{len(text)}")
    return text


def read_governance(directory: Path) -> list[dict]:
    """Identity-only ordering record: cell_id -> split/target/component, no pK."""
    rows = [json.loads(line) for line
            in (Path(directory) / "governance.jsonl").read_text(
                encoding="utf-8").splitlines() if line.strip()]
    for row in rows:
        if set(row) != {"cell_id", "split", "target_id", "protein_group_40"}:
            raise ValueError("governance.jsonl carries unexpected fields; "
                             "it must contain no labels")
    return rows


def write_cells(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 so the artifact hash depends on content only and a rebuild is
    # verifiable against the recorded digest.
    with gzip.GzipFile(path, "wb", mtime=0) as raw:
        for row in rows:
            raw.write((json.dumps(row, sort_keys=True) + "\n").encode("utf-8"))


class GovernedSplitView:
    """Mount exactly the label artifacts named, and verify their bindings.

    A development mount cannot reach `meta_test`: the artifact is not in the
    directory, and the manifest asserts it was never emitted there. This is the
    property `QPSMPData(include_meta_test=False)` cannot provide, because the
    sealed rows live in the single file it must open anyway.
    """

    def __init__(self, directory: Path, visible: tuple[str, ...],
                 authorization: str | None = None,
                 sealed_directory: Path | None = None):
        self.directory = Path(directory).resolve()
        self.manifest = self._read_manifest(self.directory)
        visible = tuple(visible)
        if not visible:
            raise ValueError("a view must make at least one split visible")
        unknown = [split for split in visible if split not in ALL_SPLITS]
        if unknown:
            raise ValueError(f"unknown split(s) requested: {unknown}")
        if SEALED_SPLIT in visible:
            authorization = check_authorization(authorization)
            if sealed_directory is None:
                raise ValueError(
                    "mounting meta_test requires an explicit out-of-tree "
                    "sealed_directory; it is not part of the development "
                    "surface")
            sealed_directory = Path(sealed_directory).resolve()
            if (sealed_directory == self.directory
                    or self.directory in sealed_directory.parents):
                raise ValueError(
                    "the sealed_directory must be outside the development "
                    "surface; a path inside it is not an out-of-tree artifact")
            if set(visible) & set(DEVELOPMENT_SPLITS):
                raise ValueError(
                    "refusing to mount meta_test alongside a development "
                    "split in one process")
        elif authorization is not None:
            raise ValueError(
                "authorization supplied without requesting meta_test; "
                "refusing an ambiguous mount")
        if (SEALED_SPLIT not in visible
                and self.manifest.get("meta_test_label_artifact_emitted") is not False):
            raise ValueError(
                "development manifest does not assert that the meta_test "
                "label artifact was withheld")

        self.visible = visible
        self.authorization = authorization
        self.sealed_directory = sealed_directory

        # Identity-only ordering record. Reading it restores the row order the
        # all-label loader produces, so mounting the view renumbers nothing and
        # every recorded episode index stays valid. It carries no label.
        governance = read_governance(self.directory)
        recorded_governance = self.manifest["governance_jsonl_sha256"]
        actual_governance = sha256_file(self.directory / "governance.jsonl")
        if recorded_governance != actual_governance:
            raise ValueError(
                f"governance.jsonl hash mismatch: {actual_governance} != "
                f"{recorded_governance}")
        corpus_order = {row["cell_id"]: index
                        for index, row in enumerate(governance)}
        governed_split = {row["cell_id"]: row["split"] for row in governance}

        self.cells: list[dict] = []
        for split in visible:
            root = (sealed_directory if split == SEALED_SPLIT
                    else self.directory)
            path = root / split / "cells.jsonl.gz"
            recorded = self.manifest["artifacts"].get(f"{split}/cells.jsonl.gz")
            actual = sha256_file(path)
            if recorded != actual:
                raise ValueError(
                    f"{split} label artifact hash mismatch: {actual} != {recorded}")
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                rows = [json.loads(line) for line in handle if line.strip()]
            if any(row["split"] != split for row in rows):
                raise ValueError(f"{split} artifact contains another split")
            if any(governed_split.get(row["cell_id"]) != split for row in rows):
                raise ValueError(
                    f"{split} artifact disagrees with governance.jsonl on "
                    f"cell membership")
            self._check_bindings(split, rows)
            self.cells.extend(rows)
        self.cells.sort(key=lambda row: corpus_order[row["cell_id"]])

    @staticmethod
    def _read_manifest(directory: Path) -> dict:
        """Fail closed on an absent, unparseable or incomplete manifest."""
        path = directory / "manifest.json"
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError as error:
            raise ValueError(
                f"no governed-split-view manifest at {path}") from error
        except json.JSONDecodeError as error:
            raise ValueError(
                f"malformed governed-split-view manifest at {path}: "
                f"{error}") from error
        if not isinstance(manifest, dict):
            raise ValueError("malformed governed-split-view manifest: "
                             "expected an object")
        if manifest.get("schema") != SCHEMA:
            raise ValueError("unsupported governed-split-view schema")
        missing = [field for field in REQUIRED_MANIFEST_FIELDS
                   if field not in manifest]
        if missing:
            raise ValueError(
                f"governed-split-view manifest is missing binding(s): {missing}")
        return manifest

    def _check_bindings(self, split: str, rows: list[dict]) -> None:
        """Count- and index-equivalence against the recorded governed source."""
        joined = "\n".join(sorted(row["cell_id"] for row in rows))
        recorded_index = self.manifest["cell_id_index_sha256"].get(split)
        actual_index = hashlib.sha256(joined.encode("utf-8")).hexdigest()
        if recorded_index != actual_index:
            raise ValueError(
                f"{split} cell_id index mismatch: {actual_index} != "
                f"{recorded_index}")
        recorded_counts = self.manifest["counts"].get(split, {})
        actual_counts = {
            "rows": len(rows),
            "targets": len({row["target_id"] for row in rows}),
            "components": len({row["protein_group_40"] for row in rows}),
        }
        if any(recorded_counts.get(key) != value
               for key, value in actual_counts.items()):
            raise ValueError(
                f"{split} count mismatch: {actual_counts} != {recorded_counts}")

    @property
    def sealed_cell_count(self) -> int:
        """Rows withheld from this mount, read from the manifest count only.

        This is the row count the builder recorded; no sealed label is opened
        to produce it, and it is the only meta_test quantity the development
        surface carries.
        """
        if SEALED_SPLIT in self.visible:
            return 0
        return int(self.manifest["counts"][SEALED_SPLIT]["rows"])

    def isolation_record(self) -> dict:
        return {
            "level": ("open" if SEALED_SPLIT in self.visible
                      else "physically_isolated"),
            "visible_splits": list(self.visible),
            "labels_parsed_in_process": list(self.visible),
            "physically_isolated": SEALED_SPLIT not in self.visible,
            "meta_test_label_artifact_present_in_view": (
                SEALED_SPLIT in self.visible),
            "authorization": self.authorization,
            "source_manifest_sha256": self.manifest["source_manifest_sha256"],
            "split_assignment_sha256": self.manifest["split_assignment_sha256"],
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sealed-output", type=Path, required=True,
                        help="path for the meta_test label artifact; must be "
                             "outside --output")
    parser.add_argument("--authorization", required=True,
                        help="written reason, recorded in the sealing record")
    arguments = parser.parse_args()

    try:
        check_authorization(arguments.authorization)
    except ValueError as error:
        parser.error(str(error))
    output = arguments.output.resolve()
    sealed = arguments.sealed_output.resolve()
    if sealed == output or output in sealed.parents:
        parser.error("--sealed-output must be outside --output")
    if output.exists() or sealed.exists():
        raise FileExistsError("refusing to overwrite an existing view")

    manifest_path = arguments.corpus / "manifest.json"
    assignment_path = arguments.split_directory / "assignment.json"
    split_manifest = json.loads(
        (arguments.split_directory / "manifest.json").read_text(encoding="utf-8"))
    assignment = json.loads(assignment_path.read_text(encoding="utf-8"))
    recorded = split_manifest.get("assignment_sha256")
    actual = hashlib.sha256(
        json.dumps(assignment, sort_keys=True).encode("utf-8")).hexdigest()
    if recorded != actual:
        raise ValueError(f"split assignment hash mismatch: {actual} != {recorded}")

    # The single read of the all-label corpus. Rows are partitioned as they are
    # streamed; no structure ever holds two splits at once beyond the three
    # output buffers, and no label is inspected.
    buckets: dict[str, list[dict]] = {"meta_train": [], "meta_val": [],
                                      SEALED_SPLIT: []}
    governance: list[dict] = []
    with gzip.open(arguments.corpus / "cells.jsonl.gz", "rt",
                   encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            cell = json.loads(line)
            if cell["cell_id"] not in assignment:
                continue
            cell["split"] = assignment[cell["cell_id"]]
            buckets[cell["split"]].append(cell)
            governance.append({"cell_id": cell["cell_id"],
                               "split": cell["split"],
                               "target_id": cell["target_id"],
                               "protein_group_40": cell["protein_group_40"]})

    for split in DEVELOPMENT_SPLITS:
        write_cells(output / split / "cells.jsonl.gz", buckets[split])
    write_cells(sealed / SEALED_SPLIT / "cells.jsonl.gz", buckets[SEALED_SPLIT])
    (output / "governance.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in governance),
        encoding="utf-8")

    artifacts = {f"{split}/cells.jsonl.gz":
                 sha256_file(output / split / "cells.jsonl.gz")
                 for split in DEVELOPMENT_SPLITS}
    artifacts[f"{SEALED_SPLIT}/cells.jsonl.gz"] = sha256_file(
        sealed / SEALED_SPLIT / "cells.jsonl.gz")

    def cell_index_sha256(rows: list[dict]) -> str:
        joined = "\n".join(sorted(row["cell_id"] for row in rows))
        return hashlib.sha256(joined.encode("utf-8")).hexdigest()

    manifest = {
        "schema": SCHEMA,
        "built": "2026-08-16",
        "source_corpus": arguments.corpus.name,
        "source_manifest_sha256": sha256_file(manifest_path),
        "split_assignment_sha256": actual,
        "artifacts": artifacts,
        "cell_id_index_sha256": {
            split: cell_index_sha256(rows) for split, rows in buckets.items()},
        "counts": {split: {"rows": len(rows),
                           "targets": len({r["target_id"] for r in rows}),
                           "components": len({r["protein_group_40"] for r in rows})}
                   for split, rows in buckets.items()},
        # The development surface withholds the sealed label artifact. Its hash
        # is recorded so a confirmation run can verify what it mounts; its path
        # is not.
        "meta_test_label_artifact_emitted": False,
        "meta_test_label_artifact_location": "out of tree, supplied at mount time",
        "governance_jsonl_sha256": sha256_file(output / "governance.jsonl"),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    (sealed / "SEALING_RECORD.json").write_text(json.dumps({
        "schema": "MetaSieve.SealingRecord.v1",
        "built": "2026-08-16",
        "authorization": arguments.authorization,
        "development_surface": str(output),
        "source_manifest_sha256": manifest["source_manifest_sha256"],
        "split_assignment_sha256": actual,
        "artifacts": artifacts,
        "counts": manifest["counts"],
        "metrics_computed": None,
        "note": ("this process read meta_test labels once to partition them "
                 "and computed no statistic of them beyond a row count and a "
                 "file hash"),
    }, indent=1, sort_keys=True), encoding="utf-8")

    print(f"development surface: {output}")
    for split in DEVELOPMENT_SPLITS:
        print(f"  {split:<11} {manifest['counts'][split]}")
    print(f"sealed artifact:     {sealed} "
          f"({manifest['counts'][SEALED_SPLIT]['rows']} rows, withheld)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
