"""Freeze exact-assay choices for X1A-R without reading affinity values."""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RECT_ROOT = ROOT / "research" / "crossed_interaction" / "artifacts" / "x1a_r_direct_dd"
INPUT = RECT_ROOT / "rectangles.jsonl"
INPUT_SHA256 = "22f3e738f4dbc7b53ca9ef23e995e2a398cbca280a9cdde12c546be21500d0a5"
OUTPUT = ROOT / "research" / "crossed_interaction" / "artifacts" / "x1a_r_assays"
REQUIRED_N = 245
MAX_SHARE = 0.25


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def choose_common_assay(left: dict, right: dict) -> tuple[str | None, list[int], list[int]]:
    left_map = left["assay_activity_ids"]
    right_map = right["assay_activity_ids"]
    common = sorted(set(left_map) & set(right_map))
    if not common:
        return None, [], []
    assay = min(common, key=lambda key: (-min(len(left_map[key]), len(right_map[key])), key))
    return assay, sorted(map(int, left_map[assay])), sorted(map(int, right_map[assay]))


def select_rectangle(row: dict) -> dict:
    cells = row["cells"]
    assay_a, ids_aa, ids_ab = choose_common_assay(cells[0], cells[1])
    assay_b, ids_ba, ids_bb = choose_common_assay(cells[2], cells[3])
    reason = None
    if assay_a is None:
        reason = "protein_a_has_no_common_exact_assay"
    elif assay_b is None:
        reason = "protein_b_has_no_common_exact_assay"
    return {
        "rectangle_id": row["rectangle_id"],
        "endpoint": row["endpoint"],
        "dependency_cluster": row["dependency_cluster"],
        "panel_id": row["panel_id"],
        "selected_at_frozen_cap": row["selected_at_frozen_cap"],
        "eligible_primary": reason is None,
        "exclusion_reason": reason,
        "protein_a_assay_id": assay_a,
        "protein_b_assay_id": assay_b,
        "cell_activity_ids": [ids_aa, ids_ab, ids_ba, ids_bb],
    }


def run(output: Path = OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"no-clobber output exists: {output}")
    if sha256_file(INPUT) != INPUT_SHA256:
        raise RuntimeError("rectangle manifest hash mismatch")
    with INPUT.open(encoding="utf-8") as handle:
        selected = [select_rectangle(json.loads(line)) for line in handle if line.strip()]

    endpoint_report = {}
    for endpoint in ("Ki", "Kd"):
        rows = [row for row in selected if row["endpoint"] == endpoint]
        capped = [row for row in rows if row["selected_at_frozen_cap"]]
        primary = [row for row in capped if row["eligible_primary"]]
        by_cluster = Counter(row["dependency_cluster"] for row in primary)
        total = len(primary)
        largest_share = max(by_cluster.values(), default=0) / total if total else 1.0
        endpoint_report[endpoint] = {
            "rectangles": len(rows),
            "frozen_cap_rectangles": len(capped),
            "primary_exact_assay_rectangles": total,
            "primary_clusters": len(by_cluster),
            "largest_primary_cluster_share": largest_share,
            "exclusions": dict(Counter(row["exclusion_reason"] for row in capped
                                       if not row["eligible_primary"])),
            "nominal_design_pass": total >= REQUIRED_N and largest_share <= MAX_SHARE,
        }

    output.mkdir(parents=True)
    rows_path = output / "assay_selection.jsonl"
    with rows_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(selected, key=lambda item: (item["endpoint"], item["rectangle_id"])):
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    report = {
        "stage": "E-AFF-X1A-R_EXACT_ASSAY_SELECTION",
        "label_blind": True,
        "input_rectangles_sha256": INPUT_SHA256,
        "assay_selection_sha256": sha256_file(rows_path),
        "endpoints": endpoint_report,
        "required_nominal_units": REQUIRED_N,
        "maximum_cluster_share": MAX_SHARE,
        "affinity_value_fields_selected": 0,
        "training_performed": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
