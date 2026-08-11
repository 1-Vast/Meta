import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "dataset" / "raw" / "adambind_public_01a169a6"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def test_adambind_public_snapshot_matches_pinned_manifest():
    manifest = json.loads((SNAPSHOT / "acquisition_manifest.json").read_text())
    assert manifest["source_commit"] == "01a169a6d62fba0d6c003f47bfba539e55f5b344"
    assert manifest["exact_reproduction_authorized"] is False
    assert manifest["training_authorized"] is False
    for name, expected in manifest["files"].items():
        path = SNAPSHOT / name
        assert path.stat().st_size == expected["bytes"]
        assert sha256(path) == expected["sha256"]
