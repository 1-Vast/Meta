"""Download immutable external inputs directly on the remote machine.

Run from the project root. Dataset downloads are accepted only when their
SHA-256 matches the locally validated source copy. Hugging Face snapshots use
the standard cache layout consumed by preprocessing/protein/embedding.py.
"""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import sys
import tempfile
import urllib.request


ROOT = Path(__file__).resolve().parents[2]
DATASETS = {
    "davis.tab": {
        "url": "https://dataverse.harvard.edu/api/access/datafile/5219748",
        "sha256": "6d4c6809dcb7c5da2b91a32d594d6935b75484940bde4d18055eb5e1059262f4",
    },
    "kiba.tab": {
        "url": "https://dataverse.harvard.edu/api/access/datafile/5255037",
        "sha256": "eb77bed3ba64cc0cd07e05ff7f4b94f94a866a54cf673bb2c798cfeb0ebe2322",
    },
}
MODELS = (
    "facebook/esm2_t6_8M_UR50D",
    "facebook/esm2_t30_150M_UR50D",
    "facebook/esm2_t33_650M_UR50D",
)
MODEL_FILES = (
    "config.json",
    "model.safetensors",
    "special_tokens_map.json",
    "tokenizer_config.json",
    "vocab.txt",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_dataset(name: str, spec: dict[str, str]) -> None:
    destination = ROOT / "dataset" / "raw" / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and sha256(destination) == spec["sha256"]:
        print(f"dataset ok: {name}")
        return

    fd, temporary_name = tempfile.mkstemp(prefix=f"{name}.", dir=destination.parent)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(
            spec["url"], headers={"User-Agent": "MetaSieve-DTA/phase3"})
        with urllib.request.urlopen(request, timeout=120) as response:
            with temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1 << 20)
        actual = sha256(temporary)
        if actual != spec["sha256"]:
            raise RuntimeError(
                f"{name}: upstream SHA-256 {actual} != expected {spec['sha256']}")
        temporary.replace(destination)
        print(f"dataset downloaded: {name} ({destination.stat().st_size} bytes)")
    finally:
        temporary.unlink(missing_ok=True)


def download_weights(endpoint: str | None = None) -> None:
    if endpoint:
        # Mirror endpoints generally proxy regular repository files, not Xet CAS.
        os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
    try:
        from huggingface_hub import snapshot_download
    except ImportError as exc:
        raise RuntimeError(
            "huggingface-hub is required; install config/remote/requirements.txt") from exc

    cache = ROOT / "weights" / "hub"
    cache.mkdir(parents=True, exist_ok=True)
    for model in MODELS:
        snapshot = snapshot_download(
            repo_id=model,
            cache_dir=cache,
            endpoint=endpoint,
            allow_patterns=MODEL_FILES,
        )
        print(f"weights ready: {model} -> {snapshot}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets-only", action="store_true")
    parser.add_argument("--weights-only", action="store_true")
    parser.add_argument(
        "--hf-endpoint",
        help="Hugging Face-compatible endpoint used when huggingface.co is unreachable",
    )
    args = parser.parse_args()
    if args.datasets_only and args.weights_only:
        parser.error("choose at most one of --datasets-only and --weights-only")
    if not args.weights_only:
        for name, spec in DATASETS.items():
            download_dataset(name, spec)
    if not args.datasets_only:
        download_weights(args.hf_endpoint)
    return 0


if __name__ == "__main__":
    sys.exit(main())
