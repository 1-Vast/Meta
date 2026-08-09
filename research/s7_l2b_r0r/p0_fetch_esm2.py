"""Phase 0 item 5 — range-resuming acquisition of esm2_t33_650M_UR50D.

The HuggingFace client restarted rather than range-resumed on this proxy, so the
transfer never completed. This uses curl with `-C -` (byte-range resume) against
a stable local path, retries until the declared size is reached, then verifies
the git-LFS SHA-256 published by the HuggingFace API.

Afterwards the model is loaded from the local directory, so no network access is
needed at extraction or training time.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(r"D:\MetaSieve")
DEST = ROOT / "dataset" / "raw" / "esm2_t33_650M_UR50D"
OUT = ROOT / "report" / "s7_l2b_r0r"
REPO = "facebook/esm2_t33_650M_UR50D"
BASE = f"https://huggingface.co/{REPO}/resolve/main"
FILES = ["config.json", "vocab.txt", "tokenizer_config.json",
         "special_tokens_map.json", "pytorch_model.bin"]
MAX_ATTEMPTS = 40


def api_metadata():
    url = f"https://huggingface.co/api/models/{REPO}?blobs=true"
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode())


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def fetch(name: str, expect_size: int | None) -> Path:
    dest = DEST / name
    for attempt in range(1, MAX_ATTEMPTS + 1):
        have = dest.stat().st_size if dest.exists() else 0
        if expect_size and have >= expect_size:
            return dest
        print(f"  {name}: attempt {attempt}, have {have} bytes"
              f"{'' if not expect_size else f' of {expect_size}'}", flush=True)
        cmd = ["curl", "-sS", "-L", "-C", "-", "--fail-with-body",
               "--connect-timeout", "30", "--max-time", "900",
               "--retry", "0", "-o", str(dest), f"{BASE}/{name}"]
        subprocess.run(cmd, capture_output=True, text=True)
        new = dest.stat().st_size if dest.exists() else 0
        if expect_size and new >= expect_size:
            return dest
        if not expect_size and new > 0:
            return dest
        if new == have:
            time.sleep(3)          # no progress; brief backoff then resume
    return dest


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    meta = api_metadata()
    sizes, oids = {}, {}
    for s in meta.get("siblings", []):
        fn = s.get("rfilename")
        if fn in FILES:
            sizes[fn] = s.get("size")
            lfs = s.get("lfs") or {}
            oids[fn] = lfs.get("oid") or lfs.get("sha256")
    print(f"declared sizes: { {k: sizes.get(k) for k in FILES} }", flush=True)

    results = {}
    for fn in FILES:
        p = fetch(fn, sizes.get(fn))
        got = p.stat().st_size if p.exists() else 0
        h = sha256(p) if p.exists() else None
        ok_size = (sizes.get(fn) is None) or (got == sizes[fn])
        ok_hash = (oids.get(fn) is None) or (h == oids[fn])
        results[fn] = {"bytes": got, "declared_bytes": sizes.get(fn),
                       "sha256": h, "declared_lfs_sha256": oids.get(fn),
                       "size_ok": bool(ok_size), "hash_ok": bool(ok_hash),
                       "verdict": "PASS" if (ok_size and ok_hash and got > 0)
                                  else "FAIL"}
        print(f"  {fn}: {results[fn]['verdict']} ({got} bytes)", flush=True)

    allok = all(v["verdict"] == "PASS" for v in results.values())
    manifest = {
        "schema": "MetaSieve.S7L2B.P0.ESM2AcquisitionManifest.v1",
        "model_id": REPO,
        "revision_sha": meta.get("sha"),
        "local_dir": str(DEST.relative_to(ROOT)).replace("\\", "/"),
        "method": "curl -C - byte-range resume, repeated until declared size",
        "files": results,
        "all_files_verified": allok,
        "licence": "MIT (facebookresearch/esm)",
        "elapsed_sec": round(time.time() - t0, 1),
        "verdict": "ESM2_WEIGHTS_ACQUIRED_AND_VERIFIED" if allok
                   else "ESM2_WEIGHT_ACQUISITION_INCOMPLETE",
    }
    (OUT / "P0_ESM2_ACQUISITION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in
                      ("revision_sha", "all_files_verified", "verdict",
                       "elapsed_sec")}, indent=2))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
