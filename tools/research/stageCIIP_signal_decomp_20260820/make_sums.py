"""Build SHA256SUMS for all frozen stage artifacts (deliverable)."""
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKIP = {"SHA256SUMS", "SMOKE_RESULT.json", "seed1.log", "gen_erased_s1.log",
        "__pycache__"}
lines = []
for p in sorted(HERE.rglob("*")):
    if p.is_dir() or p.name in SKIP or "__pycache__" in str(p) or p.name.startswith("_"):
        continue
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    lines.append(f"{h} *{p.relative_to(HERE).as_posix()}")
(HERE / "SHA256SUMS").write_text("
".join(lines) + "
", encoding="utf-8")
print(f"wrote SHA256SUMS with {len(lines)} entries")
