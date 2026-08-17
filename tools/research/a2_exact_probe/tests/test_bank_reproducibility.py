"""Episode banks must be bitwise reproducible across independent processes.

Stage R's first run seeded its episode draws with Python's built-in `hash()`
over a tuple containing a string. `hash()` on `str` is salted per process by
`PYTHONHASHSEED`, so the "fixed" banks differed between runs:

    $ python -c "print(abs(hash(('a2-exact', 20260818, 7, 0))) % (2**32))"
    1872535333
    $ python -c "print(abs(hash(('a2-exact', 20260818, 7, 0))) % (2**32))"
    3725077289

Within a single process every arm still drew the same episodes, so the paired
contrasts in the superseded run were internally valid. The artifact was not
reproducible, which is a governance defect on its own, and it is the kind that
silently invalidates any later attempt to re-derive a number.

The repair uses `scripts.qpsmp_data.stable_seed` (sha256) keyed on the target's
*name* rather than its positional index, so the bank is also immune to a change
in the ligand cap or the eligibility filter renumbering the targets.

The first test below is the one that matters: it launches genuinely separate
interpreters with hostile `PYTHONHASHSEED` values and compares the emitted
episode identities.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import stable_seed                       # noqa: E402

FEATURES = ROOT / "tools/research/a2_exact_probe/features"

EMIT = """
import json, sys
sys.path.insert(0, r"{root}")
from tools.research.a2_exact_probe.run_probe import Corpus
corpus = Corpus(r"{path}", "cpu")
banks = corpus.nested_banks(4, 73101)
out = {{str(k): corpus.episode_identity(v) for k, v in banks.items()}}
print(json.dumps(out))
"""


def emit(hash_seed: str) -> dict:
    """Run the bank construction in a fresh interpreter with a forced salt."""
    environment = dict(os.environ, PYTHONHASHSEED=hash_seed)
    script = EMIT.format(root=str(ROOT), path=str(FEATURES / "meta_val.npz"))
    finished = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True,
        env=environment, cwd=str(ROOT), timeout=600)
    if finished.returncode != 0:
        raise AssertionError(finished.stderr[-3000:])
    return json.loads(finished.stdout)


# --- the load-bearing test ------------------------------------------------

@pytest.mark.skipif(not (FEATURES / "meta_val.npz").exists(),
                    reason="feature cache absent; regenerate with "
                           "extract_ligand_features.py")
def test_episode_identities_are_bitwise_equal_across_processes():
    """Two interpreters with different hash salts must emit identical banks.

    `PYTHONHASHSEED=0` disables randomisation; `=1` and `=4242` force two
    different salts. Under the old `hash()` seeding these three disagree.
    """
    baseline = emit("0")
    assert baseline and all(baseline.values()), "the bank is empty"
    for salt in ("1", "4242"):
        assert emit(salt) == baseline, (
            f"episode identities changed under PYTHONHASHSEED={salt}; "
            "the bank is not reproducible across processes")


@pytest.mark.skipif(not (FEATURES / "meta_val.npz").exists(),
                    reason="feature cache absent")
def test_the_bank_is_nested_and_disjoint():
    """Reproducibility is worthless if the bank is malformed.

    Support at size k must be the k-prefix of the size-5 support, the query
    panel must be shared across k, and support and query ligands must not
    overlap — the same contract `fixed_nested_episode_banks` holds.
    """
    banks = emit("0")
    sizes = ("0", "1", "2", "3", "5")
    assert {k: len(v) for k, v in banks.items()} == {
        k: len(banks["5"]) for k in sizes}
    for episodes in zip(*(banks[k] for k in sizes)):
        largest = episodes[-1]
        for size, item in zip((0, 1, 2, 3, 5), episodes):
            assert item[0] == largest[0]                  # same target
            assert item[1] == largest[1][:size]           # prefix support
            assert item[2] == largest[2]                  # shared query panel
            assert not set(item[1]) & set(item[2])        # disjoint ligands


# --- the seeding primitive -------------------------------------------------

def test_stable_seed_is_deterministic_across_processes():
    script = ("import sys; sys.path.insert(0, r'%s');"
              "from scripts.qpsmp_data import stable_seed;"
              "print(stable_seed('a2-exact', 73101, 'P12345', 3))" % ROOT)
    values = set()
    for salt in ("0", "1", "4242"):
        finished = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True,
            env=dict(os.environ, PYTHONHASHSEED=salt), cwd=str(ROOT),
            timeout=120)
        assert finished.returncode == 0, finished.stderr[-2000:]
        values.add(finished.stdout.strip())
    assert len(values) == 1, f"stable_seed is not stable: {values}"


def test_builtin_hash_is_not_stable_across_processes():
    """Documents the defect, so a future refactor cannot reintroduce it quietly."""
    script = "print(abs(hash(('a2-exact', 20260818, 7, 0))) % (2**32))"
    values = set()
    for salt in ("1", "4242"):
        finished = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True,
            env=dict(os.environ, PYTHONHASHSEED=salt), timeout=120)
        values.add(finished.stdout.strip())
    assert len(values) > 1, (
        "builtin hash() appears stable here; if PYTHONHASHSEED is pinned "
        "globally this test cannot detect the defect it documents")


def test_seeding_keys_on_the_target_name_not_its_index():
    """Renumbering targets must not change any episode."""
    assert stable_seed("a2-exact", 73101, "P00001", 0) != stable_seed(
        "a2-exact", 73101, "P00002", 0)
    assert stable_seed("a2-exact", 73101, "P00001", 0) == stable_seed(
        "a2-exact", 73101, "P00001", 0)


def test_run_probe_does_not_call_builtin_hash():
    """Parse the AST, so prose about the defect does not trip the guard."""
    import ast
    path = ROOT / "tools/research/a2_exact_probe/run_probe.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    offending = [node.lineno for node in ast.walk(tree)
                 if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Name)
                 and node.func.id == "hash"]
    assert offending == [], f"builtin hash() called at lines {offending}"
