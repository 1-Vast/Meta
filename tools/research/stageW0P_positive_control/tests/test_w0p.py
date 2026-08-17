"""Stage W0-P panel tests."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import sys
REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
STAGE = Path(__file__).resolve().parents[1]
PREREG_SHA = "ba0b51ec419b0275a129e69e4cb45db1bccbdd138000893ee7daf881e7bacbf1"

def test_prereg_and_panel():
    assert hashlib.sha256((STAGE/'PREREGISTRATION.md').read_bytes()).hexdigest() == PREREG_SHA
    d=json.loads((STAGE/'W0P_PANEL.json').read_text(encoding='utf-8'))
    assert d['preregistration_sha256'] == PREREG_SHA
    assert d['summary']['pairs'] == 6
    assert d['summary']['total_rows'] == 32
    assert d['summary']['sufficient'] is True
    for p in d['pairs']:
        assert 1 <= len(p['mutation_positions']) <= 5
        assert p['n_shared_ligands'] >= 3
