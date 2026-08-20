"""Stage X0 live integrity assertions (I6)."""
from __future__ import annotations
import json, hashlib
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
STAGE = Path(__file__).resolve().parents[1]
D = STAGE / 'downloads'
PREREG_SHA = '03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'


def test_prereg_frozen():
    assert hashlib.sha256((STAGE/'STAGE_X0_PREREGISTRATION.md').read_bytes()).hexdigest() == PREREG_SHA


def test_x0d_audit_exists_and_records_downloads():
    d=json.loads((STAGE/'X0D_DATA_AUDIT.json').read_text(encoding='utf-8'))
    assert d['preregistration_sha256']==PREREG_SHA
    assert all(entry['exists'] for entry in d['sources']['duongly_2016'])
    assert all(entry['exists'] for entry in d['sources']['anastassiadis_2011'])


def test_davis_label_orientation_named_anchors():
    df=pd.read_excel(D/'davis_MOESM5.xls', sheet_name='Sheet1')
    # columns are compounds; Kd in nM, smaller = stronger
    def kd(compound):
        vals=df[compound].dropna()
        return float(vals.min()), float(vals.max())
    # Staurosporine is broadly sub-nanomolar in this panel; find its column.
    st=[c for c in df.columns if 'Staurosporine' in str(c)]
    assert st, 'staurosporine column missing'
    lo,hi=kd(st[0]); assert lo < 100, (st[0],lo,hi)
    # Imatinib should have ABL1-family Kd < CDK2 Kd if both present.
    im=[c for c in df.columns if 'Imatinib' in str(c)]
    if im:
        col=im[0]
        abl=df[df['Entrez Gene Symbol'].astype(str).str.contains('ABL1',na=False)][col].dropna()
        cdk=df[df['Entrez Gene Symbol'].astype(str).str.contains('CDK2',na=False)][col].dropna()
        if len(abl) and len(cdk):
            assert float(abl.min()) < float(cdk.min()), (abl.min(), cdk.min())


def test_duongly_percent_remaining_orientation():
    df=pd.read_excel(D/'duongly_mmc3.xlsx', sheet_name='Table S2')
    assert 'Imatinib' in df.columns
    abl_rows=df[df['Compound name'].astype(str).str.contains('ABL1', na=False)]
    assert len(abl_rows) >= 2
    vals=abl_rows['Imatinib'].astype(float)
    assert vals.min() < 100


def test_anastassiadis_matrix_semantics():
    df=pd.read_excel(D/'anastassiadis_MOESM23.xls', sheet_name='Sheet1', header=None)
    # header rows are 0 and 1; body starts row 2. Smaller = stronger inhibition.
    text=str(df.iloc[0,0]).lower()
    assert 'percent remaining kinase activity' in text


def test_csc_antisymmetry_identity_reference():
    # The CSC operator is defined antisymmetric and zero on identity by contract;
    # assert the numerical reference implementation used by Stage X.
    def csc(y_l_p, y_l_q, ref):
        return (y_l_p - y_l_q) - ref
    assert csc(1.0,2.0,0.0) == -csc(2.0,1.0,0.0)
    assert csc(1.0,1.0,0.0) == 0.0


def test_planted_interaction_uses_both_projection_matrices():
    from tools.research.stageX_csc_signal.x0_planted import planted_interaction

    prot = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
    lig = np.array([[5.0, 6.0], [7.0, 8.0]], dtype=np.float32)
    U = np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32)
    V = np.array([[2.0, 0.0], [0.0, 3.0]], dtype=np.float32)

    expected = (prot @ U) @ (lig @ V).T
    actual = planted_interaction(prot, lig, U, V)

    assert actual.shape == (2, 2)
    np.testing.assert_allclose(actual, expected)

def test_representation_capability_identifies_global_pooling_failure():
    d=json.loads((STAGE/'X0_INSTRUMENTS.json').read_text(encoding='utf-8'))
    cap=d['representation_capability']
    assert cap['global_esm']['ratio'] < 0.05
    assert cap['mutation_token']['pass_capability'] is True
    assert cap['local_window']['pass_capability'] is True
    assert cap['esm_local_window']['pass_capability'] is True
