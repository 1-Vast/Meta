import numpy as np

from .sample_complexity import _fit_budget_certificate


def test_budget_certificate_rejects_an_uninformative_predictor():
    rng = np.random.default_rng(1)
    features = rng.normal(size=(90, 14))
    margins = rng.normal(scale=0.01, size=90)
    folds = np.repeat(np.arange(3), 30)
    targets = np.repeat(np.arange(45), 2)
    _, _, report = _fit_budget_certificate(features, margins, folds, targets)
    assert not report["source_gate_passed"]


def test_budget_certificate_report_is_deterministic():
    rng = np.random.default_rng(2)
    features = rng.normal(size=(90, 14))
    margins = 0.01 * features[:, 0] + 0.02
    folds = np.repeat(np.arange(3), 30)
    targets = np.arange(90)
    _, q1, report1 = _fit_budget_certificate(features, margins, folds, targets)
    _, q2, report2 = _fit_budget_certificate(features, margins, folds, targets)
    assert q1 == q2
    assert report1 == report2
