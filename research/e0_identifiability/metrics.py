"""Metrics shared by the synthetic E0 identifiability audits."""

import numpy as np


def concordance(labels: np.ndarray, predictions: np.ndarray) -> float:
    left, right = np.triu_indices(len(labels), 1)
    delta_y = labels[left] - labels[right]
    keep = delta_y != 0
    if not np.any(keep):
        return 0.5
    delta_p = predictions[left][keep] - predictions[right][keep]
    return float(((np.sign(delta_y[keep]) == np.sign(delta_p)).sum()
                  + 0.5 * (delta_p == 0).sum()) / len(delta_p))
