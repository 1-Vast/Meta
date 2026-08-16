import pytest

from research.meta_fewshot.select_v1_population_width import (
    choose_population_width,
)


def test_population_width_selection_uses_seed_mean():
    rows = [
        {"population_hidden_dim": 32, "best_combined_val_score": 3.0},
        {"population_hidden_dim": 32, "best_combined_val_score": 4.0},
        {"population_hidden_dim": 64, "best_combined_val_score": 3.1},
        {"population_hidden_dim": 64, "best_combined_val_score": 3.2},
    ]
    selected = choose_population_width(rows)
    assert selected["population_hidden_dim"] == 64
    assert selected["mean_best_combined_val_score"] == pytest.approx(3.15)
