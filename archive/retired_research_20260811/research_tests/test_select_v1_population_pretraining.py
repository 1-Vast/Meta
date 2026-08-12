from research.meta_fewshot.select_v1_population_pretraining import (
    choose_pretraining,
)


def test_pretraining_selection_uses_mean_and_prefers_shorter_tie():
    rows = [
        {"population_pretrain_steps": 0, "best_combined_val_score": 3.0},
        {"population_pretrain_steps": 0, "best_combined_val_score": 4.0},
        {"population_pretrain_steps": 100, "best_combined_val_score": 3.0},
        {"population_pretrain_steps": 100, "best_combined_val_score": 3.0},
    ]
    selected = choose_pretraining(rows)
    assert selected["population_pretrain_steps"] == 100
    assert selected["mean_best_combined_val_score"] == 3.0
