from research.meta_fewshot.select_v1_section_hyperparameters import (
    choose_by_validation,
)


def test_section_selection_uses_seed_mean_and_deterministic_tie_break():
    rows = [
        {"section_dim": 1, "ridge": 1.0, "best_combined_val_score": 2.0},
        {"section_dim": 1, "ridge": 1.0, "best_combined_val_score": 4.0},
        {"section_dim": 2, "ridge": 0.1, "best_combined_val_score": 2.5},
        {"section_dim": 2, "ridge": 0.1, "best_combined_val_score": 2.5},
    ]
    selected = choose_by_validation(rows)
    assert selected == {
        "section_dim": 2,
        "ridge": 0.1,
        "mean_best_combined_val_score": 2.5,
    }
