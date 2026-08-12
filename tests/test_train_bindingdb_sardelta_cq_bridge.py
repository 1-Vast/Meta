import numpy as np

from research.crossed_interaction.train_bindingdb_sardelta_cq_bridge import (
    component_contrast,
)


def test_component_contrast_passes_only_when_correct_beats_zero():
    rows = []
    for component in ("a", "b", "c"):
        rows.append({
            "dependency_component": component,
            "squared_error": 0.25,
            "zero_squared_error": 1.0,
        })

    result = component_contrast(rows, draws=999, seed=1)

    assert result["components"] == 3
    assert result["component_macro_reduction"] == 0.75
    assert result["one_sided_95_lcb"] > 0.0
    assert result["pass"] is True


def test_component_contrast_fails_when_correct_is_worse():
    rows = []
    for component in ("a", "b", "c"):
        rows.append({
            "dependency_component": component,
            "squared_error": 2.0,
            "zero_squared_error": 1.0,
        })

    result = component_contrast(rows, draws=999, seed=1)

    assert np.isclose(result["component_macro_reduction"], -1.0)
    assert result["pass"] is False
