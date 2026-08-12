import numpy as np

from research.crossed_interaction.train_bindingdb_sardelta_attribution import (
    antisymmetry_audit,
    component_contrast,
    fit_positive_ridge_no_intercept,
    interaction_feature,
    predict,
)


def test_interaction_feature_is_outer_product_flattened():
    protein = np.asarray([2.0, 3.0], dtype=np.float64)
    delta = np.asarray([5.0, 7.0], dtype=np.float64)

    feature = interaction_feature(protein, delta)

    assert np.allclose(feature, [10.0, 14.0, 15.0, 21.0])


def test_no_intercept_ligand_delta_model_is_antisymmetric():
    x = np.asarray([[1.0, 2.0], [-2.0, 1.0], [0.5, -1.0]], dtype=np.float64)
    y = np.asarray([1.0, -0.5, 0.25], dtype=np.float64)
    model = fit_positive_ridge_no_intercept(x, y, ridge=1.0)

    audit = antisymmetry_audit(model, x, -x)

    assert audit["max_abs_sum"] < 1e-12
    assert np.allclose(predict(model, -x), -predict(model, x))


def test_component_contrast_passes_when_correct_has_lower_mse():
    rows = []
    for component in ("a", "b", "c"):
        rows.append({
            "dependency_component": component,
            "I_squared_error": 0.25,
            "L_squared_error": 1.0,
        })

    result = component_contrast(rows, "I", "L", draws=999, seed=1)

    assert result["component_macro_reduction"] == 0.75
    assert result["pass"] is True
