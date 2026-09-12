import numpy as np

from gvpy.ice_cauldron import IceCauldron
from gvpy.plume import predict_forest, predict_tree


class TestGetPlumeFluxes:
    """get_plume_fluxes, translated from @IceCauldron/getPlumeFluxes.m.

    Only the disable_plume_flux (default) path is exercised here - the
    random-forest emulator path is an in-development/test mode that can't be
    run end-to-end without the (missing) emulator JSON file.
    """

    def test_disabled_plume_flux_returns_zeros(self):
        c = IceCauldron()
        assert c.disable_plume_flux is True
        h_m, r_c, qw0, qwc, qs0, qsc = c.get_plume_fluxes(Ze=10.0, Q=5e7, a=5.0, L=1000.0)
        assert (h_m, r_c, qw0, qwc, qs0, qsc) == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


class TestPredictForest:
    """predict_tree/predict_forest, translated from getPlumeFluxes.m's embedded
    local functions - exercised here with small synthetic trees since the
    real emulator JSON is not available in this repo.
    """

    def test_predict_tree_leaf(self):
        leaf = {"value": [1.0, 2.0, 3.0]}
        result = predict_tree(leaf, np.array([0.0]))
        assert np.allclose(result, [1.0, 2.0, 3.0])

    def test_predict_tree_single_split(self):
        tree = {
            "feature": 0,
            "threshold": 0.5,
            "left": {"value": [1.0, 1.0]},
            "right": {"value": [9.0, 9.0]},
        }
        assert np.allclose(predict_tree(tree, np.array([0.1])), [1.0, 1.0])
        assert np.allclose(predict_tree(tree, np.array([0.9])), [9.0, 9.0])

    def test_predict_forest_averages_trees(self):
        tree_a = {"value": [1.0, 2.0]}
        tree_b = {"value": [3.0, 4.0]}
        forest = [tree_a, tree_b]
        X = np.array([[0.0]])
        result = predict_forest(forest, X)
        assert result.shape == (1, 2)
        assert np.allclose(result[0], [2.0, 3.0])
