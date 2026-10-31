import numpy as np

from gvpy.ice_cauldron import IceCauldron
from gvpy.plume import predict_forest, predict_tree


class TestGetPlumeFlows:
    """get_plume_flows, translated from @IceCauldron/getPlumeFluxes.m.

    The "off" (default) and "test-fixed" paths are exercised here - the
    random-forest emulator path is an in-development/test mode that can't be
    run end-to-end without the (missing) emulator JSON file.
    """

    def test_disabled_plume_flow_returns_zeros(self):
        c = IceCauldron(n_cauldrons=2)
        assert c.plume_flow_mode == "off"
        for out in c.get_plume_flows(t=1000.0, Ze_n=10.0, Q_n=c.params["Q_n"].values, a_n=5.0, L_n=1000.0):
            assert np.array_equal(out, np.zeros(2))

    def test_fixed_flows_zero_while_closed(self):
        c = IceCauldron(n_cauldrons=2, plume_flow_mode="test-fixed")
        for out in c.get_plume_flows(t=1000.0, Ze_n=10.0, Q_n=c.params["Q_n"].values, a_n=5.0, L_n=1000.0):
            assert np.array_equal(out, np.zeros(2))

    def test_fixed_flows_profile_follows_each_open_time(self):
        """Cauldron 0 opens at 100 s, cauldron 1 at 1000 s: each profile is zero before its
        open time, at the vent-input plateau from +600 s to +2400 s, and zero after +6000 s."""
        c = IceCauldron(n_cauldrons=2, plume_flow_mode="test-fixed", Q_0=4e7)
        c.open_cauldron[:] = True
        c.open_time[:] = [100.0, 1000.0]
        Q_n = c.params["Q_n"].values
        q_s_peak = Q_n * (1 - c.n_w) / c.rho_p
        q_w_peak = Q_n * c.n_w / c.CONSTANTS.RHO_W

        def flows(t):
            return c.get_plume_flows(t=t, Ze_n=10.0, Q_n=Q_n, a_n=5.0, L_n=1000.0)

        h_m, q_w_p_n, q_s_p_n = flows(500.0)  # cauldron 1 not started yet
        assert q_s_p_n[0] > 0 and q_s_p_n[1] == 0.0 and h_m[1] == 0.0

        _, q_w_p_n, q_s_p_n = flows(1700.0)  # both on their plateau
        assert np.allclose(q_s_p_n, q_s_peak) and np.allclose(q_w_p_n, q_w_peak)

        h_m, q_w_p_n, q_s_p_n = flows(6200.0)  # cauldron 0 ended, cauldron 1 in its tail
        assert q_s_p_n[0] == 0.0 and 0 < q_s_p_n[1] < q_s_peak[1]

        h_m, _, _ = flows(1700.0)
        assert np.allclose(h_m, 345 * (c.rho_p * q_s_peak) ** 0.226)


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
