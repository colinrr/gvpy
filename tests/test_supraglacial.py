import numpy as np
import pytest

from gvpy.ice_cauldron import IceCauldron


class TestGetQS:
    """get_q_s (supraglacial drainage), translated from @IceCauldron/get_q_s.m."""

    def test_off_mode_returns_zeros(self):
        c = IceCauldron(supraglacial_drainage_mode="off")
        V_fluid_n = np.array([50.0])
        V_cavity_n = np.array([100.0])
        dV_fluid_n = np.array([1.0])
        dV_cavity_n = np.array([2.0])

        q_s_n = c.get_q_s(V_fluid_n, V_cavity_n, dV_fluid_n, dV_cavity_n)

        assert np.all(q_s_n == 0)

    def test_overflow_mode_no_overflow(self):
        c = IceCauldron()
        assert c.supraglacial_drainage_mode == "overflow"

        V_fluid_n = np.array([50.0])
        V_cavity_n = np.array([100.0])
        dV_fluid_n = np.array([1.0])
        dV_cavity_n = np.array([2.0])

        q_s_n = c.get_q_s(V_fluid_n, V_cavity_n, dV_fluid_n, dV_cavity_n)

        assert np.all(q_s_n == 0)

    def test_overflow_condition_computes_expected_value(self):
        """The overflow branch's assignment is a normal numpy operation under
        well-shaped inputs (the earlier assert already guarantees matching
        shapes), so it succeeds without ever reaching the `faafo("fix",
        "it")` placeholder in its `except` clause - that placeholder (see
        get_q_s's TRANSLATION NOTE) is dead/defensive code here, not a bug
        that fires in normal use."""
        c = IceCauldron()
        V_fluid_n = np.array([99.0])
        V_cavity_n = np.array([100.0])
        dV_fluid_n = np.array([5.0])
        dV_cavity_n = np.array([2.0])

        q_s_n = c.get_q_s(V_fluid_n, V_cavity_n, dV_fluid_n, dV_cavity_n)

        assert q_s_n[0] == pytest.approx(5.0 - 2.0)

    def test_mismatched_shapes_raise(self):
        c = IceCauldron()
        with pytest.raises(AssertionError):
            c.get_q_s(np.array([1.0, 2.0]), np.array([1.0]), np.array([1.0]), np.array([1.0]))
