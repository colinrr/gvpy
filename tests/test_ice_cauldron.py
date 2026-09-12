import numpy as np

from gvpy.ice_cauldron import IceCauldron


class TestVentMeltingIceInflow:
    """Smoke tests for the vent/melting & ice inflow methods translated from
    @IceCauldron/get_u_ice.m, get_f_i.m, get_u_melt.m, get_da_dt.m,
    getMaterialHeights.m.

    These exercise only the closed-cauldron / default-state code paths that
    are reachable without the not-yet-translated ODE solver and events
    system - the open-cauldron and full-output branches remain untested.
    """

    def test_get_u_ice_off(self):
        c = IceCauldron()
        assert c.ice_inflow_mode == "off"
        u_ice_0, u_ice_bar = c.get_u_ice()
        assert u_ice_0 == 0
        assert u_ice_bar == 0

    def test_get_u_ice_fixed_glen(self):
        c = IceCauldron(ice_inflow_mode="fixed-glen")
        u_ice_0, u_ice_bar = c.get_u_ice()
        assert np.all(u_ice_0 > 0)
        assert np.all(u_ice_bar < 0)

    def test_get_f_i_default(self):
        c = IceCauldron()
        f_i_n = c.get_f_i()
        assert f_i_n.shape == (c.n_cauldrons, 1)
        assert np.all(f_i_n == c.f_i)

    def test_get_f_i_ode_mode(self):
        c = IceCauldron()
        a = np.array([[5.0]])
        cc = np.array([[10.0]])
        f_i_n = c.get_f_i(a, cc)
        assert f_i_n.shape == a.shape
        # Not ice-free by default, so f_i_n should equal the base f_i parameter
        assert np.all(f_i_n == c.f_i)

    def test_get_u_melt_closed_cauldron(self):
        c = IceCauldron()
        a = np.array([[5.0]])
        cc = np.array([[10.0]])
        f_i = c.get_f_i(a, cc)
        u_melt, v_melt = c.get_u_melt(a, cc, f_i)
        assert u_melt.shape == a.shape
        assert v_melt.shape == a.shape
        assert np.all(np.isfinite(u_melt))
        assert np.all(np.isfinite(v_melt))

    def test_get_da_dt(self):
        c = IceCauldron()
        a = np.array([[5.0]])
        cc = np.array([[10.0]])
        f_i = c.get_f_i(a, cc)
        da_dt, dc_dt = c.get_da_dt(a, cc, f_i, 0.0)
        assert da_dt.shape == a.shape
        assert dc_dt.shape == a.shape

    def test_get_material_heights_closed_cauldron(self):
        c = IceCauldron()
        shape = (2, 1)
        V_ice_n = np.full(shape, 100.0)
        V_w_n = np.full(shape, 10.0)
        V_p_n = np.full(shape, 5.0)
        V_cavity_n = np.full(shape, 200.0)
        a_n = np.full(shape, 5.0)
        c_n = np.full(shape, 10.0)

        H_i, H_w, H_p, H_cum = c.get_material_heights(V_ice_n, V_w_n, V_p_n, V_cavity_n, a_n, c_n)

        assert H_i.shape == shape
        assert H_w.shape == shape
        assert H_p.shape == shape
        assert H_cum.shape == shape
        assert np.all(H_w >= 0)
        assert np.allclose(H_cum, H_p + H_w + H_i)
