import numpy as np
import pytest

from gvpy.glaciovolcano import glaciovolcano
from gvpy.gv_main import gv_main
from gvpy.ice_cauldron import IceCauldron


@pytest.mark.parametrize("ice_inflow_mode", ["off", "fixed-glen"])
def test_par_exposes_per_cauldron_flows(ice_inflow_mode):
    """Flow terms added to glaciovolcano's par are per cauldron, so
    gv_main stacks them as (time, cauldron)."""
    gv = IceCauldron(n_cauldrons=2, G_n=[300, 400], ice_inflow_mode=ice_inflow_mode)
    _, par = glaciovolcano(0.0, gv.get_initial_conditions(), gv, return_par=True)
    for name in ("Q_n", "u_ice_bar", "q_d_n", "q_c_n"):
        assert np.shape(par[name]) == (2,), name


def test_fixed_plume_flows_start_at_open_events():
    """In a "test-fixed" run, open_time records each openCauldron event, and each cauldron's
    plume flows are zero before its event and on at the vent-input plateau after it."""
    result = gv_main(n_cauldrons=2, G_n=[300, 400], plume_flow_mode="test-fixed")
    gv, data, events = result.gv, result.data, result.events
    q_s_peak = gv.params["Q_n"].values * (1 - gv.n_w) / gv.rho_p
    is_open = events["indexName"].values == "openCauldron"
    for ci in range(gv.n_cauldrons):
        t_open = events["t"].values[is_open & (events["cauldronIndex"].values == ci)][0]
        assert gv.open_time[ci] == t_open
        q_s = data["q_s_p_n"].sel(cauldron=ci)
        assert np.all(q_s.where(data.time < t_open, drop=True) == 0)
        plateau = (data.time > t_open + 600) & (data.time < t_open + 2400)
        assert plateau.any()
        assert np.allclose(q_s.where(plateau, drop=True), q_s_peak[ci])
