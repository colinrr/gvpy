import numpy as np
import pytest

from gvpy.glaciovolcano import glaciovolcano
from gvpy.ice_cauldron import IceCauldron


@pytest.mark.parametrize("ice_inflow_mode", ["off", "fixed-glen"])
def test_par_exposes_per_cauldron_fluxes(ice_inflow_mode):
    """Flux terms added to glaciovolcano's par are per cauldron, so
    gv_main stacks them as (time, cauldron)."""
    gv = IceCauldron(n_cauldrons=2, G_n=[300, 400], ice_inflow_mode=ice_inflow_mode)
    _, par = glaciovolcano(0.0, gv.get_initial_conditions(), gv, return_par=True)
    for name in ("Q_n", "u_ice_bar", "q_d_n", "q_c_n"):
        assert np.shape(par[name]) == (2,), name
