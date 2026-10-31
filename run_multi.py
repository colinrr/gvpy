"""Sample model run - currently just based on the simple cauldron growth
benchmark, to update as model evolves.

Translated from MATLAB source: run_single.m

TRANSLATION NOTE: the model run below calls gv_main (translated from gvMain.m), which
builds its own IceCauldron from the same parameter overrides as the MATLAB script's
gvIn struct (Python field names, e.g. nCauldrons -> n_cauldrons). The `c`/`c2`/`c3`
objects are just for inspecting an IceCauldron directly.

# Handy for fiddling
attrs.asdict(c, recurse=False),
"""

import time

from gvpy.gv_main import gv_main
from gvpy.ice_cauldron import IceCauldron
from gvpy.plotting import plot_results

# ---- Set up IceCauldron control params ----

run1 = {
    'alpha' : 1,  # Symettrical vertical/radial melt rate
    'n_cauldrons' : 1,  # Single cauldron
    'G_n' : 400,  # 400 m ice thickness
    'ice_inflow_mode' : "off",  # No ice flow
    'supraglacial_drainage_mode' : "off",
    'plume_flow_mode' : "off",
}

run2 = {
    'alpha' : 1,  # Symettrical vertical/radial melt rate
    'n_cauldrons' : 1,  # Single cauldron
    'G_n' : 400,  # 400 m ice thickness
    'ice_inflow_mode' : "off",  # No ice flow
    'supraglacial_drainage_mode' : "overflow",
    'plume_flow_mode' : "test-fixed",
}

# ---- Run the model ----
tic = time.perf_counter()
dat = [
    gv_main(**run1), 
    gv_main(**run2),
    ]
print(f"Elapsed time: {time.perf_counter() - tic:.2f} s")  # toc

# ---- Results plot ----
# TODO - model run dashboard plot
plot_results(dat, show_flow_modes=True).show()  # opens in a browser; blocks until stopped (Ctrl-C)


# gv_opts = {
#     'alpha' : 1,  # Symettrical vertical/radial melt rate
#     'n_cauldrons' : 1,  # Single cauldron
#     'G_n' : 400,  # 400 m ice thickness
#     'ice_inflow_mode' : "off",  # No ice flow
#     'supraglacial_drainage_mode' : "off",
#     'T_m' : -100,
#     }
# c3 = IceCauldron(**gv_opts)