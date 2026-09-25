"""ODE solver bookkeeping helpers shared by the right-hand side (glaciovolcano.py)
and the solver driver (gv_main.py) - not tied to a specific model component.

Translated from MATLAB source: utilities/getEventsTable.m and
utilities/assignIntegratedValues.m
"""

from typing import Optional

import numpy as np
import pandas as pd


def get_events_table(n_cauldrons: int) -> tuple[pd.DataFrame, list[str]]:
    """[indexTable, eventNames] = getEventsTable(nCauldrons)
    Helper function to tell which indices in ode solver events index correspond
    to which events, since they vary with number of cauldrons.

    Output:
      index_table:  names and cauldron number for all positional events
                    in the ode events function
      event_names : list of all possible event names
    """
    # TRANSLATION NOTE: MATLAB returns eventNames only when nargout==2; always returned here.
    # cauldronIndex is 0-based (project convention; MATLAB's is 1-based) and NaN for events that
    # aren't per-cauldron, and row positions are 0-based (MATLAB's `ie` event indices are
    # 1-based). MATLAB's unused locals (n_indices, dn) and the string RowNames are dropped.
    per_cauldron_events = ["openCauldron", "iceFreeCauldron", "activeFlood", "activePlume", "inactiveVent"]

    single_events = ["toeOutburst", "endFlood"]

    n_per_cauldron = len(per_cauldron_events)
    n_single = len(single_events)

    # Assign vectors of event names and cauldron indices
    index_name = np.concatenate([np.repeat(per_cauldron_events, n_cauldrons), single_events])
    cauldron_index = np.concatenate([np.tile(np.arange(n_cauldrons), n_per_cauldron), np.full(n_single, np.nan)])

    index_table = pd.DataFrame({"indexName": index_name, "cauldronIndex": cauldron_index})

    return index_table, per_cauldron_events + single_events


def assign_integrated_values(y: np.ndarray, gv, dat: Optional[dict] = None) -> dict:
    """dat = assignIntegratedValues(y,gv,solutions,dat)
    Get integrated value variables from ode y output. Generalized format to
    allow flexible development and addition of new variables.

    y         = ode solution vector (n_states,) for a single time step, or
                matrix (n_times, n_states) for a full solution
    gv        = IceCauldron object
    dat       = optional input dict in which to assign output variables
    """
    # TRANSLATION NOTE: MATLAB's `solutions` flag (choosing y(idx) vs y(:,idx)) is dropped - with
    # the cauldron-axis-last convention y[..., idx] handles both a single (n_states,) vector and a
    # (n_times, n_states) matrix, so the flag would be redundant. Returns plain NumPy arrays, not
    # xarray (also used inside the performance-critical right-hand side).
    if dat is None:
        dat = {}

    var_indices, _ = gv.get_solution_indices()

    for name, idx in var_indices.items():
        dat[name] = y[..., idx]

    return dat
