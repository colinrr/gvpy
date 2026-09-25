"""Solver driver: builds the IceCauldron, integrates the ODEs in glaciovolcano.py through
terminal events (cauldron opening / going ice-free), and assembles the run output.

Translated from MATLAB source: gvMain.m (including its local functions gvEvents, parseEvents,
checkTerminalConditions).
"""

import copy
from typing import Mapping, Optional

import attrs
import numpy as np
import pandas as pd
import xarray as xr
from scipy.integrate import solve_ivp

from .glaciovolcano import glaciovolcano
from .ice_cauldron import IceCauldron
from .utilities.solver_helpers import assign_integrated_values, get_events_table


@attrs.define
class GVResult:
    """Output of gv_main (MATLAB's `dat` struct).

    gv     : the IceCauldron in its final state (dat.GV)
    data   : xarray Dataset with dims (time, cauldron) - integrated and derived time series
    events : xarray Dataset with dim (event) - one entry per solver event (dat.events)
    """

    gv: IceCauldron
    data: xr.Dataset
    events: xr.Dataset


def gv_main(gv_in: Optional[Mapping] = None, **kwargs) -> GVResult:
    """dat = gvMain(varargin)
    COLIN ROWELL, MEGHAN SHARP, MARK JELLINEK, 2024
    Solving coupled ODEs in galciovolcano.m
    All input as NAME/VALUE pairs as shown in INPUT DEFAULTS section below
    (OR struct of NAME/VALUE pairs)


    rebuild 06/2024 -- CR

    revised 8/2022 -- MS,CR. Updated to correct parameters for final sims.
           Packaged code for final simulation runs.
                   - see original submission version (BVsubmit branch) for
                     previous versions
    """
    # TRANSLATION NOTE: MATLAB's pathConfig(PROJ_DIR) (MATLAB path setup), tic/toc and the
    # unused `global step_ct` are dropped. Inputs are a dict and/or keyword arguments using the
    # Python IceCauldron field names (MATLAB nCauldrons -> n_cauldrons, etc.).
    # SOLVER NOTE: ode15s -> scipy.integrate.solve_ivp(method="BDF") (closest stiff variable-order
    # method); rtol = 1e-7 as in MATLAB, atol left at SciPy's default (1e-6, same as MATLAB's
    # AbsTol default). MATLAB's odeset 'NonNegative' (all states) has no SciPy equivalent and is
    # NOT replicated.

    # %% INPUT DEFAULTS
    # default = getGVdefaults;

    # %% PARSE INPUT
    # Initialize cauldron class
    gv = IceCauldron(**{**(gv_in or {}), **kwargs})

    # %% GET CONSTANTS

    # ---- Calc constants ----
    # gv.chi = const.c_m .* (C.T_m - C.T_ice) ./ (C.rho_ice .* C.l_ice);
    #  --> gv = gv.get_chi

    # ---- Get Dimensional scales ----
    # tau_melt, V_melt = getDimScales(gv);
    # --> gv = gv.getDimScales

    # %% BUILD GEOMETRY
    # Calculate gradients, thicknesses, etc
    # Setup drainage conductances?

    # %% LOAD PLUME EMULATOR

    # %% SETUP INITIAL CONDITIONS and ODE options

    y0 = gv.get_initial_conditions()
    timespan = np.linspace(0, gv.t_final, gv.timesteps)
    n_states = len(y0)

    refine = 4
    event_funcs = _make_event_functions(gv, y0)

    t_out = [timespan[0]]
    y_out = [y0]
    te_out = []
    ye_out = []
    ie_out = []
    # %% RUN ODE SOLVER
    end_sim = False
    current_time_vector = timespan
    first_step = None
    while not end_sim:
        sol = solve_ivp(
            lambda t, y: glaciovolcano(t, y, gv),
            (current_time_vector[0], current_time_vector[-1]),
            y0,
            method="BDF",
            t_eval=current_time_vector,
            events=event_funcs,
            rtol=1e-7,
            first_step=first_step,
        )

        t = sol.t
        y = sol.y.T

        # Gather events in order of occurrence (MATLAB's te, ye, ie; ie is 0-based here)
        found = sorted(
            (te_i, ye_i, i) for i, (tes, yes) in enumerate(zip(sol.t_events, sol.y_events)) for te_i, ye_i in zip(tes, yes)
        )
        te = np.array([f[0] for f in found])
        ye = np.array([f[1] for f in found]).reshape(-1, n_states)
        ie = np.array([f[2] for f in found], dtype=int)

        # TRANSLATION NOTE: like ode15s, the terminal event point is appended as the last row of
        # the segment's (t, y) (solve_ivp with t_eval only returns the t_eval grid points).
        if len(te):
            t = np.append(t, te[-1])
            y = np.vstack([y, ye[-1]])

        nt = len(t)
        t_out += list(t[1:nt])
        y_out += list(y[1:nt])
        te_out += list(te)
        ye_out += list(ye)
        ie_out += list(ie)

        if np.all(t < gv.t_final):
            try:
                assert np.all(ye - y[-1, :] == 0), "Event solutions do not match ode solution end points"  # A check for now
            except:  # noqa: E722 - TRANSLATION NOTE: MATLAB's bare `catch ME` (unused ME) catches everything; preserved as a bare except.
                # TRANSLATION NOTE: MATLAB's catch body is literally `faafo` - undefined-function
                # placeholder dev text, preserved as-is per CLAUDE.md's in-development-components
                # policy (would raise NameError if the assert ever failed).
                faafo()
        gv, y0 = parse_events(te, ye, ie, gv)

        # Check for simulation end
        end_sim = check_terminal_conditions(t[-1], gv)

        # Reset for continuing
        if not end_sim:
            #         y0 = y(end,:)';
            #         options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
            #           'MaxStep',t(nt)-t(1), 'Events',@(t,y) gvEvents(t,y,gv));
            # TRANSLATION NOTE: MATLAB's t(nt-refine) errors when fewer than refine+1 points precede
            # the event; Python negative indexing would silently wrap, so guard explicitly. The
            # step is also capped at the remaining span (solve_ivp raises if first_step exceeds it;
            # ode15s just clips it). MATLAB's `options` Events re-creation is unnecessary here:
            # event functions read `gv` (mutated in place by parse_events) at call time.
            if nt - 1 - refine < 0:
                raise IndexError(f"Fewer than {refine + 1} output points before the event; cannot set InitialStep.")
            # TRANSLATION NOTE (deliberate change from MATLAB): MATLAB continues on
            # timespan(timespan>=t(end)), i.e. from the next output-grid point, while y0 is the
            # state at the event time - that freezes the state for up to one output step and delays
            # every later event (probe: 3e-3 RMSE vs 3e-11 when restarted at the event time). Restart
            # exactly at the event time instead; the grid points after it are kept.
            current_time_vector = np.concatenate([[t[-1]], timespan[timespan > t[-1]]])
            first_step = min(t[nt - 1] - t[nt - 1 - refine], current_time_vector[-1] - current_time_vector[0])

    # %%  MAIN SOLVER OUTPUT
    t_out = np.array(t_out)
    y_out = np.array(y_out)
    te_out = np.array(te_out)
    ye_out = np.array(ye_out).reshape(-1, n_states)
    ie_out = np.array(ie_out, dtype=int)

    dat = assign_integrated_values(y_out, gv)
    events, _ = parse_events(te_out, ye_out, ie_out, gv, True)

    # dat.te = teout;
    # dat.ye = yeout;
    # dat.ie = ieout;
    # %% ADDITIONAL CALCS FOR RECONSTRUCTED VARS

    dat["V_w_plus_p"] = dat["V_w_n"] + dat["V_p_n"]
    dat["V_cum"] = dat["V_w_n"] + dat["V_p_n"] + dat["V_ice_n"]

    # if gv.nCauldrons>1                                              % Melting velocity
    #     [~,dat.horizontalMeltingRate_n] = gradient(dat.a_n,1,dat.t);
    #     [~,dat.verticalMeltingRate_n] = gradient(dat.c_n,1,dat.t);
    # else
    #     dat.horizontalMeltingRate_n = gradient(dat.a_n,dat.t);
    #     dat.verticalMeltingRate_n = gradient(dat.c_n,dat.t);
    # end
    dat["f_i_n"] = gv.get_f_i(dat["a_n"], dat["c_n"], t_out, events)

    # Reconstruct event times for da_dt
    # TRANSLATION NOTE: MATLAB grows openTime/iceFreeTime by assignment, so a cauldron without the
    # event is silently zero-filled ("open/ice free from t = 0") - an accident, translated
    # literally as zeros. Not reachable in practice: get_f_i above already fails unless every
    # cauldron went ice-free, and a cauldron must open before it can go ice-free.
    ev_name = events["indexName"].values
    ev_ci = events["cauldronIndex"].values
    ev_t = events["t"].values
    open_time = np.zeros(gv.n_cauldrons)
    evoc = ev_name == "openCauldron"
    open_time[ev_ci[evoc].astype(int)] = ev_t[evoc]
    ice_free_time = np.zeros(gv.n_cauldrons)
    evif = ev_name == "iceFreeCauldron"
    ice_free_time[ev_ci[evif].astype(int)] = ev_t[evif]

    open_idx = np.zeros((len(t_out), gv.n_cauldrons), dtype=bool)
    ice_free_idx = np.zeros((len(t_out), gv.n_cauldrons), dtype=bool)

    for ci in range(gv.n_cauldrons):
        open_idx[:, ci] = t_out >= open_time[ci]
        ice_free_idx[:, ci] = t_out >= ice_free_time[ci]

    dat["horizontalMeltingRate_n"], dat["verticalMeltingRate_n"] = gv.get_u_melt(
        dat["a_n"], dat["c_n"], dat["f_i_n"], open_idx, ice_free_idx
    )

    # V_cavity is identically the CV volume MINUS ice volume
    dat["V_CV_n"] = gv.get_cv_control_volume(dat["a_n"])
    # dat.V_cavity_n = pi/2 .* dat.a_n .* dat.c_n .* gv.L_n;  % Cavity volume - not open
    # V_cav_open = 2*dat.a_n .* dat.c_n .* gv.L_n - dat.V_ice_n;
    # dat.V_cavity_n(openIdx) = V_cav_open(openIdx);
    dat["V_cavity_n"] = dat["V_CV_n"] - dat["V_ice_n"]
    dat["H_i_n"], dat["H_w_n"], dat["H_p_n"], dat["H_cum_n"] = gv.get_material_heights(
        dat["V_ice_n"], dat["V_w_n"], dat["V_p_n"], dat["V_cavity_n"], dat["a_n"], dat["c_n"]
    )

    # Gradient fields best calculated inside the solver function
    # TRANSLATION NOTE: MATLAB's IceCauldron is a value class (this copy is independent); the
    # attrs IceCauldron is mutable, so an explicit shallow copy is made - only open_cauldron /
    # ice_free_cauldron are reassigned on it, never mutated in place.
    gv_temp = copy.copy(gv)
    par_series = {}
    for ii in range(len(t_out)):
        # TODO: fix this to accommodate Open/IceFree condition changes
        gv_temp.open_cauldron = open_idx[ii, :]
        gv_temp.ice_free_cauldron = ice_free_idx[ii, :]
        _, gv_out = glaciovolcano(t_out[ii], y_out[ii, :], gv_temp, return_par=True)
        for name, value in gv_out.items():
            if ii == 0:
                par_series[name] = np.zeros((len(t_out),) + np.shape(value))
            par_series[name][ii] = value
    dat.update(par_series)

    # Quick internal check
    # error('TODO: YOU NEED TO PLACE LIMITS ON c==G_n')
    assert set(gv.derived_vars()) <= set(dat), "Not all derived values were calculated."

    data = xr.Dataset(
        {name: (("time", "cauldron")[: np.ndim(arr)], np.asarray(arr)) for name, arr in dat.items()},
        coords={"time": t_out, "cauldron": np.arange(gv.n_cauldrons)},
    )

    return GVResult(gv=gv, data=data, events=events)


def _make_event_functions(gv: IceCauldron, y0: np.ndarray) -> list:
    """Wrap the vector-valued gvEvents as one scalar callable per event, as solve_ivp requires.
    isterm/direction are read from a probe call at y0 (they are constant in gv_events).
    """
    index_table, _ = get_events_table(gv.n_cauldrons)
    _, is_term, direction = gv_events(0.0, y0, gv, index_table)

    funcs = []
    for k in range(len(is_term)):

        def event(t, y, k=k):
            return gv_events(t, y, gv, index_table)[0][k]

        event.terminal = bool(is_term[k])
        event.direction = float(direction[k])
        funcs.append(event)
    return funcs


def gv_events(t: float, y: np.ndarray, gv: IceCauldron, index_table: pd.DataFrame) -> tuple:
    """[value,isterm,direction] = gvEvents(t,y,gv)"""
    # TRANSLATION NOTE: the events index table is built once by the caller (MATLAB rebuilds it on
    # every call) and passed in, since building a DataFrame per solver step is needlessly slow.
    sol_indices, _ = gv.get_solution_indices()

    # hi = gv.G_n' - y(indices.c_n);

    # Events: (TEMPORARY, IN-DEV)
    # To track:
    # cauldron opening
    # flood initiates
    # plume initiates
    # eruption ends
    # glacier toe outburst
    # flood ends (in cauldron & from toe)

    # TRANSLATION NOTE: MATLAB grows value/isterm/direction implicitly to the highest assigned
    # index; only the first two per-cauldron event types (open, ice-free) are implemented, so the
    # arrays have length 2 * n_cauldrons, matching the event index table rows 0..2n-1.
    n_events = 2 * gv.n_cauldrons
    value = np.zeros(n_events)
    isterm = np.zeros(n_events)
    direction = np.zeros(n_events)

    # ----- Non-terminal conditions ----

    # ---- CURRENT TERMINAL CONDITIONS ----
    # -> Check ice-free cauldron conditions
    cauldron_idx = np.flatnonzero(index_table["indexName"] == "iceFreeCauldron")
    value[cauldron_idx], isterm[cauldron_idx], direction[cauldron_idx] = gv.check_ice_free_cauldron_conditions(
        y[sol_indices["V_ice_n"]], y[sol_indices["a_n"]]
    )

    # -> slowest cauldron melts through: c_n = G_n
    cauldron_idx = np.flatnonzero(index_table["indexName"] == "openCauldron")
    value[cauldron_idx], isterm[cauldron_idx], direction[cauldron_idx] = gv.check_open_cauldron_conditions(
        y[sol_indices["c_n"]]
    )

    # value(3)        = value(gv.slowestCauldronIndex);
    # isterm(3)       = 1;
    # direction(3)    = -1;

    # Events: (TEMPORARY, IN-DEV)
    #   (1) "open cauldron": c_n = G_n
    #   (2) V_ice = 0 (terminal for now)
    #   (3) stop draining if t>1 and hw=0
    # Can include the first two later but functionality is fine for now

    return value, isterm, direction


def parse_events(te: np.ndarray, ye: np.ndarray, ie: np.ndarray, gv: IceCauldron, as_output: bool = False) -> tuple:
    """[out, new_y0] = parseEvents(te, ye, ie, gv, asOutput)
    asOutput = true | false [default]
      -> When 'true', will parse the output a structure to parse all events for
      the final simulation output.
      -> When false, updates simulation parameters as necessary for each
      event type.

    Returns (events Dataset, None) when as_output, else (gv, new_y0).
    """
    # Parse ode events to determine changes

    # Get event index reference table
    index_table, event_names = get_events_table(gv.n_cauldrons)
    sol_ind, _ = gv.get_solution_indices()

    if as_output:
        # Outputs after simulation finish - build events table
        # TRANSLATION NOTE: an xarray Dataset (dim `event`) instead of a MATLAB table; the
        # per-cauldron solution values of ALL cauldrons are saved for each event (dims
        # event x cauldron) - use `cauldronIndex` to pull the value pertaining to a specific one.
        ev_table = index_table.iloc[ie]
        data = {
            "indexName": ("event", ev_table["indexName"].to_numpy()),
            "cauldronIndex": ("event", ev_table["cauldronIndex"].to_numpy()),
            "t": ("event", np.asarray(te, dtype=float)),
        }
        for name, idx in sol_ind.items():
            # This approach saves all per_cauldron values for each event.
            # Must use the 'cauldronIndex' field to pull the value
            # pertaining to a specific cauldron.
            data[name] = (("event", "cauldron")[: ye[:, idx].ndim], ye[:, idx])

            # An alternative approach saves only the value of cauldron
            # associated with the event, but this makes non-per-cauldron
            # events a bit more cumbersome.
            #             for ii = 1:length(ie)
            #                 evTable.(sol_vars{fi}) = ye( solInd.(sol_vars{fi})(evTable.cauldronIndex) );
            #             end

        return xr.Dataset(data, coords={"cauldron": np.arange(gv.n_cauldrons)}), None

    # TRANSLATION NOTE: simultaneous events make y0 a matrix in MATLAB (an error on the next
    # solve); kept as an error, raised here explicitly. Todo: decide whether to handle properly.
    if len(ie) > 1:
        raise NotImplementedError("Simultaneous solver events are not supported.")

    new_y0 = ye[0].copy() if len(ie) == 1 else ye.copy()
    # Updates during solution
    for ii in range(len(ie)):
        event_name = index_table["indexName"].iloc[ie[ii]]
        which_caul = int(index_table["cauldronIndex"].iloc[ie[ii]])
        if event_name == "openCauldron":
            # Specify which cauldron has opened
            gv.open_cauldron[which_caul] = True
            # Recalculate new c_n for rectangular geometry (assuming step change)
            # TRANSLATION NOTE: MATLAB reads gv.L_n(whichCaul)/gv.G_n(whichCaul); here the
            # broadcast per-cauldron values in gv.params are used (raw gv.L_n/gv.G_n may be scalars).
            new_y0[sol_ind["c_n"].start + which_caul] = min(
                (ye[0, sol_ind["V_w_n"].start + which_caul] + ye[0, sol_ind["V_p_n"].start + which_caul])
                / (ye[0, sol_ind["a_n"].start + which_caul] * gv.params["L_n"].values[which_caul]),
                gv.params["G_n"].values[which_caul],
            )

            #                     new_y0(solInd.c_n) = min( [(ye(solInd.V_w_n) + ye(solInd.V_p_n)) ./ (ye(solInd.a_n) .* gv.L_n) ; gv.G_n],[],1);

        elif event_name == "iceFreeCauldron":
            gv.ice_free_cauldron[which_caul] = True
        #             case "activeFlood"
        #             case "activePlume"
        #             case "inactiveVent"
        #             case "toeOutburst"
        #             case "endFlood"
    return gv, new_y0


def check_terminal_conditions(t: float, gv: IceCauldron) -> bool:
    """endSim = checkTerminalConditions(t,gv)
    Check conditions to fully end the simulation
    """
    # temporary
    #     all_open = all(gv.openCauldron); % open cauldrons
    all_ice_free = np.all(gv.ice_free_cauldron)  # Ice free cauldrons
    # todo: eruption ends + flood ends (target)

    # End if any
    end_if_any = [
        t >= gv.t_final,  # Max time elapsed
        all_ice_free,
    ]

    return bool(np.any(end_if_any))
