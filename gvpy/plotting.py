"""Plotting tools for model results (GVResult), using HoloViews/Panel with the Bokeh backend
(native zoom/pan; axes with the same dimension are linked).

Translated/adapted from MATLAB sources:
  plot_tools/getEventsCmap.m                             -> EVENT_COLORS
  plot_tools/plotVariableTimeSeries.m,
  plot_tools/plotSimulationTimeSeries.m                  -> plot_time_series
  plot_tools/plotCauldronGeometry.m                      -> plot_cauldron_geometry
  plot_tools/plot_layout.m                               -> plot_results
Variable labels/units (MATLAB's getVarLabels.m) are NOT defined here: they are read from each
variable's attrs, set from IceCauldron.VAR_INFO. Only display unit conversions live here.

Functions accept one GVResult, a list of them, or a {label: GVResult} dict, so runs made in the
same session (or loaded from a save) can be overlaid. In a notebook, call `pn.extension()` once and
display the returned object; from a script, call `.show()` on it.
"""

from collections.abc import Mapping, Sequence
import warnings

import holoviews as hv
from matplotlib.colors import to_hex
import numpy as np
import panel as pn
from scipy.integrate import cumulative_trapezoid
import xarray as xr

from .gv_main import GVResult
from .utilities.geometry import s_spheroid
from .utilities.solver_helpers import get_events_table

hv.extension("bokeh", logo=False)

Results = GVResult | Sequence[GVResult] | Mapping[str, GVResult]

# Display conversions, keyed by model units: {units: (display units, scale factor)}. Every unit
# plotted must be listed (identity entries included) - an unlisted unit warns and plots natively,
# flagging drift between model metadata and plotting. Override per call with display_units=.
DISPLAY_UNITS = {
    "s": ("hours", 1 / 3600),
    "m": ("m", 1),
    "m^3": ("km^3", 1e-9),
    "m/s": ("m/s", 1),
    "m^3/s": ("m^3/s", 1),
    "kg/s": ("kg/s", 1),
    "1": ("-", 1),
}

# ---------------------------------------------------------------------------
# Translated from plot_tools/getEventsCmap.m
# TRANSLATION NOTE: the first three rows are cubehelix(3,3,-0.45,1.56,0.93,[0.17 0.67],[0.14 0.73])
# evaluated once and hard-coded (cubehelix.m is not ported). Colours map to events by position in
# get_events_table's event list, as in MATLAB.
# ---------------------------------------------------------------------------
EVENTS_CMAP = [
    (0.0833, 0.1784, 0.3615),
    (0.0971, 0.5791, 0.5603),
    (0.4485, 0.8084, 0.5320),
    (0.85, 0.45, 0.1),
    (0.55, 0.3, 0.05),
    (0.1, 0.55, 0.85),
    (1, 1, 1),
]
EVENT_COLORS = {name: to_hex(rgb) for name, rgb in zip(get_events_table(1)[1], EVENTS_CMAP)}

# Line styling: colour steps through (variable, cauldron) within a panel; dash distinguishes runs
LINE_COLORS = [  # MATLAB default ColorOrder
    (0, 0.447, 0.741), (0.85, 0.325, 0.098), (0.929, 0.694, 0.125), (0.494, 0.184, 0.556),
    (0.466, 0.674, 0.188), (0.301, 0.745, 0.933), (0.635, 0.078, 0.184),
]
RUN_DASHES = ["solid", "dashed", "dashdot", "dotted"]


def _as_labelled_results(results: Results) -> dict:
    """Normalize a GVResult, list, or {label: GVResult} dict to a {label: GVResult} dict."""
    if isinstance(results, GVResult):
        return {"": results}
    if isinstance(results, Mapping):
        return dict(results)
    return {f"run {i}": r for i, r in enumerate(results)}


def _display(da: xr.DataArray, display_units: Mapping) -> tuple:
    """(scale factor, long_name, display units) for a variable, from its attrs and the display conversions."""
    if "long_name" not in da.attrs:
        raise ValueError(f"'{da.name}' has no long_name attr - add it to IceCauldron.VAR_INFO.")
    units = da.attrs.get("units")
    if units is None:
        return 1, da.attrs["long_name"], None
    if units not in display_units:
        warnings.warn(f"No display conversion for units '{units}' ('{da.name}') - plotting in model units.")
        return 1, da.attrs["long_name"], units
    display_unit, scale = display_units[units]
    return scale, da.attrs["long_name"], display_unit


def _axis_label(long_names: Sequence[str], display_unit: str | None) -> str:
    label = ", ".join(long_names)
    return f"{label} ({display_unit})" if display_unit else label


def _legend_text(da: xr.DataArray) -> str:
    return (da.attrs.get("symbol") or str(da.name)).replace("{", "").replace("}", "")


def _plot_data(result: GVResult) -> xr.Dataset:
    """result.data plus plotting-only derived variables."""
    data = result.data
    if all(q in data for q in ("q_s_n", "q_d_n", "q_c_n")):
        # Cumulative discharge: time integral of all outflow volume fluxes (display-only, not model output)
        q_out = data["q_s_n"] + data["q_d_n"] + data["q_c_n"]
        V_discharge_cum = cumulative_trapezoid(q_out.values, data["time"].values, axis=0, initial=0)
        attrs = {"units": "m^3", "long_name": "Cum. Discharge", "symbol": "V_{out}"}
        data = data.assign(V_discharge_cum=(q_out.dims, V_discharge_cum, attrs))
    return data


# Translated from plot_tools/plotVariableTimeSeries.m (with plotSimulationTimeSeries.m's labelling)
def plot_time_series(
    results: Results,
    var: str | Sequence[str],
    y_dim: hv.Dimension | None = None,
    display_units: Mapping | None = None,
) -> hv.Overlay:
    """Plot time series with events for output variable(s).

    One line per (run, cauldron, variable), plus event markers coloured by event type.
    var           = variable name, or list of names to share one axis
    y_dim         = optional shared y Dimension - plots with the same y Dimension name link y axes
    display_units = optional overrides/additions to DISPLAY_UNITS
    """
    # TRANSLATION NOTE: cauldron labels in legends are 0-based (the Python `cauldron` coordinate);
    # MATLAB's were 1-based. Event values: solution variables read the stored event state, derived
    # variables are interpolated (np.interp for interp1), as in MATLAB.
    display_units = {**DISPLAY_UNITS, **(display_units or {})}
    runs = _as_labelled_results(results)
    variables = [var] if isinstance(var, str) else list(var)
    if y_dim is None:
        data0 = _plot_data(next(iter(runs.values())))
        parts = [_display(data0[v], display_units) for v in variables]
        y_dim = hv.Dimension("+".join(variables), label=_axis_label([p[1] for p in parts], parts[0][2]))

    curves = []
    ev_t, ev_y, ev_name, ev_cauldron, ev_run = [], [], [], [], []
    for ri, (run, result) in enumerate(runs.items()):
        data = _plot_data(result)
        t_scale, t_name, t_unit = _display(data["time"], display_units)
        t_dim = hv.Dimension("time", label=_axis_label([t_name], t_unit))
        t = data["time"].values
        for vi, v in enumerate(variables):
            y_scale = _display(data[v], display_units)[0]
            y = data[v].values * y_scale
            per_cauldron = y.ndim == 2
            y = y if per_cauldron else y[:, None]
            for ci in range(y.shape[1]):
                label = f"{_legend_text(data[v])}_{ci}" if per_cauldron else _legend_text(data[v])
                label = f"{run}: {label}" if run else label
                color = to_hex(LINE_COLORS[(vi * y.shape[1] + ci) % len(LINE_COLORS)])
                curves.append(
                    hv.Curve((t * t_scale, y[:, ci]), t_dim, y_dim, label=label).opts(
                        color=color, line_dash=RUN_DASHES[ri % len(RUN_DASHES)]
                    )
                )
            if not per_cauldron:
                continue

            # Parse and plot values at events (per-cauldron events only)
            events = result.events
            for ei in range(events.sizes["event"]):
                ci = events["cauldronIndex"].values[ei]
                if np.isnan(ci):
                    continue
                ci = int(ci)
                te = events["t"].values[ei]
                if v in events:
                    ye = events[v].values[ei, ci] * y_scale
                else:
                    ye = np.interp(te, t, y[:, ci])
                ev_t.append(te * t_scale)
                ev_y.append(ye)
                ev_name.append(str(events["indexName"].values[ei]))
                ev_cauldron.append(ci)
                ev_run.append(run)

    markers = hv.Scatter(
        (ev_t, ev_y, ev_name, ev_cauldron, ev_run), t_dim, [y_dim, "event", "cauldron", "run"]
    ).opts(color="event", cmap=EVENT_COLORS, size=8, line_color="black", tools=["hover"], show_legend=False)

    return hv.Overlay(curves + [markers]).opts(
        legend_position="top_left", legend_opts={"click_policy": "hide", "label_text_font_size": "8pt"}
    )


def _rgba2rgb(rgb: Sequence[float], alphas: np.ndarray) -> list:
    """Colour at each alpha, blended over a white background."""
    # TRANSLATION NOTE: rgba2rgb.m is not in the katlaGV repo; implemented as alpha blending over white.
    rgb = np.asarray(rgb)
    return [to_hex(rgb * alpha + (1 - alpha)) for alpha in alphas]


# Translated from plot_tools/plotCauldronGeometry.m
def plot_cauldron_geometry(
    result: GVResult, c_idx: int = 0, steps: int | Sequence[int] = 4, display_units: Mapping | None = None
) -> hv.Overlay:
    """plotCauldronGeometry

    Cross-section snapshots of one cauldron: spheroidal roof while closed, rectangular outline
    (with remaining ice lid) while open, and walls only once ice free.
    steps = number of log-spaced snapshots, or explicit (0-based) time indices
    """
    display_units = {**DISPLAY_UNITS, **(display_units or {})}
    data = result.data
    t = data["time"].values
    G = result.gv.params["G_n"].values[c_idx]
    L = result.gv.params["L_n"].values[c_idx]
    z_b = result.gv.params["z_b_n"].values[c_idx]
    a_n = data["a_n"].values[:, c_idx]
    c_n = data["c_n"].values[:, c_idx]
    V_ice_n = data["V_ice_n"].values[:, c_idx]

    if np.isscalar(steps):
        # Up to length of t, but can stop before for shorter cauldron evolutions
        steps = np.round(np.logspace(0, np.log10(len(t)), steps + 1)).astype(int)
        steps = steps[1:] - 1  # TRANSLATION NOTE: MATLAB 1-based indices -> 0-based
    steps = np.asarray(steps)
    t_stops = t[steps]

    # Find event times
    # TRANSLATION NOTE: MATLAB errors if the cauldron never opens / goes ice free (empty event
    # time); a missing event is treated here as never happening (t = inf).
    events = result.events
    names = events["indexName"].values.astype(str)
    cauldron = events["cauldronIndex"].values
    t_open = events["t"].values[(names == "openCauldron") & (cauldron == c_idx)]
    t_if = events["t"].values[(names == "iceFreeCauldron") & (cauldron == c_idx)]
    t_open = t_open[0] if t_open.size else np.inf
    t_if = t_if[0] if t_if.size else np.inf

    closed_steps = steps[t_stops < t_open]
    open_steps = steps[(t_stops >= t_open) & (t_stops < t_if)]
    if_steps = steps[t_stops >= t_if]

    # Construct colormap
    alpha_range = [0.5, 1]
    closed_cmap = _rgba2rgb(LINE_COLORS[0], np.linspace(*alpha_range, len(closed_steps)))
    open_cmap = _rgba2rgb(EVENTS_CMAP[0], np.linspace(*alpha_range, len(open_steps)))
    if_cmap = _rgba2rgb(EVENTS_CMAP[1], np.linspace(*alpha_range, len(if_steps)))
    cmap = closed_cmap + open_cmap + if_cmap
    lw = 1.5

    r_dim = hv.Dimension("radial_distance", label="Radial distance (m)")
    z_dim = hv.Dimension("elevation", label="Elevation (m)")

    # Plot bedrock and ice surfaces
    xmax = np.max(a_n)
    x_surf = xmax * np.array([-1.2, 1.2])
    h = [
        hv.Curve((x_surf, [z_b, z_b]), r_dim, z_dim).opts(color="black", line_width=lw),
        hv.Curve((x_surf, [G + z_b, G + z_b]), r_dim, z_dim).opts(color="#808080", line_width=lw),
    ]

    outlines = []
    # ---- CLOSED cauldron ----
    if closed_steps.size:
        # Get cavity roof curve x,y
        _, _, x_sphd, z_sphd = s_spheroid(a_n[closed_steps], c_n[closed_steps], return_xz=True)
        z_sphd = z_sphd + z_b
        outlines += [(x_sphd[:, i], z_sphd[:, i]) for i in range(len(closed_steps))]
    # ---- OPEN cauldron ----
    for si in open_steps:
        x_rect = np.array([-1, -1, np.nan, 1, 1, np.nan, -1, 1, 1, -1, -1]) * a_n[si]

        c_over_G = c_n[si] / G
        hi_over_G = V_ice_n[si] / (2 * L * a_n[si]) / G
        y_rect = np.array([0, 1, np.nan, 1, 0, np.nan, 1, 1, 1, 1, 1], dtype=float)
        y_rect[[6, 7, 10]] = c_over_G
        y_rect[8:10] = c_over_G + hi_over_G
        y_rect = y_rect * G + z_b
        outlines.append((x_rect, y_rect))
    # ---- ICE FREE Cauldron ----
    for si in if_steps:
        x_rect = np.array([-1, -1, np.nan, 1, 1]) * a_n[si]
        y_rect = np.array([0, 1, np.nan, 1, 0]) * G + z_b
        outlines.append((x_rect, y_rect))

    # Wrapup
    t_scale, t_name, t_unit = _display(data["time"], display_units)
    labels = [f"{t_stop * t_scale:.2f}" for t_stop in t_stops]
    h += [
        hv.Curve(xy, r_dim, z_dim, label=label).opts(color=color, line_width=lw)
        for xy, label, color in zip(outlines, labels, cmap)
    ]
    return hv.Overlay(h).opts(
        data_aspect=1,
        ylim=(-0.05 * G + z_b, 1.05 * G + z_b),
        show_grid=True,
        legend_position="right",
        legend_opts={"title": _axis_label([t_name], t_unit), "label_text_font_size": "8pt"},
    )


# ---------------------------------------------------------------------------
# Translated from plot_tools/plot_layout.m
# ---------------------------------------------------------------------------

# Axes Names
NAMES_LEFT = [
    "Cumulative Volume",
    "Component Volume",
    "Sheet Volume",
    "Cumulative Discharge",
    "Cauldron Dimensions",
    "Cauldron Geometry",
    "Melting Rates",
    "Heating Efficiency",
]

NAMES_RIGHT = [
    "Conduit Flux",
    "Ice Flux",
    "Plume Flux",
    "Intracauldron Flux",
    "Supraglacial Flux",
    "Subglacial Flux",
]

# Variables per panel (MATLAB's var2ax, completed for every panel). None = no model output yet.
PANEL_VARS = {
    "Cumulative Volume": ["V_cavity_n", "V_cum"],
    "Component Volume": ["V_ice_n", "V_w_n", "V_p_n"],
    "Sheet Volume": None,  # no subglacial sheet variable yet
    "Cumulative Discharge": ["V_discharge_cum"],
    "Cauldron Dimensions": ["a_n", "c_n"],
    "Cauldron Geometry": "geometry",
    "Melting Rates": ["horizontalMeltingRate_n", "verticalMeltingRate_n"],
    "Heating Efficiency": ["f_i_n"],
    "Conduit Flux": ["Q_n"],
    "Ice Flux": ["u_ice_bar"],
    "Plume Flux": None,  # plume fluxes are not computed in the ODE right-hand side
    "Intracauldron Flux": ["q_c_n"],
    "Supraglacial Flux": ["q_s_n"],
    "Subglacial Flux": ["q_d_n"],
}

# Linking Groups: the volume panels share one y dimension (linked x and y); every time-series
# panel shares the time dimension (linked x).
# TRANSLATION NOTE: MATLAB linked the volume panels' x separately from the other time panels
# (two linkaxes calls); here all time axes are linked together.
LINK_GROUP_XY = NAMES_LEFT[0:4]

# Formatting Groups (Which axes should have X-labels hidden?)
HIDE_X_LABELS = [NAMES_LEFT[i] for i in (0, 2, 4, 6)] + NAMES_RIGHT[:-1]


def _placeholder(name: str, reason: str) -> pn.pane.Markdown:
    return pn.pane.Markdown(f"**{name}**\n\n*{reason}*", styles={"border": "1px solid #ccc"})


def _panel(name: str, runs: dict, hide_x: bool, display_units: Mapping) -> pn.viewable.Viewable:
    """One named dashboard panel."""
    panel_vars = PANEL_VARS[name]
    if panel_vars is None:
        return _placeholder(name, "not yet output by model")
    if panel_vars == "geometry":
        # TRANSLATION NOTE: geometry is drawn for the first run and cauldron 0 only
        plot = plot_cauldron_geometry(next(iter(runs.values())), display_units=display_units)
    else:
        missing = [v for v in panel_vars for r in runs.values() if v not in _plot_data(r)]
        if missing:
            warnings.warn(f"Panel '{name}': variables {sorted(set(missing))} not in results - showing placeholder.")
            return _placeholder(name, f"missing from results: {', '.join(sorted(set(missing)))}")
        y_dim = None
        if name in LINK_GROUP_XY:
            display_unit = _display(_plot_data(next(iter(runs.values())))[panel_vars[0]], display_units)[2]
            y_dim = hv.Dimension("volume", label=_axis_label(["Volume"], display_unit))
        plot = plot_time_series(runs, panel_vars, y_dim=y_dim, display_units=display_units)
        if hide_x:
            plot = plot.opts(xaxis="bare")
    return pn.pane.HoloViews(plot.opts(title=name, responsive=True, show_grid=True), sizing_mode="stretch_both")


def plot_results(
    results: Results,
    panels: Sequence[str] | None = None,
    display_units: Mapping | None = None,
    ncols: int = 3,
    height: int = 1000,
) -> pn.viewable.Viewable:
    """Results dashboard for one or more model runs.

    panels        = None for the full layout (from plot_layout.m: a 2x2 grid of stacked panel
                    pairs, beside a column of 6 flux panels), or a list of panel names (see
                    PANEL_VARS) to show only those, in order, in an ncols-wide grid.
    display_units = optional overrides/additions to DISPLAY_UNITS, e.g. {"s": ("min", 1/60)}
    Returns a Panel layout - display it in a notebook, or call .show() from a script.
    """
    display_units = {**DISPLAY_UNITS, **(display_units or {})}
    runs = _as_labelled_results(results)
    event_legend = pn.pane.HTML(
        "<b>Events:</b> "
        + " ".join(
            f"<span style='color:{c}; -webkit-text-stroke: 0.5px black'>&#9679;</span> {n}"
            for n, c in EVENT_COLORS.items()
        )
    )

    if panels is not None:
        items = [_panel(p, runs, hide_x=False, display_units=display_units) for p in panels]
        for item in items:
            item.height = height // 3
        grid = pn.GridBox(*items, ncols=ncols, sizing_mode="stretch_width")
        return pn.Column(event_legend, grid, sizing_mode="stretch_width")

    # Left: pairs of panels stacked in a 2x2 grid (4 x 2 panels, each 3 rows tall).
    # Right: 6 stacked panels (each 2 rows tall). Three equal-width columns.
    grid = pn.GridSpec(sizing_mode="stretch_width", height=height)
    for i, name in enumerate(NAMES_LEFT):
        pair, j = divmod(i, 2)
        row = 6 * (pair // 2) + 3 * j
        col = pair % 2
        grid[row : row + 3, col] = _panel(name, runs, name in HIDE_X_LABELS, display_units)
    for i, name in enumerate(NAMES_RIGHT):
        grid[2 * i : 2 * i + 2, 2] = _panel(name, runs, name in HIDE_X_LABELS, display_units)

    return pn.Column(event_legend, grid, sizing_mode="stretch_width")
