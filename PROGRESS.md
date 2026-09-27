# Translation Progress

Snapshot of MATLAB -> Python translation status. This is a manually-maintained
snapshot, not derived automatically - update it as translation work happens,
and treat it as advisory (check current file contents/`git log` if in doubt).

See also `TRANSLATION_NOTES.md` for known bugs, behavioral changes, and
anticipated design decisions arising from this translation work.

Last updated: 2026-09-27

## Translated

- `gvpy/utilities/math_helpers.py` - general math/geometry helpers
- `gvpy/utilities/melting.py`
- `gvpy/utilities/constants.py`
- `gvpy/utilities/geometry.py` - now also `s_spheroid` (from `S_spheroid.m`, with
  user-approved edge-case corrections - see `TRANSLATION_NOTES.md`) and
  `ellipse_perimeter` (Ramanujan approximation shared by `s_spheroid` and
  `get_elliptical_cylinder_sa`).
- `gvpy/utilities/validators.py` - property-value validators (translated from
  IceCauldron.m's local validators plus MATLAB's built-in mustBeNonnegative/
  mustBePositive), wired into `IceCauldron` via `attrs.field(validator=...)`.
- `gvpy/plume.py` - plume model component (standalone functions), translated
  from `getPlumeFluxes.m`'s embedded local functions (`predictForest`,
  `predictTree`).
- `gvpy/subglacial.py` - subglacial drainage model component (standalone
  functions), translated from `get_L_lambda.m`'s and `solveDeltaP.m`'s
  embedded local functions. **Entire component is in-development/test-only**
  per project direction - see `TRANSLATION_NOTES.md`.
- `gvpy/ice_cauldron.py` - `IceCauldron` class, now including:
  - Constructor, dimensional scales, initial conditions, ODE event
    conditions, control-volume geometry, solution-variable bookkeeping,
    `get_elevation_profiles` (from `IceCauldron.m` itself).
  - Vent/melting & ice inflow: `get_u_ice`, `get_f_i`, `get_u_melt`,
    `get_da_dt`, `get_material_heights`.
  - Plume: `get_plume_fluxes`.
  - Supraglacial drainage: `get_q_s`. (No standalone helpers needed, so no
    `gvpy/supraglacial.py` was created - that pattern is reserved for
    components with local-function dependencies, per `plume.py`/`subglacial.py`.)
  - Subglacial drainage (in development): `get_fluxural_rigidity`,
    `get_l_lambda`, `get_u_tip`, `get_del_h`, `drainage_density`,
    `get_drainage_fluxes` (dummy stub, unimplemented in MATLAB too),
    `solve_delta_p`.
- `gvpy/glaciovolcano.py` - the ODE right-hand side (`glaciovolcano.m`);
  `glaciovolcano(t, y, gv, return_par=False)`.
- `gvpy/gv_main.py` - the solver driver (`gvMain.m` and its local functions):
  `gv_main(gv_in=None, **kwargs) -> GVResult`, `gv_events`, `parse_events`,
  `check_terminal_conditions`. Uses `scipy.integrate.solve_ivp` (BDF) with
  terminal events and restarts. `GVResult` holds `.gv`, `.data` (xarray
  Dataset, dims time x cauldron) and `.events` (xarray Dataset, dim event).
- `gvpy/utilities/solver_helpers.py` - `get_events_table`,
  `assign_integrated_values` (from `getEventsTable.m` /
  `assignIntegratedValues.m`).
- `gvpy/plotting.py` - results dashboard and plot functions (Panel + HoloViews,
  Bokeh backend), from `plot_tools/plot_layout.m`, `plotVariableTimeSeries.m`/
  `plotSimulationTimeSeries.m`, `plotCauldronGeometry.m` and `getEventsCmap.m`.
  `plot_results(results, panels=None, display_units=None)` takes one `GVResult`,
  a list, or a `{label: GVResult}` dict and returns a Panel layout: the 14
  panels from `plot_layout.m` (placeholders for Sheet Volume and Plume Flux,
  which have no model output yet), with all time axes linked and the volume
  panels also linked in y. Works inline in notebooks (`pn.extension()`) or via
  `.show()`. Holds only display-unit conversions (`DISPLAY_UNITS`, keyed by
  unit); labels and units come from the variables' attrs.
- **`getVarLabels.m` -> `IceCauldron.VAR_INFO`** (decided this session,
  superseding the plan to put it in `plotting.py`): one registry of units,
  `long_name` and `symbol` for every `IceCauldron` field, computed `params`
  entry, solver output and event variable. `add_var_attrs` attaches these as
  xarray attrs on `gv.params` and on `gv_main`'s `data`/`events`; a missing
  entry raises `MissingVarInfoWarning` at runtime and fails tests (pytest
  `filterwarnings` in `pyproject.toml`). Units are physical and pint-parseable;
  `add_var_attrs` is the single place to change for non-dimensionalization.
- `run_single.py` - now calls `gv_main` with the same overrides as
  `run_single.m` (plus scratch `IceCauldron` objects `c`, `c2`, `c3` for
  inspection; `c3` deliberately/accidentally fails the `T_m` validator), then
  `plot_results(dat).show()` (blocks until Ctrl-C).
- Other `ice_cauldron.py` / solver changes this session:
  - New field `glen_n` (Glen's flow-law exponent, default 3); `A`'s units in
    `VAR_INFO` follow it (`Pa^(-{glen_n}) s^-1`).
  - `glaciovolcano`'s `par` also returns `Q_n`, `u_ice_bar`, `q_d_n`, `q_c_n`
    (for the flux panels); `get_u_ice` returns per-cauldron zeros in "off" mode.
  - `get_material_heights` returns zero fill heights when there is no cavity
    (`c_n` or `V_cavity_n` ~ 0 via `np.isclose`), warning if volumes are
    inconsistent beyond `1e-6 * V_CV_n`; `get_u_melt` evaluates its ice-free
    division only on ice-free elements. Full runs are now warning-free.
  - Computed (`init=False`) fields that `__attrs_post_init__` always sets are
    non-optional with no default (`params`, `chi`, `tau_ND`, `E_prime`, ...);
    `t_Qdecay` and `plume_emulator` stay optional. `allow_none` removed.
- Type hints throughout `gvpy/` (every function; `X | None` syntax).
- Dependencies: `panel`, `hvplot`, `holoviews` added to `pyproject.toml`
  (lock file updated by the user).
- Shape convention adopted with the solver: the cauldron axis is always LAST -
  `(n,)` for a single ODE state, `(n_times, n_cauldrons)` for a series.
  Accordingly `get_u_melt`, `get_material_heights` and `get_f_i`'s full-output
  branch in `ice_cauldron.py` were edited (see `TRANSLATION_NOTES.md`).

## Not yet translated

- Inter-cauldron flow physics (`getDrainageFluxes`'s inter-cauldron term is
  currently part of its dummy-zeros stub, not real physics).
- `plot_tools/cubehelix.m` (its 3 event colours are hard-coded in
  `plotting.py`) and `tightSubplot.m` (not needed - Panel handles layout).

## Next steps

- **Look at the dashboard.** It has been verified structurally (renders,
  axes link, no warnings) but not yet viewed - no browser was available to
  screenshot it. Run `python run_single.py` and review layout/legends/styling.
- **Test designs for the solver.** Design new benchmarks/tests for
  `glaciovolcano.py` / `gv_main.py` together (they conceptually cover the same
  physics/use cases as the MATLAB benchmarks - basic cauldron growth,
  multi-cauldron, supraglacial drainage, drainage density - but more robustly;
  the MATLAB RMSE thresholds are deliberately not being used).
- Work through the Todo items below individually.
- Dashboard follow-ups (see `TRANSLATION_NOTES.md`): saving/loading
  `GVResult` (netCDF keeps the attrs), run/cauldron selectors for the
  geometry panel, tabs/alternate views.

## Todo

Translation details flagged while planning the `glaciovolcano.m` / `gvMain.m`
translation. Each is to be revisited individually (file - issue). They are all
implemented as described, with inline `TRANSLATION NOTE`s, but still await
review:

- `gvpy/glaciovolcano.py` (RHS) - per-cauldron values must come from
  `gv.params[...]`, not the raw `gv.G_n` / `gv.L_n` fields, which are stored
  un-broadcast (scalar or list as passed). Also decide how to carry over the
  MATLAB oddities: the dead `phi_w = 0.4` assignment (overwritten by
  `drainage_density`), the `try/assert/except faafo` length check, and the
  `par.dV_cavity_dt = dV_i_melt_dt` naming. The unused `global step_ct`,
  `pathConfig` and `tic/toc` are to be dropped and noted.
- `gvpy/gv_main.py` (post-processing loop) - `IceCauldron` is a mutable attrs
  class but MATLAB's is a value class, so the per-step reconstruction loop
  needs `copy.copy(gv)` for `gv_temp`; `parse_events` mutates in place and
  returns the object.
- `gvpy/gv_main.py` (`openTime` / `iceFreeTime`) - MATLAB zero-fills these
  when a cauldron never had the event, which looks accidental rather than
  designed. Decide whether to translate literally or fail loudly.
- `gvpy/gv_main.py` (segment restart) - MATLAB's
  `InitialStep = t(nt) - t(nt-refine)` errors when fewer than `refine + 1`
  points precede the event; Python negative indexing would silently wrap
  instead, so it needs an explicit guard that raises.
- `gvpy/gv_main.py` (`parse_events`) - simultaneous events make `y0` a matrix
  in MATLAB (an error on the next solve). Keep as an error and flag, or
  handle properly.
- `gvpy/ice_cauldron.py` `get_material_heights`, called from the RHS - its
  results only feed `get_drainage_fluxes` (dummy zeros) yet it runs two
  `minimize_scalar` calls per cauldron per RHS call, which is removable
  overhead to measure. Its `assert H_w >= 0` may also trip on optimizer
  tolerance (SciPy vs `fminbnd`); if it does, stop and discuss rather than
  loosening it. (Its divide-by-zero warnings at `c_n = 0` are resolved.)
  - Check on the squeeze usage at the end of this function
- `gvpy/gv_main.py` and `run_single.py` (inputs) - MATLAB input names map to
  the Python field names (`nCauldrons` -> `n_cauldrons`, etc.) for the
  `gv_in` dict / keyword arguments.

## Tests

- `tests/test_data.py` - pre-existing placeholder stub (`assert False`),
  untouched.
- `tests/test_utilities.py` - `IceCauldron` property validators (attrs), and
  geometry (`ellipse_perimeter`, `get_elliptical_cylinder_sa`, `s_spheroid`:
  sphere/disc values, limit continuity, mixed vectors, roof curve).
- `tests/test_ice_cauldron.py` - vent/melting & ice inflow methods (closed-
  cauldron / default-state paths only; open-cauldron and full-output paths
  remain untested), `get_u_ice` per-cauldron shape, and
  `get_material_heights`' no-cavity handling (no warnings on solver-sized
  noise; warns on real volume inconsistency).
- `tests/test_glaciovolcano.py` - the flux terms added to `par` are per cauldron.
- `tests/test_var_info.py` - every field has complete `VAR_INFO`; `params`,
  `result.data` and `result.events` carry units/long_name after 1- and
  2-cauldron `gv_main` runs; `A`'s units follow `glen_n`.
- `tests/test_plotting.py` - full dashboard renders for 1 cauldron, 2
  cauldrons and 2 overlaid runs with no warnings and the expected axis
  linking; panel subset; event markers; missing-long_name error;
  unlisted-units warning; missing-variable placeholder; single closed
  geometry snapshot.
- Suite: 47 passing, plus the `test_data.py` stub.
- `tests/test_plume.py` - `get_plume_fluxes` (disabled-flux path only) and
  `predict_tree`/`predict_forest` against small synthetic trees.
- `tests/test_supraglacial.py` - `get_q_s` (off mode, overflow mode with and
  without the overflow condition triggered, shape-mismatch assertion). Also
  caught and corrected a mistaken assumption that its `faafo("fix", "it")`
  placeholder was reachable under normal inputs - it isn't (dead/defensive
  code, confirmed by this test).
- No test file for the subglacial drainage component yet (deliberately
  skipped - see `TRANSLATION_NOTES.md` for its known-broken state).
- No tests yet for `glaciovolcano.py` / `gv_main.py`. Per project direction the
  MATLAB benchmarks' RMSE thresholds are not being ported; new, more robust
  benchmarks are to be designed (see Next steps). Only informal checks so
  far: the RHS at `y0` reproduces the exact cavity mass balance (n=1, n=2, both
  supraglacial modes), and one `gv_main` run (`G_n=400`, no inflow, no
  supraglacial drainage) completed in ~1 s with the events at 606 s (open) and
  2108 s (ice-free). Since then, 2-cauldron runs (overflow mode, and
  `fixed-glen` ice inflow) run cleanly and are exercised by the metadata and
  plotting tests - but their physics outputs are not yet checked by any test.

## Notes

- The MATLAB source itself is mid-development (see `../katlaGV/README.md`):
  some files are marked `DEPRECATED` or `NOT CURRENTLY USED`, and
  `solveDeltaP.m` (subglacial pressure-wave solver) is explicitly known to be
  incomplete/in-progress. Translate what exists as-is per CLAUDE.md's
  in-development-components policy - don't fix it while porting.
- Two pre-existing bugs unrelated to translation work were found and fixed
  while testing this session: `gvpy/__init__.py` importing a deleted
  `gvpy.config` module, and `ice_cauldron.py` importing `ThermoConstants`
  from the wrong path (`.constants` instead of `.utilities.constants`).
- Project structure: `gvpy/plume.py` and `gvpy/subglacial.py` were placed as
  top-level `gvpy/` modules (not `gvpy/utilities/`), per explicit project
  direction that major model components' standalone functions warrant their
  own top-level module rather than being grouped under `utilities/` (the
  same will apply to `gvpy/supraglacial.py`). `getVarLabels.m`'s content now
  lives in `IceCauldron.VAR_INFO` (this session's decision), with only display
  conversions in `gvpy/plotting.py` - CLAUDE.md's Project Structure and Data
  Structures sections still describe the older plan and need updating.
- Deliberate change from MATLAB in `gv_main.py`: after a terminal event the
  next segment restarts at the event time, not at the next output-grid point
  (MATLAB's `timespan(timespan>=t(end))` seeds a later time with the event-time
  state, freezing the state for up to one output step and delaying later
  events). See `TRANSLATION_NOTES.md`.
