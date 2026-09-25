# Translation Progress

Snapshot of MATLAB -> Python translation status. This is a manually-maintained
snapshot, not derived automatically - update it as translation work happens,
and treat it as advisory (check current file contents/`git log` if in doubt).

See also `TRANSLATION_NOTES.md` for known bugs, behavioral changes, and
anticipated design decisions arising from this translation work.

Last updated: 2026-09-20

## Translated

- `gvpy/utilities/math_helpers.py` - general math/geometry helpers
- `gvpy/utilities/melting.py`
- `gvpy/utilities/constants.py`
- `gvpy/utilities/geometry.py`
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
- `run_single.py` - now calls `gv_main` with the same overrides as
  `run_single.m` (plus scratch `IceCauldron` objects `c`, `c2`, `c3` for
  inspection; `c3` deliberately/accidentally fails the `T_m` validator).
- Shape convention adopted with the solver: the cauldron axis is always LAST -
  `(n,)` for a single ODE state, `(n_times, n_cauldrons)` for a series.
  Accordingly `get_u_melt`, `get_material_heights` and `get_f_i`'s full-output
  branch in `ice_cauldron.py` were edited (see `TRANSLATION_NOTES.md`).

## Not yet translated

- **`getVarLabels.m`** (plotting variable labels/metadata) - per project
  direction this is plotting tooling, not model physics/state, so it belongs
  in a `gvpy/plotting.py` module (not yet created), not as an `IceCauldron`
  method. Whether to further fold its content into `xarray` variable `attrs`
  (per CLAUDE.md's flagged architecture question) is still open.
- Inter-cauldron flow physics (`getDrainageFluxes`'s inter-cauldron term is
  currently part of its dummy-zeros stub, not real physics).
- Plotting/output utilities (`plot_tools/*.m`).

## Next steps

- **Test designs and/or a results plotting dashboard, next session.** Design
  new benchmarks/tests for `glaciovolcano.py` / `gv_main.py` together (they
  conceptually cover the same physics/use cases as the MATLAB benchmarks -
  basic cauldron growth, multi-cauldron, supraglacial drainage, drainage
  density - but more robustly; the MATLAB RMSE thresholds are deliberately not
  being used). And/or build the run-results dashboard plot (the `TODO` in
  `run_single.py`), which also ties in `getVarLabels.m` / `gvpy/plotting.py`.
- Work through the Todo items below individually.

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
  loosening it.
  - Check on the squeeze usage at the end of this function
- `gvpy/gv_main.py` and `run_single.py` (inputs) - MATLAB input names map to
  the Python field names (`nCauldrons` -> `n_cauldrons`, etc.) for the
  `gv_in` dict / keyword arguments.

## Tests

- `tests/test_data.py` - pre-existing placeholder stub (`assert False`),
  untouched.
- `tests/test_utilities.py` - `IceCauldron` property validators (attrs).
- `tests/test_ice_cauldron.py` - vent/melting & ice inflow methods (closed-
  cauldron / default-state paths only; open-cauldron and full-output paths
  remain untested).
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
  2108 s (ice-free). Multi-cauldron and supraglacial-overflow runs through
  `gv_main` have NOT been run yet.

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
  same will apply to `gvpy/supraglacial.py`). Plotting/output metadata
  methods (e.g. `getVarLabels.m`) are a further exception, belonging in a
  future `gvpy/plotting.py` rather than as `IceCauldron` methods. CLAUDE.md's
  Project Structure section has been updated to reflect both.
- Deliberate change from MATLAB in `gv_main.py`: after a terminal event the
  next segment restarts at the event time, not at the next output-grid point
  (MATLAB's `timespan(timespan>=t(end))` seeds a later time with the event-time
  state, freezing the state for up to one output step and delaying later
  events). See `TRANSLATION_NOTES.md`.
