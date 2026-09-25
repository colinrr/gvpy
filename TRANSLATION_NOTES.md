# Translation Notes: Flags, Caveats, and Open Decisions

Consolidated view of what came up during MATLAB -> Python translation that's
worth a human's attention, beyond what's inline as `# TRANSLATION NOTE:`
comments in the code (this file summarizes/indexes those, it doesn't replace
them - see the source files for full context on each). Grouped by relevance:
what changed or was dropped, what's a known bug or in-development gap, and
what decisions are still ahead.

Last updated: 2026-09-20

## Changed Functionality / Removed Components

Things that behave differently from the MATLAB source, or MATLAB constructs
that were deliberately dropped/simplified during translation (not bugs -
intentional, flagged decisions).

- **`checkVectorOrientation` removed entirely.** MATLAB's defensive
  row/column-vector orientation checking is dropped throughout (e.g.
  `get_u_melt`, `get_material_heights`, `get_l_lambda`, `validate_vector_length`).
  This project resolves vector-shape ambiguity structurally (consistent array
  shapes, named `xarray` dimensions) rather than with runtime orientation
  checks/transposes. Callers are expected to pass already-aligned shapes.
- **1-based -> 0-based indexing conventions**, beyond the mechanical
  index-arithmetic conversions: `fastest_cauldron_index`/`slowest_cauldron_index`
  are 0-based (`np.argmin`/`np.argmax`) throughout the Python codebase, a
  deliberate consistent convention shift. `plume.py`'s `predict_tree` drops
  MATLAB's `tree.feature + 1` adjustment since Python/JSON are already
  0-based.
- **Property validation**: MATLAB's declarative `{mustBeX}` property
  constraints (~28 fields) are now `attrs.field(validator=...)`, using
  `attrs` (added as a project dependency) rather than hand-written
  `__post_init__` checks - see `gvpy/utilities/validators.py`.
- **MATLAB optimizer/root-finder -> SciPy substitutions** (flagged
  individually inline, listed here for visibility - defaults/tolerances are
  not guaranteed equivalent):
  - `fminbnd` -> `scipy.optimize.minimize_scalar(method="bounded")`
    (`get_material_heights`)
  - `fminsearch` -> `scipy.optimize.minimize(method="Nelder-Mead")`
    (`solve_delta_p`)
  - `fzero` -> `scipy.optimize.fsolve` (`get_x_equals_l_lambda`)
  - `integral` -> `scipy.integrate.quad` (`get_elliptic_heat_intensity`)
  - MATLAB's `FunValCheck` optimizer option (error on NaN/Inf objective
    values) has no direct SciPy equivalent and was dropped (`solve_delta_p`).
- **Dropped a MATLAB debug/QC plotting block** inside `solve_delta_p`
  (figure/plot/`get(gca,...)` calls, plus a source typo referencing `hi`
  instead of `h_i`) - left as a flagged comment rather than ported to
  matplotlib, since this project's plotting utilities (`plot_tools/*.m`)
  haven't been translated at all yet.
- **Project structure**: `gvpy/plume.py` and `gvpy/subglacial.py` (and,
  once translated, `gvpy/supraglacial.py`) are top-level `gvpy/` modules, not
  `gvpy/utilities/<component>.py` - per explicit project direction that
  major model components warrant their own top-level module. Plotting/output
  metadata methods (e.g. `getVarLabels.m`) are a further exception to the
  "`@IceCauldron` -> single `ice_cauldron.py`" rule: they belong with
  plotting tooling (`gvpy/plotting.py`, not yet created), not as `IceCauldron`
  methods. CLAUDE.md's Project Structure section has been updated to reflect
  both.
- **Two unrelated pre-existing bugs fixed** (not translation decisions, just
  noting they were touched): `gvpy/__init__.py` imported a deleted
  `gvpy.config` module (removed); `ice_cauldron.py` imported `ThermoConstants`
  from the wrong path, `.constants` instead of `.utilities.constants` (fixed).

### Solver (`glaciovolcano.py`, `gv_main.py`, `utilities/solver_helpers.py`)

- **Shape convention: the cauldron axis is always LAST** - `(n,)` for a single
  ODE state, `(n_times, n_cauldrons)` for a series (matches xarray
  `(time, cauldron)`). Adopting it required user-approved edits to
  already-translated `IceCauldron` methods, because the RHS passes 1-D state
  slices and two methods could not handle them:
  - `get_u_melt`: the fallback open/ice-free masks are now
    `np.broadcast_to(self.open_cauldron, a.shape)` (my earlier `np.tile`
    substituted `n_cauldrons` for MATLAB's `size(a,2)` and failed for every
    n - MATLAB's `repmat(obj.openCauldron', size(a,2))` is only right for a
    single time column, the RHS case). The per-cauldron start-index branch is
    laid out `(time, cauldron)` and 0-based.
  - `get_material_heights`: accepts `(n,)` (promoted to one time step,
    outputs squeezed back). The squeeze logic is on the Todo list for review.
  - `get_f_i` full-output branch: reads the events columns as plain numpy so it
    works for any events container; `cauldronIndex` is 0-based.
- **`ode15s` -> `scipy.integrate.solve_ivp(method="BDF")`**, `rtol=1e-7` as in
  MATLAB, `atol` left at SciPy's default (1e-6, the same as MATLAB's default
  `AbsTol`). Not guaranteed equivalent step control. MATLAB's `odeset`
  `NonNegative` (all states) has no SciPy equivalent and is NOT replicated.
- **Restart at the event time (deliberate change from MATLAB).** MATLAB
  restarts each post-event segment on `timespan(timespan>=t(end))`, the next
  output-grid point, but seeds it with the state at the event time. That
  freezes the state for up to one output step (~5 s here) and delays every
  later event. On a throwaway replica of a single-cauldron run this gave a
  relative RMSE of ~3e-3 on cavity volume with MATLAB-style restarts versus
  ~3e-11 restarting at the event time. It is a pure time relabeling while the
  ODE is autonomous, but would become a real error once explicit-time terms
  (e.g. `t_Qdecay`) exist. Python restarts at the event time and keeps the
  grid points after it.
- **The MATLAB benchmark RMSE thresholds are not being used** as test targets
  (too weak/rough, per project direction); new benchmarks are to be designed.
- **Terminal event point handling**: like `ode15s`, the terminal event point is
  appended as the last row of each segment (SciPy with `t_eval` only returns the
  grid points).
- **Events**: MATLAB's vector-valued `gvEvents` becomes one scalar callable per
  event (`solve_ivp` requirement), with `terminal`/`direction` read from a probe
  call at `y0`; the events index table is built once and passed in (MATLAB
  rebuilds it every call). The events table (`get_events_table`) is a pandas
  DataFrame with 0-based `cauldronIndex` (NaN for non-per-cauldron events).
- **Result container**: `gv_main` returns a `GVResult` (`.gv`, `.data` =
  xarray Dataset dims time x cauldron, `.events` = xarray Dataset dim `event`)
  instead of MATLAB's single `dat` struct. Variable names keep the MATLAB keys
  (e.g. `horizontalMeltingRate_n`) so `derived_vars()` still lines up.
- **Inputs**: a dict and/or keyword arguments using the Python `IceCauldron`
  field names (`nCauldrons` -> `n_cauldrons`, ...), instead of a struct /
  name-value pairs.
- **Dropped or replaced**: `pathConfig` (MATLAB path setup), `tic`/`toc` (kept
  only in `run_single.py`), the unused `global step_ct`, re-creating the ODE
  `options` after each event (event functions read `gv`, mutated in place, at
  call time), and `assign_integrated_values`'s `solutions` flag (redundant with
  the axis-last convention).
- **Value vs. reference semantics**: MATLAB's `IceCauldron` is a value class;
  the attrs one is mutable. `parse_events` mutates `gv` in place, and the
  post-processing reconstruction loop works on a `copy.copy(gv)`.
- **`InitialStep` guard**: MATLAB's `t(nt)-t(nt-refine)` errors when fewer than
  `refine + 1` points precede the event; Python negative indexing would silently
  wrap, so it raises `IndexError` explicitly. The step is also capped at the
  remaining span (SciPy raises if `first_step` exceeds it; `ode15s` just clips).

## Known Bugs and In-Development Issues

Preserved from the MATLAB source as-is, per CLAUDE.md's in-development-
components policy - not fixed during translation. Flagged here for
visibility since some of these make entire code paths currently unusable.

- **`solve_delta_p` cannot currently complete its one implemented regime.**
  `objective_delta_p_1` (`gvpy/subglacial.py`) references `h_i` and `l`,
  neither of which is a parameter or locally defined - a genuine bug in the
  MATLAB source (local functions there don't share the parent function's
  workspace). Calling `solve_delta_p` with `x <= l_lambda` (the only
  implemented regime) raises `NameError` as soon as the optimizer evaluates
  the objective. Confirmed by hand; not covered by a test per project
  direction (see PROGRESS.md).
- **`solve_delta_p`'s `x > l_lambda` regime is unimplemented** - raises
  `NotImplementedError`, matching MATLAB's `error('Regime of x > l_lambda
  not yet implemented.')`.
- **`drainage_density`'s `flood_density_model == 6`** is only a comment stub
  in MATLAB (`% case 6 % Based on pressure/discharge/proportions in
  cauldron?`), never an implemented branch - selecting it (a value the attrs
  validator otherwise allows, range 1-6) raises `UnboundLocalError`.
- **`get_drainage_fluxes` is a dummy stub** returning zero fluxes - not
  implemented in the MATLAB source either. This is the actual
  subglacial-drainage AND inter-cauldron-flow flux calculation, so both
  components are effectively unimplemented pending this.
- **`get_u_melt`'s and `get_q_s`'s `faafo`/`faafo fix it` placeholder bugs**:
  MATLAB's `catch` blocks contain literal undefined-function placeholder
  text (dev scratch, not real error handling). Translated as `faafo()`/
  `faafo("fix", "it")` calls, which would raise `NameError` if ever
  triggered - preserving the "this isn't handled" intent rather than
  writing invalid Python syntax or silently fixing it. For `get_q_s`
  specifically, this is confirmed (via test) to be unreachable dead code
  under normal, already-shape-validated inputs - the `try` body is a plain
  numpy assignment that doesn't fail in ordinary use.
- **`get_f_i`'s full-output/`events` branch** is now exercised by `gv_main`'s
  post-processing (one informal single-cauldron run), but has no pytest
  coverage yet. Like MATLAB, a run needs every cauldron to reach ice-free before
  `t_final`: otherwise `ev_t[event_idx]` is empty and numpy raises a broadcast
  error (MATLAB: `t >= []`), after the solve has finished.
- **Simultaneous solver events** raise `NotImplementedError` in `parse_events`
  (MATLAB turns `y0` into a matrix and errors on the next solve).
- **`openTime`/`iceFreeTime` zero-fill**: MATLAB silently zero-fills these for a
  cauldron without the event; translated literally (zeros). Unreachable in
  practice because `get_f_i` already fails unless every cauldron went ice-free,
  and a cauldron opens before it can go ice-free.
- **`get_material_heights` in the RHS**: its results only feed
  `get_drainage_fluxes` (dummy zeros), yet it runs two `minimize_scalar` calls
  per cauldron per RHS call (removable overhead, MATLAB does the same). At the
  initial state (`c_n = 0`, `V_cavity = 0`) it also divides by zero and emits
  `RuntimeWarning`s, and its `assert H_w >= 0` could in principle trip on
  optimizer tolerance (SciPy vs `fminbnd`) - not seen so far.
- **Only one `gv_main` configuration has been run so far** (single cauldron,
  `G_n=400`, no ice inflow, no supraglacial drainage). Multi-cauldron and
  supraglacial-overflow runs are untested.
- **`t_Qdecay`** is declared as an `IceCauldron` property in MATLAB but never
  assigned anywhere in the source - left as `None` in Python, not guessed at.
- **Plume emulator JSON is intentionally absent.** The default
  `plume_emulator_json = "hydroplume_emulator_2025-04-22.json"` was removed
  from the repo (disk space, per project direction) - the `get_plume_fluxes`
  random-forest path can't be exercised end-to-end until it's restored (with
  corrections) as a future task. The `disable_plume_flux` zeros-path is the
  reliable path for now.
- **`get_u_tip` recomputes the Nikuradse roughness height locally** rather
  than reusing `self.k_nikuradse` (computed identically in
  `__attrs_post_init__`) - preserved as redundant, not deduplicated, since
  that's how the MATLAB source does it too.

## Anticipated Future Decisions

Design forks flagged during this pass that aren't resolved yet - things to
expect needing a call on in upcoming translation/development work.

- **Fixing `solve_delta_p`'s `objective_delta_p_1` bug** (undefined
  `h_i`/`l`) will require figuring out what those were actually supposed to
  be - likely re-deriving them from `x` via `get_elevation_profiles`/
  `get_l_lambda`, but that's a physics judgment call, not a translation one.
- **Test designs for the solver** (next session): new, more robust
  benchmarks/tests for `glaciovolcano.py` / `gv_main.py` that conceptually cover
  the same physics and use cases as the MATLAB benchmarks (basic cauldron
  growth, multi-cauldron, supraglacial drainage, drainage density) without their
  RMSE thresholds - to be worked out together.
- **Results dashboard plot** (the `TODO` in `run_single.py`), together with
  `getVarLabels.m` / a future `gvpy/plotting.py`.
- **`NonNegative` handling**: not replicated; whether any state can go
  negative in practice (and if so, clipping in the RHS vs. leaving it to events)
  is open.
- **`getVarLabels.m`**: per project direction, this is plotting tooling, not
  model physics/state - it belongs with plotting output (a `gvpy/plotting.py`
  module, not yet created), not as an `IceCauldron` method. Whether its
  content should further fold into `xarray` variable `attrs` (per CLAUDE.md's
  flagged architecture question) is still open.
- **`open_cauldron`/`ice_free_cauldron` migrating out of `IceCauldron`**
  into a separate mutable solver-state structure (already flagged inline in
  `ice_cauldron.py`) - now unblocked, since `gv_main` exists: it currently
  mutates them in place (`parse_events`) and reassigns them on a copy for the
  post-processing loop.
- **Plume emulator revamp**: both the emulator JSON itself (needs
  regeneration/correction) and `predictForest`/`predictTree`'s home -
  currently sourced from `physics_sandbox_2023/predictForest.m`, a
  dev-scratch file, not a settled utility - will need re-examination
  together.
- **Inter-cauldron flow and subglacial drainage physics**: `getDrainageFluxes`
  currently stubs both; once real physics is designed, it may make sense to
  split them into separate functions/components rather than one combined
  (currently dummy) call.
