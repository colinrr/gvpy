# Translation Notes: Flags, Caveats, and Open Decisions

Consolidated view of what came up during MATLAB -> Python translation that's
worth a human's attention, beyond what's inline as `# TRANSLATION NOTE:`
comments in the code (this file summarizes/indexes those, it doesn't replace
them - see the source files for full context on each). Grouped by relevance:
what changed or was dropped, what's a known bug or in-development gap, and
what decisions are still ahead.

Last updated: 2026-09-12

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
- **`get_u_melt`'s fallback `i_open`/`i_ice_free` construction** (when
  `open_idx`/`ice_free_idx` aren't supplied) translates MATLAB's
  `repmat(obj.openCauldron', size(a,2))` literally via `np.tile`, but it's
  genuinely unclear whether the resulting shape was intended to generalize
  to multi-cauldron cases, or only happens to work for `n_cauldrons == 1`.
  Not guessed at further - flagged inline.
- **`get_f_i`'s full-output/`events` branch is untested** - it depends on an
  `events` table structure from the not-yet-translated ODE solver/events
  system, whose exact shape (DataFrame? dict of arrays?) isn't settled yet.
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
- **`get_q_s.m` translation** (supraglacial drainage) is still ahead. Per
  project direction, this is its own distinct model component (like plume
  and subglacial drainage), so it gets the same treatment: a top-level
  `gvpy/supraglacial.py` module for any standalone helpers, and its
  `IceCauldron` method grouped under its own flagged section in
  `ice_cauldron.py` - plus the same `checkVectorOrientation`-removal and
  `faafo`-bug-preservation treatment applied elsewhere.
- **`getVarLabels.m`**: per project direction, this is plotting tooling, not
  model physics/state - it belongs with plotting output (a `gvpy/plotting.py`
  module, not yet created), not as an `IceCauldron` method. Whether its
  content should further fold into `xarray` variable `attrs` (per CLAUDE.md's
  flagged architecture question) is still open.
- **`open_cauldron`/`ice_free_cauldron` migrating out of `IceCauldron`**
  into a separate mutable solver-state structure (already flagged inline in
  `ice_cauldron.py`) - deferred until the ODE solver/events system exists
  and the methods that read/write this state are all translated.
- **Plume emulator revamp**: both the emulator JSON itself (needs
  regeneration/correction) and `predictForest`/`predictTree`'s home -
  currently sourced from `physics_sandbox_2023/predictForest.m`, a
  dev-scratch file, not a settled utility - will need re-examination
  together.
- **Inter-cauldron flow and subglacial drainage physics**: `getDrainageFluxes`
  currently stubs both; once real physics is designed, it may make sense to
  split them into separate functions/components rather than one combined
  (currently dummy) call.
- **ODE solver + events system translation** (`glaciovolcano.m`, `gvMain.m`)
  is the biggest remaining unblock - several already-translated branches
  (`get_f_i`'s full-output mode, the eventual open-cauldron/ice-free state
  transitions) can't be exercised or tested until it exists.
