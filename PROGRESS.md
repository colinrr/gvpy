# Translation Progress

Snapshot of MATLAB -> Python translation status. This is a manually-maintained
snapshot, not derived automatically - update it as translation work happens,
and treat it as advisory (check current file contents/`git log` if in doubt).

See also `TRANSLATION_NOTES.md` for known bugs, behavioral changes, and
anticipated design decisions arising from this translation work.

Last updated: 2026-09-12

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

## Not yet translated

- The next thing to do is build a simple test script to load and Ice Cauldron 
  class and fiddle with methods, that will eventually become a tool for benchmark testing.
- **`getVarLabels.m`** (plotting variable labels/metadata) - per project
  direction this is plotting tooling, not model physics/state, so it belongs
  in a `gvpy/plotting.py` module (not yet created), not as an `IceCauldron`
  method. Whether to further fold its content into `xarray` variable `attrs`
  (per CLAUDE.md's flagged architecture question) is still open.
- The ODE solver / model driver (`glaciovolcano.m`, `gvMain.m`), and the
  events system it depends on (`get_f_i`'s full-output/events branch is
  translated but untested and unable to run until this exists).
- Inter-cauldron flow physics (`getDrainageFluxes`'s inter-cauldron term is
  currently part of its dummy-zeros stub, not real physics).
- Plotting/output utilities (`plot_tools/*.m`).

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
