# Project Instructions for Claude

## Project Context
This is a scientific forward model, originally written in MATLAB, being translated to Python. Many physics components are not yet written - the MATLAB source itself is a work in progress. Existing "benchmarks" are coarse tests confirming fundamental behaviours in a couple of run modes; the complete model is not yet built or fully running. The immediate task is translating what exists, preserving equations and logic, then continuing development (including new physics) in Python. Coarse existing benchmarks will be turned into more rigorous tests over time - it's fine for those to be finalized in Python rather than matched exactly against MATLAB.

## Behavioural Rules

### 1. Think Before Coding
Don't assume. Don't hide confusion. Surface tradeoffs.

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First
Minimum code that solves the problem. Nothing speculative.

- No features beyond what was asked.
- No abstractions for single-use code.
- Do not create abstractions unless they reduce real, measurable complexity, or facilitate important elements of code modularity.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.
- Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.
- Three similar lines of code is better than a premature abstraction.
- This is scientific model code - interpretability of an equation or algorithm is more important than "clean abstraction" using typical software development principles.

### 3. Surgical Changes
Touch only what you must. Clean up only your own mess.

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently. Do make recommendations for what you'd change to improve structure or use alternative tools better suited to a task.
- If you notice unrelated dead code, mention it - don't delete it.


When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.
- The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution
Define success criteria. Loop until verified. Every changed line should trace directly to the task.

Transform tasks into verifiable goals:

"Add validation" → "Write tests for invalid inputs, then make them pass"
"Fix the bug" → "Write a test that reproduces it, then make it pass"
"Refactor X" → "Ensure tests pass before and after"
For multi-step tasks, state a brief plan:

1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

## Strict Rules for Modifying Code and translating code to python from matlab
1. DO NOT modify equations, or function logic inputs/outputs unless strictly necessary and then only after explicit approval. Variable and function names should be preserved in meaning/identity, but may be reformatted to PEP8 convention during translation (e.g. MATLAB `camelCase` -> Python `snake_case`) - see the PEP8 Formatting section below for what's in and out of scope for this.
2. DO NOT change or remove any existing code comments - replicate those that are already present EXACTLY. Code comments often capture specific logical, mathematical, or physics reasoning from when the code was written, and are treated as historical record, not just description - preserve verbatim, even if terse, outdated-looking, or oddly placed. Docstrings have more flexibility: syntax/format may be converted to Python-appropriate conventions (e.g. numpy/Google style), but the underlying descriptive content/text must be preserved, not rewritten or reworded. In both cases, keep any new additions minimal and lightweight - the goal is that what's new stays easy to distinguish from what's preserved, and added commentary doesn't bury or dilute the original.
3. DO feel free to flag necessary changes or things that are now out of place as a result of the above strict rules. User will make judgement calls on what to incorporate.

## Tests
Never modify existing tests unless the change intentionally updates the behaviour they cover, or we are explicitly working to fix the test itself.
If a test fails after your change, that means YOUR change is wrong, not the test.
Fix your code to make the existing test pass.

## Comments
Default to writing only brief (single line or less) comments per logical grouping of lines.
Reference existing code for style.
Add a comment when the WHY is non-obvious, or to briefly state the purpose of a block of code.

## Git
Never run git operations without a user's explicit request.

## Dependencies
Don't add new dependencies for something that can be done in a few lines of code.
Only recommend new dependencies, explaining their intended purpose and advantages.
Environment/dependency management uses conda and pyproject.toml.

## Communication
When done, explain what you changed and why in 2-3 sentences.
If you found something unexpected, mention it.
If you made a judgment call, explain your reasoning.

## Numerical Equivalence Validation
Strict adherence to MATLAB numerical benchmarks is neither necessary nor, in most cases, strictly possible - the model isn't yet advanced enough to require that level of exactness. Existing coarse benchmarks will be turned into more rigorous tests over time; those can be finalized in Python rather than matched precisely against MATLAB output.

## ODE Solver Translation Notes
No need to explicitly favour or replicate MATLAB's default solver settings. Do flag, however, where MATLAB and Python solver choices/defaults are likely to differ (e.g. `ode45`/`ode15s` vs. `solve_ivp` method/tolerance defaults), so it's a known and visible decision rather than a silent one.

## Units & Physical Constants
- Note the units of key physical quantities near their first use/definition if not already documented in the MATLAB source.
- Do not silently "fix" unit inconsistencies found during translation - flag them instead.
- If physical constants are hardcoded in MATLAB, preserve the same values/precision in Python; don't substitute a different-precision constant from a library without flagging it.

## Project Structure / Environment
- Keep translated Python files organized in a way that's traceable back to their MATLAB source (e.g. mirrored filenames/structure, or a comment noting source file).
- The new Python-based project structure broadly follows the Cookiecutter Data Science framework - not the existing MATLAB structure. For translation work specifically, use Cookiecutter Data Science as a guide, but adhere to the logical structure of the existing MATLAB project wherever possible; flag cases where the two conventions conflict so a judgment call can be made.
- Environment/dependency management uses conda and pyproject.toml as the single source of truth for dependencies.
- Flag when a MATLAB built-in (e.g. `interp1`, `fmincon`) doesn't have a drop-in Python equivalent and requires a judgment call.
- Standalone utility functions (i.e. not IceCauldron class methods) are consolidated into topic-based modules rather than mirrored 1:1 with MATLAB source files. Per the MATLAB README's "General Structure" section, the model is organized around inflow/outflow components of the IceCauldron control volume (vent/melting, plume, ice inflow, inter-cauldron flow, supraglacial drainage, subglacial drainage). When a utility function is used exclusively by one such component, it belongs in a component-named module (e.g. `utilities/subglacial.py`, `utilities/plume.py`). Generic, reusable math/geometry helpers not tied to a specific component (e.g. a smooth-min blending function) go in a general module (e.g. `utilities/math_helpers.py`, `utilities/geometry.py`) instead.
- The `@IceCauldron` classdef folder translates to a single `ice_cauldron.py` module containing the `IceCauldron` class with all its methods (this does not change based on the component-based utility organization above - it's specifically about standalone utility functions outside the class).

## In-Development Physics Components
Per the MATLAB README's "General Structure" section, each inflow/outflow component of the IceCauldron control volume (vent, plume, ice inflow, inter-cauldron flow, supraglacial drainage, subglacial drainage) is intended to have "off"/"test"/final modes, controlled by IceCauldron switches (e.g. `ice_inflow_mode`, `supraglacial_drainage_mode`, `disable_plume_flux`). Some components' physics (e.g. subglacial drainage/pressure-wave propagation in `solveDeltaP` and its dependents) are still actively in development in the MATLAB source, and may contain incomplete, inconsistent, or broken code (undefined variables, placeholder tokens, dead branches) reflecting that in-progress state.

- When translating an in-development/broken component: preserve the logic AS IS, bugs included, rather than fixing it. Do not "helpfully" correct an undefined variable, complete a stub, or resolve an inconsistency - that risks losing the actual state of in-progress thinking.
- Flag anything broken, incomplete, or confusing with a `# TRANSLATION NOTE:` (per the matlab-to-python-translation skill's flagging protocol) rather than silently working around it or silently reproducing it without comment.
- If genuinely unclear what a piece of broken code was attempting to do, say so explicitly in the flag rather than guessing.
- This preservation-over-correction approach applies specifically to physics/logic. Mechanical syntax translation (MATLAB -> Python conversion itself) still must produce valid, running Python - "preserve the bug" means preserve the logical/mathematical error, not intentionally write broken Python syntax.

## Type Hints
- DO add type hints to Python translations. (See Strict Rule #2 for docstring/comment handling.)

## PEP8 Formatting During Translation
When translating MATLAB to Python, apply PEP8 formatting mechanically where it's purely cosmetic; treat anything that changes a name's meaning or signals new structure as a judgment call to flag, not to auto-apply.

**Apply directly (no need to flag):**
- `camelCase`/`PascalCase` variable and function names -> `snake_case`. The underlying name/meaning must be preserved (e.g. `windSpeed` -> `wind_speed`, not renamed to something else) - this does not count as changing a variable name under Strict Rule #1.
- Indentation -> exactly 4 spaces.
- General whitespace/operator spacing per PEP8.

**Use judgment, flag the choice:**
- Line length: don't rigidly enforce a PEP8 line-length limit at the expense of breaking up an equation in a way that obscures it or its correspondence to the MATLAB source. Prefer readability/traceability over strict line length.
- `UPPER_SNAKE_CASE` for true constants: only apply when a value is clearly a fixed constant (not reassigned, not a function input) - flag which variables you treated as constants.
- Leading-underscore (`_helper`) for internal-only helpers: this signals something structural (internal API), not just cosmetic - flag rather than apply automatically.
- Class naming (`PascalCase`): only relevant if/when a MATLAB struct becomes a Python class - ties to the Abstractions guidance above.
- Boolean semantic prefixes (`is_`/`has_`): don't add automatically, but flag boolean variables that lack one, since we may want to adopt this convention going forward - a call to make deliberately, not silently in either direction.

**Don't do, even though it's common Python style:**
- Don't expand or contract abbreviations in names during translation.

## Reproducibility
- Preserve random seeds / RNG usage exactly as in MATLAB where present.
- RNG use should be rare, if present at all - flag any place a random seed or RNG is used during translation, even if handled correctly, so it gets explicit attention.

## Data Structures: dataclasses and xarray
The Python translation uses `dataclasses` for model objects (e.g. `IceCauldron`) rather than replicating MATLAB's `classdef`, and uses `xarray` for data with meaningful dimensions (e.g. per-cauldron, or time x cauldron), rather than plain NumPy arrays with implicit dimension order. General conventions established during translation of `IceCauldron`:
- Per-cauldron model parameters (both user-supplied init values and computed secondary values derived from them, e.g. `Q_n`, `tau_melt`) live together in one `xarray.Dataset` (e.g. `self.params`), built in `__post_init__`, with a `cauldron` dimension/coordinate. This resolves MATLAB's row/column vector-orientation ambiguity structurally (dimensions are named, not positional) rather than needing defensive orientation-checking code.
- Mutable per-cauldron simulation state (e.g. `openCauldron`, `iceFreeCauldron`, which change during the ODE solve via events) lives in a separate structure from the init/computed params - don't mix state that changes during simulation with parameters that don't.
- ODE solver state/output (time series of solution and derived variables across time and cauldron) uses `xarray` with `(time, cauldron)` dims once assembled.
- Performance note: the ODE right-hand-side function itself is called every solver step and is performance-critical - prefer plain NumPy arrays/operations inside that inner loop rather than xarray, to avoid per-step overhead. Assemble results into xarray structures after the solve completes (or at designated output/checkpoint points), not inside the hot loop. Flag any place this tradeoff is ambiguous.
- Metadata (labels, units, plotting scales - see `getVarLabels.m`) is a good candidate to eventually fold into `xarray` variable `attrs` rather than a separate parallel lookup structure, but this is a structural change to how downstream (e.g. plotting) code consumes that metadata - flag rather than doing it automatically during core model translation.
- When translating a function/method whose inputs are per-cauldron vectors, check whether it should now accept/return `xarray.DataArray`/`Dataset` objects (dimension-aware) vs. plain NumPy arrays; flag the choice made and reasoning, especially for functions called inside performance-critical loops (see above).
