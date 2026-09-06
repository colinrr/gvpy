---
name: matlab-to-python-translation
description: Use this skill whenever translating MATLAB (.m) code to Python as part of this forward-model project - including translating a single function, a file, or reasoning about how a MATLAB construct should map to Python/NumPy/SciPy. Also use it when reviewing or debugging previously-translated code for correctness against the MATLAB source, or when asked to check whether a translation preserved equations, logic, comments, or docstrings correctly. Trigger on requests like "translate this MATLAB function", "convert this .m file to Python", "port this to Python", or "does this Python match the MATLAB logic". This skill defers to CLAUDE.md for project-wide constraints (comment/docstring preservation, abstraction philosophy, scope discipline) and focuses specifically on the mechanics of MATLAB -> Python translation.
---

# MATLAB to Python Translation

This skill governs how MATLAB source code for this forward model gets translated to Python. It assumes familiarity with the project's `CLAUDE.md`, which takes precedence on anything not specific to translation mechanics (comment/docstring preservation, abstraction restraint, scope discipline, communication style).

## Before translating anything

1. Read the full MATLAB file/function first, don't translate line-by-line blind. Understand what the function does physically/mathematically before converting syntax.
2. Identify which parts are: core equations/physics (translate exactly), boilerplate/utility (e.g. plotting, I/O - flag if out of scope for this pass), and MATLAB-specific scaffolding (e.g. `nargin` checks) that may need a Python-appropriate equivalent.
3. Check whether the function is already partially translated elsewhere in the Python project (per CLAUDE.md's traceability-to-MATLAB-source structure) - don't duplicate.
4. If a MATLAB benchmark/test result exists for this function, note it, but remember (per CLAUDE.md) exact numerical matching is not required at this stage - approximate sanity-checking is enough unless told otherwise.

## Core translation principle

**Translate logic and equations as-is. Translate syntax and naming style as needed.**

MATLAB syntax constructs (indexing, loops, matrix operations) must be converted to correct, idiomatic Python/NumPy - that's necessary, not optional, and doesn't violate the "don't modify logic" rule. Variable/function names may be reformatted to PEP8 (`camelCase` -> `snake_case`) as long as the underlying name/meaning is preserved - see CLAUDE.md's "PEP8 Formatting During Translation" section for exactly what's mechanical vs. what needs flagging. What must NOT change without approval: the mathematical operations themselves, the identity/meaning of variable and function names, function signatures (inputs/outputs), and comments/docstring content (see CLAUDE.md).

When a MATLAB idiom has no clean Python equivalent, or when multiple reasonable translations exist, flag it rather than silently picking one (see "Flagging" below).

## Comments and docstrings during translation

Per CLAUDE.md: existing code comments are preserved EXACTLY, verbatim, in their new location in the Python file - even if they reference MATLAB-specific syntax now (e.g. a comment mentioning `end` or `~`). Do not "fix" or "modernize" a comment's wording during translation. If a comment is rendered confusing or inaccurate purely because of the syntax change (not the logic), flag it - don't rewrite it.

Docstrings may be reformatted into Python-appropriate syntax (e.g. numpy/Google style docstring format) but their descriptive content must be preserved, not reworded.

## Common MATLAB -> Python translation points

These are the recurring places translation bugs happen. Check each deliberately rather than relying on pattern-matching.

**Indexing**
- MATLAB is 1-indexed; Python/NumPy is 0-indexed. Every index and loop bound needs conversion, not just the obvious ones.
- MATLAB ranges `a:b` are inclusive of `b`; Python `range(a, b)` and NumPy slices are exclusive of `b`. Off-by-one is the single most common translation bug - double check every converted range.
- `end` (as in `x(end)`, `x(2:end)`) becomes `-1` / `x[1:]` etc. - watch for `end-1`, `end-2` type expressions.
- MATLAB's `:` for "all elements along a dimension" maps to `:` in NumPy too, but confirm dimension order (see below).

**Array/matrix semantics**
- MATLAB arrays are column-major (Fortran order); NumPy defaults to row-major (C order). This usually doesn't matter for correctness but can matter for performance-sensitive reshaping/memory layout - flag if reshape/memory order looks performance-relevant.
- MATLAB matrix operations (`*`, `/`, `^`) are matrix operations by default; elementwise requires `.*`, `./`, `.^`. NumPy `*` is elementwise by default; matrix multiplication requires `@` or `np.matmul`. This is a frequent silent-bug source - check every arithmetic operator on arrays individually, don't assume.
- MATLAB's backslash `\` (mldivide, linear system solve) maps to `numpy.linalg.solve` (square systems) or `numpy.linalg.lstsq` (non-square/least-squares) - confirm which applies and flag the choice.
- `size()`, `numel()`, `length()` map to `.shape`, `.size`, and `len()`/`max(.shape)` respectively - these are not interchangeable, confirm which MATLAB meant in context.

**Data structures**
- MATLAB structs -> Python: prefer a plain approach consistent with existing project patterns (e.g. `dataclass`, `dict`, or `SimpleNamespace`) rather than introducing a new struct-like abstraction; check CLAUDE.md's abstraction guidance if unsure which fits.
- MATLAB cell arrays -> Python `list` (heterogeneous) or `dict` (if keyed).
- MATLAB struct arrays -> typically a `list` of structs/dataclasses/dicts in Python; flag if a NumPy structured array or pandas DataFrame would be a better fit given how it's used downstream, but don't switch without asking.

**Control flow & functions**
- `switch`/`case` -> `if`/`elif` (or `match`/`case` if the codebase already uses Python 3.10+ patterns - check first).
- `nargin`/`nargout`/`varargin`/`varargout` -> Python default arguments, `*args`, `**kwargs`, or explicit `None` checks. Flag the chosen approach since it changes the function's call signature semantics.
- MATLAB function handles (`@myfun`, anonymous `@(x) ...`) -> Python functions or `lambda`.
- `persistent` variables -> flag explicitly; typically need a class attribute, closure, or module-level state in Python, which is a real structural change worth approval.
- `global` variables -> flag explicitly, same reasoning.

**Errors, NaN, and numeric edge cases**
- `try`/`catch` -> `try`/`except`; note that MATLAB's `catch` without a variable catches everything, Python's bare `except:` is generally discouraged - flag rather than silently narrowing/widening the exception scope.
- NaN handling: MATLAB and NumPy both support NaN, but comparison/propagation behavior in functions like `max`/`min`/`sum` can differ (e.g. `nanmax` vs `max` with NaN-aware flags) - check explicitly when NaNs are plausible in the data.
- Integer vs float division: MATLAB is float by default; confirm Python isn't doing integer division unexpectedly (relevant in Python 2 style code, less so in Python 3, but still worth a glance with `//` vs `/`).

**ODE solvers** (see also CLAUDE.md's ODE Solver Translation Notes)
- Map solver choice deliberately: `ode45` ~ `solve_ivp(method='RK45')`, `ode23` ~ `RK23`, `ode15s`/`ode23s` (stiff) ~ `BDF`/`Radau`/`LSODA`. Don't default to `RK45` without checking what the MATLAB source used.
- MATLAB's `odeset` options (`RelTol`, `AbsTol`, `MaxStep`, event functions) need explicit equivalents in `solve_ivp`'s `rtol`, `atol`, `max_step`, `events` - don't rely on SciPy defaults matching MATLAB defaults (they don't).
- MATLAB event functions (`[value, isterminal, direction]`) map to SciPy event functions with `.terminal`/`.direction` attributes set on the function object - syntax differs meaningfully, translate carefully.
- Flag solver/tolerance choices explicitly per CLAUDE.md, since exact benchmark matching isn't required but the choice should still be visible and deliberate.

**I/O, printing, plotting**
- `fprintf`/`disp` -> Python `print`/f-strings. Low-risk, but preserve any formatting precision (e.g. `%0.4f`) that reflects a deliberate choice about displayed precision.
- Plotting (if present) is likely out of scope for the core forward-model translation pass - flag rather than translate unless asked.

## Flagging protocol

When something needs a judgment call, a heads-up, or can't be translated as a mechanical 1:1 mapping, flag it clearly rather than silently deciding. Use a consistent, greppable marker distinct from preserved original comments, e.g.:

```python
# TRANSLATION NOTE: MATLAB used ode45 (default tol); mapped to solve_ivp(method='RK45'), default rtol/atol used - confirm if MATLAB odeset specified tighter tolerances.
```

Keep these tight and specific - one line where possible. Summarize all flags for a given file in your end-of-task communication (per CLAUDE.md) rather than only leaving them scattered in code, so they're easy to review in one place.

## After translating

1. Read through the translated function once as Python, independent of the MATLAB source, to sanity-check it's coherent Python (not just syntax-converted MATLAB).
2. Re-diff mentally against the MATLAB source for the specific things that must not have changed: equations, variable names, comment text, docstring content, function signature.
3. Summarize per CLAUDE.md's Communication rule: what was translated, any flags raised, any judgment calls made and why.
4. Don't add tests, docstring additions, or refactors beyond what's needed for a faithful translation, unless asked - that's separate follow-on work per CLAUDE.md's Scope rule.
