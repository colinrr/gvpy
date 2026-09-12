# Claude Code permissions (`settings.json`)

This project grants Claude Code access to the parent `katla` directory (which
holds sibling projects, including the original MATLAB source this repo is
being translated from), with different rules depending on the folder:

| Folder | Read (Read/Grep/Glob) | Write (Edit/Write) |
|---|---|---|
| `gvpy` (this project) | Auto-allowed, no prompt | **Ask** — prompts for confirmation every time |
| `katlaGV` (MATLAB source) | Auto-allowed, no prompt | **Denied** — blocked outright, no override |
| `hydroplume_emulator` (sibling project) | Auto-allowed, no prompt | **Denied** — blocked outright, no override |

Additionally:
- `git *` commands always **ask** for confirmation, reinforcing this
  project's CLAUDE.md rule that git operations require explicit user request.
- `rm -rf`, `git push --force`, `git reset --hard`, and `git clean -f` are
  **denied** outright.

## Why

The intent is to let Claude read across the whole `katla` tree for context
(e.g. comparing against the original MATLAB source during translation) while
only ever writing inside `gvpy`, and never touching the MATLAB source or
other sibling projects.

## Caveat

The `katlaGV`/`hydroplume_emulator` deny rules use a bare directory path
(`Edit(/path/to/katlaGV)`), matching the one documented example for
folder-scoped rules. This has not been independently verified to block
edits to files nested arbitrarily deep inside those folders (vs. only
direct children) — worth testing if this protection is safety-critical.

This file is documentation only; Claude Code does not read or act on it.
See `settings.json` in this folder for the actual enforced config.
