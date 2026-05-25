# AGENTS — Guidance for AI coding agents

Purpose: Help AI coding agents be productive in this repository with minimal friction.

What to know immediately
- Repository root: see [README.md](README.md) for project overview and installation.
- Package: defined in [pyproject.toml](pyproject.toml) — `allomorph` CLI entrypoint is `allomorph`.
- Tests: pytest configuration is in `pyproject.toml` (`tests/`, pattern `test_*.py`).
- Docs: user docs are under [docs/](docs/index.md).

Quick developer commands
- Create dev environment and install dev deps: `uv venv && uv pip install -e ".[dev]"` (see README.md).
- Run tests: `pytest` (tests live in `tests/`).
- Lint/format: `ruff` is configured in `pyproject.toml`.
- Build package (local wheel/sdist): use standard Python build tooling; `python -m build` after activating environment.

Important places to look
- Core package: `src/allomorph` — generation logic and constants.
- Tests: `tests/` — unit tests and pytest entry points.
- Extras & research tools: `extras/`, `EAM/`, `InitStruct/`, `MDsim/` — these contain scripts and resources used in research workflows.
- Installation & external requirements: [DEPENDENCIES.md](DEPENDENCIES.md).
- Contribution guidelines: [CONTRIBUTING.md](CONTRIBUTING.md).

Agent behavior guidelines (concise)
- Link, don’t duplicate: prefer linking to existing docs instead of copying large sections.
- Minimal edits: make the smallest focused change required; avoid refactors unless requested.
- Tests first: when adding or changing behaviour, run `pytest` and ensure new tests accompany changes.
- Tooling: respect `pyproject.toml` and dev tools (`ruff`, `pytest`, `mkdocs`).

Agent customization files
- This repository also provides root-level AI customization assets:
  - `.prompt.md` for focused, reusable prompts.
  - `.skill.md` for multi-step workflows.
  - `.agent.md` for discoverable custom agents.

If you need more context
- Read the README, docs, and `pyproject.toml` before making assumptions about packaging or test configuration.
- Ask the repo owner before changing research-oriented directories (`extras/`, `EAM/`, `InitStruct/`, `MDsim/`).
