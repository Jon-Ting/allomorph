# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-05-25

### Added

- Initial release of the package.
- Core functionality for generating monometallic, bimetallic, and trimetallic nanoparticle structural datasets for machine learning.
- MD simulation orchestration for LAMMPS.
- Feature extraction using NCPac and MD outputs.
- Comprehensive test suite and CI/CD workflow configuration.
- Basic documentation and usage examples in `README.md`.

### Changed

- Migrated documentation system from MkDocs to Sphinx.
- Updated project dependencies in `pyproject.toml` (`sphinx`, `sphinx-rtd-theme`, `myst-parser`).
- Configured Sphinx with `sphinx_rtd_theme` and `myst-parser` for Markdown support.
- Initialized documentation structure in `docs/` with `toctree`.
- Removed obsolete `mkdocs.yml`.