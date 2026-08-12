# Changelog

All notable changes to ForeFire are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Version numbers are the tags in this repository and the values reported by
`forefire -v`; they are also what `pip install forefire==<version>` resolves.

Entries below v2.5.0 were reconstructed from the release notes and the commit
history, so they summarise each release rather than list every change. The
`Full Changelog` link on each version gives the complete commit range.

## [Unreleased]

Merged since v2.5.0. Work still on the `dev` branch is listed in that branch's
copy of this file.

### Changed

- Linux wheels are built against a NetCDF without DAP and HDF4, roughly halving
  the wheel. ([#154])
- Free-threaded (`cp314t`) wheels are no longer published, because the core is
  not yet thread-safe; CPython 3.14 is declared supported. ([#155])

## [v2.5.0] — 2026-08-11

### Added

- **Pip-installable wheels** for Linux (x86_64, aarch64) and macOS (Apple
  Silicon and Intel), CPython 3.9 and newer. `pip install forefire` gives both
  the `forefire` command-line interpreter and the `pyforefire` module, with
  NetCDF bundled inside the wheel. Wheels are built without MPI and without
  `-march=native`. ([#151])

### Changed

- The Dockerfile is a multi-stage build, and only `libforefireL` is built in
  the builder stage. ([#147])
- The `install-forefire.sh` build is driven by an option-based `CMakeLists.txt`
  (`FOREFIRE_ENABLE_MPI`, `FOREFIRE_NATIVE_ARCH`, `FOREFIRE_BUILD_PYTHON`,
  `FOREFIRE_STATIC_CORE`, `FOREFIRE_BUILD_TOOLS`, `FOREFIRE_CHECK_LFS`).

### Fixed

- A double free in `FireFront` cleanup. ([#145])
- `CMakeLists.txt` was matched by `.gitignore`, so it was missing from the
  source distribution. ([#152])
- Bugs in `tests/run.bash`, including a `runANN` typo and an unreachable
  `./clean.bash`. ([#146])

**Full Changelog**: <https://github.com/forefireAPI/forefire/compare/v2.4.2...v2.5.0>

## [v2.4.2] — 2025-11-27

- Simplified the installation process, following review comments on the JOSS
  submission, and corrected the installation instructions.
- Aligned the coupling with
  [PACK-MNH-V5-7-2](https://src.koda.cnrs.fr/mesonh/mesonh-code/-/releases/PACK-MNH-V5-7-2).
- Aligned the repository with the accepted JOSS paper
  ([10.21105/joss.08680](https://doi.org/10.21105/joss.08680)).

**Full Changelog**: <https://github.com/forefireAPI/forefire/compare/v2.1.122...v2.4.2>

## [v2.1.122] — 2025-09-16

Tagged, but never published as a GitHub release.

- Docker images published to the GitHub container registry, and tested in CI.
- Python bindings reworked, with the example from [#103] applied.
- macOS CI runs the test script.

**Full Changelog**: <https://github.com/forefireAPI/forefire/compare/v2.0...v2.1.122>

## [v2.0] — 2025-06-05

The V2 release: the HTTP interface, and alignment with Meso-NH V4.7.2.

### Added

- The built-in **HTTP command server and web UI** (`listenHTTP[]`, or
  `forefire -l`), serving a map view of the simulation.
- The **Read the Docs documentation site**, built with Sphinx, Breathe and
  Doxygen.
- `install-forefire.sh`, with `-y` to add `forefire` to `PATH` and set
  `FOREFIREHOME`.
- A working Dockerfile, and GitHub Actions for Linux and macOS.
- The `RothermelAndrews2018` propagation model. ([#34])

### Changed

- CMake replaces SCons throughout; references to SCons were removed.

### Fixed

- The `-lnetcdf_c++4` link failure. ([#26])
- A division by zero in `setArrivalTime`/`getArrivalTime`. ([#36])

**Full Changelog**: <https://github.com/forefireAPI/forefire/compare/v1.2...v2.0>

## [v1.2] — 2024-01-16

### Added

- The `geojson` dump mode.

## [v1.1.10] — 2022-10-25

### Changed

- CMake became the default build system. ([#9])

## [v1.1.0] — 2022-09-28

Tagged before the sources were moved into `src/`, so that this point in the
repository stays easy to return to.

[Unreleased]: https://github.com/forefireAPI/forefire/compare/v2.5.0...dev
[v2.5.0]: https://github.com/forefireAPI/forefire/releases/tag/v2.5.0
[v2.4.2]: https://github.com/forefireAPI/forefire/releases/tag/v2.4.2
[v2.1.122]: https://github.com/forefireAPI/forefire/releases/tag/v2.1.122
[v2.0]: https://github.com/forefireAPI/forefire/releases/tag/v2.0
[v1.2]: https://github.com/forefireAPI/forefire/releases/tag/v1.2
[v1.1.10]: https://github.com/forefireAPI/forefire/releases/tag/v1.1.10
[v1.1.0]: https://github.com/forefireAPI/forefire/releases/tag/v1.1.0

[#9]: https://github.com/forefireAPI/forefire/issues/9
[#26]: https://github.com/forefireAPI/forefire/pull/26
[#34]: https://github.com/forefireAPI/forefire/pull/34
[#36]: https://github.com/forefireAPI/forefire/pull/36
[#103]: https://github.com/forefireAPI/forefire/issues/103
[#145]: https://github.com/forefireAPI/forefire/pull/145
[#146]: https://github.com/forefireAPI/forefire/pull/146
[#147]: https://github.com/forefireAPI/forefire/pull/147
[#151]: https://github.com/forefireAPI/forefire/pull/151
[#152]: https://github.com/forefireAPI/forefire/pull/152
[#154]: https://github.com/forefireAPI/forefire/pull/154
[#155]: https://github.com/forefireAPI/forefire/pull/155
