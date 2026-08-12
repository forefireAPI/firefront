# Testing ForeFire

This document describes how to run the automated tests for the ForeFire wildfire simulator. Tests are located within the `tests/` directory.

## Test Dependencies

Running the test verification scripts requires:

*   **Python 3**
*   **Python Libraries:** `lxml`, `xarray`, `netCDF4` (and their C dependencies like `libnetcdf-dev`).

Install Python libraries via pip:
```bash
pip3 install lxml xarray netCDF4
```

## Running the Unit Tests

The C++ unit tests exercise propagation and flux models one call at a time,
without running a simulation. They are built alongside everything else and run
through CTest:

```bash
cmake -S . -B build && cmake --build build -j
ctest --test-dir build --output-on-failure
```

They need no Python and no test data. Configure with
`-DFOREFIRE_BUILD_TESTS=OFF` to skip building them; wheel builds already do.

`tests/unit/README.md` describes what they cover and how to add one.

## Running the Core Test (`runff`)

The primary automated test, validated in our CI pipeline, is located in `tests/runff/`. This test verifies core simulation, save/reload functionality, and NetCDF/KML output generation against reference files.

**To run this test manually:**

1.  Ensure ForeFire is compiled (e.g., via `install-forefire.sh`).
2.  Navigate to the test directory: `cd tests/runff`
3.  Execute the test script: `bash ff-run.bash`

**Test Logic:**

The `ff-run.bash` script:
1.  Runs an initial simulation (`real_case.ff`) generating NetCDF output (`ForeFire.0.nc`) and a reload file (`to_reload.ff`).
2.  Runs a second simulation (`reload_case.ff`) using the reload file, which generates KML output (`real_case.kml`).
3.  Uses Python scripts (`compare_kml.py`, `compare_nc.py`) to compare the generated KML and NetCDF files against reference files (`*.ref`) with numerical tolerance, accounting for minor floating-point variations.
4.  Exits with status 0 on success, non-zero on failure.

## Other Tests

The `tests/` directory contains other subdirectories (`mnh_*`, `python`, `runANN`) for potentially testing specific features like coupled simulations or Python bindings. A main `tests/run.bash` script exists but is not currently fully validated in CI. Refer to specific subdirectories for details if needed.

## Compiler Warnings

ForeFire's own sources compile with `-Wall -Wextra` by default. The warnings
are not yet clean, so they are informational rather than fatal; two options
control this:

*   `-DFOREFIRE_ENABLE_WARNINGS=OFF` — build quietly.
*   `-DFOREFIRE_WARNINGS_AS_ERRORS=ON` — fail the build on any warning. Useful
    on a subset of files while clearing them; not yet usable repository-wide.

The flags apply to `libforefireL`, the `forefire` executable and the unit
tests. NetCDF's headers are included as system headers so their warnings do
not appear.

## Contributing

Please see `CONTRIBUTING.md` for guidelines on contributing to ForeFire, including adding new tests. Report any issues via the repository's issue tracker.