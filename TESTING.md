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

## Running the Model Invariants (`test_moisture_invariants.py`)

The second test validated in CI, by the `invariants.yml` workflow. Where
`runff` compares ForeFire against frozen ForeFire output — and so cannot tell a
physics fix from a physics regression — this suite holds no reference data.
Every assertion follows from the published spread equations, so it stays valid
across recalibration.

It asserts, for each propagation model that consumes dead fuel moisture, that
the rate of spread stays finite for any moisture, decreases as moisture rises,
reaches zero at the moisture of extinction, responds to a dynamic
dead-moisture layer, and that `DataBroker` resolves every property the model
registers.

**To run it manually:**

1.  Install the Python package, which builds the `pyforefire` extension:
    ```bash
    python3 -m venv .venv
    ./.venv/bin/python -m pip install .
    ```
2.  Run the suite (add `-v` to print every probe's spread rate):
    ```bash
    ./.venv/bin/python tests/python/test_moisture_invariants.py
    ```

Restrict it while iterating with `--model NAME` and `--test NAME`, both
repeatable. It needs no fixtures — fuel, wind, temperature and moisture layers
are built in memory — and takes well under a minute.

Note that each probe runs in its own interpreter. The C++ core keeps mutable
global state, so a second `ForeFire()` in one process inherits the first one's
parameters and a parameter sweep silently returns one identical result. Keep
that in mind when writing any new Python test that varies parameters.

## Other Tests

The `tests/` directory contains other subdirectories (`mnh_*`, `runANN`) for testing specific features like coupled simulations. A main `tests/run.bash` script exists but is not currently fully validated in CI. Refer to specific subdirectories for details if needed.

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