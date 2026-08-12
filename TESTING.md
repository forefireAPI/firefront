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

## Running the Concurrency Stress Test (`test_threading.py`)

Not run in CI, deliberately. On an ordinary interpreter the GIL serialises
every call into the extension, so the test skips and would report a green
result that proves nothing. It is a reproduction tool, run by hand when
working on the shared-state problem in
[#175](https://github.com/forefireAPI/forefire/issues/175).

It runs eight simulations in threads and requires each to reproduce the fire
node count it produces when run alone, plus a case that hammers object
construction to exercise the shared id counter.

**To run it**, you need a free-threaded build of CPython (`python3.14t` or
later) and `PYTHON_GIL=0`. The variable is required: `_pyforefire` does not
declare `py::mod_gil_not_used()`, so importing it switches the GIL back on and
the test skips.

```bash
python3.14t -m venv .venv
./.venv/bin/python -m pip install .
PYTHON_GIL=0 ./.venv/bin/python tests/python/test_threading.py
```

**It does not pass today.** On `dev` it segfaults or returns wrong node counts,
which is the point — it is the failing test the work in #175 has to make pass.
Wiring it into CI belongs with the last step of that issue, once it can pass
for the right reason.

## Other Tests

The `tests/` directory contains other subdirectories (`mnh_*`, `runANN`) for testing specific features like coupled simulations. A main `tests/run.bash` script exists but is not currently fully validated in CI. Refer to specific subdirectories for details if needed.

## Sanitizers

`-DFOREFIRE_SANITIZE=<list>` builds with `-fsanitize=<list>`, applied to the
compile line, the executables and the shared library. It also switches the
optimisation flags to `-g -O1 -fno-omit-frame-pointer`, since the default
release set (`-O3 -flto -fomit-frame-pointer`) makes sanitizer reports hard to
read.

```bash
cmake -S . -B build-asan -DFOREFIRE_SANITIZE=address
cmake --build build-asan -j
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build-asan --output-on-failure
```

`address` is what CI runs, on both the unit suite and `runff`, as a blocking
check. Other values are passed straight through — `undefined`, or
`address,undefined` for both — but only `address` is currently verified clean.

**`detect_leaks=0` is deliberate, not a workaround.** ForeFire reports zero
ASan *errors* — no use-after-free, no overflow, no double free — on either test
path, which is what makes a blocking job possible. It does leak: nothing owns a
`PropagationModel` (#159), so every one is reported. Leaving leak detection on
would produce a permanently failing job that everyone learns to ignore. The CI
workflow runs the leak check anyway as an informational step, so the number
stays visible, and it can be made blocking once #159 lands.

Note that the sanitizer build writes `bin/forefire` and `lib/libforefireL.so`
like any other build, so it replaces a release build in the source tree.

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