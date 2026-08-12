# Vendored third-party code

Dependencies that are checked in rather than fetched, so that a build and a CI
run need no network access beyond cloning the repository.

Nothing here is compiled into `libforefireL`: `CMakeLists.txt` globs `src/*.cpp`
only, and this directory is on the include path of the test target alone.

| Directory | Version | License |
| --- | --- | --- |
| `doctest/` | 2.5.3 (2026-07-06) | MIT |

## doctest

The unit-test framework, used by `tests/unit/`. A single header, taken
unmodified from
<https://github.com/doctest/doctest/blob/v2.5.3/doctest/doctest.h>.

The 2.5 series gates its C++17 features behind `DOCTEST_CPLUSPLUS` checks, so
it still builds under the project's `CMAKE_CXX_STANDARD 11`.

To update it, replace the file with a newer tagged release and run
`ctest --test-dir build`. Nothing in `tests/unit/` uses anything beyond
`TEST_SUITE`, `TEST_CASE`, `CHECK`, `REQUIRE`, `CAPTURE` and `doctest::Approx`.
