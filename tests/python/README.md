# Python examples and tests

Everything here needs the `pyforefire` module: either a build configured with
`-DFOREFIRE_BUILD_PYTHON=ON`, or `pip install forefire`.

`../run.bash` runs only the first two, through `$PYTHONEXE`.

## `percolation.py`

Four fires side by side, each in an 80×60 band of fuel filled at random to a
different density (`k_coeffs = [.10, .3, .4, .5]`), so the bands sit either
side of the percolation threshold. Writes `percolation.nc`. A demonstration of
building a heterogeneous fuel map in numpy and handing it to ForeFire with
`addIndexLayer`.

## `idealizedwind.py`

One fire under a wind rotating from 0° to 360°, which should trace a circular
front. Writes `360wind.png` through ForeFire's own `plot[]` command — the one
artefact `../run.bash` checks for.

## `farsite_flat.py`

Reproduces a FARSITE benchmark case: flat terrain, a 3 mph north wind, seven
hours of spread, compared against the FARSITE result.

**It cannot run as checked out.** It reads `flatland.lcp`, a landscape file
that is not in this repository. Download it first:

```bash
curl -LO https://github.com/mbedward/farsite/raw/refs/heads/master/examples/flatland/Inputs/a_lcpFiles/flatland.lcp
```

The weather it uses, `flatland_3mph0deg7hr.raws`, *is* in this directory.

## `test_wheel.py`

Not part of this suite. It is the smoke test `cibuildwheel` runs against a
built wheel — it checks that the module imports, that its vendored NetCDF
resolves, that the propagation models registered, and that a trivial
simulation advances. Run it against an *installed* `forefire`, never from the
source tree:

```bash
python tests/python/test_wheel.py
```
