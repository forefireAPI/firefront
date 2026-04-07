# *Python-based* Test Suite

## Running in Docker

Build the image from the repository root:

```bash
docker build -t forefire .
```

Run any script:

```bash
# Interactive shell
docker run --rm -it forefire bash

# Run a script directly
docker run --rm forefire python3 tests/python/farsite_flat.py
```

## Scripts

> Looking for the **fire-spread predictor** and its Streamlit UI? They moved
> to [`apps/spread_predictor/`](../../apps/spread_predictor/).

### farsite_flat.py

Runs a simulation based on Farsite software for a north-wind of 3 mph.

Required input `.lcp` file can be downloaded from:
https://github.com/mbedward/farsite/raw/refs/heads/master/examples/flatland/Inputs/a_lcpFiles/flatland.lcp
