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
docker run --rm forefire python3 tests/python/predict_spread.py

# Copy outputs (HTML map + GeoJSON) to your host
docker run --rm -v $(pwd)/output:/output forefire \
  bash -c "cd tests/python && python3 predict_spread.py && cp fire_prediction.html fire_prediction.geojson /output/"
```

## Scripts

### predict_spread.py

Quick fire spread prediction for a given lat/lon ignition point using Rothermel on flat terrain with uniform fuel and constant wind. Outputs perimeter stats at 30 min and 2 h, writes `fire_prediction.geojson`, and generates `fire_prediction.html` — an interactive CartoDB map with the fire spread polygons drawn on it.

```bash
docker run --rm forefire python3 tests/python/predict_spread.py
```

Edit the constants at the top of the file to customize:
- `IGNITION_LAT` / `IGNITION_LON` — ignition coordinates
- `WIND_SPEED_MS` — wind speed in m/s
- `WIND_DIRECTION_DEG` — meteorological wind direction (wind FROM)

> ⚠️ The simulation is currently location-agnostic: fuel is uniform (Rothermel type 6) and terrain is flat. Lat/lon only positions the result on the map. To make location matter, plug a real fuel map and DEM into the `fuel_map` and `altitude` arrays.

To run with a one-off override (without editing the file in place), mount a patched copy on top of the one inside the image:

```bash
docker run --rm \
  -v $(pwd)/output:/output \
  -v $(pwd)/tests/python/predict_spread.py:/forefire/tests/python/predict_spread.py \
  forefire \
  bash -c "cd tests/python && python3 predict_spread.py && cp fire_prediction.html fire_prediction.geojson /output/"
```

### streamlit_app.py

Interactive Streamlit UI that wraps `predict_spread.py`. The sidebar exposes ignition coordinates, wind speed and wind direction; on submit it patches a temporary copy of `predict_spread.py`, runs it inside the `forefire` Docker image, and renders the resulting GeoJSON perimeters on a CartoDB map.

Prerequisites:
- The `forefire` Docker image must be built (`docker build -t forefire .` from the repo root)
- `streamlit` and `folium` installed on the host (the app shells out to `docker`, so it must run on the host — not inside a container)

```bash
pip install streamlit folium
streamlit run tests/python/streamlit_app.py
```

Then open http://localhost:8501.

### farsite_flat.py

Runs a simulation based on Farsite software for a north-wind of 3 mph.

Required input `.lcp` file can be downloaded from:
https://github.com/mbedward/farsite/raw/refs/heads/master/examples/flatland/Inputs/a_lcpFiles/flatland.lcp
