# ForeFire Spread Predictor

Interactive fire-spread prediction app built on top of the ForeFire C++
simulation engine. Two entry points:

- **`predict_spread.py`** — standalone Python script that runs a Rothermel
  simulation for a configured lat/lon, wind, and fuel, then dumps a GeoJSON
  perimeter (and optionally an HTML map if `folium` is available).
- **`streamlit_app.py`** — Streamlit UI that exposes ignition coordinates,
  wind speed and direction in the sidebar, patches a temporary copy of
  `predict_spread.py`, runs it inside the `forefire` Docker image, and renders
  the resulting perimeters on a CartoDB map.

> ⚠️ The simulation is currently location-agnostic: fuel is uniform
> (Rothermel type 6) and terrain is flat. Lat/lon only positions the result on
> the map. To make location matter, plug a real fuel map and DEM into the
> `fuel_map` and `altitude` arrays in `predict_spread.py`.

## Quick start

```bash
# 1. Build the forefire Docker image (only needed once)
make build

# 2. Install host dependencies
pip install streamlit folium

# 3. Launch the app
make app
```

Then open http://localhost:8501.

Each click of **Run Simulation** spawns a fresh `docker run --rm forefire ...`,
which keeps the C++ engine state clean between runs.

## Running `predict_spread.py` directly

Edit the constants at the top of `predict_spread.py`
(`IGNITION_LAT`, `IGNITION_LON`, `WIND_SPEED_MS`, `WIND_DIRECTION_DEG`) and
run it inside the Docker image:

```bash
docker run --rm \
  -v $(pwd)/output:/output \
  -v $(pwd)/apps/spread_predictor/predict_spread.py:/work/predict_spread.py \
  -w /work \
  forefire \
  bash -c "python3 predict_spread.py && cp fire_prediction.geojson /output/"
```

The script writes `fire_prediction.geojson` (always) and `fire_prediction.html`
(only when `folium` is installed — it isn't in the slim `forefire` image).

## Files

| File | Purpose |
|------|---------|
| `predict_spread.py` | Standalone ForeFire Rothermel simulation script |
| `streamlit_app.py`  | Streamlit UI wrapper |
| `Makefile`          | `make build` / `make app` shortcuts |
