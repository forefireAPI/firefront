#!/usr/bin/env python3
"""
Streamlit app to run fire spread predictions using ForeFire.
Runs predict_spread.py inside the `forefire` Docker image, with lat/lon/wind
substituted into a temporary copy of the script.
"""

import base64
import json
import os
import re
import subprocess
import tempfile

import folium
import streamlit as st

st.set_page_config(page_title="ForeFire Spread Predictor", layout="wide")
st.title("ForeFire Fire Spread Predictor")

APP_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_SCRIPT = os.path.join(APP_DIR, "predict_spread.py")

# --- Sidebar controls ---
with st.sidebar.form("simulation_form"):
    st.header("Ignition Point")
    latlon_input = st.text_input("Lat, Lon", value="44.58836, 4.22524")

    st.header("Wind")
    wind_speed = st.slider("Wind speed (m/s)", 0.0, 30.0, 10.0, 0.5)
    wind_direction = st.slider(
        "Wind FROM direction (°)", 0, 360, 25, 5,
        help="Meteorological convention: 270° = wind from the west",
    )

    submitted = st.form_submit_button("Run Simulation", type="primary")

# Parse lat/lon
try:
    _parts = [p.strip() for p in latlon_input.split(",")]
    ignition_lat, ignition_lon = float(_parts[0]), float(_parts[1])
except (ValueError, IndexError):
    st.error("Enter coordinates as: lat, lon")
    st.stop()


def patch_script(template_path: str, lat: float, lon: float, wspeed: float, wdir: float) -> str:
    """Read predict_spread.py and substitute the configuration constants."""
    with open(template_path, "r") as f:
        src = f.read()
    src = re.sub(r"^IGNITION_LAT\s*=.*$", f"IGNITION_LAT = {lat}", src, count=1, flags=re.M)
    src = re.sub(r"^IGNITION_LON\s*=.*$", f"IGNITION_LON = {lon}", src, count=1, flags=re.M)
    src = re.sub(r"^WIND_SPEED_MS\s*=.*$", f"WIND_SPEED_MS = {wspeed}", src, count=1, flags=re.M)
    src = re.sub(r"^WIND_DIRECTION_DEG\s*=.*$", f"WIND_DIRECTION_DEG = {wdir}", src, count=1, flags=re.M)
    return src


# --- Run simulation via docker ---
if submitted:
    with st.spinner("Running ForeFire simulation in Docker..."):
        try:
            patched = patch_script(TEMPLATE_SCRIPT, ignition_lat, ignition_lon, wind_speed, wind_direction)

            with tempfile.TemporaryDirectory() as tmpdir:
                script_path = os.path.join(tmpdir, "predict_spread.py")
                output_dir = os.path.join(tmpdir, "output")
                os.makedirs(output_dir, exist_ok=True)
                with open(script_path, "w") as f:
                    f.write(patched)

                cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{output_dir}:/output",
                    "-v", f"{script_path}:/work/predict_spread.py",
                    "-w", "/work",
                    "forefire",
                    "bash", "-c",
                    "python3 predict_spread.py && cp fire_prediction.geojson /output/",
                ]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

                if result.returncode != 0:
                    st.error(f"Simulation failed:\n{result.stderr[-800:]}")
                else:
                    geojson_path = os.path.join(output_dir, "fire_prediction.geojson")
                    with open(geojson_path, "r") as f:
                        geojson = json.load(f)
                    st.session_state["results"] = {
                        "geojson": geojson,
                        "stdout": result.stdout,
                        "params": {
                            "ignition_lat": ignition_lat,
                            "ignition_lon": ignition_lon,
                            "wind_speed_ms": wind_speed,
                            "wind_direction_deg": wind_direction,
                        },
                    }
        except Exception as e:
            st.error(f"Error: {e}")


# --- Display results ---
if "results" in st.session_state:
    results = st.session_state["results"]
    geojson = results["geojson"]
    params = results["params"]

    st.subheader("Results")
    with st.expander("Simulation log"):
        st.code(results["stdout"])

    # Build map
    st.subheader("Fire Spread Map")
    snapshot_colors = {"30min": "#FFA500", "2h": "#CC0000"}
    lat = params["ignition_lat"]
    lon = params["ignition_lon"]
    m = folium.Map(
        location=[lat, lon], zoom_start=13,
        tiles="https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
        attr='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/">CARTO</a>',
    )
    folium.Marker(
        [lat, lon], popup="Ignition point",
        icon=folium.Icon(color="red", icon="fire", prefix="fa"),
    ).add_to(m)
    for feat in geojson.get("features", []):
        label = feat["properties"]["label"]
        coords = feat["geometry"]["coordinates"][0]
        latlngs = [[la, lo] for lo, la in coords]
        color = snapshot_colors.get(label, "#FF6600")
        folium.Polygon(
            locations=latlngs, color=color, fill=True, fill_opacity=0.25,
            popup=f"Fire perimeter at {label}", tooltip=label,
        ).add_to(m)
    map_html = m.get_root().render()
    map_data_uri = "data:text/html;base64," + base64.b64encode(map_html.encode("utf-8")).decode("ascii")
    st.iframe(map_data_uri, height=500)

    st.download_button(
        "Download GeoJSON",
        data=json.dumps(geojson, indent=2),
        file_name="fire_prediction.geojson",
        mime="application/geo+json",
    )

    with st.expander("Simulation parameters"):
        st.json(params)
