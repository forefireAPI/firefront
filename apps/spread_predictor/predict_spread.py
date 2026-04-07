#!/usr/bin/env python3
"""
Quick fire spread prediction for a given lat/lon ignition point.
Uses idealized flat terrain with uniform fuel and constant wind.

WARNING: This is a rough estimate. For real predictions you need:
  - Actual fuel map (vegetation type) from land cover data
  - Digital elevation model (DEM)
  - Real-time wind data (speed + direction)
  - Fuel moisture conditions
"""

import math
import numpy as np
import json

import pyforefire as pyff

try:
    import folium
except ImportError:
    folium = None

# --- Configuration ---
IGNITION_LAT = 44.58836445769429
IGNITION_LON = 4.225245582315042
WIND_SPEED_MS = 10.0       # m/s (~18 km/h)
WIND_DIRECTION_DEG = 25   # degrees, meteorological (wind FROM west)

# Domain: ~20km x 20km around ignition (large enough for 2h spread)
DOMAIN_HALF_KM = 10.0

# --- Coordinate math ---
M_PER_DEG_LAT = 111195.0
M_PER_DEG_LON = M_PER_DEG_LAT * math.cos(math.radians(IGNITION_LAT))

lat_offset = DOMAIN_HALF_KM * 1000 / M_PER_DEG_LAT
lon_offset = DOMAIN_HALF_KM * 1000 / M_PER_DEG_LON

sw_lon = IGNITION_LON - lon_offset
sw_lat = IGNITION_LAT - lat_offset
ne_lon = IGNITION_LON + lon_offset
ne_lat = IGNITION_LAT + lat_offset

domain_width = (ne_lon - sw_lon) * M_PER_DEG_LON
domain_height = (ne_lat - sw_lat) * M_PER_DEG_LAT
center_x = domain_width / 2
center_y = domain_height / 2

# Wind vector (meteorological: "from" direction, so flip 180 for propagation)
wind_rad = math.radians(WIND_DIRECTION_DEG + 180)
wind_u = WIND_SPEED_MS * math.sin(wind_rad)
wind_v = WIND_SPEED_MS * math.cos(wind_rad)

# --- ForeFire setup ---
ff = pyff.ForeFire()

ff["propagationModel"] = "Rothermel"
ff["fuelsTable"] = pyff.helpers.extendedRothermelFuelTable()
ff["spatialIncrement"] = 1.0
ff["perimeterResolution"] = 10
ff["minimalPropagativeFrontDepth"] = 20
ff["relax"] = 0.5
ff["bmapLayer"] = 1
ff["windReductionFactor"] = 0.4

grid_n = 200

ff.execute(
    f"FireDomain[sw=(0,0,0);ne=({domain_width:.0f},{domain_height:.0f},0);"
    f"t=0;BBoxWSEN=({sw_lon},{sw_lat},{ne_lon},{ne_lat})]"
)

ff.addLayer("propagation", "Rothermel", "propagationModel")

# Uniform fuel type 6 (dormant brush/hardwood slash — moderate spread)
fuel_map = np.full((1, 1, grid_n, grid_n), 6, dtype=np.float64)
ff.addIndexLayer("table", "fuel", 0, 0, 0, domain_width, domain_height, 0, fuel_map)

# Flat terrain
altitude = np.zeros((1, 1, grid_n, grid_n))
ff.addScalarLayer("data", "altitude", 0, 0, 0, domain_width, domain_height, 0, altitude)

# Wind layers
windU = np.zeros((1, 2, grid_n, grid_n))
windU[0, 0, :, :] = 1.0
windV = np.zeros((1, 2, grid_n, grid_n))
windV[0, 0, :, :] = 0.0
ff.addScalarLayer("windScalDir", "windU", 0, 0, 0, domain_width, domain_height, 0, windU)
ff.addScalarLayer("windScalDir", "windV", 0, 0, 0, domain_width, domain_height, 0, windV)

# Ignite at center
ff.execute(f"startFire[loc=({center_x:.1f},{center_y:.1f},0);t=0]")
ff.execute(f"trigger[wind;loc=(0.,0.,0.);vel=({wind_u:.2f},{wind_v:.2f},0);t=0]")


def extract_perimeter(output):
    """Parse print[] output into list of (lon, lat) coordinates."""
    nodes = []
    for line in output.splitlines():
        if "FireNode" not in line:
            continue
        loc_start = line.find("loc=(")
        if loc_start == -1:
            continue
        loc_start += 5
        loc_end = line.find(")", loc_start)
        parts = line[loc_start:loc_end].split(",")
        if len(parts) >= 2:
            x, y = float(parts[0]), float(parts[1])
            if math.isnan(x) or math.isnan(y):
                continue
            lon = sw_lon + x / M_PER_DEG_LON
            lat = sw_lat + y / M_PER_DEG_LAT
            nodes.append((lon, lat))
    return nodes


def perimeter_stats(nodes):
    """Compute extent of fire perimeter."""
    if not nodes:
        return {}
    lons = [n[0] for n in nodes]
    lats = [n[1] for n in nodes]
    width_m = (max(lons) - min(lons)) * M_PER_DEG_LON
    height_m = (max(lats) - min(lats)) * M_PER_DEG_LAT
    area_ha = width_m * height_m / 10000  # rough bbox area in hectares
    return {
        "nodes": len(nodes),
        "width_m": round(width_m),
        "height_m": round(height_m),
        "bbox_area_ha": round(area_ha, 1),
        "lon_range": (round(min(lons), 5), round(max(lons), 5)),
        "lat_range": (round(min(lats), 5), round(max(lats), 5)),
    }


def to_geojson(nodes, label):
    """Convert perimeter nodes to a GeoJSON Polygon."""
    if len(nodes) < 3:
        return None
    coords = [[lon, lat] for lon, lat in nodes]
    coords.append(coords[0])  # close the polygon
    return {
        "type": "Feature",
        "properties": {"label": label},
        "geometry": {"type": "Polygon", "coordinates": [coords]},
    }


# --- Run predictions ---
print(f"Ignition: {IGNITION_LAT}, {IGNITION_LON}")
print(f"Wind: {WIND_SPEED_MS} m/s from {WIND_DIRECTION_DEG}° (u={wind_u:.1f}, v={wind_v:.1f})")
print(f"Fuel: Rothermel type 6 (dormant brush), flat terrain")
print(f"Domain: {domain_width:.0f} x {domain_height:.0f} m")
print()

features = []

# 30 minutes
ff.execute("step[dt=1800]")
out_30m = ff.execute("print[]")
nodes_30m = extract_perimeter(out_30m)
stats_30m = perimeter_stats(nodes_30m)
print(f"=== 30 minutes ===")
print(f"  Fire perimeter: {stats_30m.get('nodes', 0)} nodes")
print(f"  Extent: {stats_30m.get('width_m', 0)} x {stats_30m.get('height_m', 0)} m")
print(f"  Bounding box area: ~{stats_30m.get('bbox_area_ha', 0)} ha")
print(f"  Lon range: {stats_30m.get('lon_range')}")
print(f"  Lat range: {stats_30m.get('lat_range')}")
feat = to_geojson(nodes_30m, "30min")
if feat:
    features.append(feat)

# Continue to 2 hours (need 90 more minutes)
ff.execute("step[dt=5400]")
out_2h = ff.execute("print[]")
nodes_2h = extract_perimeter(out_2h)
stats_2h = perimeter_stats(nodes_2h)
print(f"\n=== 2 hours ===")
print(f"  Fire perimeter: {stats_2h.get('nodes', 0)} nodes")
print(f"  Extent: {stats_2h.get('width_m', 0)} x {stats_2h.get('height_m', 0)} m")
print(f"  Bounding box area: ~{stats_2h.get('bbox_area_ha', 0)} ha")
print(f"  Lon range: {stats_2h.get('lon_range')}")
print(f"  Lat range: {stats_2h.get('lat_range')}")
feat = to_geojson(nodes_2h, "2h")
if feat:
    features.append(feat)

# Save GeoJSON
if features:
    geojson = {
        "type": "FeatureCollection",
        "features": features,
    }
    with open("fire_prediction.geojson", "w") as f:
        json.dump(geojson, f, indent=2)
    print(f"\nGeoJSON saved to fire_prediction.geojson")

# Build CartoDB HTML map (only if folium is available)
if folium is not None:
    m = folium.Map(
        location=[IGNITION_LAT, IGNITION_LON],
        zoom_start=13,
        tiles="https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
        attr='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/">CARTO</a>',
    )

    folium.Marker(
        [IGNITION_LAT, IGNITION_LON],
        popup="Ignition point",
        icon=folium.Icon(color="red", icon="fire", prefix="fa"),
    ).add_to(m)

    colors = {"30min": "orange", "2h": "red"}
    for feat in features:
        label = feat["properties"]["label"]
        coords = feat["geometry"]["coordinates"][0]
        # folium expects (lat, lon)
        latlngs = [[lat, lon] for lon, lat in coords]
        folium.Polygon(
            locations=latlngs,
            color=colors.get(label, "red"),
            fill=True,
            fill_opacity=0.3,
            popup=f"Fire perimeter at {label}",
        ).add_to(m)

    map_file = "fire_prediction.html"
    m.save(map_file)
    print(f"Map saved to {map_file}")
