# Data Sources for Fire Spread Prediction

ForeFire requires three input layers: terrain elevation (DEM), fuel/vegetation type, and wind. All can be sourced globally for free.

## DEM (Digital Elevation Model)

| Source | Resolution | Coverage | Auth | Format | URL |
|--------|-----------|----------|------|--------|-----|
| Copernicus GLO-30 | 30m | Global | None | GeoTIFF | https://registry.opendata.aws/copernicus-dem/ |
| SRTM | 30m | 60N-56S | NASA Earthdata login | GeoTIFF/HGT | https://earthexplorer.usgs.gov/ |
| IGN BD ALTI | 25m | France | IGN key | ASC/GeoTIFF | https://geoservices.ign.fr/ |

**Usage in ForeFire**: load into the `altitude` scalar layer as a 2D array (meters above sea level).

## Fuel / Vegetation

| Source | Resolution | Coverage | Auth | Format | URL |
|--------|-----------|----------|------|--------|-----|
| ESA WorldCover | 10m | Global | None | GeoTIFF | https://esa-worldcover.org/ |
| Corine Land Cover | 100m | Europe | None | GeoTIFF | https://land.copernicus.eu/pan-european/corine-land-cover |
| OSO (Theia) | 10m | France | None | GeoTIFF | https://theia.cnes.fr/ |

**Usage in ForeFire**: land cover classes must be remapped to Rothermel fuel model numbers, then loaded into the `fuel` index layer.

### Land cover to Rothermel fuel model mapping (approximate)

| Land cover class | Rothermel fuel model | Description |
|-----------------|---------------------|-------------|
| Tree cover (broadleaf) | 8 | Compact timber litter |
| Tree cover (needleleaf) | 9 | Hardwood long-needle litter |
| Shrubland | 6 | Dormant brush / hardwood slash |
| Grassland | 1 | Short grass (< 0.3m) |
| Cropland | 3 | Tall grass (1m) |
| Built-up / urban | 0 | Non-burnable |
| Water / wetland | 0 | Non-burnable |
| Bare / rock | 0 | Non-burnable |

> This mapping is approximate and region-dependent. Mediterranean shrublands behave differently from boreal shrublands. For production use, a locally-validated fuel map is strongly recommended.

## Wind

| Source | Resolution | Coverage | Auth | Frequency | URL |
|--------|-----------|----------|------|-----------|-----|
| Open-Meteo | ~10km | Global | None | Hourly | https://open-meteo.com/ |
| GFS (NOAA) | 25km | Global | None | 3-hourly | https://nomads.ncep.noaa.gov/ |
| AROME (Meteo-France) | 1.3km | France | API key | Hourly | https://portail-api.meteofrance.fr/ |

**Usage in ForeFire**: wind is provided as U (east-west) and V (north-south) components in m/s, loaded into `windU` / `windV` scalar layers. Can be constant or time-varying via `trigger[wind;...]` commands.

### Open-Meteo example query

```
https://api.open-meteo.com/v1/forecast?latitude=48.26&longitude=2.70&hourly=wind_speed_10m,wind_direction_10m&forecast_days=1
```

Returns hourly wind speed (m/s) and direction (degrees) at 10m height. Convert to U/V:

```python
import math
wind_u = -wind_speed * math.sin(math.radians(wind_direction))
wind_v = -wind_speed * math.cos(math.radians(wind_direction))
```

## What the current script uses

`predict_spread.py` does **not** use any real data. All layers are faked:

- **DEM**: `np.zeros(...)` (flat terrain)
- **Fuel**: `np.full(..., 6)` (uniform fuel type 6 everywhere)
- **Wind**: hardcoded constant (5 m/s from west)
