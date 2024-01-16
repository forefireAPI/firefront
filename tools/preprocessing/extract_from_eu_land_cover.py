#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:38:26 2023

@author: filippi_j
"""

import rasterio
from rasterio.windows import from_bounds
from pyproj import Transformer
import matplotlib.pyplot as plt
import numpy as np
from pyproj import CRS
import osmnx as ox  
from rasterio.features import geometry_mask
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import box
import json
from PIL import Image
import rasterio
from rasterio.warp import calculate_default_transform, reproject, transform_bounds, Resampling
import numpy as np
from .ffToGeoJson import create_kml,generate_indexed_png_and_legend


attribute_widths_road_edge_full = {
    'secondary': 0.8,
    'track': 0.8,
    'tertiary': 0.8,
    'path': 0.8,
    'trunk': 0.8,
    'unclassified': 0.8,
    'service': 0.8,
    'living_street': 0.8,
    "['living_street', 'pedestrian']": 0.8,
    "['living_street', 'path']": 0.8,
    'residential': 0.8,
    "['track', 'residential']": 0.8,
    "['track', 'unclassified']": 0.8,
    'primary': 0.8,
    'trunk_link': 0.8,
    'primary_link': 0.8,
    "['steps', 'footway']": 0.8,
    "['unclassified', 'residential']": 0.8,
    'steps': 0.8,
    'pedestrian': 0.8,
    'secondary_link': 0.8,
    "['path', 'service']": 0.8,
    "['unclassified', 'path']": 0.8,
    'footway': 0.8,
    'tertiary_link': 0.8,
    "['pedestrian', 'footway']": 0.8,
    "['track', 'service']": 0.8,
    "['secondary', 'tertiary']": 0.8,
    "['residential', 'tertiary']": 0.8,
    "['pedestrian', 'residential']": 0.8,
    "['track', 'path']": 0.8,
    "['footway', 'residential']": 0.8,
    "['path', 'residential']": 0.8,
    "['track', 'path', 'residential']": 0.8,
    "['service', 'residential']": 0.8,
    "['path', 'pedestrian']": 0.8,
    "['service', 'pedestrian']": 0.8,
    "['track', 'unclassified', 'residential']": 0.8,
    "['residential', 'footway']": 0.8,
    "['unclassified', 'service']": 0.8,
    'cycleway': 0.8,
    "['steps', 'path', 'footway']": 0.8,
    "['service', 'footway']": 0.8,
    "['track', 'service', 'residential']": 0.8,
    "['steps', 'pedestrian', 'footway']": 0.8,
    "['service', 'steps', 'footway']": 0.8,
    "['path', 'footway']": 0.8,
}
attribute_widths_road_Edge_heavu = {
    'secondary': 3,
    'track': 1,
    'tertiary': 2,
    'path': 1,
    'trunk': 1,
    'service': 1,
    'living_street': 1,
    'residential': 1,
    'primary': 5,
    'trunk_link': 1,
    'primary_link': 3,
    'secondary_link': 2,
    'tertiary_link': 1,
    'cycleway': 1,
}
attribute_widths_road_Edge = {
    'secondary': 1,
    'track': 0.5,
    'tertiary': 0.7,
    'primary': 1.5,
}

def extract_subregion(input_tif, output_tif, westI, southI, eastI, northI):
    print(f"clipping {input_tif} into {output_tif} bounds W{westI}, S{southI}, E{eastI}, N{northI}")
    print(f"NW({northI},{westI}), SW({southI},{westI}), SE({southI},{eastI}), NE({northI},{eastI})")
    
    with rasterio.open(input_tif) as src:
        # Transformer pour convertir les coordonnées lat/lon (WGS84) en coordonnées du système du raster
        transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        sw = transformer.transform( westI,southI)
        nw = transformer.transform( westI,northI)
        ne = transformer.transform( eastI,northI)
        se = transformer.transform( eastI,southI)
       
    #    print(sw,nw,se,ne)
        # Convertir les coordonnées
        west  = np.min([sw[0],nw[0]])
        north =  np.max([nw[1],ne[1]])
        east = np.max([se[0],ne[0]])
        south = np.min([sw[1],se[1]])
        
        print(west, south, east, north) 
        
        # Vérifications supplémentaires
        if west >= east or south >= north:
            print("Invalid bounds.",  west, south, east, north)
            return
        
        if west < src.bounds.left or east > src.bounds.right or \
           south < src.bounds.bottom or north > src.bounds.top:
            print("Bounds are outside the raster.")
            return
        
        # Créer une fenêtre (window) pour la sous-région
        window = from_bounds(west, south, east, north, src.transform)
        
        # Lire la sous-région
        data = src.read(1, window=window)
        
        # Calculer la nouvelle transformation affine pour la sous-région
        new_transform = src.window_transform(window)
        
        
        
        
        # Écrire la sous-région dans un nouveau fichier TIF
        with rasterio.open(output_tif, 'w', driver='GTiff',
                           height=data.shape[0], width=data.shape[1],
                           count=1, dtype=data.dtype,
                           crs=src.crs, transform=new_transform) as dst:
            dst.write(data, 1)



def warp_and_clip_raster(input_file_path: str, warped_file_path: str, output_file_path: str, west, south, east, north , target_width=2000, target_height = 2000):
    """
    Warps and clips a raster file based on the given bounds.

    Parameters:
        input_file_path (str): The path to the input raster file.
        output_file_path (str): The path to save the output raster file.
        bounds (dict): Dictionary containing the west, south, east, and north bounds.

    Returns:
        None: The function saves the output to a file.
    """
    # Define target CRS
    target_crs = 'EPSG:4326'
    
    # Warp raster to target CRS
    with rasterio.open(input_file_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, target_crs, src.width, src.height, *src.bounds)
        
        out_meta = src.meta.copy()
        out_meta.update({
            'crs': target_crs,
            'transform': transform,
            'width': int(width),
            'height': int(height)
        })

        temp_file_path = warped_file_path
        with rasterio.open(temp_file_path, 'w', **out_meta) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.nearest)
    
    # Create bounding box and GeoDataFrame 
    bbox = box(west, south, east, north)
    geo = gpd.GeoDataFrame({'geometry': [bbox]}, index=[0], crs="EPSG:4326")
    
    # Clip the raster
    with rasterio.open(temp_file_path) as src:
        coords = [json.loads(geo.to_json())['features'][0]['geometry']]
        out_img, out_transform = mask(src, shapes=coords, crop=True)
        out_meta = src.meta.copy()
        
        out_meta.update({
            "height": target_height,
            "width": target_width,
            "transform": out_transform,
        })

    # Save the clipped raster
    with rasterio.open(output_file_path, "w", **out_meta) as dest:
        dest.write(out_img)
    


def extract_roads_from_geotiff(west,south,east,north, road_shape_file):
    print(f"Getting roads with bounds from {west,south,east,north} saving {road_shape_file}")
    try:
        G = ox.graph_from_bbox(north, south, east, west, network_type='all')
    except Exception as e:
        return str(e)
    ox.save_graph_shapefile(G, filepath=road_shape_file )

def rasterize_shapefile(shapefile_path, ref_tif, output_path, attribute_widths=attribute_widths_road_Edge, default_width = 0.3, imbounds =None,code=62):
    print(f"Rasterizing roads {shapefile_path} to {output_path}  with {attribute_widths}")
    reference_properties = {}
    width=2000
    height = 2000
    
    gdf = gpd.read_file(shapefile_path)
    bounds = gdf.total_bounds  # minx, miny, maxx, maxy
   
 

    if imbounds is not None:
        bounds=imbounds
        
        
    x_min, y_min, x_max, y_max = bounds
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    refT =  rasterio.open(ref_tif) 
    ref_crs = refT.crs
    width = refT.width
    height = refT.height  

    # Calculate the resolution
    resolutionx = x_range / width
    resolutiony = y_range / height 
    # Create a transform
    transform = from_origin(x_min, y_max, resolutionx, resolutiony)
    
    # Create an empty raster
    raster = refT.read(1)
    refT.close()
    
    # Rasterize each feature based on the "highway" attribute
    for _, feature in gdf.iterrows():
        highway_type = feature['highway']
        line_width = attribute_widths.get(highway_type, default_width)  # Default line width is 1 pixel
        mask = geometry_mask([feature['geometry'].buffer(line_width * resolutionx / 2)],
                             transform=transform, invert=True, out_shape=(height, width))
        raster[mask] = code  # Set pixel value to 255 (white) where there is a feature


    # Save the reprojected data to a new GeoTIFF file
    with rasterio.open(output_path, 'w', driver='GTiff',
                       width=width, height=height,
                       count=1, dtype=raster.dtype,
                       crs=ref_crs, transform=transform) as dst:
        dst.write(raster, 1)    

def fake_fuel(ref_tif, output_path,fuel_test,imbounds =None):
    print(f"faking fuels to {output_path}")
        
    bounds=imbounds
    x_min, y_min, x_max, y_max = bounds
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    refT =  rasterio.open(ref_tif) 
    ref_crs = refT.crs
    width = refT.width
    height = refT.height  

    # Calculate the resolution
    resolutionx = x_range / width
    resolutiony = y_range / height 
    # Create a transform
    transform = from_origin(x_min, y_max, resolutionx, resolutiony)
    
    # Create an empty raster
    raster = refT.read(1)
    n_bandes = len(fuel_test)
    limites = np.linspace(0, width, n_bandes + 1).astype(int)
    # Remplissage de fuelMap
    for i in range(n_bandes):
        raster[:, limites[i]:limites[i + 1]] = fuel_test[i]
        
    

    
    refT.close()

    # Save the reprojected data to a new GeoTIFF file
    with rasterio.open(output_path, 'w', driver='GTiff',
                       width=width, height=height,
                       count=1, dtype=raster.dtype,
                       crs=ref_crs, transform=transform) as dst:
        dst.write(raster, 1)    

def fake_fuel(legend_file_path, WSEN, LBRT,output_dir,fuel_resolution = 10):
    fuel_indices_warped_clipped = f"{output_dir}/fuel_indices_warped_clipped.tif"
    fuel_road_indices = f'{output_dir}/fuel.tif'
    fuel_kml = f'{output_dir}/fuel.kml'
    fuel_png =f'{output_dir}/fuel.png'
    fuel_cbar_png =f'{output_dir}/fuel_bar.png'
      
    west, south, east, north =WSEN
    l,b,r,t = LBRT
    
    fake_fuel(fuel_indices_warped_clipped, fuel_road_indices,fuel_test,imbounds=(west, south,east ,north))   
    generate_indexed_png_and_legend(legend_file_path,fuel_road_indices, fuel_png, fuel_cbar_png)
    create_kml(west, south, east, north, "FUEL", fuel_png,fuel_kml , pngcbarfile=fuel_cbar_png) 
    
    

def landcover_roads_to_fuel(S2GLC_tif,legend_file_path, WSEN, LBRT,output_dir,fuel_resolution = 10):

    fuel_indices_origin = f"{output_dir}/fuel_indices_S2GLC.tif"
    fuel_indices_warped = f"{output_dir}/fuel_indices_warped.tif"
    fuel_indices_warped_clipped = f"{output_dir}/fuel_indices_warped_clipped.tif"
    
    roads_shape = f"{output_dir}/roads_shape.shp"
    
    fuel_road_indices = f'{output_dir}/fuel.tif'
     
    
    fuel_kml = f'{output_dir}/fuel.kml'
    fuel_png =f'{output_dir}/fuel.png'
    
    fuel_cbar_png =f'{output_dir}/fuel_bar.png'
    
    
      
    west, south, east, north =WSEN
    l,b,r,t = LBRT
    
        

    
    
    extract_subregion(S2GLC_tif, fuel_indices_origin, west, south, east, north)
    warp_and_clip_raster(fuel_indices_origin,fuel_indices_warped,fuel_indices_warped_clipped ,west, south, east, north,target_width = int((r-l)/fuel_resolution),target_height = int((t-b)/fuel_resolution))
 
    extract_roads_from_geotiff(west, south, east, north, roads_shape)
    rasterize_shapefile(roads_shape,fuel_indices_warped_clipped,fuel_road_indices,imbounds=(west, south,east ,north),code=62)

    generate_indexed_png_and_legend(legend_file_path,fuel_road_indices, fuel_png, fuel_cbar_png)
    create_kml(west, south, east, north, "FUEL", fuel_png,fuel_kml , pngcbarfile=fuel_cbar_png)
