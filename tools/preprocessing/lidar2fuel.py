#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 17:45:24 2023

@author: filippi_j
"""
import numpy as np
import matplotlib.pyplot as plt
 
import os.path
import struct
import laspy 
import geopandas as gpd
from fiona.crs import from_epsg
from shapely.geometry import box
 
import rasterio
import requests

from rasterio.warp import reproject, Resampling
 
import vtk
from PIL import Image
# Convertir le point SW
import contextily as cx
import json
from rasterio.mask import mask
import pycrs 
import zipfile
import io


def gen_density_altitude_and_filter(in_points):
    filtered_points['heightStat'] = np.array([]).reshape(0, 3)
    
    # Dictionnaire pour stocker les points par carré
    squares = defaultdict(list)
    
    # Partitionnement des points dans les carrés
    for i, (x, y, z) in enumerate(filtered_points['coords']):
        square_x = int(x // square_size)
        square_y = int(y // square_size)
        squares[(square_x, square_y)].append(i)
    
    # Listes pour stocker les résultats et les indices des points centraux
    height_stats = []
    central_indices = []
    
    # Calculs pour chaque carré
    for (square_x, square_y), indices in squares.items():
        altitudes = [filtered_points['coords'][i, 2] for i in indices]
        
        # 1) Différence de hauteur entre le point le plus haut et le plus bas
        max_z = max(altitudes)
        min_z = min(altitudes)
        height_diff = max_z - min_z
        
        # 2) Moyenne des 70% des points les plus hauts en retranchant l'altitude du point le plus bas
        sorted_altitudes = sorted(altitudes, reverse=True)
        n_top_70 = int(0.7 * len(sorted_altitudes))
        top_70_mean = np.mean(sorted_altitudes[:n_top_70]) - min_z
        
        # 3) Altitude du point médian - l'altitude du point le plus bas
        median_z = np.median(sorted_altitudes)
        median_diff = median_z - min_z
        
        # Stockage des résultats
        height_stats.append([height_diff, top_70_mean, median_diff])
        
        # Sélection du point central (point avec l'altitude médiane)
        central_idx = indices[altitudes.index(median_z)]
        central_indices.append(central_idx)
    
    # Mettre à jour 'heightStat' avec les statistiques de hauteur calculées
    filtered_points['heightStat'] = np.array(height_stats)
    
    # Filtrer les autres champs pour ne garder que les points centraux
    filtered_points['coords'] = filtered_points['coords'][central_indices]
    filtered_points['intensity'] = np.array(filtered_points['intensity'])[central_indices]
    filtered_points['color'] = filtered_points['color'][central_indices]
    
    return filtered_points

def filter_las_files_with_attributesO(file_paths, left, bottom, right, top):
    # Initialize the dictionary with specific structures for coordinates and any known attributes
    all_filtered_points = {
        'coords': np.array([]).reshape(0, 3),
        'color': np.array([]).reshape(0, 3)
    }

    for file_path in file_paths:
        inFile = laspy.read(file_path)

        # Calculate the mask for the bounding box
        points = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
        mask = (points[:, 0] >= left) & (points[:, 0] <= right) & (points[:, 1] >= bottom) & (points[:, 1] <= top)
        print("Got ", np.sum(mask), " from file ", file_path)

        # Handle coordinates, intensity, and color specifically
        filtered_points = points[mask]
        all_filtered_points['coords'] = np.vstack((all_filtered_points['coords'], filtered_points))

        filtered_intensity = inFile.intensity[mask]
        all_filtered_points['intensity'].extend(filtered_intensity)

        if hasattr(inFile, 'red') and hasattr(inFile, 'green') and hasattr(inFile, 'blue'):
            # Assuming color data is available and stored in 'red', 'green', 'blue'
            filtered_color = np.vstack((inFile.red[mask], inFile.green[mask], inFile.blue[mask])).transpose()
            all_filtered_points['color'] = np.vstack((all_filtered_points['color'], filtered_color))

        # Dynamically handle other attributes
        for attribute in inFile.point_format.dimension_names:
            if attribute not in ['x', 'y', 'z', 'intensity', 'red', 'green', 'blue']:  # Exclude already handled
                if attribute not in all_filtered_points:
                    all_filtered_points[attribute] = []

                filtered_attribute = getattr(inFile, attribute)[mask]
                all_filtered_points[attribute].extend(filtered_attribute)

    # Convert lists to numpy arrays for numeric data types
    for key in all_filtered_points:
        print(key, " attribute found in filtered laz file")
        if isinstance(all_filtered_points[key], list):  # Convert lists to numpy arrays
            all_filtered_points[key] = np.array(all_filtered_points[key])

    print("Filtered bounding box ", left, bottom, right, top, " from files ", file_paths)
    return all_filtered_points


def filter_las_files_with_attributes(file_paths, left, bottom, right, top, remove_ratio=0):
    # Initialize the dictionary with specific structures for coordinates and any known attributes
    all_filtered_points = {
        'coords': np.array([]).reshape(0, 3),  # Pre-initialized as an empty array for coords
        'color': np.array([]).reshape(0, 3),   # Pre-initialized as an empty array for color
        'intensity': []                        # Other attributes initialized as lists
    }

    for file_path in file_paths:
        inFile = laspy.read(file_path)

        # Calculate the mask for the bounding box
        points = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
        mask = (points[:, 0] >= left) & (points[:, 0] <= right) & (points[:, 1] >= bottom) & (points[:, 1] <= top)
        print("Got ", np.sum(mask), " points from file ", file_path)

        # Append filtered coordinates
        filtered_points = points[mask]
        all_filtered_points['coords'] = np.vstack((all_filtered_points['coords'], filtered_points))

        # Append filtered intensity, check and initialize if first access
        if 'intensity' in all_filtered_points:
            all_filtered_points['intensity'].extend(inFile.intensity[mask])
        else:
            all_filtered_points['intensity'] = list(inFile.intensity[mask])

        # Append filtered color, assuming RGB channels are available
        if hasattr(inFile, 'red') and hasattr(inFile, 'green') and hasattr(inFile, 'blue'):
            filtered_color = np.vstack((inFile.red[mask], inFile.green[mask], inFile.blue[mask])).transpose()
            all_filtered_points['color'] = np.vstack((all_filtered_points['color'], filtered_color))

        # Dynamically handle other attributes
        for attribute in inFile.point_format.dimension_names:
            if attribute not in ['x', 'y', 'z', 'intensity', 'red', 'green', 'blue']:
                if attribute not in all_filtered_points:
                    all_filtered_points[attribute] = []
                filtered_attribute = getattr(inFile, attribute)[mask]
                all_filtered_points[attribute].extend(filtered_attribute)

    # Convert lists to numpy arrays for all attributes collected as lists
    for key in all_filtered_points:
        if isinstance(all_filtered_points[key], list):
            all_filtered_points[key] = np.array(all_filtered_points[key])

    # Apply removal ratio
    if remove_ratio > 0 and all_filtered_points['coords'].size > 0:
        total_points = len(all_filtered_points['coords'])
        keep_indices = np.random.choice(total_points, int(total_points * (1 - remove_ratio)), replace=False)
        for key in all_filtered_points:
            if all_filtered_points[key].size > 0:  # Ensure the array is not empty
                all_filtered_points[key] = all_filtered_points[key][keep_indices]
            else:
                print(f"Warning: {key} attribute has no data to filter.")


    # Final print to confirm attributes and their filtered status
    for key in all_filtered_points:
        print(key, " attribute found in filtered laz file, total points after filtering: ", len(all_filtered_points[key]))

    print("Filtered bounding box ", left, bottom, right, top, " from files ", file_paths)
    return all_filtered_points


def filter_las_files_with_single_attributes(file_paths, left, bottom, right, top):
    all_filtered_points = {'coords': np.array([]).reshape(0, 3), 'intensity': [], 'color': np.array([]).reshape(0, 3)}

    for file_path in file_paths:
        inFile = laspy.read(file_path)

        points = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
        tr = np.max(points[:,0])
        tl = np.min(points[:,0])
        tt = np.max(points[:,1])
        tb = np.min(points[:,1]) 

        intensity = inFile.intensity
  

        mask = (points[:, 0] >= left) & (points[:, 0] <= right) & (points[:, 1] >= bottom) & (points[:, 1] <= top)
        print("Got ", np.sum(mask)," from file ",file_path )
        filtered_points = points[mask]
        filtered_intensity = intensity[mask]
        

        all_filtered_points['coords'] = np.vstack((all_filtered_points['coords'], filtered_points))
        
        all_filtered_points['intensity'].extend(filtered_intensity)


    all_filtered_points['intensity'] = np.array(all_filtered_points['intensity'])
 
    print("filtered box ", left, bottom, right, top," from files ", file_paths)
    return all_filtered_points

def save_filtered_las_color_classification_and_coords(all_filtered_points, output_file_path):
    # Set up LAS file with a point format that supports RGB colors and extended classification
    header = laspy.LasHeader(version="1.4", point_format=7)  # Point format 7 supports RGB and extended classification
    outFile = laspy.LasData(header)
    
    outFile.x = all_filtered_points['coords'][:, 0]
    outFile.y = all_filtered_points['coords'][:, 1]
    outFile.z = all_filtered_points['coords'][:, 2]
    outFile.intensity = all_filtered_points['intensity']
    
    # Set classification; make sure it is compatible with extended values
    outFile.classification = all_filtered_points['classification'].astype(np.uint8)
    
    # Set RGB colors; scale colors to 16-bit if necessary
    if np.max(all_filtered_points['color']) <= 255:
        scaleFactor = 65535 / 255
        outFile.red = (all_filtered_points['color'][:, 0] * scaleFactor).astype(np.uint16)
        outFile.green = (all_filtered_points['color'][:, 1] * scaleFactor).astype(np.uint16)
        outFile.blue = (all_filtered_points['color'][:, 2] * scaleFactor).astype(np.uint16)
    else:
        outFile.red = all_filtered_points['color'][:, 0].astype(np.uint16)
        outFile.green = all_filtered_points['color'][:, 1].astype(np.uint16)
        outFile.blue = all_filtered_points['color'][:, 2].astype(np.uint16)

    outFile.write(output_file_path)
    print("Saved LAS as ", output_file_path)

def save_filtered_las(all_filtered_points, output_file_path):
    # Create a new LAS file with appropriate header setup
    header = laspy.LasHeader(point_format=2, version="1.2")  # Assuming point format 2 for RGB color
    outFile = laspy.LasData(header)
    
    # Set x, y, z coordinates
    outFile.x = all_filtered_points['coords'][:, 0]
    outFile.y = all_filtered_points['coords'][:, 1]
    outFile.z = all_filtered_points['coords'][:, 2]
    
    # Remove 'coords' key as it's already processed
    del all_filtered_points['coords']
    
    # Set colors if available
    if 'color' in all_filtered_points:
        outFile.red = all_filtered_points['color'][:, 0]
        outFile.green = all_filtered_points['color'][:, 1]
        outFile.blue = all_filtered_points['color'][:, 2]
        del all_filtered_points['color']  # Remove 'color' key as it's already processed

    # Dynamically set other attributes
    for attribute, values in all_filtered_points.items():
        print("handling ",attribute)
        if hasattr(outFile, attribute):  # Check if the attribute exists in the LAS file structure
            setattr(outFile, attribute, np.array(values))

    # Write the file to disk
    outFile.write(output_file_path)
    print("Saved LAS as ", output_file_path)

    
def convert_LASCRS_to_latlon(filtered_points, lidarCRS = 2154):
    # Initialiser le transformateur pour convertir de Lambert-93 (EPSG:2154) à WGS84 (EPSG:4326)
    transformer = Transformer.from_crs(lidarCRS, 4326)
    
    # Extraire les coordonnées x, y, et z
    x_lambert = filtered_points['coords'][:, 0]
    y_lambert = filtered_points['coords'][:, 1]
    z_lambert = filtered_points['coords'][:, 2]
    
    # Effectuer la transformation
    lat, lon = transformer.transform(x_lambert, y_lambert)
    
    # Recréer le dictionnaire des points avec les nouvelles coordonnées
    filtered_points_latlon = {
        'coords': np.vstack((lat, lon, z_lambert)).T,
        'intensity': filtered_points['intensity'],
        'color': filtered_points['color']
    }
    print("Projected points to latlon from CRS ", lidarCRS)
    return filtered_points_latlon    
    

def toImage(red,green,blue, outF):
 
    
    # Convert the normalized bands to 8-bit (0-255)
    red_8bit = (red * 255).astype('uint8')
    green_8bit = (green * 255).astype('uint8')
    blue_8bit = (blue * 255).astype('uint8')
    
    # Create an RGB image
    rgb_array = np.stack((red_8bit, green_8bit, blue_8bit), axis=2)
    rgb_image = Image.fromarray(rgb_array, 'RGB')
    
    # Save the image
    print("saved RVB bands to ", outF)
    rgb_image.save(outF)

def assign_colors_from_tif(filtered_points, tif_path, bands=[1,2,3],scale=[255.0,255.0,255.0]):
    # Lire l'image TIF avec rasterio
    with rasterio.open(tif_path) as src:
 
        transform = src.transform

        inv_transform = ~transform
        
        red  = src.read(bands[0]).astype('float64')/scale[0]
        green  = src.read(bands[1]).astype('float64')/scale[0]
        blue  = src.read(bands[2]).astype('float64')/scale[0]
        
        # Normalisation des valeurs de pixel
    
        toImage(red,green,blue, tif_path+"testIMG.jpg")
        
        # Pour chaque point, trouver la valeur de pixel correspondante dans l'image TIF
        new_colors = []
        for point in filtered_points['coords']:
            # Convertir les coordonnées du point en coordonnées de pixel
            lamx, lamby = point[0], point[1]
            px, py = inv_transform * (lamx, lamby)
            col,row = int(px), int(py)
    
            # Récupérer et normaliser la valeur de pixel
            try:
                # Utiliser les bandes NIR (3), R (0) et V (1) pour créer un point RVB
                vr = (red[row, col]*255.0).astype(np.uint8)
                vv = (green[row, col]*255.0).astype(np.uint8)
                vb = (blue[row, col]*255.0).astype(np.uint8)
                
                new_colors.append([vr,vv,vb])
            except IndexError:
                # Si le point est en dehors des limites de l'image, attribuer une couleur par défaut (par exemple, noir)
                print("error assigning tif color")
                new_colors.append([0, 0, 0])
    
        # Mettre à jour l'attribut de couleur des points
        filtered_points['color'] = np.array(new_colors)
        print("Assigned colors to points using ", tif_path)
        return filtered_points
        


# Ouvrir le fichier TIF source
def reprojectTif(inTif, outTif, destCRS=2154):
    with rasterio.open(inTif) as src:
        transform = src.transform
        crs = src.crs
        array = src.read()
    
        # Définir les métadonnées pour le fichier de destination
        dst_crs = destCRS  # WGS 84
        dst_transform, dst_width, dst_height = rasterio.warp.calculate_default_transform(
            crs, dst_crs, src.width, src.height, *src.bounds
        )
        dst_kwargs = src.meta.copy()
        dst_kwargs.update({
            'crs': dst_crs,
            'transform': dst_transform,
            'width': dst_width,
            'height': dst_height
        })
    
        # Créer le fichier TIF de destination
        with rasterio.open(outTif, 'w', **dst_kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=transform,
                    src_crs=crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest
                )
        print(inTif, " converted to crs ", destCRS, " saved to ", outTif)


def getLeftRightBottomTop(inTif):
    with rasterio.open(inTif) as src:
        
        wsen = src.bounds.left,  src.bounds.right, src.bounds.bottom, src.bounds.top
        print("got bounds of ", inTif, " as ", wsen)
        return wsen

def plot_colored_points_3d(points, dpi=300, filename="3D_colored_points.png", z_scale=1):
    coords = points['coords']
    colors = points['color'] / 255.0  # Assurez-vous que les couleurs sont normalisées entre 0 et 1
    
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')
    
    # Appliquer un facteur d'échelle à l'axe Z pour le rendre moins grand
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2] * z_scale, c=colors, s=0.05)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # Sauvegarder la figure en PNG avec une résolution spécifique (dpi)
    plt.savefig(filename, dpi=dpi)
    print("Plotted points in ", filename)
import vtk
import numpy as np

def save_points_to_vtk(points, filename='points.vtk'):
    # Create a vtkPoints object to store the coordinates
    vtk_points = vtk.vtkPoints()
    tb = np.min(points['coords'][:, 1])
    tl = np.min(points['coords'][:, 0])

    for coord in points['coords']:
        ncoord = coord - [tl, tb, 0]
        vtk_points.InsertNextPoint(ncoord)

    # Create a vtkCellArray to store the connectivity of the points
    vertices = vtk.vtkCellArray()
    vertices.InsertNextCell(len(points['coords']))
    for i in range(len(points['coords'])):
        vertices.InsertCellPoint(i)

    # Create a vtkPolyData object to store the points and their connectivity
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetVerts(vertices)

    # Create a vtkUnsignedCharArray to store the colors
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetName("Colors")
    for color in points['color']:
        colors.InsertNextTuple3(*color)
    
    # Add the colors to the polydata
    polydata.GetPointData().SetScalars(colors)

    # Create a vtkUnsignedCharArray for the classification
    classif = vtk.vtkUnsignedCharArray()
    classif.SetNumberOfComponents(1)
    classif.SetName("Classification")
    for kindOf in points['classification']:
        classif.InsertNextValue(kindOf)
    
    # Add the classification to the polydata as an additional array
    polydata.GetPointData().AddArray(classif)

    # Write the data to a VTK file
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()
    print("Saved vtk file points in ", filename)


 

def webMapsToTif(west, south, east, north, outF, providerSRC=cx.providers.GeoportailFrance.orthos, zoomLevel=19):
    tempOUT = outF+"_temp.tif"
 
 
    
    cx.bounds2raster(west, south, east, north ,
                         ll=True,
                         path=tempOUT,
                                         zoom=zoomLevel,
                                         source=cx.providers.GeoportailFrance.orthos
     
                        )
    
    
    data = rasterio.open(tempOUT) 
    
    bbox = box(west, south, east, north)
    
    geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=from_epsg(4326))
    geo = geo.to_crs(crs=data.crs.data)

    coords = [json.loads(geo.to_json())['features'][0]['geometry']]
      
    out_img, out_transform = mask(data, shapes=coords, crop=True)
    epsg_code = int(data.crs.data['init'][5:])
    out_meta = data.meta.copy()
    print(out_meta)
    
    
    out_meta.update({"driver": "GTiff",
                        "height": out_img.shape[1],
                        "width": out_img.shape[2],
                         "transform": out_transform,
                        "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4()}
                               )
    
    with rasterio.open(outF, "w", **out_meta) as dest:
        dest.write(out_img)
    
    print("Extracted image from contextily bounds:",west, south, east, north," zoom ", zoomLevel, " out files ",outF," and temporary ",tempOUT)

def findfileToDownload(west, south, east, north,laslist = [], lidarCRS=2154,step=1000,prefix="https://wxs.ign.fr/2s53j8r4gxyr0bfcounndiea/telechargement/prepackage/LIDARHD_PACK_IP_2021$LIDARHD_1-0_LAZ_IP-"):
    points = []
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:%d"%lidarCRS)
 
    
    # Convertir NE et SW en Lambert 93
    x_ne, y_ne = transformer.transform( north,east)
    x_sw, y_sw = transformer.transform(south,west)
    
    for x in np.arange(x_sw-step, x_ne+step, step):
        for y in np.arange(y_sw, y_ne, step):
            points.append((int(x/1000), int(y/1000)))
    filenames = []
    for x, y in points:
        x_str = str(int(x))
        y_str = str(int(y))
        filename = f"{prefix}{x_str}_{y_str}-2022/file/LIDARHD_1-0_LAZ_IP-{x_str}_{y_str}-2022.7z"
        filenames.append(filename)
        
    print("\n".join(filenames))

def list_laz_files(directory_path):

    laz_files = []
    
    # Walk through the directory
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.laz'):
                full_path = os.path.join(root, file)
                laz_files.append(full_path)
                
    return laz_files

import shapefile

def rectangles_intersect(A, B):
    A_x1, A_y1 = A['pbs']
    A_x2, A_y2 = A['phd']
    B_x1, B_y1 = B['pbs']
    B_x2, B_y2 = B['phd']

    return (A_x1 < B_x2 and B_x1 < A_x2 and
            A_y1 < B_y2 and B_y1 < A_y2)

    
def filter_shapefile_by_coordinates(west, south, east, north, filename,lidarCRS=2154, buffer = 500):
    points = []
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:%d"%lidarCRS)
    bs = transformer.transform(south,west)
    hd = transformer.transform( north,east)
    rect_A = {'pbs': (bs[0]-buffer,bs[1]-buffer), 'phd': (hd[0]+buffer,hd[1]+buffer)}
    print("goffu ",rect_A)
    sf = shapefile.Reader(filename)
    
    filtered_records = []
    for shape_record in sf.iterShapeRecords():
        rect_B = {'pbs': shape_record.shape.points[0], 'phd': shape_record.shape.points[2]}
        if "1220" in shape_record.record['nom_pkk']:
            print("goffu ",rect_A, rect_B, shape_record.record['nom_pkk'])

        #)
        if rectangles_intersect(rect_A, rect_B):
            print(rect_A, rect_B)
            filtered_records.append({
                'nom_pkk': shape_record.record['nom_pkk'],
                'url_telech': shape_record.record['url_telech']
            })

    # Filter records
    
    # Print filtered records
    for record in filtered_records:
        print(f"Nom PKK: {record['nom_pkk']}, URL: {record['url_telech']}")
 
    
def analyze_las_file(file_path):
    # Define classification descriptions
    classifications_description = {
        0: "Jamais classé",
        1: "Non attribuée",
        2: "Sol",
        3: "Végétation basse",
        4: "Moyenne végétation",
        5: "Haute végétation",
        6: "Bâtiment",
        7: "Point bas",
        8: "Réservé",
        9: "Eau",
        10: "Ferroviaire",
        11: "Surface routière",
        12: "Réservé",
        13: "Fil métallique (protection)",
        14: "Conducteur métallique (Phase)",
        15: "Tour de transmission",
        16: "Connecteur de structure métallique (Isolant)",
        17: "Tablier de pont",
        18: "Niveau sonore élevé"
    }

    # Open the LAS file
    inFile = laspy.read(file_path)
    
    # Total number of points
    total_points = inFile.header.point_count
    print("Total number of points:", total_points)
    
    # Get classifications and count occurrences
    classifications = inFile.classification
    unique, counts = np.unique(classifications, return_counts=True)
    
    # Display counts for each classification
    classification_counts = dict(zip(unique, counts))
    print("Counts by classification:")
    for classification, count in classification_counts.items():
        description = classifications_description.get(classification, "Unknown classification")
        print(f"Classification {classification} ({description}): {count} points")

    return total_points, classification_counts


def las_to_compressed_bin2(file_path, output_filename='output.bin'):
    # Open the LAS file
    inFile = laspy.read(file_path)
    
    # Total number of points
    total_points = inFile.header.point_count
    print("Total number of points:", total_points)
    
    # Get the minimum values for x, y, z to normalize coordinates
    min_x, min_y, min_z = np.min(inFile.x), np.min(inFile.y), np.min(inFile.z)
    
    # Prepare data for binary packing
    coords = np.vstack((inFile.x - min_x, inFile.y - min_y, inFile.z - min_z)).T
    scaleFactor = 255 / 65535
    colors = np.vstack((inFile.red * scaleFactor, inFile.green * scaleFactor, inFile.blue * scaleFactor)).T
    classifications = inFile.classification
    
    # Open a binary file to write the compressed data
    with open(output_filename, 'wb') as file:
        # Writing header with minimal coordinates and total points
        header = struct.pack('<3dI', min_x, min_y, min_z, total_points)
        file.write(header)
        
        # Prepare binary data
        binary_data = bytearray()
        for coord, color, classification in zip(coords, colors, classifications):
            # Pack data: 3 floats for coordinates, 3 uint8 for colors, 1 uint8 for classification
            packed_data = struct.pack('<3f3B1B', *coord, *np.round(color).astype(np.uint8), classification)
            binary_data.extend(packed_data)
        
        if output_filename != "None":    
            with zipfile.ZipFile(output_filename, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                # Create an in-memory file-like object for the binary data
                data_bin = io.BytesIO(binary_data)
                # Add this object to the zip, naming it 'data.bin'
                zip_file.writestr('data.bin', data_bin.getvalue())
      
        print(f"Data saved to {output_filename}")

def las_to_compressed_bin(file_path, output_filename='output.zip'):
    inFile = laspy.read(file_path)
    total_points = inFile.header.point_count
    min_x, min_y, min_z = np.min(inFile.x), np.min(inFile.y), np.min(inFile.z)
    
    coords = np.vstack((inFile.x - min_x, inFile.y - min_y, inFile.z - min_z)).T
    scaleFactor = 255 / 65535
    colors = np.vstack((inFile.red * scaleFactor, inFile.green * scaleFactor, inFile.blue * scaleFactor)).T
    classifications = inFile.classification
    
    binary_data = bytearray()
    binary_data.extend(struct.pack('<3dI', min_x, min_y, min_z, total_points))
    for coord, color, classification in zip(coords, colors, classifications):
        packed_data = struct.pack('<3f3B1B', *coord, *np.round(color).astype(np.uint8), classification.astype(np.uint8))
        binary_data.extend(packed_data)
    
    with zipfile.ZipFile(output_filename, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        zip_file.writestr('data.bin', binary_data)

    print(f"Data saved to {output_filename}")

# Function to download a file
def download_file(url,basepath=""):
    local_filename = basepath+url.split('/')[-1]
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            print("downloading ",local_filename)
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

def download_urs_set(urls,basepath):
    for url in urls:
        print(f"Downloading {url}")
        download_file(url,basepath)
        print(f"Finished downloading {url.split('/')[-1]}")


las1 = "/Users/filippi_j/data/2023/prunelli/LIDARHD_1-0_LAZ_VT-1224_6123-2021/Semis_2021_1224_6122_LA93_IGN78.laz"
las2 = "/Users/filippi_j/data/2023/prunelli/LIDARHD_1-0_LAZ_VT-1224_6123-2021/Semis_2021_1224_6123_LA93_IGN78.laz"
las3 = "/Users/filippi_j/data/2023/prunelli/LIDARHD_1-0_LAZ_VT-1224_6123-2021/Semis_2021_1225_6122_LA93_IGN78.laz"
las4 = "/Users/filippi_j/data/2023/prunelli/LIDARHD_1-0_LAZ_VT-1224_6123-2021/Semis_2021_1225_6123_LA93_IGN78.laz"
las5 = "/Users/filippi_j/data/2023/prunelli/LIDARHD_1-0_LAZ_VT-1226_6123-2021/Semis_2021_1226_6122_LA93_IGN78.laz"
las6 = "/Users/filippi_j/data/2023/prunelli/LIDARHD_1-0_LAZ_VT-1226_6123-2021/Semis_2021_1226_6123_LA93_IGN78.laz"
lidarHDIndex = "/Users/filippi_j/data/2023/lidargrid/TA_diff_pkk_lidarhd.shp"

# List of URLs to download you must select that first
prunellieasturls = [
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1224_6122_PTS_O_LAMB93_IGN69.copc.laz",
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1224_6123_PTS_O_LAMB93_IGN69.copc.laz",
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1225_6122_PTS_O_LAMB93_IGN69.copc.laz",
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1225_6123_PTS_O_LAMB93_IGN69.copc.laz"
]
pieraggiurls = [
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1226_6132_PTS_O_LAMB93_IGN69.copc.laz",
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1227_6132_PTS_O_LAMB93_IGN69.copc.laz",
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1227_6133_PTS_O_LAMB93_IGN69.copc.laz",
    "https://storage.sbg.cloud.ovh.net/v1/AUTH_63234f509d6048bca3c9fd7928720ca1/ppk-lidar/VT/LHD_FXX_1226_6133_PTS_O_LAMB93_IGN69.copc.laz"
]


#bbox de Yolanda
lat_sw, lon_sw = 42.0063067 , 9.3263082
lat_ne, lon_ne = 42.0104249 ,  9.3327945


#bbox réduite
#lat_sw, lon_sw = 42.007884 , 9.327366
#lat_ne, lon_ne = 42.008428 ,  9.328064

#bbox de corbara
#lat_sw, lon_sw = 42.6016, 8.9103
#lat_ne, lon_ne = 42.6049 ,8.9216

#bbox de monze
#lat_sw, lon_sw = 43.1289 , 2.3986
#lat_ne, lon_ne = 43.187127, 2.5797

#bbox de Prunelli smaller
#bbox de Prunelli tiny
lat_sw, lon_sw = 42.0075 , 9.3270
lat_ne, lon_ne = 42.0085 ,  9.3280

llPoint = [42.094471,   9.362104]
#bbox de Prunelli tiny
lat_sw, lon_sw = llPoint[0]-0.0005 , llPoint[1]-0.0005
lat_ne, lon_ne = llPoint[0]+0.0005 ,  llPoint[1]+0.0005


lat_sw, lon_sw =  42.091867 , 9.359683
lat_ne, lon_ne =  42.098067,    9.365923

dirout = "/Users/filippi_j/data/2024/prunelli/pieraggi/"
lidarDir = "/Users/filippi_j/data/2024/prunelli/pieraggi/lidarHD/"
dirout = "/Users/filippi_j/data/2024/prunelli/tmp/"

subsetFDS = "%ssubsetFDS.TIF"%dirout
subsetFDSLamb = "%ssubsetFDSLAMBERT.TIF"%dirout
lidplot = "%sslidarColor.png"%dirout
outLAS = "%sslidarAreaS.laz"%dirout
lidarVTKout = "%sslidarArea.vtk"%dirout
orthoTIF = "%sortho.tif"%dirout
orthoTIFLamb = "%sorthoLambert.tif"%dirout

compressedBin = "%spieraggiS.zip"%dirout

selectedLidarHDURLS = pieraggiurls

west, south, east, north = lon_sw,lat_sw, lon_ne,lat_ne

download_lidar_files = False
# downloading ultraHD zoom18 image of orthophoto
if download_lidar_files:
    download_urs_set(selectedLidarHDURLS,lidarDir)

get_orthophoto = False
# downloading ultraHD zoom18 image of orthophoto
if get_orthophoto:
    webMapsToTif(west, south, east, north, orthoTIF, zoomLevel=19)


# reprojection using lambert 2154 to match crs in corsica of lidar
project_orthophoto = False
if project_orthophoto:
    reprojectTif(orthoTIF,orthoTIFLamb,destCRS=2154)


# filter las points from a BBOX that corresponds to those found in the ortho tiff file and from files that are contained in a directory
# save it as vtk and las file
filter_las_points = True
if filter_las_points:
    left,right,bottom,top = getLeftRightBottomTop(orthoTIFLamb)  
    all_filtered_points = filter_las_files_with_attributes(  list_laz_files(lidarDir), left, bottom, right, top,remove_ratio=0.9 )

    all_filtered_points_colored = assign_colors_from_tif(all_filtered_points, orthoTIFLamb)
    save_filtered_las_color_classification_and_coords(all_filtered_points_colored, outLAS)
    save_points_to_vtk(all_filtered_points_colored, filename=lidarVTKout)


# gives some info
info_colored_las= True
if info_colored_las:
    analyze_las_file(outLAS)


# gives some info
save_compressed_bin = True
if save_compressed_bin:
    las_to_compressed_bin(outLAS,output_filename=compressedBin)



#plot_colored_points_3d(all_filtered_points_colored, filename=lidplot)
#print("plotted")
#all_filtered_points_latlon = convert_LASCRS_to_latlon(all_filtered_points_colored)
#save_filtered_las(all_filtered_points_latlon, outLAS+"latlon.las")


