#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 19:25:06 2024

@author: filippi_j
"""
# A file to upgrade DIR/HDR files in order to stop having Halo problems in MesoNH

import numpy as np
import matplotlib.pyplot as plt

def load_and_plot_data(file_path, rows, cols):
    # Charger les données du fichier
    data = np.fromfile(file_path, dtype=np.int16)  # Assurez-vous que le dtype correspond au type de données
    data = data.reshape((rows, cols))  # Remodeler les données selon les dimensions fournies

    # Gérer les valeurs 'nodata'
    data = np.where(data == -9999, np.nan, data)  # Remplacer -9999 par NaN

    # Afficher les données
    plt.imshow(data, cmap='terrain')  # Choisir une colormap appropriée
    plt.colorbar(label='Altitude')
    plt.title('Visualisation du Fichier .dir')
    plt.show()

# Utilisation de la fonction


import numpy as np
import xarray as xr
import scipy.ndimage

def geotiff_to_xarray(filename):
    import rasterio
    """ Load a GeoTIFF file as an xarray.DataArray, assuming WGS84 coordinates """
    with rasterio.open(filename) as src:
        data = src.read(1)  # Read the first band; adjust if necessary for multiple bands
        # Get the bounds from the file
        left, bottom, right, top = src.bounds.left, src.bounds.bottom, src.bounds.right, src.bounds.top
        # Calculate the resolution of the raster
        res_lat = (top - bottom) / src.height
        res_lon = (right - left) / src.width
        
        # Generate coordinates based on the bounds and resolution
        lats = np.linspace(bottom + res_lat / 2, top - res_lat / 2, src.height)
        lons = np.linspace(left + res_lon / 2, right - res_lon / 2, src.width)
        
        # Create the xarray DataArray
        da = xr.DataArray(
            np.flipud(data),
            dims=['latitude', 'longitude'],
            coords={
                'latitude': lats,
                'longitude': lons
            },
            attrs={
                'description': 'Raster data from {}'.format(filename),
                'units': 'unknown',  # Modify according to the actual units
                'nodata': src.nodatavals[0]
            }
        )
    return da

def read_hdr(filenamehdr):
    """ Lire le fichier .hdr et extraire les métadonnées """
    metadata = {}
    with open(filenamehdr, 'r') as file:
        for line in file:
            parts = line.strip().split(': ')
            if len(parts) == 2:
                key, value = parts
                metadata[key] = value
    return metadata

def dirHDRtoXarray(filenamedir, filenamehdr):
    """ Charger un fichier .dir en tant que xarray.DataArray """
    metadata = read_hdr(filenamehdr)
    print(metadata)
    rows = int(metadata['rows'])
    cols = int(metadata['cols'])
    data_type = np.int16  # Modifiez cela en fonction du type de données spécifié dans le .hdr

    data = np.fromfile(filenamedir, dtype=data_type).reshape((rows, cols))
    data=np.flipud(data)
    # Créer des coordonnées (exemple simple, ajustez selon vos besoins réels)
    lats = np.linspace(float(metadata['south']), float(metadata['north']), rows)
    lons = np.linspace(float(metadata['west']), float(metadata['east']), cols)

    return xr.DataArray(data, coords=[lats, lons], dims=['latitude', 'longitude'])

def xarrayToDirHdr(xrarray, filenamedir, filenamehdr):
    """ Enregistrer un xarray.DataArray dans un fichier .dir en tant que données 16 bits et créer un fichier .hdr """
    # Convertir les données en int16
    nparray = np.flipud(xrarray.values.astype(np.int16))

    # Enregistrer les données dans un fichier .dir
    nparray.tofile(filenamedir)

    # Créer les métadonnées pour le fichier .hdr
    metadata = {
        'nodata': '-9999',
        'north': str(xrarray.latitude.max().item()),
        'south': str(xrarray.latitude.min().item()),
        'east': str(xrarray.longitude.max().item()),
        'west': str(xrarray.longitude.min().item()),
        'rows': str(xrarray.shape[0]),
        'cols': str(xrarray.shape[1]),
        'recordtype': 'integer 16 bit'
    }

    with open(filenamehdr, 'w') as file:
        file.write("PROCESSED SRTM DATA VERSION 4.1, orography model\n")
        for key, value in metadata.items():
            file.write(f"{key}: {value}\n")
            


def makeSubset():
    dirin = '/Users/filippi_j/soft/Meso-NH/PGD/srtm_ne_250.dir'
    hdrin = '/Users/filippi_j/soft/Meso-NH/PGD/srtm_ne_250.hdr'
    dirout = '/Users/filippi_j/soft/Meso-NH/PGD/ffDEM.dir'
    hdrout = '/Users/filippi_j/soft/Meso-NH/PGD/ffDEM.hdr'

    # Charger les données
    data_xarray = dirHDRtoXarray(dirin, hdrin)
    print("data loaded")
    # Sélectionner le sous-ensemble
    subset = data_xarray.sel(latitude=slice(30, 60.02), longitude=slice(-15, 35))
    #subset.plot()
    # Augmenter la résolution
    # Le facteur 4 signifie que chaque dimension sera 4 fois plus grande
    zoom_factor = 4
    #resampled_data = scipy.ndimage.zoom(subset, (zoom_factor, zoom_factor), order=3)  # order=3 pour une interpolation cubique
    
    # Convert subset to float
    
    # Create a mask for the -9999 values
    mask = subset == -9999
    
    # Replace -9999 values with NaN using xarray's where method
    data_with_0 = subset.where(~mask, 0)
    
    
    
    resampled_data = scipy.ndimage.zoom(data_with_0, (zoom_factor, zoom_factor), order=3)
    print("resampled")
    resampled_mask = scipy.ndimage.zoom(mask, (zoom_factor, zoom_factor), order=0)
    # Reapply the mask to the resampled data
    resampled_data_masked = np.where(resampled_mask, -9999, resampled_data)
    print("masked")
    
    
    # Créer un nouvel xarray.DataArray avec les nouvelles coordonnées
    new_lats = np.linspace(subset.latitude.min(), subset.latitude.max(), resampled_data_masked.shape[0])
    new_lons = np.linspace(subset.longitude.min(), subset.longitude.max(), resampled_data_masked.shape[1])
    resampled_xarray = xr.DataArray(resampled_data_masked, coords=[new_lats, new_lons], dims=['latitude', 'longitude'])
    
    # Sauvegarder les données
    xarrayToDirHdr(resampled_xarray, dirout, hdrout)

