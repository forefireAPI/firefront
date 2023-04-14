from rasterio.crs import CRS
from rasterio.enums import Resampling
# from rasterio.transform import Affine
from rasterio.vrt import WarpedVRT
import rasterio as rio
from math import floor
import numpy as np
from datetime import *
import netCDF4 as netcdf
from pyproj import Transformer
import affine

# functions may need to be adapted depending on the specific use case

def default_wind_generator(elevation_array):
    '''generates 8 bands in 8 directions of default 10 m/s wind'''
    wind_dict = {}
    wind_dict['wind_u'] = []
    wind_dict['wind_v'] = []
    wind_dict['wind_shape'] = (floor(elevation_array.shape[0]/10), floor(elevation_array.shape[1]/10))
    theta = np.linspace(0, 315, 8)
    wind_speed = 10

    for i in range(8):
        wind_direction = np.radians(theta[i])

        a = -wind_speed * np.cos(wind_direction)
        b = -wind_speed * np.sin(wind_direction)
        
        u = np.ones(wind_dict['wind_shape'])*a
        v = np.ones(wind_dict['wind_shape'])*b
        
        wind_dict['wind_u'].append(u)
        wind_dict['wind_v'].append(v)

def prop_vrt_Warp(src_path, epsg):         
    src = rio.open(src_path)
    bbox = src.bounds

    min_long, min_lat = (bbox[0], bbox[1])
    max_long, max_lat = (bbox[2], bbox[3])

    transformer = Transformer.from_crs("epsg:4326", f'epsg:{epsg}')
    min_x, min_y = transformer.transform( min_long, min_lat)
    max_x, max_y = transformer.transform( max_long, max_lat)
    
    right = max_x
    bottom = min_y
    left = min_x
    top = max_y
    
    dst_crs = CRS.from_epsg(epsg)
    resolution = 30

    dst_width = floor((right - left)/resolution)
    dst_height = floor((top - bottom)/resolution)

    xres = (right - left) / dst_width
    yres =  (top - bottom) / dst_height

    dst_transform = affine.Affine(xres, 0.0, left,
                                  0.0, -yres, top)
    
    vrt_options = {
        'resampling': Resampling.nearest,
        'crs': dst_crs,
        'transform': dst_transform,
        'height': dst_height,
        'width': dst_width,
    }

    with rio.open(src_path) as src:

        with WarpedVRT(src, **vrt_options) as vrt:
            
            data = vrt.read(1)

    return vrt, data

def fuel_model_map_generator(fuel_filepath, epsg):
    '''fuel_filepath is path to TIF fuel file where fuel classes correspond to the fuels.ff (sample or custom)'''
    fuel_ds = prop_vrt_Warp(fuel_filepath, epsg)
    fuel_model_map = fuel_ds[1]
    return fuel_model_map

def elevation_generator(elevation_filepath, epsg):
    '''elevation_filepath is path to TIF DEM file'''
    elevation_ds = prop_vrt_Warp(elevation_filepath, epsg)
    elevation_map = elevation_ds[1]
    return elevation_map

def domainGenerator(fuelModelMap):
    
    meta = fuelModelMap[0].meta['transform']
    res_y = meta[0]
    res_x = meta[4]

    shp_y = fuelModelMap[0].shape[0]
    shp_x = fuelModelMap[0].shape[1]
    
    x = meta[2] #+ res_x*shp_x
    y = meta[5] - res_y*shp_y
    
    Lx = shp_x*res_x
    Ly = shp_y*res_y
   
    domainProperties= {}
    domainProperties['SWx']  = x
    domainProperties['SWy']  = y
    domainProperties['SWz']  = 0
    domainProperties['Lx']   = -Lx
    domainProperties['Ly']   = Ly
    domainProperties['Lz']   = 0
    domainProperties['t0']   = 0
    domainProperties['Lt']   = 0
    
    print(domainProperties)
    return domainProperties

def parameter_generator(projection):
    today = datetime.now()
    year = today.year
    day = today.day
    iso_today = today.isoformat(sep ='_', timespec="seconds")
    
    parameters_properties= {}
    parameters_properties['date']  = f"{iso_today}"
    parameters_properties['duration']  = 3600
    parameters_properties['projection']  = f"{projection}"
    parameters_properties['refYear']  = year
    parameters_properties['refDay']  = day
    return parameters_properties

def landscape_generator(filename, domain_properties, parameters_properties, projection, fuel_model_map, wind_dict, elevation=None):
    ncfile =  netcdf.netcdf_file(filename, 'w')
    ncfile.createDimension('wind_dimensions', 2)
    ncfile.createDimension('wind_directions', 8)
    ncfile.createDimension('wind_rows', wind_dict['wind_shape'][0])
    ncfile.createDimension('wind_columns', wind_dict['wind_shape'][1])
    
    #FUEL  
    ncfile.createDimension('ft', 1)
    ncfile.createDimension('fz', 1)
    ncfile.createDimension('fy', fuel_model_map.shape[0])
    ncfile.createDimension('fx', fuel_model_map.shape[1])

    #ELEVATION
    ncfile.createDimension('nt', 1)
    ncfile.createDimension('nz', 1)
    ncfile.createDimension('ny', elevation.shape[0])
    ncfile.createDimension('nx', elevation.shape[1])
    
    ncfile.projection = projection
    
    domain = ncfile.createVariable('domain', 'S1', ())
    domain.type = "domain" 
    domain.SWx = domain_properties['SWx'] 
    domain.SWy = domain_properties['SWy'] 
    domain.SWz = domain_properties['SWz']  
    domain.Lx =  domain_properties['Lx']  
    domain.Ly =  domain_properties['Ly']  
    domain.Lz =  domain_properties['Lz']  
    domain.t0 =  domain_properties['t0']  
    domain.Lt =  domain_properties['Lt'] 
    parameters = ncfile.createVariable('parameters', 'S1', ())
    parameters.type = "parameters"       

    if (parameters_properties is not None):
        parameters.projection = parameters_properties['projection']        
    
    fuel = ncfile.createVariable('fuel', 'i4', ('ft', 'fz', 'fy', 'fx'))
    fuel[0,0,:,:] = np.flip(fuel_model_map, axis=0)
    fuel.type = "fuel"
    
    if (elevation is not None):
        altitude = ncfile.createVariable('altitude', 'f8', ('nt', 'nz', 'ny', 'nx'))
        altitude[0,0,:,:] = np.flip(elevation, axis=0)
        altitude.type = "data" 
    
    wind = ncfile.createVariable('wind', 'f8', ('wind_dimensions', 'wind_directions', 'wind_rows', 'wind_columns'))
    wind.type = "data"          
        
    for i in range(8):           
        wind[0,i,:,:] = np.flip(wind_dict['wind_u'][i], axis=0)
        wind[1,i,:,:] = np.flip(wind_dict['wind_v'][i], axis=0)
    
    print("writing ", filename)
    ncfile.sync()
    ncfile.close()
    return


### WRAPPING EVERYTHING IN AN EXAMPLE FUNCTION TO GENERATE LANDSCAPE

def example_landscape_gen():
    fuel_filepath = 'path/to/fuel.tif'
    elevation_filepath = 'path/to/elevation.tif'
    epsg = 5880

    fuel_model_map = fuel_model_map_generator(fuel_filepath, epsg)
    elevation_map = elevation_generator(elevation_filepath, epsg)
    wind_dict = default_wind_generator(elevation_map)

    src = prop_vrt_Warp(fuel_filepath, epsg)
    domain = domainGenerator(src)
    parameters = parameter_generator(epsg)
    
    landscape_generator('custom_landscape.nc', domain, parameters, epsg, fuel_model_map, wind_dict, elevation_map)
