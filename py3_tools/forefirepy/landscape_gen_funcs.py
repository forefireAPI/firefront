import rasterio as rio
from math import floor
import numpy as np
from datetime import *
import netCDF4 as netcdf

def fuel_model_map_generator(fuel_filepath):
    '''fuel_filepath is path to TIF fuel file where fuel classes correspond to the fuels.ff (sample or custom)'''
    fuel_ds = rio.open_rasterio(fuel_filepath)
    fuel_model_map = fuel_ds[0]
    return fuel_model_map

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

def domain_generator(fuel_model_map):
    meta = fuel_model_map[1]
    res_y = meta[0]
    res_x = meta[4]

    shp_y = fuel_model_map[0].shape[0]
    shp_x = fuel_model_map[0].shape[1]
    
    x = meta[2] #+ res_x*shp_x
    y = meta[5] - res_y*shp_y
    
    Lx = shp_x*res_x
    Ly = shp_y*res_y
   
    domain_properties= {}
    domain_properties['SWx']  = x
    domain_properties['SWy']  = y
    domain_properties['SWz']  = 0
    domain_properties['Lx']   = -Lx
    domain_properties['Ly']   = Ly
    domain_properties['Lz']   = 0
    domain_properties['t0']   = 0
    domain_properties['Lt']   = 0
    # print(domain_properties)
    return domain_properties

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
    
    ncfile.projection = projection # e.g. 5880
    
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