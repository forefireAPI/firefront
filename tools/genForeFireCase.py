# Copyright (C) 2012
# Author(s): Jean Baptiste Filippi, Vivien  Mallet
#
# This file is part of pyFireScore, a tool for scoring wildfire simulation
#
# pyFireScore is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# pyFireScore is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

# USAGE :
# This script is to be used in order to generate a packed data landscape data in cdf format for ForeFire.
# Usage: FiretoNC(filename, file name
#         domainProperties, the domain extension (map matching forefire parameters SWx, SWy, SWz, Lx, Ly, Lz, t0, Lt)
#     parametersProperties, the other optional properties you may want to put in the list
#             fuelModelMap, a numpy integer array containing the indexes of fuel type
#           elevation=None, a numby real array with the elevation
#                wind=None, a map with a ["zonal"] and ["meridian"] numpy real array values
#       fluxModelMap=None): a map with a (["table"] and ["name"] ) and fMap["data"] numpy int array values containing indices to the corresponding flux model
# 

import numpy as np
from scipy.io import netcdf
import time


def FiretoNC(filename, domainProperties, parametersProperties, fuelModelMap, elevation=None, wind=None, fluxModelMap =None):
 
        ncfile =  netcdf.netcdf_file(filename, 'w')   
        ncfile.createDimension('DIMX', fuelModelMap.shape[1])
        ncfile.createDimension('DIMY', fuelModelMap.shape[0])
        ncfile.createDimension('DIMZ', 1)
        ncfile.createDimension('DIMT', 1)
        ncfile.version = "FF.1.0"
        domain = ncfile.createVariable('domain', 'S1', ())
        domain.type = "domain" 
        domain.SWx = domainProperties['SWx'] 
        domain.SWy = domainProperties['SWy'] 
        domain.SWz = domainProperties['SWz']  
        domain.Lx =  domainProperties['Lx']  
        domain.Ly =  domainProperties['Ly']  
        domain.Lz =  domainProperties['Lz']  
        domain.t0 =  domainProperties['t0']  
        domain.Lt =  domainProperties['Lt'] 
        
        parameters = ncfile.createVariable('parameters', 'S1', ())
        parameters.type = "parameters"       

        parameters.projectionproperties = parametersProperties['projectionproperties'] 
        parameters.date = parametersProperties['date'] 
        parameters.duration = parametersProperties['duration'] 
        parameters.projection = parametersProperties['projection']
        parameters.refYear = parametersProperties['refYear'] 
        parameters.refDay = parametersProperties['refDay']
        
        
        fuel = ncfile.createVariable('fuel', 'i4', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
        fuel[0,0,:,:] = fuelModelMap 
        fuel.type = "fuel" ;
        
        if (elevation != None):
            altitude = ncfile.createVariable('altitude', 'f8', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
            altitude[0,0,:,:] = elevation
            altitude.type = "data" 
        
        if (wind != None):
            windU = ncfile.createVariable('windU', 'f8', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
            windU.type = "data" 
            windU[0,0,:,:] = wind["zonal"]
            windV = ncfile.createVariable('windV', 'f8', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
            windV.type = "data" 
            windV[0,0,:,:] = wind["meridian"]
            
        if (fluxModelMap != None):
            numOfModels = 0
            for fMap in fluxModelMap:
                numOfModels += len(fMap["table"])
            
            for fMap in fluxModelMap:  
                fVar = ncfile.createVariable(fMap["name"], 'i4', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
                fVar.type = "flux" ;
                for entry in fMap["table"].keys():
                    setattr(fVar, "model%dname"%fMap["table"][entry], entry)
                fVar.indices = np.array(fMap["table"].values(),dtype=('i4'))
                fVar[0,0,:,:] = fMap["data"]
        print "writing ", filename
        ncfile.sync()
        ncfile.close()






