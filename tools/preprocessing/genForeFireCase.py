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

def addFieldToNcFile(ncfile, field, fieldname, typeName, dvartype):
        sp = np.shape(field)
        ncfile.createDimension('%sNX'%fieldname, sp[1])
        ncfile.createDimension('%sNY'%fieldname, sp[0])
        if (len(sp) > 2):
            ncfile.createDimension('%sNZ'%fieldname, sp[2])
        else:
            ncfile.createDimension('%sNZ'%fieldname, 1)
            ncfile.createDimension('%sNT'%fieldname, 1)
        if (len(sp) > 3):
            ncfile.createDimension('%sNT'%fieldname, sp[3])
        
        variable = ncfile.createVariable(fieldname, dvartype, ('%sNT'%fieldname, '%sNZ'%fieldname,'%sNY'%fieldname, '%sNX'%fieldname))
        if (len(sp) == 4):
            variable[:,:,:,:] = field 
        if (len(sp) == 3):
            variable[0,:,:,:] = field 
        if (len(sp) == 2):
            variable[0,0,:,:] = field 
            
        variable.type = typeName;
        
        return variable
    

def FiretoNC(filename, domainProperties, parametersProperties, fuelModelMap, elevation=None, wind=None, fluxModelMap =None, bmap=None, cellMap=None):
 
        ncfile =  netcdf.netcdf_file(filename, 'w')   
     #   ncfile.createDimension('DIMX', elevation.shape[1])
      #  ncfile.createDimension('DIMY', elevation.shape[0])
      #  ncfile.createDimension('DIMZ', 1)
      #  ncfile.createDimension('DIMT', 1)
        ncfile.version = "FF.1.0"
        domain = ncfile.createVariable('domain', 'S1', ())
        domain.type = "domain" 
        domain.SWx = float(domainProperties['SWx'])
        domain.SWy = float(domainProperties['SWy'] )
        domain.SWz = float(domainProperties['SWz']  )
        domain.Lx  = float(domainProperties['Lx']  )
        domain.Ly  = float(domainProperties['Ly']  )
        domain.Lz  = float(domainProperties['Lz']  )
        domain.t0  = float(domainProperties['t0']  )
        domain.Lt  = float(domainProperties['Lt'] )

        parameters = ncfile.createVariable('parameters', 'S1', ())
        parameters.type = "parameters"       

        if (parametersProperties is not None):
            parameters.date = parametersProperties['date'] 
            parameters.duration = parametersProperties['duration'] 
            parameters.refYear = parametersProperties['refYear'] 
            parameters.refDay = parametersProperties['refDay']
            parameters.year = parametersProperties['year']  
            parameters.month = parametersProperties['month']  
            parameters.day = parametersProperties['day'] 
        
    
        
        if (fuelModelMap is not None):
            addFieldToNcFile(ncfile, fuelModelMap, 'fuel', 'fuel', 'i4')
            
        if (elevation is not None):
            addFieldToNcFile(ncfile, elevation, 'altitude', 'data', 'f8')
        
        if (wind is not None):
            addFieldToNcFile(ncfile, wind["zonal"], 'windU', 'data', 'f8')
            addFieldToNcFile(ncfile, wind["meridian"], 'windV', 'data', 'f8')
                
      #      windU = ncfile.createVariable('windU', 'f8', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
       #     windU.type = "data" 
       #     windU[0,0,:,:] = wind["zonal"]
       #     windV = ncfile.createVariable('windV', 'f8', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
       #     windV.type = "data" 
       #     windV[0,0,:,:] = wind["meridian"]
            
        if (fluxModelMap is not None):
            numOfModels = 0
            for fMap in fluxModelMap:
                numOfModels += len(fMap["table"])
            
            for fMap in fluxModelMap:  
                fVar = addFieldToNcFile(ncfile, fMap["data"], fMap["name"], 'flux', 'i4')
                #ncfile.createVariable(fMap["name"], 'i4', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
                #fVar.type = "flux" ;
                for entry in fMap["table"].keys():
                    setattr(fVar, "model%dname"%fMap["table"][entry], entry)
                fVar.indices = np.array(list(fMap["table"].values()),dtype=('i4'))
                #fVar[0,0,:,:] = fMap["data"]

        
        print("writing ", filename)
        ncfile.sync()
        ncfile.close()






