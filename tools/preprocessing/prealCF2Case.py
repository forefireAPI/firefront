#!/usr/bin/python
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

import numpy as np
import netCDF4 as nc4
from PIL import Image
from .ffToGeoJson import get_WSEN_LBRT_ZS_From_Pgd
 
def addFieldToNcFile(ncfile, field, fieldname, typeName, dvartype):
        print("adding ", fieldname)
        sp = np.shape(field)
        
        
        if (len(sp) == 2):
            ncfile.createDimension('%sNX'%fieldname, sp[1])
            ncfile.createDimension('%sNY'%fieldname, sp[0])
            ncfile.createDimension('%sNZ'%fieldname, 1)
            ncfile.createDimension('%sNT'%fieldname, 1)
        if (len(sp) == 3):
            ncfile.createDimension('%sNX'%fieldname, sp[2])
            ncfile.createDimension('%sNY'%fieldname, sp[1])
            ncfile.createDimension('%sNZ'%fieldname, sp[0])
            ncfile.createDimension('%sNT'%fieldname, 1)
        if (len(sp) == 4):
            ncfile.createDimension('%sNX'%fieldname, sp[3])
            ncfile.createDimension('%sNY'%fieldname, sp[2])
            ncfile.createDimension('%sNZ'%fieldname, sp[1])
            ncfile.createDimension('%sNT'%fieldname, sp[0])
        
        
        
        # ncfile.createDimension('%sNX'%fieldname, sp[1])
        # ncfile.createDimension('%sNY'%fieldname, sp[0])
        # if (len(sp) > 2):
        #     ncfile.createDimension('%sNZ'%fieldname, sp[2])
        # else:
        #     ncfile.createDimension('%sNZ'%fieldname, 1)
        #     ncfile.createDimension('%sNT'%fieldname, 1)
        # if (len(sp) > 3):
        #     ncfile.createDimension('%sNT'%fieldname, sp[3])
        
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
 
        ncfile =  nc4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
        ncfile.version = "FF.1.0"
        domain = ncfile.createVariable('domain', 'S1', ())
        domain.type = "domain"
        for key, value in domainProperties.items():
            setattr(domain, key, value)

        parameters = ncfile.createVariable('parameters', 'S1', ())
        parameters.type = "parameters"       

        if (parametersProperties is not None):
            for key, value in parametersProperties.items():
                setattr(parameters, key, value)

 
        if (fuelModelMap is not None):
            addFieldToNcFile(ncfile, fuelModelMap, 'fuel', 'fuel', 'i2')
            
        if (elevation is not None):
            addFieldToNcFile(ncfile, elevation, 'altitude', 'data', 'f4')
        
        if (wind is not None):
            addFieldToNcFile(ncfile, wind, 'wind', 'data', 'f4')
            # 2 types de vents, 
            # soit windU/windV statique (2 arrays, windU et WindV de forme:  1,1,NI,NJ)
            # soit variable dans le temps, dynamique 2 arrays, windU et WindV de forme NT, 1, NI,NJ)
            # soit champ potentiel nom wind de forme NAXES,NDIRCOEFFS,NI,NJ 
#            addFieldToNcFile(ncfile, wind["zonal"], 'windU', 'data', 'f4')
 #           addFieldToNcFile(ncfile, wind["meridian"], 'windV', 'data', 'f4')
                
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
                fVar = addFieldToNcFile(ncfile, fMap["data"], fMap["name"], 'flux', 'i1')
                #ncfile.createVariable(fMap["name"], 'i4', ('DIMT', 'DIMZ', 'DIMY', 'DIMX'))
                #fVar.type = "flux" ;
                for entry in fMap["table"].keys():
                    setattr(fVar, "model%dname"%fMap["table"][entry], entry)
                fVar.indices = np.array(list(fMap["table"].values()),dtype=('i1'))
                #fVar[0,0,:,:] = fMap["data"]

        
        print("writing ", filename)
        ncfile.sync()
        ncfile.close()

def PGD2Case(pgd_path, png_path, out_path, dateStartDom, fuel_test=None, gen_wind=None):
    WSEN, LBRT, ZS = get_WSEN_LBRT_ZS_From_Pgd(pgd_path)
    image2Case(WSEN, LBRT, ZS, png_path, out_path, dateStartDom, fuel_test=fuel_test, gen_wind=gen_wind)

    
def image2Case(WSEN, LBRT, ZS, png_path, out_path, dateStartDom, fuel_test=None, gen_wind=None):
    
 
    fout = out_path 
    imgPath = png_path
    W,S,E,N=WSEN
    L, B, R, T = LBRT
    out_size = (int((LBRT[2]-LBRT[0])),int((LBRT[3]-LBRT[1])))
    wind=None
     
 
    
    domainProperties= {}
    domainProperties['SWx']  = np.float32(L)
    domainProperties['SWy']  = np.float32(B)
    domainProperties['SWz']  = np.float32(0)
    domainProperties['Lx']   = np.float32(out_size[0])
    domainProperties['Ly']   = np.float32(out_size[1])
    domainProperties['Lz']   = np.float32(0)
    domainProperties['t0']   = np.float32(0)
    domainProperties['Lt']   = np.Inf
    domainProperties['BBoxWSEN']  = f"{W},{S},{E},{N}"
    
    
  #  secondsSinceMidnight = (dateStartDom - dateStartDom.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    
    parametersProperties= {}
    parametersProperties['date']  = dateStartDom.strftime("%Y-%m-%d_%H:%M:%S")
    parametersProperties['duration']  = np.int32(3600*24)
    parametersProperties['refYear']  =  np.int32(dateStartDom.year) 
    parametersProperties['refDay']  =  np.int32(dateStartDom.timetuple().tm_yday) 
    parametersProperties['year']  =  np.int32(dateStartDom.year) 
    parametersProperties['month']  =  np.int32(dateStartDom.month) 
    parametersProperties['day']  =  np.int32(dateStartDom.day) 
    
   # dom= "FireDomain[sw=(%d,%f,0);ne=(%f,%f,0);t=%d]"%(domainProperties['SWx'],domainProperties['SWy'],domainProperties['SWx']+domainProperties['Lx'],domainProperties['SWy']+domainProperties['Ly'],secondsSinceMidnight)
   # print( dom )
    
    elevation =  ZS
   
    fuelMap   =None
    if imgPath is not None:
        img = Image.open(imgPath)
        if img.mode != "P":
            raise ValueError("L'image n'est pas en mode palette (P).")
            
        fuelMap = np.flipud(np.asarray(img).astype('i2'))
        #fuelMap = np.zeros(np.shape(fuelData),dtype=('i4'))
        #fuelMap[fuelData == 0] = 1
        #fuelMap[fuelData == 1] = 0
        #fuelMap[fuelData == 2] = 1
        #fuelMap[fuelData == 3] = 0
    else:
        fuelMap   = 1*np.ones(np.shape(elevation),dtype=('u1'))
    
    if fuel_test is not None :
        fuelData = np.flipud(np.asarray(Image.open(imgPath)))
        fuelMap = np.zeros(np.shape(fuelData),dtype=('i1'))
        n_bandes = len(fuel_test)
        limites = np.linspace(0, fuelMap.shape[1], n_bandes + 1).astype(int)
        elevation = elevation*0
        # Remplissage de fuelMap
        for i in range(n_bandes):
            fuelMap[:, limites[i]:limites[i + 1]] = fuel_test[i]
            
    if gen_wind is not None:
        #elevation = elevation*0
        ni=3
        nj=3
        wind=np.full((2, 4, ni, nj), 0)
        #wind params examples : 
        #Southerly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(0.0,10.0,0.);t=0.]
        #Westerly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(10.0,0.0,0.);t=0.]
        #Northerly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(0.0,-10.0,0.);t=0.]
        #Easterly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(-10.0,0.0,0.);t=0.]
        N = 8  # Nombre de directions
        angle_step = 2 * np.pi / N
        
        wind = np.zeros((2, N, ni, nj))
        
        for i in range(N):
            angle = i * angle_step
            u = 10 * np.cos(angle)  # Composante U
            v = 10 * np.sin(angle)  # Composante V
            wind[0][i] = np.full((ni, nj), u)
            wind[1][i] = np.full((ni, nj), v)
        # #for U
        # #stherly
        # wind[0][0] = np.full((ni, nj), 0)
        # #westerly
        # wind[0][1] = np.full((ni, nj), 10)
        # #northrly
        # wind[0][2] = np.full((ni, nj), 0)
        # #eastrly
        # wind[0][3] = np.full((ni, nj), -10)
        
        # #for V
        # wind[1][0] = np.full((ni, nj), 10)
        # #westerly
        # wind[1][1] = np.full((ni, nj), 0)
        # #northrly
        # wind[1][2] = np.full((ni, nj), -10)
        # #eastrly
        # wind[1][3] = np.full((ni, nj), 0)   
        
    
    heatFluxModelMap = {}
    heatFluxModelMap["name"] =  "heatFlux"
    heatFluxModelMap["data"] =  np.zeros(np.shape(elevation),dtype=('u1'))
    heatFluxModelMap["table"] =  {"heatFluxBasic":0,}
    
    vaporFluxModelMap = {}
    vaporFluxModelMap["name"] =  "vaporFlux"
    vaporFluxModelMap["data"] =  np.ones(np.shape(elevation),dtype=('u1'))
    vaporFluxModelMap["table"] =  {"vaporFluxBasic":1,}
    
     
    
    FiretoNC(fout, domainProperties,parametersProperties,fuelMap,elevation=elevation, wind=wind, fluxModelMap = (heatFluxModelMap,vaporFluxModelMap))
