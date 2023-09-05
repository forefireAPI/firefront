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
import sys
import xarray as xr
import netCDF4 as nc4
import datetime
import matplotlib.pyplot as plt
from PIL import Image



def forceDim(arrIn, newD=None, kind='edge'):
    if(newD==None):
        return arrIn
   
    if (newD[0][1] < 0):
        newr = arrIn[newD[0][0]:newD[0][1],newD[1][0]:newD[1][1]]
        print("new shape ",np.shape(newr))
        return(newr)
    

    return np.pad(arrIn, newD, kind)

 
exec(compile(open("genForeFireCase.py", "rb").read(), "genForeFireCase.py", 'exec'))


fname = None
fout = None 
imgPath = None

#### configuration Pigna
fname = "/Users/filippi_j/Volumes/orsu/firecaster/2023/nest150Ref/001_pgd/PGD_D80mA.nested.nc"
fout = "/Users/filippi_j/data/2023/corbara20230727/mnhCase/lpsimplef.nc"
imgPath = "/Users/filippi_j/data/2023/corbara20230727/pignaFUEL.png"
dateString = "202307252000"


if  (len(sys.argv)<2):
   print("Usage prealCF2Case inMNH_PGD_File.nc outFFDCaseFile.nc YYYYMMDDHHMM opt:pngFuelFile")
   if(fname == None) :
    sys.exit(0)
else :
    fname = sys.argv[1]
    fout = sys.argv[2]
    dateString = sys.argv[4]


if(len(sys.argv)>4):
    imgPath=sys.argv[4]

 
 
dateStartDom = datetime.datetime.strptime(dateString, "%Y%m%d%H%M")

 
print("reading %s generationg %s with %s as fuel"%(fname,fout,imgPath))
#.png"


#geometry 
geomCDFBin = fname 
 
nc = nc4.Dataset(geomCDFBin, 'r') 

# recuperer l'altitude
# decider d'une résolution de combustible, faire carte de combustible
# pas mettre de vent
# faire les maps de flux et de modèles de propa en se basant sur 

DeltaY = nc.variables['YHAT'][1]-nc.variables['YHAT'][0]
DeltaX = nc.variables['XHAT'][1]-nc.variables['XHAT'][0]

domainProperties= {}
domainProperties['SWx']  = nc.variables['XHAT'][0]
domainProperties['SWy']  = nc.variables['YHAT'][0]
domainProperties['SWz']  = 0
domainProperties['Lx']   = nc.variables['XHAT'][-1]+DeltaX-nc.variables['XHAT'][0]
domainProperties['Ly']   = nc.variables['YHAT'][-1]+DeltaY-nc.variables['YHAT'][0]
domainProperties['Lz']   = 0
domainProperties['t0']   = 0
domainProperties['Lt']   = np.Inf

secondsSinceMidnight = (dateStartDom - dateStartDom.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()

parametersProperties= {}
parametersProperties['date']  = dateStartDom.strftime("%Y-%m-%d_%H:%M:%S");
parametersProperties['duration']  = 3600*24;
parametersProperties['refYear']  = dateStartDom.year ;
parametersProperties['refDay']  = dateStartDom.timetuple().tm_yday ;
parametersProperties['year']  = dateStartDom.year ;
parametersProperties['month']  = dateStartDom.month ;
parametersProperties['day']  = dateStartDom.day ;

dom= "FireDomain[sw=(%d,%f,0);ne=(%f,%f,0);t=%d]"%(domainProperties['SWx'],domainProperties['SWy'],domainProperties['SWx']+domainProperties['Lx'],domainProperties['SWy']+domainProperties['Ly'],secondsSinceMidnight)
print( dom, parametersProperties['date']  )
print ("   ori ", nc.variables['LAT0'].getValue(), nc.variables['LON0'].getValue())
  
SVM = nc.variables['ZS'][:]
print(np.shape(SVM))
 
NJ = len(SVM)-2
NI = len(SVM[0])-2
 
NJT = len(SVM) 
NIT = len(SVM[0]) 
 
padDim = ((4,4),(4,4))

DeltaY = nc.variables['YHAT'][1]-nc.variables['YHAT'][0]
DeltaX = nc.variables['XHAT'][1]-nc.variables['XHAT'][0]

dt = 0
lastTime = 0

lastTime = 0#nc.variables['DTCUR'].getValue()
firstTime = 0#nc.variables['DTCUR'].getValue()

elevation =  nc.variables['ZS'][:,:]



# Tout le combustible a 1
fuelMap   =None
if imgPath is not None:
    fuelMap = np.flipud(np.asarray(Image.open(imgPath))+10)
    fuelMap[fuelMap == 14] = 4
    fuelMap[fuelMap == 13] = 3
    fuelMap[fuelMap == 12] = 0
    fuelMap[fuelMap == 11] = 1
    fuelMap[fuelMap == 10] = 2
    print("Image size", np.shape(fuelMap))
    
#    fuelMap   =  np.zeros(np.shape(np.flipud(np.asarray(Image.open(imgPath)))))   
#    fuelMap[:1200,:] = 1
#    fuelMap[:900,:] = 2
#    fuelMap[:600,:] = 3
#    fuelMap[:300,:] = 4
    
else:
    fuelMap   = 1*np.ones(np.shape(elevation),dtype=('i4'))
    
fig, ax = plt.subplots()
plt.imshow(fuelMap,origin="lower")
plt.plot()
plt.colorbar()



heatFluxModelMap = {}
heatFluxModelMap["name"] =  "heatFlux"
heatFluxModelMap["data"] =  np.zeros(np.shape(elevation),dtype=('i4'))
heatFluxModelMap["table"] =  {"heatFluxBasic":0,}

vaporFluxModelMap = {}
vaporFluxModelMap["name"] =  "vaporFlux"
vaporFluxModelMap["data"] =  np.ones(np.shape(elevation),dtype=('i4'))
vaporFluxModelMap["table"] =  {"vaporFluxBasic":1,}

nc.close()

FiretoNC(fout, domainProperties,parametersProperties,fuelMap,elevation=elevation, wind=None, fluxModelMap = (heatFluxModelMap,vaporFluxModelMap))
