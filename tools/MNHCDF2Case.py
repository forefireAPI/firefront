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

import netCDF4 as nc4

    

execfile("genForeFireCase.py")

#if(len(sys.argv)==1):
#    print "Usage MNHCDF2Case MNHSrcFile.nc FFDestFile.nc a.x,a.y b.x,b.y ...."
#    exit(0)
#fname = "/Volumes/brando/cases/paugham/02_mesonh/ff2ideal.nc" #sys.argv[1]
fname="/Volumes/brando/cases/couleeLaveLES/STRAP.1.SEG01.001.nc4"


fout = "../swig/lava/caseTule.nc" #sys.argv[2]
d2dvar=["ZS",]
d2dnames=["altitude",]
scalarNames=["THT","UM","VM"]
names=["temperature","windU","windV"]

#geometry 
geomCDFBin = fname 
 
nc = nc4.Dataset(geomCDFBin, 'r') 

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

dom= "FireDomain[sw=(%d,%f,0);ne=(%f,%f,0);]"%(domainProperties['SWx'],domainProperties['SWy'],domainProperties['SWx']+domainProperties['Lx'],domainProperties['SWy']+domainProperties['Ly'])
print domainProperties['Lx'] ,domainProperties['Ly'] ,len(nc.variables['XHAT'][:])
print "   ori ", nc.variables['LAT0'].getValue(), nc.variables['LON0'].getValue()
 
 
SVM = nc.variables['ZS'][:]

 
NJ = len(SVM)-2
NI = len(SVM[0])-2

 
NJT = len(SVM) 
NIT = len(SVM[0]) 
 

DeltaY = nc.variables['YHAT'][1]-nc.variables['YHAT'][0]
DeltaX = nc.variables['XHAT'][1]-nc.variables['XHAT'][0]

dt = 0
lastTime = 0
lastTime = nc.variables['DTCUR__TIME'].getValue()
firstTime = nc.variables['DTCUR__TIME'].getValue()

domainProperties= {}
domainProperties['SWx']  = nc.variables['XHAT'][0]
domainProperties['SWy']  = nc.variables['YHAT'][0]
domainProperties['SWz']  = 0
domainProperties['Lx']   = nc.variables['XHAT'][-1]+DeltaX-nc.variables['XHAT'][0]
domainProperties['Ly']   = nc.variables['YHAT'][-1]+DeltaY-nc.variables['YHAT'][0]
domainProperties['Lz']   = 0
domainProperties['t0']   = 0
domainProperties['Lt']   = np.Inf

dom= "FireDomain[sw=(%d,%f,0);ne=(%f,%f,0);t=%f]"%(domainProperties['SWx'],domainProperties['SWy'],domainProperties['SWx']+domainProperties['Lx'],domainProperties['SWy']+domainProperties['Ly'],lastTime)
print dom
print "\n".join(nc.variables)

parametersProperties= {}
parametersProperties['projectionproperties']  = "41.551998,8.828396,41.551998,8.828396" ;
parametersProperties['date']  = "2009-07-23_12:00:00" ;
parametersProperties['duration']  = 360000;
parametersProperties['projection']  = "OPENMAP" ;
parametersProperties['refYear']  = 2013 ;
parametersProperties['refDay']  = 30 ;


elevation =  nc.variables['ZS'][:,:]

# Tout le combustible a 1
fuelMap   =  1*np.ones(np.shape(elevation),dtype=('i4'))

wind =  {}

wind["zonal"]    = np.array(nc.variables['UT'][0,:,:])
wind["meridian"] = np.array(nc.variables['VT'][0,:,:])

heatFluxModelMap = {}
heatFluxModelMap["name"] =  "heatFlux"
heatFluxModelMap["data"] =  np.zeros(np.shape(fuelMap),dtype=('i4'))
heatFluxModelMap["table"] =  {"heatFluxBasic":0,}

vaporFluxModelMap = {}
vaporFluxModelMap["name"] =  "vaporFlux"
vaporFluxModelMap["data"] =  np.ones(np.shape(fuelMap),dtype=('i4'))
vaporFluxModelMap["table"] =  {"vaporFluxBasic":1,}

nc.close()

FiretoNC(fout, domainProperties,parametersProperties,fuelMap,elevation=elevation, wind=wind, fluxModelMap = (heatFluxModelMap,vaporFluxModelMap))
