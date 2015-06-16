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
from scipy.io import netcdf
import sys

execfile("genForeFireCase.py")

#if(len(sys.argv)==1):
#    print "Usage MNHCDF2Case MNHSrcFile.nc FFDestFile.nc a.x,a.y b.x,b.y ...."
#    exit(0)
fname = "/Users/filippi/workspace/fireflux2/ffcase/ff2ideal.nc" #sys.argv[1]
fout = "/Users/filippi/workspace/fireflux2/ffcase/case.nc" #sys.argv[2]
d2dvar=["ZS",]
d2dnames=["altitude",]
scalarNames=["THT","UM","VM"]
names=["temperature","windU","windV"]

#geometry 
geomCDFBin = fname 
 
nc = netcdf.netcdf_file(geomCDFBin, 'r') 

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
 
SVM = nc.variables[scalarNames[0]][:]

NK = (len(SVM))-2
NJ = len(SVM[0])-2
NI = len(SVM[0][0])-2

NKT = (len(SVM)) 
NJT = len(SVM[0]) 
NIT = len(SVM[0][0]) 
 
DeltaZ = nc.variables['ZHAT'][1]-nc.variables['ZHAT'][0]
DeltaY = nc.variables['YHAT'][1]-nc.variables['YHAT'][0]
DeltaX = nc.variables['XHAT'][1]-nc.variables['XHAT'][0]

dt = 0
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

parametersProperties= {}
parametersProperties['projectionproperties']  = "41.551998,8.828396,41.551998,8.828396" ;
parametersProperties['date']  = "2009-07-23_12:00:00" ;
parametersProperties['duration']  = 360000;
parametersProperties['projection']  = "OPENMAP" ;
parametersProperties['refYear']  = 2013 ;
parametersProperties['refDay']  = 30 ;

elevation =  nc.variables['ZS'][:,:]
fuelMap   =  5*np.ones(np.shape(elevation),dtype=('i4'))

wind =  {}

wind["zonal"]    = nc.variables['UT'][0,:,:]
wind["meridian"] = nc.variables['VT'][0,:,:]

heatFluxModelMap = {}
heatFluxModelMap["name"] =  "heatFlux"
heatFluxModelMap["data"] =  np.zeros(np.shape(fuelMap),dtype=('i4'))
heatFluxModelMap["table"] =  {"heatFluxNominal":0,}

vaporFluxModelMap = {}
vaporFluxModelMap["name"] =  "vaporFlux"
vaporFluxModelMap["data"] =  np.ones(np.shape(fuelMap),dtype=('i4'))
vaporFluxModelMap["table"] =  {"vaporFluxNominal":1,}

FiretoNC(fout, domainProperties,parametersProperties,fuelMap,elevation=elevation, wind=wind, fluxModelMap = (heatFluxModelMap,vaporFluxModelMap))
