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

'''
netcdf dataROB {
dimensions:
	DIMX = 302 ;
	DIMY = 302 ;
	DIMZ = 1 ;
	DIMT = 1 ;
variables:
	int fuel(DIMT, DIMZ, DIMY, DIMX) ;
		fuel:type = "fuel" ;
	int heatFlux(DIMT, DIMZ, DIMY, DIMX) ;
		heatFlux:type = "flux" ;
		heatFlux:model0name = "IndFlux" ;
		heatFlux:indices = 0 ;
	int vaporFlux(DIMT, DIMZ, DIMY, DIMX) ;
		vaporFlux:type = "flux" ;
		vaporFlux:model1name = "IndVap" ;
		vaporFlux:indices = 1 ;
	char domain ;
		domain:type = "domain" ;
		domain:SWx = 288520. ;
		domain:SWy = 282520. ;
		domain:SWz = 0 ;
		domain:Lx = 24160. ;
		domain:Ly = 24160. ;
		domain:Lz = 0 ;
		domain:t0 = 0 ;
		domain:Lt = Infinityf ;
	char parameters ;
		parameters:type = "parameters" ;
		parameters:projectionproperties = "41.551998,8.828396,41.551998,8.828396" ;
		parameters:date = "2017-06-17_14:00:00" ;
		parameters:duration = 100000 ;
		parameters:projection = "OPENMAP" ;
		parameters:refYear = 2017 ;
		parameters:refDay = 168 ;

// global attributes:
		:version = "FF.1.0" ;
}

    '''
import numpy as np
import sys
import xarray as xr
import netCDF4 as nc4

def forceDim(arrIn, newD=None, kind='edge'):
    print(np.shape(arrIn))
    if (newD[0][1] < 0):
        newr = arrIn[newD[0][0]:newD[0][1],newD[1][0]:newD[1][1]]
        print("new shape ",np.shape(newr))
        return(newr)
    
    if(newD==None):
        return arrIn
    return np.pad(arrIn, newD, kind)


execfile("genForeFireCase.py")

#if(len(sys.argv)==1):
#    print "Usage MNHCDF2Case MNHSrcFile.nc FFDestFile.nc a.x,a.y b.x,b.y ...."
#    exit(0)
#fname = "/Volumes/brando/cases/paugham/02_mesonh/ff2ideal.nc" #sys.argv[1]
fname="/Users/filippi_j/soft/MNH-V5-5-1/MY_RUN/KTEST/016_PEDROGAO/002_mesonh/ForeFire/RELIEF3D.nc"


fout = "/Users/filippi_j/soft/MNH-V5-5-1/MY_RUN/KTEST/016_PEDROGAO/002_mesonh/ForeFire/ppsimple.nc" #sys.argv[2]
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
print( domainProperties['Lx'] ,domainProperties['Ly'] ,len(nc.variables['XHAT'][:]))
print ("   ori ", nc.variables['LAT0'].getValue(), nc.variables['LON0'].getValue())
 
 
SVM = nc.variables['ZS'][:]

 
NJ = len(SVM)-2
NI = len(SVM[0])-2

 
NJT = len(SVM) 
NIT = len(SVM[0]) 
 
padDim = ((4,4),(4,4))
#padDim = ((135,-135),(135,-135))


DeltaY = nc.variables['YHAT'][1]-nc.variables['YHAT'][0]
DeltaX = nc.variables['XHAT'][1]-nc.variables['XHAT'][0]

dt = 0
lastTime = 0

lastTime = 0#nc.variables['DTCUR'].getValue()
firstTime = 0#nc.variables['DTCUR'].getValue()


domainProperties= {}
domainProperties['SWx']  = -400#nc.variables['XHAT'][0]
domainProperties['SWy']  = -400#nc.variables['YHAT'][0]
domainProperties['SWz']  = 0
domainProperties['Lx']   = 24400#nc.variables['XHAT'][-1]+DeltaX-nc.variables['XHAT'][0]
domainProperties['Ly']   = 24400#nc.variables['YHAT'][-1]+DeltaY-nc.variables['YHAT'][0]
domainProperties['Lz']   = 0
domainProperties['t0']   = 0
domainProperties['Lt']   = np.Inf

dom= "FireDomain[sw=(%d,%f,0);ne=(%f,%f,0);t=%f]"%(domainProperties['SWx'],domainProperties['SWy'],domainProperties['SWx']+domainProperties['Lx'],domainProperties['SWy']+domainProperties['Ly'],lastTime)
 

parametersProperties= {}
parametersProperties['projectionproperties']  = "41.551998,8.828396,41.551998,8.828396" ;
parametersProperties['date']  = "2017-06-17_14:00:00" ;
parametersProperties['duration']  = 32400;
parametersProperties['projection']  = "OPENMAP" ;
parametersProperties['refYear']  = 2017 ;
parametersProperties['refDay']  = 168 ;


elevation =  nc.variables['ZS'][:,:]


# Tout le combustible a 1
fuelMap   =  1*np.ones(np.shape(elevation),dtype=('i4'))

wind =  {}
 

wind["zonal"]    = np.array(nc.variables['UT'][0,1,:,:])
wind["meridian"] = np.array(nc.variables['VT'][0,1,:,:])

print("max wind is",np.max(wind["zonal"] ),np.shape(forceDim(wind["zonal"],padDim)))

heatFluxModelMap = {}
heatFluxModelMap["name"] =  "heatFlux"
heatFluxModelMap["data"] =  forceDim(np.zeros(np.shape(fuelMap),dtype=('i4')),padDim)
heatFluxModelMap["table"] =  {"heatFluxBasic":0,}

vaporFluxModelMap = {}
vaporFluxModelMap["name"] =  "vaporFlux"
vaporFluxModelMap["data"] =  forceDim(np.ones(np.shape(fuelMap),dtype=('i4')),padDim)
vaporFluxModelMap["table"] =  {"vaporFluxBasic":1,}

nc.close()

FiretoNC(fout, domainProperties,parametersProperties,forceDim(fuelMap,padDim),elevation=None, wind=None, fluxModelMap = (heatFluxModelMap,vaporFluxModelMap))
