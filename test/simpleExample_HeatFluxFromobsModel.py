import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import pdb 
import imp 
from scipy.io import netcdf
import glob
import datetime
import socket 
import urllib
sys.path.append('../swig/')
import forefire

##################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


##################################################
def download_inputFile():
    url_nc = "https://www.dropbox.com/s/hxd4ivkkvoxbe8n/mydata.nc?dl=1"  # dl=1 is important
    filename_nc = "./Inputs/mydata.nc"

    url_bmap = "https://www.dropbox.com/s/q57wy04219a0m4x/mybmap.nc?dl=1"  # dl=1 is important
    file_name_bmap = "./Inputs/mybmap.nc"

    url_fuel = "https://www.dropbox.com/s/m3bri1n85c5yjn2/fuels.ff?dl=1"  # dl=1 is important
    file_name_fuel = "./Inputs/fuels.ff"

    for url, file_name in zip([url_nc,url_bmap,url_fuel],[filename_nc,file_name_bmap,file_name_fuel]):
        urllib.urlretrieve (url, file_name)

###################################################
def getDomainExtent(line):
    print line
    llv = line.split("sw=(")
    llr = llv[1].split("ne=(");
    return( float( llr[0].split(",")[0]), float(llr[1].split(",")[0]), float(llr[0].split(",")[1]), float(llr[1].split(",")[1]) )

####################################################
#start here

#set input     
if not os.path.isfile('./Inputs/mydata.nc'):
    ensure_dir('./Inputs/')
    ensure_dir('./ForeFire/')
    ensure_dir('./ForeFire/Outputs/')
    ensure_dir('./MODEL1/')
    download_inputFile()


ff = forefire.PLibForeFire()


#set param
ff.setString('ForeFireDataDirectory','Inputs')
ff.setString('fireOutputDirectory','ForeFire/Outputs')
ff.setInt('outputsUpdate',10)

ff.setString('NetCDFfile',    'mydata.nc')
ff.setString('fluxNetCDFfile','mydata.nc')
ff.setString('fuelsTableFile','fuels.ff')
ff.setString('BMapFiles', 'mybmap.nc')

ff.setDouble("spatialIncrement",.3)
ff.setDouble("perimeterResolution",1)
ff.setDouble("minimalPropagativeFrontDepth",1)

ff.setDouble("nominalHeatFlux",1.e6)
ff.setDouble("nominalVaporFlux",0.2)
ff.setDouble("burningDuration",7200)

ff.setDouble("bmapOutputUpdate",300)

#ff.setInt("defaultHeatType",0)
#ff.setInt("defaultFRPType",2)
#ff.setInt("defaultVaporType",1)

#ff.setInt('bmapLayer',1)
ff.setDouble("InitTime",1.e6)

#set domain
ff.setInt("atmoNX",62)
ff.setInt("atmoNY",62)
ff.execute("FireDomain[sw=(%f,%f,0.);ne=(%f,%f,0.);t=%f]"%(-50,-50,\
                                                            3050,3050,  \
                                                            0))
atmo_dx=62
atmo_dy=62

extentLocal= getDomainExtent(ff.execute("print[]").split("\n")[0]);


#set propagation model
ff.addLayer("propagation","TroisPourcent","propagationModel")
#ff.addLayer("flux",ronan_param['heatFluxModel'],"defaultHeatType")
#ff.addLayer("flux","FRP","defaultFRPType")
#ff.addLayer("flux","vaporFluxBasic","defaultVaporType")

#residenceTime = np.zeros((nx,ny,1)) + 60
#ff.addScalarLayer("double","residenceTime",0 , 0, 0,extentLocal[1]-extentLocal[0], extentLocal[3]-extentLocal[2] , 0, residenceTime)

fuelmap = ff.getDoubleArray("fuel").astype("int32")
ff.addIndexLayer("table","fuel", extentLocal[0], extentLocal[2],0,  extentLocal[1]-extentLocal[0], extentLocal[3]-extentLocal[2], 0, fuelmap)

print "resolution of bmap is ", ff.getString("bmapResolution")


#---------------------------------------------
# run ForeFire simulation
#---------------------------------------------
pathes = []
step = 10
N_step = 200/step
flux_out_ff_history = []
for i in np.arange(1,N_step):
    
    ff_time = i*step

    print "goTo[t=%f]"%(i*step)
    ff.execute("goTo[t=%f]"%(i*step))
    
    FRP = ff.getDoubleArray("FRP")[:,:,0]
    flux2d = ff.getDoubleArray('heatFlux')[:,:,0]
    
    ax = plt.subplot(121)
    ax.imshow(FRP.T,origin='lower',interpolation='nearest')
    bx = plt.subplot(122)
    bx.imshow(flux2d.T,origin='lower',interpolation='nearest')
    plt.show()
    

