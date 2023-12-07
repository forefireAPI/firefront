import numpy as np
from scipy.io import netcdf
import sys
from scipy.interpolate import RectBivariateSpline

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cbook as cbook
from scipy.interpolate import CloughTocher2DInterpolator
from scipy import interpolate
import random

fnameShape = "/Users/filippi/data/pedrograo/ForeFireMOD.1.nc"
fnameout = "/Users/filippi/data/pedrograo/ForeFire.interp.nc"

dataTimed = []
dataFile = '/Users/filippi/Desktop/image2567b.png' # pedrogao
dataFile = '/Users/filippi/Desktop/imageBMapqueaius.png' # pedrogao

def readCoordinatesInKmlFile(filename,value):
    data=""
    x=[]
    y=[]
    z=[]
    with open(filename, 'r') as file:
        data = file.read()
    lls = data.split("<coordinates>")[1].split("</coordinates>")[0].replace('\t', '').replace('\n', '').split(" ")
    for ll in lls:
        llvals = ll.split(",")
        if len(llvals) > 2:
            x.append(float(llvals[0]))
            y.append(float(llvals[1]))
            z.append(value)
            
    return x,y,z

contoursKML = {14.5:"/Users/filippi/data/pedrograo/kml2bmap/14v5.kml",
               15:"/Users/filippi/data/pedrograo/kml2bmap/15.kml",
               16:"/Users/filippi/data/pedrograo/kml2bmap/16.kml",
               17:"/Users/filippi/data/pedrograo/kml2bmap/17.kml",
               18:"/Users/filippi/data/pedrograo/kml2bmap/18.kml",
               18:"/Users/filippi/data/pedrograo/kml2bmap/18.kml",
               19:"/Users/filippi/data/pedrograo/kml2bmap/19.kml",
               20:"/Users/filippi/data/pedrograo/kml2bmap/20.kml",
               21:"/Users/filippi/data/pedrograo/kml2bmap/21.kml",
               22:"/Users/filippi/data/pedrograo/kml2bmap/22.kml",
               23:"/Users/filippi/data/pedrograo/kml2bmap/23.kml",
               24:"/Users/filippi/data/pedrograo/kml2bmap/24.kml"
              }

bx,by,bz= readCoordinatesInKmlFile("/Users/filippi/data/pedrograo/kml2bmap/bounds.kml",10)


vx=[]
vy=[]
vz=[]

for key in contoursKML.keys():
    ax,ay,az= readCoordinatesInKmlFile(contoursKML[key],key)
    vx = vx + ax
    vy = vy + ay
    vz = vz + az
    print az

    
x = np.array(vx)
y = np.array(vy)
z = np.array(vz)

#for i in range(dim):
#    x[i] = x[i] * random.random() * 10000
#    y[i] = y[i] * random.random() * 10000
#    z[i] = z[i] * random.random() * 55000
    
#y = radom.random() - 0.5

print np.shape(z), np.shape(x), np.shape(y) 

X = np.linspace(min(bx), max(bx),num=1208)
Y = np.linspace(min(by), max(by),num=1208)

X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation
#interp = CloughTocher2DInterpolator(list(zip(x, y)), z)
interp = interpolate.griddata((x, y), z, (X, Y), method='linear')

Z = interp#(X, Y)
Z= np.where(Z < 24, Z, 24)

Z= np.where(Z > 14, Z, 14)
print np.min(Z),np.max(Z),np.min(z),np.max(z)

#plt.pcolormesh(X, Y, Z, shading='auto')
plt.imshow(Z,origin='lower')

#plt.plot(x, y, "ok", label="input point")
plt.legend()
plt.colorbar()
plt.axis("equal")
plt.show()
# je charge les donnees
exit(0)
# parameters
# origin time
# otime = 48600 #pedrogao 3fires reference time 
otime = 11*3600 # queaius
# Shape 
FLx = 24160.0#300*300 case reference
FLy = 24160.0#300*300 case  reference
#FSWx = 288520.0+5*400#pedrogao 3fires reference
#FSWy = 282520.0+6*400#pedrogao 3fires reference
FSWx = 428520.0#queaius  reference
FSWy = 190520.00#queaius  reference
#date 
#FrefYear = 2017 #pedrogao 
# FrefDay = 168 #pedrogao 
FrefYear = 2017 # queaius
FrefDay = 288 # queaius
# Duration .. how long spread between min and max in bmap
#NHours = 11 # pedrogao
NHours = 15

image_file = cbook.get_sample_data(dataFile)
raw = (plt.imread(image_file))
print np.max(raw[:,:,1]),np.min(raw[:,:,1])
sumRaw =  raw[:,:,1]
im = sumRaw*3600*NHours + otime
print np.shape(im)
image = im 
image = np.where(image == np.max(image),-9999, image)
im=image.astype(int)


print np.shape(image), np.max(image),np.min(image)

nc = netcdf.netcdf_file(fname, 'r')
nco = netcdf.netcdf_file(fnameout, 'w')
values = image.astype(int)

tvar={}
tvar["arrival_time_of_front"]=[["DIMY", "DIMX"] ,values]
tvar["cell_active"]=[["C_DIMY", "C_DIMX"] ,np.ones([302,302])]
tvar["domain"]=[["domdim"] , [1]]
for key in tvar.keys():
    varn = tvar[key][0]
    vsize = np.shape(tvar[key][1])
    for i in range(len(varn)):
        nco.createDimension(varn[i], vsize[i])
        print varn[i], vsize[i]


for key in ["arrival_time_of_front","cell_active"]:
    print "copying" , key, nc.variables[key].dimensions
    nco.createVariable(key, nc.variables[key][:].dtype, tvar[key][0])
    nco.variables[key][:] = tvar[key][1]
nco.createVariable("domain", nc.variables["domain"][:].dtype, tvar["domain"][0])   
for attvar in nc.variables["domain"]._attributes:
    print "   copying attribute ."+attvar ," = ", nc.variables["domain"].__dict__[attvar]
    nco.variables["domain"]._attributes[attvar] = nc.variables["domain"]._attributes[attvar]

nco.variables["domain"]._attributes["Lx"] = FLx
nco.variables["domain"]._attributes["Ly"] = FLy
nco.variables["domain"]._attributes["SWx"] = FSWx
nco.variables["domain"]._attributes["SWy"] = FSWy
nco.variables["domain"]._attributes["refYear"] = FrefYear 
nco.variables["domain"]._attributes["refDay"] = FrefDay 
nc.close()


ca2 = nco.variables['cell_active']
ca2[:].fill(1)


destAT = nco.variables['arrival_time_of_front']
print np.shape(image), np.max(image),np.min(image), np.shape(destAT[:])
destAT[:] = np.flipud(values)
fig, ax = plt.subplots()
ax.contour(destAT[:], levels=[ 3600*13.5,3600*15, 3600*17, 3600*19, 3600*21],
                       colors=['orange', 'red', 'yellow', 'green', '#A0A0A0', '#A0B0A0'], extend='both',origin='lower')
ax.imshow(destAT[:],origin='lower')
 

ax.axis('off')
plt.show()



# display verification
print np.max(destAT[:] ), np.min(destAT[:] ), np.shape(values), np.shape(destAT)
nco.close()
