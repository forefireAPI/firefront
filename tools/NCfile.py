import numpy as np
from scipy.io import netcdf
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cbook as cbook

fname = "/Users/filippi/workspace/fireflux2/ffcase/ForeFire/ForeFire.0.nc"
fnameout = "/Users/filippi/workspace/fireflux2/ffcase/ForeFire/bmapout.nc"

nc = netcdf.netcdf_file(fname, 'r')
nco = netcdf.netcdf_file(fnameout, 'w')

print nc.variables
print nc.dimensions 

for dim in nc.dimensions:
    print dim
    nco.createDimension(dim,nc.dimensions[dim] )

domain = {}
domain["SWx"] = -10. 
domain["SWy"] = -10. 
domain["Lx"] = 1020. 
domain["Ly"] = 2020. 
domain["Lz"] = 0 
domain["refYear"] = 2013 
domain["refDay"] = 30 

nco.createVariable("arrival_time_of_front", nc.variables["arrival_time_of_front"][:].dtype, nc.variables["arrival_time_of_front"].dimensions)
nco.createVariable("cell_active", nc.variables["cell_active"][:].dtype, nc.variables["cell_active"].dimensions)
nco.createVariable("domain", nc.variables["domain"][:].dtype, nc.variables["domain"].dimensions)

toto = nc.variables['arrival_time_of_front']
ca = nc.variables['cell_active']
ca2 = nco.variables['cell_active']

do2 = nco.variables['domain']
do2.SWx = -10. 
do2.SWy = -10. 
do2.Lx = 1020. 
do2.Ly = 2020. 
do2.Lz = 0 
do2.refYear = 2013 
do2.refDay = 30 



ca2[:].fill(1)

image_file = cbook.get_sample_data('/Users/filippi/workspace/fireflux2/ffcase/bmapForced.png')
image = plt.imread(image_file)
values = (image[::-1,:,2] * (8 * 60) ) + 75600

v2 = np.where(values == np.max(values),-9999, values)
 

toto2 = nco.variables['arrival_time_of_front']
toto2[:] = v2

nco.close()

print np.max(values), np.min(values), np.shape(values), np.shape(toto)
#toto = values
#nc.close()
fig, ax = plt.subplots()
im = ax.imshow(values)
#patch = patches.Circle((260, 200), radius=200, transform=ax.transData)
#im.set_clip_path(patch)

plt.axis('off')
plt.show()