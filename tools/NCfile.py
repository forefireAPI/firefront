import numpy as np
from scipy.io import netcdf
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cbook as cbook

#
# The purpose of this script is to copy a AT map to another one and impose some data (here a png file)
# It may be used as an example to look at how to write AT map for any purposes 
#

fname = "/Users/filippi/workspace/fireflux2/ffcase/ForeFire/ForeFire.0.nc"
fnameout = "/Users/filippi/workspace/fireflux2/ffcase/ForeFire/bmapout.nc"

nc = netcdf.netcdf_file(fname, 'r')
nco = netcdf.netcdf_file(fnameout, 'w')


for dim in nc.dimensions:
    nco.createDimension(dim,nc.dimensions[dim] )
for key in nc.variables.keys():
    print "copying" , key, nc.variables[key].dimensions
    nco.createVariable(key, nc.variables[key][:].dtype, nc.variables[key].dimensions)
    nco.variables[key][:] = np.array(nc.variables[key][:])
    for attvar in nc.variables[key]._attributes:
        print "   copying attribute " , key+ "."+attvar ," = ", nc.variables[key].__dict__[attvar]
        nco.variables[key]._attributes[attvar] = nc.variables[key]._attributes[attvar]


nc.close()

ca2 = nco.variables['cell_active']
ca2[:].fill(1)

image_file = cbook.get_sample_data('/Users/filippi/workspace/fireflux2/ffcase/bmapForced.png')
image = plt.imread(image_file)
values = (image[::-1,:,2] * (8 * 60) ) + 75600
destAT = nco.variables['arrival_time_of_front']
destAT[:] = np.where(values == np.max(values),-9999, values)

nco.close()

# display verification
print np.max(values), np.min(values), np.shape(values), np.shape(destAT)
fig, ax = plt.subplots()
im = ax.imshow(values,origin='lower')
plt.axis('off')
plt.show()
