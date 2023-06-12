 #!python
import sys
import pickle as pkl
import netCDF4
import numpy as np
import sys

list_lat = []
list_lng = []
list_alt = []

dbo = {}


if len(sys.argv) != 2:
     print("usage kmlDomain PGDFile")
     exit(0)
f = netCDF4.Dataset(sys.argv[1],'r')

dbo['LAT'] = f.variables["LAT"][:]
dbo['LON'] = f.variables["LON"][:]
dbo['ALT'] = f.variables["ZSMT"][:]


with open('.coordinates.pickle', 'wb') as dbfile:
	pkl.dump(dbo, dbfile)

#Recuperation de la taille de la matrice
max_x = f.variables['LAT'].shape[0]-1
max_y = f.variables['LAT'].shape[1]-1
#Recuperation des latitudes/longitudes dans un tableau trie par ordre croissant
tab_lat = f.variables['LAT'][:,0]
tab_lon = f.variables['LON'][0,:]

print("""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>kmldomain.kml</name>
<Style id="doms"><PolyStyle><fill>0</fill></PolyStyle></Style>
	<Placemark>
		<name>MNH Bound in PGD</name>
		<styleUrl>#doms</styleUrl>
			<Polygon>
			<tessellate>0</tessellate>
			<outerBoundaryIs>
				<LinearRing>
					<coordinates>
""")
print("%f,%f,0 %f,%f,0 %f,%f,0 %f,%f,0  %f,%f,0"%(f.variables['LON'][0][0],f.variables['LAT'][0][0],f.variables['LON'][0][max_y],f.variables['LAT'][0][max_y],f.variables['LON'][max_x][max_y],f.variables['LAT'][max_x][max_y],f.variables['LON'][max_x][0],f.variables['LAT'][max_x][0],f.variables['LON'][0][0],f.variables['LAT'][0][0]))

print("""
					</coordinates>
				</LinearRing>
			</outerBoundaryIs>
		</Polygon>
	</Placemark>
</Document>
</kml>""")