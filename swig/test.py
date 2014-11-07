import forefire, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches


def getLocationFromLine(line):

    llv = line.split("loc=(")
    if len(llv) < 2: 
        return None
    llr = llv[1].split(",");
    if len(llr) < 3: 
        return None
    return (float(llr[0]),float(llr[1]))

def dist(a,b):
    return np.sqrt(np.power(a[0]-b[0],2)+np.power(a[1]-b[1],2))

def printToPathe(linePrinted):
    fronts = linePrinted.split("FireFront")
    pathes = []
    for front in fronts[1:]:
        nodes =front.split("FireNode")[1:]
        if len(nodes) > 0: 
            Path = mpath.Path
           
            codes = []
            verts = []
            lastNode = getLocationFromLine(nodes[0])
            firstNode = getLocationFromLine(nodes[0])
            
 
            codes.append(Path.MOVETO)
            verts.append(firstNode)      

            for node in nodes[:]:
                newNode = getLocationFromLine(node)
                codes.append(Path.LINETO)
                verts.append(newNode)         
                lastNode = newNode
                
    
            codes.append(Path.LINETO)
            verts.append(firstNode)          
           
            pathes.append(mpath.Path(verts, codes))

    return pathes;




ff = forefire.PLibForeFire()

sizeX = 300
sizeY = 200


ff.setString("fuelsTableFile","fuels.ff")
ff.setString("ForeFireDataDirectory","test")

ff.setDouble("spatialIncrement",.3)
ff.setDouble("minimalPropagativeFrontDepth",0.1)
ff.setDouble("perimeterResolution",1)
ff.setInt("atmoNX",sizeX)
ff.setInt("atmoNY",sizeY)



ff.setDouble("initialFrontDepth",1)
ff.setDouble("relax",.5)
ff.setDouble("smoothing",4)
ff.setDouble("z0",0.)
ff.setDouble("windU",0.)
ff.setDouble("windV",12.)
ff.setInt("defaultFuelType",1)
ff.setInt("bmapLayer",1)

ff.setInt("defaultHeatType",0)
ff.setDouble("nominalHeatFlux",100000)
ff.setDouble("burningDuration",80)
ff.setDouble("maxFrontDepth",50)
ff.setDouble("minSpeed",0.0001)


ff.execute("FireDomain[sw=(0.,0.,0.);ne=(%f,%f,0.);t=0.]"%(sizeX,sizeY))

print "resolution of bmap is ", ff.getString("bmapResolution")

ff.addLayer("data","altitude","z0")
ff.addLayer("data","windU","windU")
ff.addLayer("data","windV","windV")
ff.addLayer("BRatio","BRatio","BRatio")
ff.addLayer("flux","heatFluxBasic","defaultHeatType")
ff.addLayer("propagation","BalbiUnsteady","propagationModel")

fuelmap = np.zeros((sizeX,sizeY,1), dtype=np.int32)

fuelmap[:,90:110,:] = 1
fuelmap[150:,90:110,:] = 3
ff.addIndexLayer("table","fuel",0 , 0, 0, sizeX, sizeY, 0, fuelmap)

print(np.shape(ff.getDoubleArray("fuel")))

ff.execute("\tFireFront[t=0.]")
ff.execute("\t\tFireNode[loc=(40,60,0.);vel=(-1.,-1.,0.);t=0.]")
ff.execute("\t\tFireNode[loc=(40,65,0.);vel=(-1.,1.,0.);t=0.]")
ff.execute("\t\tFireNode[loc=(260,65,0.);vel=(1.,1.,0.);t=0.]")
ff.execute("\t\tFireNode[loc=(260,60,0.);vel=(1.,-1.,0.);t=0.]")

pathes = []
step = 20
for i in range(1,20):
    print "goTo[t=%f]"%(i*step)
    ff.execute("goTo[t=%f]"%(i*step))
    pathes += printToPathe( ff.execute("print[]"))


fig, ax = plt.subplots()
 
tab = np.transpose(ff.getDoubleArray("BMap"))[0]
CS = ax.imshow(tab, cmap=plt.cm.gray, interpolation='nearest')
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('v')

for path in pathes:
    patch = mpatches.PathPatch(path,edgecolor='red', facecolor='none', alpha=1)
    ax.add_patch(patch)

ax.grid()
ax.axis('equal')
plt.show()






