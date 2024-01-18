#!/usr/bin/python
'''
Created on 4 avr. 2011

@author: filippi
'''
 
from evtk.vtk import *  
import os 
import numpy as np 
import struct
import csv
import shutil
from copy import copy


def pointsToVTK(path, x, y, z, data):
    """
        Export points and associated data as an unstructured grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D arrays with coordinates of the points.
            data: dictionary with variables associated to each point.
                  Keys should be the names of the variable stored in each array.
                  All arrays must have the same number of elements.

        RETURNS:
            Full path to saved file.

    """
    assert (x.size == y.size == z.size)
    npoints = x.size
 
    # create some temporary arrays to write grid topology
    offsets = np.arange(start = 1, stop = npoints + 1, dtype = 'int32')   # index of last node in each cell
    connectivity = np.arange(npoints, dtype = 'int32')                   # each point is only connected to itself
    cell_types = np.empty(npoints, dtype = 'uint8') 
   
    cell_types[:] = VtkVertex.tid

    w = VtkFile(path, VtkUnstructuredGrid)
    w.openGrid()
    w.openPiece(ncells = npoints, npoints = npoints)
    
    w.openElement("Points")
    w.addData("points", (x,y,z))
    w.closeElement("Points")
    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")
    
    _addDataToFile(w, cellData = None, pointData = data)

    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, cellData = None, pointData = data)

    w.save()
    return w.getFileName()

def imageToVTK(path, origin = (0.0,0.0,0.0), spacing = (1.0,1.0,1.0), cellData = None, pointData = None, start =   (0,0,0), end = None):

    assert (cellData <> None or pointData <> None)
    
    # Extract dimensions
 
    if cellData <> None:
        keys = cellData.keys()
        data = cellData[keys[0]]
        if end ==None:
            end = data.shape
    elif pointData <> None:
        keys = pointData.keys()
        data = pointData[keys[0]]
        if end ==None:
            end = data.shape
            end = (end[0] - 1, end[1] - 1, end[1] - 1)

    # Write data to file
    w = VtkFile(path, VtkImageData)
    w.openGrid(start = start, end = end, origin = origin, spacing = spacing)
    w.openPiece(start = start, end = end)
    _addDataToFile(w, cellData, pointData)
    w.closePiece()
    w.closeGrid()
    _appendDataToFile(w, cellData, pointData)
    w.save()
    return w.getFileName()
 
   
def _addDataToFile(vtkFile, cellData, pointData):
    # Point data
    if pointData <> None:
        keys = pointData.keys()

        vtkFile.openData("Point", scalars = keys[0])
        for key in keys:
            data = pointData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Point")

    # Cell data
    if cellData <> None:
        keys = cellData.keys()
        vtkFile.openData("Cell", scalars = keys[0])
        for key in keys:
            data = cellData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Cell")

def _appendDataToFile(vtkFile, cellData, pointData):
    # Append data to binary section
    if pointData <> None:
        keys = pointData.keys()
        for key in keys:
            data = pointData[key]
            vtkFile.appendData(data)

    if cellData <> None:
        keys = cellData.keys()
        for key in keys:
            data = cellData[key]
            vtkFile.appendData(data)

def iround(x):
    """iround(number) -> integer
    Round a number to the nearest integer."""
    return int(round(x) - .5) + (x > 0)
 


def gridToVTK(path, x, y, z, cellData = None, pointData = None, vectData = None,start = None, end = None):

    if(start == None):
        start = (0,0,0)
    

    ftype = VtkStructuredGrid
 
    s = x.shape
    nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1
  

    if(end == None):
        end = (nx, ny, nz)
  

    w = VtkFile(path, ftype)
    w.openGrid(start = start, end = end)
    w.openPiece(start = start, end = end)

    w.openElement("Points")
    w.addData("points", (x,y,z))
    w.closeElement("Points")
    

    _addDataToFile(w, cellData, pointData)
    w.closePiece()
    w.closeGrid()
    # Write coordinates
    w.appendData( (x,y,z) )
    # Write data
    _appendDataToFile(w, cellData, pointData)
    w.save()
    return w.getFileName()



def genATMapImage(fname,fnameZ):
    from scipy.io import netcdf
    from scipy import interpolate
    
    #fname = "/Volumes/stooges data/jonathan2014/observed.2.nc" #sys.argv[1]
    nc = netcdf.netcdf_file(fname, 'r') 
    at = nc.variables['arrival_time_of_front'][:]
    shpAT =  at.shape
    nx, ny, nz = shpAT[1], shpAT[0], 1 
    
    #fnameZ = "/Volumes/stooges data/jonathan2014/3DOM1.3.TES21.003.cdf"
    ncz = netcdf.netcdf_file(fnameZ, 'r') 
    zgrid = np.transpose(ncz.variables['ZSBIS'][:,:])
    
    ox = 167900+150
    oy = 85400+150
    
    ex = 169100-100
    ey = 86100-100
    
    res= 9#(ex-ox)/(nx+1)
    print res, nx, ny, ex-ox
    shpZ =  zgrid.shape
    znx, zny = shpZ[0], shpZ[1] 
    
    kernelIn = zgrid
    inKSize = len(kernelIn)
    
    kernelOut = np.zeros((znx,zny),np.float64)
    
    projx = np.array(range(znx))
    projy = np.array(range(zny))
    
    projz = kernelIn
    xx = np.linspace(projx.min(),projx.max(),nx+1)
    yy = np.linspace(projy.min(),projy.max(),ny+1)
    
    newKernel = interpolate.RectBivariateSpline(projx,projy,projz, kx=2,ky=2)
    
    kernelOut = newKernel(xx,yy)
    
    print kernelOut.shape, shpAT, shpZ, ex, ey
    
    
    from evtk.hl import imageToVTK 
    # Dimensions 
    
    X = np.arange(ox-50, 50+ox+(res*(nx+1)), res, dtype='float64') 
    Y = np.arange(oy-50, 50+oy+(res*(ny+1)), res, dtype='float64') 
    
    x = np.zeros(( nx + 1,ny + 1, nz + 1)) 
    y = np.zeros(( nx + 1,ny + 1, nz + 1)) 
    z = np.zeros(( nx + 1,ny + 1, nz + 1)) 
    
    # We add some random fluctuation to make the grid more interesting 
    for k in range(nz + 1): 
        for j in range(ny + 1):
            for i in range(nx + 1): 
                x[i,j,k] = X[i]
                y[i,j,k] = Y[j]
                z[i,j,k] = kernelOut[i,j]
    # Variables 
    
    atime = np.zeros((nx , ny , nz )) 
    
    print x.shape, y.shape, z.shape
    
    atime[:,:,0] = np.transpose(at)
    
    gridToVTK("./image", x, y, z, cellData = {"atime" : atime}, pointData = None)
    
    
    exit(1)

genATMapImage("/Users/filippi/observed.2.nc","/Users/filippi/FINAL.3.GAR07.001dKCL.nc")

inpattern = ""
outPath = ""
inFFpattern = inpattern

if(len(sys.argv)==1):
    print "Usage PARALLELMNH2SINGLEVTK.py infilePattern [in Domain File Pattern] [out path] +(one option \"-norecompute\" \"-clean\" \"-hdf\" \"-lidar lidarfilein lidarfileout\" \"-steps startgenstep endgenstep\" \"-gendomain swx;swy width;height\")\n example \"PARALLELMNH2SINGLEVTK.py  MODEL2/output ForeFire/output vtkOutputs/ "
    
   # exit(1)




cleanFile = False

lidarIn = None
lidarOut = None
startStep = -1
endStep = -1
genDomainOrigin = None
genDomainExtent = None
norecompute=False

if(len(sys.argv)>1):
    inpattern = sys.argv[1]
    inFFpattern = sys.argv[1]
    outPath = "."
    
if(len(sys.argv)>1):
    inFFpattern =  sys.argv[2]
    
if(len(sys.argv)>2):
    outPath =  sys.argv[3]

if(len(sys.argv)>4):
    if(sys.argv[4] == "-clean"):
        print "saving in hdf"
        cleanFile = False
    if(sys.argv[4] == "-hdf"):
        print "saving in hdf"
    if(sys.argv[4] == "-lidar"):
        if(len(sys.argv)>5):
            lidarIn  =  sys.argv[5]
            lidarOut =  sys.argv[6]    
            print "Simulating lidar scan"
    if(sys.argv[4] == "-steps"):
        if(len(sys.argv)>=5):
            startStep =  int(sys.argv[5])
            endStep =  int(sys.argv[6])    
            print "Gen from step%d to %d"%(startStep,endStep)
    if(sys.argv[4] == "-norecompute"):
        norecompute=True
    if(sys.argv[4] == "-gendomain"):
        if(len(sys.argv)>=5):
            genDomainOrigin =  sys.argv[5].split(",")
            genDomainExtent =  sys.argv[6].split(",")  
            print "Regenerates minimum domain shape"
                
                
#inpattern = "/labos/lacy/forefire/gaz2/MODEL3/output"
#inFFpattern = "/labos/lacy/forefire/gaz2/ffout3/output"
#outPath = "/labos/lacy/forefire/gaz2/vtk3/"
  
fprefix=inpattern.split("/")[-1]

gridKey  ="ZGRID"
bmapKey="bmap"

scals={}
scals["cells"] = ()
#scals["points"] = ("SO2","TKE")
scals["points"] = ("HCL","CO2","SO2","TKE","P","T","moist")
inputcdfFile= "/Volumes/stooges data/jonathan2014/3DOM1.3.TES21.003.cdf"
inputcdfvars = ()

#scals["points"] = ("ZGRID",)
Vects={}
Vects["Wind"] = ("U","V","W")

appendedSteps=None;

domains=range(1,2)
vtpShapes3Dout=False
lidarList={}
smokeLocList=[]
probeFile = None 

if lidarIn != None:
    for filename in os.listdir("%s"%(lidarIn)):
        if filename.endswith("txt"):
            shot = LaserShot("%s%s"%(lidarIn,filename),origin=(8150,2180,1440,8030-300))
            lidarList[shot.time] = (filename,shot)
    probeFile = open("%s/windInSmoke.csv"%lidarOut,'w')
    probeFile.write("time;u;v;w;average of speeds;average of modules\n") 


tsteps = []
stepzgrid = ""
varsDataIn = Vects["Wind"] + scals["points"] + scals["cells"]

 



domFFFileNames = sorted(os.listdir("%s"%(inFFpattern[0:len(inFFpattern)-(len(fprefix)+1)])))  

domFFFileName = domFFFileNames[0]  

for firstDOmEnc in domFFFileNames:
    fvars = firstDOmEnc.split(".")
    print "looking for fully composed initial domain file step",fvars, firstDOmEnc
    
    if fvars[0]== fprefix:
        domFFFileName = firstDOmEnc  
         
    break
  
domFFBaseNumber = int(domFFFileName.split(".")[-1])
print "Reading ", inpattern, " data ", varsDataIn, " FF Domain Base step ", domFFBaseNumber

for filename in os.listdir("%s"%(inpattern[0:len(inpattern)-len(fprefix)])):    
    fvars = filename.split(".")
    print filename
    if fvars[0]== fprefix:
        print fvars
        if (fvars[1] == "1" and len(fvars) == 4 and fvars[2].isdigit() and fvars[3]== varsDataIn[0]):
            tsteps.append(float(fvars[2]))
        if (fvars[1] == "1" and len(fvars) == 4 and fvars[3]== gridKey):
            stepzgrid = fvars[2]
            
        if (fvars[1] == "1" and len(fvars) == 3 and fvars[2] == varsDataIn[1]):
            print  "appending %s.%s.%s"%(inpattern,fvars[1],fvars[2])
            appendedSteps=readAllSteps("%s.%s.%s"%(inpattern,fvars[1],fvars[2]))
            print "steps ",appendedSteps
            tsteps = appendedSteps.keys()
            stepzgrid = tsteps[0]
            break

 
            
Allsteps = np.sort(tsteps)



steps = []
gf = VtkGroup("%s/%sgroupfields"%(outPath,fprefix))
for stepV in Allsteps[:]:
    outname = "%s/%s.full.%d.vts"%(outPath,fprefix,stepV) 
    fpoutGname = "%s/%sPoints%d.vtu"%(outPath,fprefix,stepV);
    if (cleanFile and os.path.isfile(outname) and os.path.isfile(fpoutGname) ):
        print "Step %d already post-processed, cleaning up"%stepV
        delCmd =  "rm %s.*.%d"%(inFFpattern,stepV)
        print delCmd
        os.system(delCmd)
        for varName in varsDataIn:
            delCmd=  "rm %s.*.%d.%s"%(inpattern,stepV,varName)
            print delCmd
            os.system(delCmd)
        
    if ((norecompute or cleanFile) and os.path.isfile(outname)):
        print "Step %d already post-processed, not recomputing"%stepV       
    else:
        steps.append(stepV)   
        
    
    gf.addFile(filepath = "%s"%outname, sim_time = stepV)
gf.save()


if len(steps) < 1:
    print "nothing more to be done"
    exit(0)
    


selectedSteps = steps[:]

MNHSteptoLidarStep={}

if len(lidarList.keys())>0:
    sortedList = np.sort(lidarList.keys())
    for val in lidarList.keys():
        for stepV in steps[:]:
            if val <= stepV:
                MNHSteptoLidarStep[stepV] = val
                break;
    selectedSteps=np.sort(MNHSteptoLidarStep.keys())


domainID = 1
numOfDomains = 0;

# looking to regenerate domains files
if genDomainExtent != None:
    print "generating domain files assuming square subdomains and square domain"
    
    while(domainID > 0):
        zBinFileName = "%s.%d.%s.%s"%(inpattern,domainID,stepzgrid,gridKey)
        if appendedSteps != None :
            zBinFileName = "%s.%d.%s"%(inpattern,domainID,gridKey)
        if os.path.isfile(zBinFileName) : 
            domainID +=1
            numOfDomains += 1
        else:
            domainID = 0
    if(numOfDomains == 0):
        print "bad path to generate domain files ","%s.%d.%s.%s"%(inpattern,domainID,stepzgrid,gridKey)
        exit(1)
    
  #  readBinS(inpattern,domnum,stepV,vkey, appendedDict, bestshape=None):
    localD = readBinS(inpattern,1,stepzgrid,gridKey,appendedSteps,zkey=gridKey)
    shZ = np.shape(localD)    
    print shZ    
    nx =  shZ[0]-1
    ny = shZ[1]-1
    nz =  shZ[2]-1        
    nlines = ncols = np.sqrt(numOfDomains)
    rG = float(genDomainExtent[0])/((ncols*nx)+2)
    if len(genDomainExtent)>2:
        rG = float(genDomainExtent[2])
    print "resolution ",rG, numOfDomains, "domains ", nx, ny,"in size" 
    inipoint = (float(genDomainOrigin[0]),float(genDomainOrigin[1]))
    dloc =[0,0]
    
    for did in range(1,numOfDomains+1):
        localD = readBinS(inpattern,did,stepzgrid,gridKey,appendedSteps,zkey=gridKey)
        shZ = np.shape(localD)
        nx =  shZ[0]-1
        ny = shZ[1]-1    
        fDomFileName = "%s.%d.%d"%(inFFpattern, did,Allsteps[0])
        domstring = "FireDomain[sw=(%.0f,%.0f,0);ne=(%.0f,%.0f,0);t=%d]\n"%(dloc[0]+inipoint[0],dloc[1]+inipoint[1],dloc[0]+inipoint[0]+rG*(nx+2),dloc[1]+inipoint[1]+rG*(ny+2),Allsteps[0])
        dloc[0] = dloc[0] + rG*nx
        if dloc[0] > (float(genDomainExtent[0])-rG*nx):
            dloc[0] = 0
            dloc[1] = dloc[1] + rG*ny
        if not os.path.isfile(fDomFileName) :
            fDomFile = open(fDomFileName, 'w')
            fDomFile.write(domstring)
            fDomFile.close;
            print "writing  ", domstring, fDomFileName
        else :
            fDomFile = open(fDomFileName, 'r')
            domS = fDomFile.readline()
            if not (domS==domstring) :
                fDomFileName = "%sM.%d.%d"%(inFFpattern, did,Allsteps[0])    
                fDomFile = open(fDomFileName, 'w')
                fDomFile.write(domstring)
                print "overwritten %s to %s"%(domstring ,fDomFileName), 
            else:
                print " %s identical domain exist"%fDomFileName   
    exit(0)
     

while(domainID > 0):
    domFFFileName = "%s.%d.%d"%(inFFpattern, domainID,domFFBaseNumber)
   
    if os.path.isfile(domFFFileName) : 
        domainID +=1
        numOfDomains += 1
    else:
        domainID = 0
if(numOfDomains == 0):
    print "Found no domains ","%s.%d.%d"%(inFFpattern, domainID,domFFBaseNumber) 
    exit(1)
domains = np.zeros(shape=(4,numOfDomains),dtype=float)

print numOfDomains , " Domains"



for domainID in range(1,numOfDomains+1):
    fDomFile = open("%s.%d.%d"%(inFFpattern, domainID,domFFBaseNumber), 'r')
    domS = fDomFile.readline()
    sw = domS.split("sw=(")[1].split(")")[0].split(',');
    ne = domS.split("ne=(")[1].split(")")[0].split(',');    
    domains[0][domainID-1] = float(sw[0])
    domains[1][domainID-1] = float(sw[1])
    domains[2][domainID-1] = float(ne[0])
    domains[3][domainID-1] = float(ne[1])
    fDomFile.close()


if startStep > -1 and endStep > -1:
    selectedSteps=selectedSteps[startStep:endStep]

print len(selectedSteps), " time steps to be computed ", selectedSteps, " "

largCell = 0;
if numOfDomains > 0:
    largCell = (domains[2][0] - domains[0][1])/2 
    # resolution is exual to half of overlap area
dx = dy = largCell

#largSubDom = (domains[2][0]) - (domains[0][0]) - largCell*2
#lonSubDom = (domains[3][0])  - (domains[1][0]) - largCell*2 
        
largDom = (max(domains[2]))  - (min(domains[0]))   - largCell*2
lonDom = (max(domains[3]))  - (min(domains[1]))  - largCell*2


tx = int(largDom/largCell)
ty =  int(lonDom/largCell)
print "domain ",tx, ty, " resolution ", dx,dy, " in  domains (%d)"%numOfDomains

domainID = 1               

paralPiecesVtkStr=""
tz = 0
curStep = 0

        
        
gf = VtkGroup("%s/%sgroupfields"%(outPath,fprefix))

gs = 0
if vtpShapes3Dout :
    gs = VtkGroup("%s/%sgroupshapes"%(outPath,fprefix))
gp = VtkGroup("%s/%sgrouppoints"%(outPath,fprefix))


xd = {}
yd = {}
zd = {}

startd_D = {}
endd_D = {}

startd_F = {}
endd_F = {}
fd = {}
xf = {}
yf = {}
zf = {}
paralPiecesVtkStr=""


nxf = 0
nyf = 0

probesLocation=[]



zzd = None 
domainshapesInPoint ={}
dLocalShape = {}



 
for inddom in range(0,numOfDomains):
    #zBin = "%s.%d.%s.%s"%(inpattern,inddom+1,stepzgrid,gridKey)
      #  readBinS(inpattern,domnum,stepV,vkey, appendedDict, bestshape=None):
    local = readBinS(inpattern,inddom+1,stepzgrid,gridKey,appendedSteps,zkey=gridKey)
    shZ = np.shape(local)  
    nx =  shZ[0]-1
    ny = shZ[1]-1
    nz =  shZ[2]-1 
       
    domainshapesInPoint[inddom] = (nx,ny,nz)
    tz = nz
    
    if (zzd == None):
        zzd = np.zeros((tx+1, ty+1, tz+1), dtype=np.float32)
    #print inddom, "res is ", nx,ny,nz,domains[0][inddom+1], (domains[0][inddom+1] - domains[2][inddom])/shZ[1], (domains[1][inddom+1] - domains[3][inddom])/shZ[0]
    
    dLocalShape[inddom] =(int((domains[1][inddom] - domains[1][0])/dy) ,int((domains[3][inddom] - domains[1][0])/dy -1), int((domains[0][inddom] - domains[0][0])/dx) ,int( (domains[2][inddom] - domains[0][0])/dx -1 ))
    
    
   # print dLocalShape[inddom]
    
    for k in range(nz+1):
        for j in range(dLocalShape[inddom][0],dLocalShape[inddom][0]+ny+1):
            for i in range(dLocalShape[inddom][2],dLocalShape[inddom][2]+nx+1):
                zzd[i,j,k] = local[i-dLocalShape[inddom][2],j-dLocalShape[inddom][0],k]  


xxd = np.zeros((tx+1, ty+1, tz+1), dtype=np.float32)
yyd = np.zeros((tx+1, ty+1, tz+1), dtype=np.float32)
X = np.arange(domains[0][0]+1.5*dx,domains[0][0]+(tx+1)*dx+1.5*dx, dx, dtype=np.float32) 
Y = np.arange(domains[1][0]+1.5*dy,domains[1][0]+(ty+1)*dy+1.5*dy, dy, dtype=np.float32)

for i in range(tx+1):
    for j in range(ty+1):
        for k in range(tz+1):
            xxd[i,j,k] = X[i]  
            yyd[i,j,k] = Y[j]       


print "geometry OK"
 
paralPiecesVtkStr = ""
# selection of smoke points
if len(MNHSteptoLidarStep.keys())>0 :
    for myshotfname, myshot in lidarList.values():
        for iShot, shotLoc in enumerate(myshot.getShotLocationAsXYZ()):
            if myshot.getValues()[iShot] > 0.1:
                smokeLocList.append(getLocationAndCoeffsInGridPoints(shotLoc,(xxd[0,0,0],yyd[0,0,0]),tx,ty,tz, largCell, zzd))
    print len(smokeLocList), "sample smoke point selected"
for stepV in selectedSteps:
     
        
   
    paralPiecesVtkStr = ""
    fronts=[]
    ntotPoints = 0
    varmapAll = {}
    
    for keycount, vkey in enumerate(varsDataIn):
        varmapAll[vkey] = np.zeros((tx+1, ty+1, tz+1), dtype=np.float32)
 #   for vkey in inputcdfvars :
 #       varmapAll[vkey] = np.zeros((tx+1, ty+1, tz+1), dtype=np.float32)
 #   
 #   for vkey in inputcdfvars:
 #       varData = readNcField(inputcdfFile, vkey)
 #       for k in range(tz):
 #           for j in range(ty):
 #               for i in range(tx):
 #                   varmapAll[vkey][i,j,k] = varData[k,j,i]
    
    for indom in range(0,numOfDomains):
        if (indom%50 == 0) : 
            print  "processing sub-domains %d to %d at step %d"%(indom,indom+50, stepV)
        ffrontFileName = "%s.%d.%d"%(inFFpattern, indom+1,stepV)
        if os.path.isfile(ffrontFileName) :
            fDomFile = open(ffrontFileName, 'r')
            fDomFile.readline()
            fDomFile.readline()
            domainLines = fDomFile.read().split("\n")
            cfront = []
            
            for dline in domainLines:
                if "FireFront" in dline:
                    fronts.append(cfront)
                    cfront = []
                if "FireNode" in dline:
                    loc = dline.split("loc=(")[1].split(")")[0].split(',');
                    vel = dline.split("vel=(")[1].split(")")[0].split(',');    
                    velMod = np.sqrt((float(vel[0])*float(vel[0]))+(float(vel[1])*float(vel[1]))+(float(vel[2])*float(vel[2])))
                    ntotPoints = ntotPoints +1
                    if np.abs(velMod) > 0.00:    
                        cfront.append(((float(loc[0]),float(loc[1]),float(loc[2])),(float(vel[0]),float(vel[1]),float(vel[2]))))
                        
            if (len(cfront)> 0):
                fronts.append(cfront)
        
        
        nx,ny,nz =         domainshapesInPoint[indom] 
        for keycount, vkey in enumerate(varsDataIn):
            fname =  "%s.%d.%d.%s"%(inpattern,indom+1,stepV,vkey)
            varData = readBinS(inpattern,indom+1,stepV,vkey,appendedSteps, bestshape=(nx+2, ny+2)) 
            if varData != None:
                try:
                    for k in range(nz+1):
                        for j in range(dLocalShape[indom][0],dLocalShape[indom][0]+ny+1):
                            for i in range(dLocalShape[indom][2],dLocalShape[indom][2]+nx+1):
                                varmapAll[vkey][i,j,k] = varData[i-dLocalShape[indom][2],j-dLocalShape[indom][0],k] 
                except IndexError:
                    print "index problem with " + fname

                            
                            
            else:
                print "problem with " + fname
            
       
       
             
    
    outname = "%s/%s.full.%d"%(outPath,fprefix,stepV) 
    ptsAll = {}
    for key in scals["points"]:
        ptsAll[key] = varmapAll[key]
        
    for key in inputcdfvars:
        ptsAll[key] = varmapAll[key]
        
    vects = {}
    for key in Vects.keys():
        ptsAll[key] =(varmapAll[Vects[key][0]],varmapAll[Vects[key][1]],varmapAll[Vects[key][2]])
    
#//    ICI diagnosticas probes &  lidarshot
    if len(MNHSteptoLidarStep.keys())>0 :
        myshotfname, myshot = lidarList[MNHSteptoLidarStep[stepV]]
        print "Start lidarsim with ", stepV, " here ",  myshot.time
        #numpointsinshot = 300
        #myshot.setNumPointsWithSameOrigin(numpointsinshot)
        shotVals = np.zeros(myshot.getNumPoints())
        
        
        for iShot, shotLoc in enumerate(myshot.getShotLocationAsXYZ()):
            shotAndCoeffLoc = getLocationAndCoeffsInGridPoints(shotLoc,(xxd[0,0,0],yyd[0,0,0]),tx,ty,tz, largCell, zzd)
            shotVals[iShot] = getProbeValue(varmapAll["BRatio"],shotAndCoeffLoc)
        
        windHere = []
        
        for smokeLoc in smokeLocList:
            if getProbeValue(varmapAll["BRatio"],smokeLoc) > 0:
                windHere.append((getProbeValue(varmapAll["U"],smokeLoc),getProbeValue(varmapAll["V"],smokeLoc),getProbeValue(varmapAll["W"],smokeLoc)))
            
        avgU = avgV = avgW = avgMM = avgM = 0
        for wu,wv,ww in windHere:
            avgU += wu
            avgV += wv
            avgW += ww
            avgMM += np.sqrt(wu*wu+wv*wv+ww*ww)
        if len(windHere) >0:
            avgU = avgU / len(windHere)
            avgV = avgV / len(windHere)
            avgW = avgW / len(windHere)
            avgMM = avgMM / len(windHere)
        avgM = np.sqrt(avgU*avgU+avgV*avgV+avgW*avgW)
        probeFile.write("%f;%f;%f;%f;%f;%f\n"%(stepV,avgU,avgV,avgW,avgMM,avgM))
        
        myshot.setValues(shotVals)
        myshot.saveAs("%s%s"%(lidarOut,myshotfname))
        print "end lidarsim" 
        
     #   print " handling probes"
     #   print getLocationAndCoeffsInGridPoints((7950,2173,1480),(xxd[0,0,0],yyd[0,0,0]),tx,ty,tz, largCell, zzd)
     #   for sde in range(0,100,10):
     #       print getProbeValue(varmapAll["BRatio"],getLocationAndCoeffsInGridPoints((7950,2173,1480+sde),(xxd[0,0,0],yyd[0,0,0]),tx,ty,tz, largCell, zzd))
     #   print "finished probes"
    else:
        gridToVTK(outname, xxd, yyd, zzd,  cellData = None, pointData = ptsAll)  
    
        gf.addFile(filepath = "%s.vts"%outname, sim_time = stepV)
    
        paralVtkStr="<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"
        paralVtkStr+="<PolyData>\n"
        
        xp = np.zeros((ntotPoints))
        yp=np.zeros((ntotPoints))
        zp=np.zeros((ntotPoints))
        pvx=np.zeros((ntotPoints))
        pvy=np.zeros((ntotPoints))
        pvz=np.zeros((ntotPoints))
        datap=np.zeros((ntotPoints))
        ptcount = 0
        
        for front in fronts:
            npoints = 0
            npoints = len(front)*2-1
            if vtpShapes3Dout : 
                paralVtkStr+="<Piece NumberOfPoints=\"%d\" NumberOfPolys=\"%d\">\n"%(npoints,npoints-2)
                paralVtkStr+="     <Points>\n"
                paralVtkStr+="        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"
            ptplus = 1
            ptleft = None
            ptright = None
            
            for point,vel in front:
                ptleft = ptright
                ptright = point
                velMod = np.sqrt((vel[0]*vel[0])+(vel[1]*vel[1])+(vel[2]*vel[2]))
                if(ptcount < ntotPoints):
                    datap[ptcount] = velMod
                    pvx[ptcount] = vel[0]
                    pvy[ptcount] = vel[1]
                    pvz[ptcount] = vel[2]
                    xp[ptcount]= point[0]
                    yp[ptcount]=point[1]
                    zp[ptcount]=point[2]
                if vtpShapes3Dout :
                    if ptplus == 2:
                        paralVtkStr+="\t\t %f %f %f\n"% ((ptleft[0]+ptright[0])/2+3*vel[0],(ptleft[1]+ptright[1])/2+3*vel[1],(ptleft[2]+ptright[2])/2+velMod*5)
                        ptplus = 1
                    paralVtkStr+="\t\t %f %f %f\n"%(point[0],point[1],point[2]-1)
                    ptplus += 1
                ptcount += 1
            if vtpShapes3Dout :
                paralVtkStr+="\n        </DataArray>\n"
                paralVtkStr+="     </Points>\n"
                paralVtkStr+="  <Polys>\n"
                paralVtkStr+="        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
                for val in range(0,npoints-2):
                    paralVtkStr+=" %d %d %d\n"%(val,val+1,val+2)
                paralVtkStr+="\n        </DataArray>\n"
                paralVtkStr+="        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
                for val in range(1,npoints-1):
                    paralVtkStr+="  %d"%(val*3)
                paralVtkStr+="\n        </DataArray>\n"
                paralVtkStr+="     </Polys>\n"
                paralVtkStr+="  </Piece>\n"
        
        paralVtkStr+="</PolyData>\n"
        paralVtkStr+="</VTKFile>\n"
        if vtpShapes3Dout :
            foutGname = "%s/%s%d%s"%(outPath,fprefix,stepV,".vtp");
            pvtkFile = open(foutGname,"w")
            pvtkFile.write(paralVtkStr)
            pvtkFile.close()
            gs.addFile(filepath = foutGname, sim_time = stepV)
            
        fpoutGname = "%s/%sPoints%d"%(outPath,fprefix,stepV);
        if(xp.size > 0):
            pointsToVTK(fpoutGname, xp, yp, zp, {"vmod":datap,"velocity":(pvx,pvy,pvz)})
            gp.addFile(filepath = "%s.vtu"%fpoutGname, sim_time = stepV)
            
if len(MNHSteptoLidarStep.keys())==0 :
    gf.save()
    gp.save()
else:
    probeFile.close()
if vtpShapes3Dout :
    gs.save()
print "saved groups , to make image :  paraview, then ffmpeg -r 5 -i \"filemame.%04d.jpg\" -qscale 2 -s 'hd720'  m.mp4"
