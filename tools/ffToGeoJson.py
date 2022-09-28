import os
import sys
from pyproj import Proj, transform
import json
from operator import isNumberType

#http://dev.firecaster.valabre.fr/forefireAPI.php?command=ignition&longitude=42.348570491488999&date_simulation=2020-01-04T09%3A37%3A09Z&date_eclosion=2020-02-04T09%3A37%3A09Z&latitude=9.017890207687699&pathdate=2020-01-07_T10_47_15Z
    
#http://dev.firecaster.valabre.fr/forefireAPI.php?command=sequence&duree_h=4&vent_direction_degre=20&vent_vitesse_kmh=50&temperature=10&hygrometrie=10&pathdate=2020-01-07_T10_47_15Z
    
    
#     $duree_h=$_REQUEST['duree_h'];
# 11 $vent_direction_degre=$_REQUEST['vent_direction_degre'];
# 12 $vent_vitesse_kmh=$_REQUEST['vent_vitesse_kmh'];
# 13 $hygrometrie=$_REQUEST['hygrometrie'];
# 14 $temperature=$_REQUEST['temperature'];
# 15 $pathdate=$_REQUEST['pathdate'];

printff = 'FireDomain[sw=(479337,4.62206e+06,0);ne=(520762,4.67757e+06,0);t=32406]\n\tFireFront[id=2;domain=0;t=32400]\n\t\tFireNode[domain=0;id=4;fdepth=100;kappa=0.0106066;loc=(511897,4.6547e+06,0);vel=(1.50339e-13,0.114793,0);t=32400;state=moving;frontId=2]\n\t\tFireNode[domain=0;id=6;fdepth=100;kappa=0.01;loc=(511977,4.65454e+06,0);vel=(0.106766,-0.0659849,0);t=32400;state=moving;frontId=2]\n\t\tFireNode[domain=0;id=8;fdepth=100;kappa=0.0124568;loc=(511817,4.65454e+06,0);vel=(-0.137972,-0.0852716,0);t=32400;state=moving;frontId=2]\n\t\tFireNode[domain=0;id=10;fdepth=100;kappa=-1.31783e-13;loc=(511843,4.65459e+06,0);vel=(-0.113095,0.0565476,0);t=32400;state=moving;frontId=2]\n\t\tFireNode[domain=0;id=12;fdepth=100;kappa=1.31783e-13;loc=(511870,4.65465e+06,0);vel=(-0.0864376,0.0432188,0);t=32400;state=moving;frontId=2]\n'
print printff
    
class ffToGeoJson:
    
    def __init__(self, finenamein,destProj='epsg:4326'):
        self.fname = finenamein
        
        simpleFname = finenamein.split("/")[-1]
        self.frontCount=int(simpleFname.split("-")[0])
        frontProj=int(simpleFname[simpleFname.find("Z")+2:simpleFname.find(".ff")].split("-")[1])
        epsgstring=('epsg:%d')%frontProj
        
        self.projSRC=Proj(epsgstring)
        self.projDest=Proj(destProj) 
        self.area=999
        laDate=simpleFname[simpleFname.find("-")+1:simpleFname.find("Z")].split("T")
         
        self.frontDate = "%s %s UTC"%(laDate[0],laDate[1].replace('-', ':'))
        print self.frontDate
       
    def parse(self, ):
        
        def proj(point,inProj = None, outProj = None):
            if inProj is None or outProj is None:
                return point
            return transform(inProj,outProj,point[0],point[1])
        
        def isPoint(element):
            if len(element) == 2:
                if isinstance(element[0],float) and isinstance(element[1],float) :
                    return True
            return False
        
        def getLocationFromLine(line,pattern="loc=("):
            llv = line.split(pattern)
            if len(llv) < 2: 
                return None
            llr = llv[1].split(",");
            if len(llr) < 3: 
                return None
            return (float(llr[0]),float(llr[1]))
    
        def printToPolygons(linePrinted, inProj = None, outProj = None, level=1):
            if level > 8:
                 return 

            fronts = linePrinted.split("\n%sFireFront"%('\t'*level))

            pointsMap = []
            if len(fronts)>0:
                nodes = fronts[0].split("FireNode")
                if len(nodes) > 1:
                    for node in nodes[1:]:
                        pointsMap.append(proj(getLocationFromLine(node),inProj, outProj))
                    pointsMap.append(proj(getLocationFromLine(nodes[1]),inProj, outProj))

                for subline in fronts[1:]:
                    pointsMap.append(printToPolygons(subline,inProj, outProj, level+1))

            return pointsMap
        
        geoJSonStyleData = []
        f = open(self.fname, 'r')

        
        for front in printToPolygons(f.read(),self.projSRC, self.projDest):
            mainfront = []
            otherfronts = []
            for element in front:
                if isPoint(element): 
                    mainfront.append(element)
                else:
                    otherfronts.append(element)

            geoJSonStyleData.append(mainfront[::-1])
            for element in otherfronts:
                geoJSonStyleData.append(element[::-1])


        data = {
        "type": "FeatureCollection",
        "name": "Front Nr:%d"%self.frontCount,
        "features": [
        { "type": "Feature", "properties": { "area": self.area, "date": self.frontDate }, 

        "geometry": { "type": "MultiPolygon", "coordinates": [geoJSonStyleData]
        } 
        }
        ]
        }
        return data

ffGSon = ffToGeoJson(sys.argv[1])
print json.dumps(ffGSon.parse())
