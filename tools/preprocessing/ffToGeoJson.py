 
from pyproj import Proj, transform
import numpy as np 
from shapely.geometry import Polygon 
from datetime import timedelta,datetime
 

import re
from PIL import Image, ImageDraw, ImageFont 
import matplotlib.cm as cm
    
import glob
import pandas as pd
import xarray as xr

 
from fastkml import kml, styles
from shapely.geometry import shape 
import math

from .geo_to_kml import *

class ffData:
    def isPoint(element):
  
        if len(element) == 3:
            if isinstance(element[0],float) and isinstance(element[1],float) and isinstance(element[2],float) :
                return True
        return False
    
    def __init__(self, filenamein=None,mode="batchWeb",destProj='epsg:4326',baseDate=datetime(2020, 1, 1, 0, 0, 0, 0),wsen=None,lbrt=None):
        
        self.fname = filenamein
 
        self.inXYpolygon = []
        self.inLLpolygon = []
        
        if mode == "batchWeb" :
            simpleFname = filenamein.split("/")[-1]
            self.frontCount=int(simpleFname.split("-")[0])
            frontProj=int(simpleFname[simpleFname.find("Z")+2:simpleFname.find(".ff")].split("-")[1])
            epsgstring=('epsg:%d')%frontProj
            self.mode = "batchWeb"
            self.projSRC=Proj(epsgstring)
            self.projDest=Proj(destProj) 
            self.area=0
            laDate=simpleFname[simpleFname.find("-")+1:simpleFname.find("Z")].split("T")
             
            self.frontDate = "%s %s UTC"%(laDate[0],laDate[1].replace('-', ':'))
     
        if mode == "mnhPGD" :
            
            self.frontCount=1
            self.projSRC=None
            self.projDest=None
            self.wsen = wsen
            self.lbrt = lbrt
            self.area=0
            localDate = baseDate 
            if filenamein is not None :
                localDate = baseDate + timedelta(seconds=int(filenamein.split(".")[-1]))
            self.mode="mnhPGD"
            self.frontDate = localDate.isoformat()
            

            
          

            
        if filenamein is not None :
            self.parse()
        self.atime = None

    def lalo2xy(self,lalo):
        if self.projSRC is None or self.projDest is None:
            wsen = self.wsen
            lbrt = self.lbrt
            if wsen is None or lbrt is None:
                return lalo
            x = lbrt[0]+((lalo[1] - wsen[0])/(wsen[2] - wsen[0]))*(lbrt[2]-lbrt[0])
            y = lbrt[1]+((lalo[0] - wsen[1])/(wsen[3] - wsen[1]))*(lbrt[3]-lbrt[1])
            return (x,y,lalo[2])
        return transform(self.projDest,self.projSRC,point[0],point[1])
    
    def xy2lalo(self,xy):
        if self.projSRC is None or self.projDest is None:
            wsen = self.wsen
            lbrt = self.lbrt
            if wsen is None or lbrt is None:
                return xy
            lo = wsen[0]+(float(xy[0] - lbrt[0])/float(lbrt[2] - lbrt[0]))*(wsen[2]-wsen[0])
            la = wsen[1]+(float(xy[1] - lbrt[1])/float(lbrt[3] - lbrt[1]))*(wsen[3]-wsen[1])
       
            return (la,lo,xy[2])       
        return transform(self.projSRC,self.projDest,point[0],point[1])
    
    def setATimeMatrix(self,at):
        self.atime = at
        
    def getATimeMatrix(self):
        return self.atime     

    def getROSMatrix(self):
        if self.atime is  None:
            return None
        
        l,b,r,t = self.lbrt
        norm_data = np.copy(self.atime )
         
        #norm_data = np.where(self.atime < -9990, np.nan, self.atime)
        norm_data[norm_data < 0] =  np.nan
        
        gradient_y, gradient_x = np.gradient(norm_data, 1)
        dspeed = np.sqrt(gradient_x**2 + gradient_y**2)
        
        BMapresolution = float( (r-l) / norm_data.shape[0] )
        print(BMapresolution, dspeed.shape, norm_data.shape)
        speed = BMapresolution / dspeed
        

        return speed
    

    
    
    
           
    def getLLPolygon(self, ):
        if len(self.inXYpolygon) > len(self.inLLpolygon):
            for front in self.inXYpolygon:
                fa = []
                for element in front:
                    if ffData.isPoint(element):
                        fa.append(self.xy2lalo(element))
                    else:
                        fb = [] 
                        for el2 in element:
                            if ffData.isPoint(el2):
                                fb.append(self.xy2lalo(el2))
                            else:
                                print("WARNING SUB-SUB-SUB front not parsed")
                        fa.append(fb)
                self.inLLpolygon.append(fa)
                    
        return self.inLLpolygon
    
    def getCartesianPolygon(self, ):
        
        return self.inXYpolygon
    
    def toKML(self, ):
        return to_kml(self.toGeoJson())
 
    
    def toVTK(self, ):
        return 0
    
    def toGeoJson(self, ):
        llPoly = self.getLLPolygon()
        
        geoJSonStyleData = []
        
        for front in llPoly:
            mainfront = []
            otherfronts = []
            for element in front:
                if ffData.isPoint(element): 
                    mainfront.append((element[1],element[0],element[2]))
                else:
                    otherfronts.append(element)
 
            geoJSonStyleData.append(mainfront[::-1])
            
            for nfront in otherfronts:
                holesFronts  = []
                for element in nfront:
                    if ffData.isPoint(element): 
                        holesFronts.append((element[1],element[0],element[2]))
                    else:
                        print("Error Json Parsing fronts, too many holes in holes") 
                geoJSonStyleData.append(holesFronts[::-1])
                
   
   
        geometryInfo={}   
        geometryInfo["numberOfPolygons"]=  len(geoJSonStyleData)
 
 
            
        self.area = self.metadata["totalArea"]
        data = {
        "type": "FeatureCollection",
        "date": self.frontDate, 
        "details":self.metadata ,
        "features": [
        { "type": "Feature", "properties": { "numberOfPolygons": len(geoJSonStyleData) }, 

        "geometry": { "type": "MultiPolygon", "coordinates": [geoJSonStyleData]
        } 
        }
        ]
        }
        return data
        return 0
    
    def toFFPrint(self, ):
        return 0        
    
    def getDomainInitString(self, secondsSinceMidnight=0):
        l,b,r,t = self.lbrt
        
        return "FireDomain[sw=(%d,%f,0);ne=(%f,%f,0);t=%d]"%(l,b,r,t,secondsSinceMidnight)
      
    def parse(self, ):            

        
        def getLocationFromLine(line):
            match = re.search(r"loc=\(([^,]+),([^,]+),([^)]+)\)", line)
            if match:
                x, y, z = map(float, match.groups())
                return (x, y, z)
            else:
                return None
        
            return (float(llr[0]),float(llr[1]),float(llr[2]))
        
        def getMarkPropFromLine(line):
      
            Mpos = getLocationFromLine(line)
           # Mvel = (1,1)#getLocationFromLine(line,pattern="vel=(")
         #   print("here ", Mvel)
            return Mpos#, Mvel #+ Mvel
    
        def printToPolygons(linePrinted, level=1):
            if level > 8:
                 return 

            fronts = linePrinted.split("\n%sFireFront"%('\t'*level))

            pointsMap = []
            if len(fronts)>0:
                nodes = fronts[0].split("FireNode")
                if len(nodes) > 1:
                    for node in nodes[1:]:
                        #ptl, sptl = getMarkPropFromLine(node)
                        pointsMap.append(getMarkPropFromLine(node))
                    #ptl, sptl =  
                    pointsMap.append(getMarkPropFromLine(nodes[1]))

                for subline in fronts[1:]:
                    pointsMap.append(printToPolygons(subline,level+1))

            return pointsMap
        
        self.inXYpolygon = []
        f = open(self.fname, 'r')
        rawdata = f.read()
        firstLine = rawdata.partition('\n')[0]
        print(firstLine)
        llv = firstLine.split("sw=(")
        llr = llv[1].split("ne=(");
        self.lbrt = [  float( llr[0].split(",")[0]), float(llr[0].split(",")[1]), float(llr[1].split(",")[0]), float(llr[1].split(",")[1]) ]
        
        self.inXYpolygon = printToPolygons(rawdata)
        
        f.close()
        
        metaHelper = []
        for front in self.inXYpolygon :
            mainfront = []
            otherfronts = []
            for element in front:
                if ffData.isPoint(element): 
                    mainfront.append(element)
                else:
                    otherfronts.append(element)
 
            metaHelper.append(mainfront[::-1])
            for element in otherfronts:
                metaHelper.append(element[::-1])
                            
        self.metadata = {}   
        self.metadata["numberOfPolygons"]=  len(metaHelper)
        totalArea = 0
        totalLength = 0
        
        for na,a in enumerate(metaHelper):
             pgon = Polygon(a)
             self.metadata["P%d"%na]={"area":pgon.area,"perimeter":pgon.length, "directWinding":pgon.exterior.is_ccw}
             totalLength = totalLength+pgon.length
             totalArea = totalArea+pgon.area
        
        self.metadata["totalLength"]=  totalLength
        self.metadata["totalArea"]=  totalArea



def get_WSEN_LBRT_ZS_From_Pgd(pgd_path):
    pgd = xr.load_dataset(pgd_path)
    dlon = float(pgd.longitude[0,1]-pgd.longitude[0,0])
    dlat = float(pgd.latitude[1,0]-pgd.latitude[0,0])
    dx = float(pgd.XHAT[1]-pgd.XHAT[0])
    dy = float(pgd.YHAT[1]-pgd.YHAT[0])
 
    return (
        (
        float(pgd.longitude.min())-dlon/2,
        float(pgd.latitude.min())-dlat/2,
        float(pgd.longitude.max())+dlon/2,
        float(pgd.latitude.max())+dlon/2
        ),
        (
        float(pgd.XHAT.min())-dx/2,
        float(pgd.YHAT.min())-dy/2,
        float(pgd.XHAT.max())+dx/2,
        float(pgd.YHAT.max())+dy/2,
        ),
        pgd.ZS
    )

 
def areaBurnt():    
    contours1  = glob.glob('/Users/filippi_j/data/2022/pedrogao/OutputsFullCoupled/output.0.*')
    contours2  = glob.glob('/Users/filippi_j/data/2022/pedrogao/OutputsNoCoupled/output.0.*')
    

        
    for contour in contours1:
        ffGSon= ffToGeoJson(contour,mode = "mnhPGD", baseDate=baseDate)
        d = ffGSon.parse()
        area = d["details"]["P0"]["area"] 
        date = datetime.fromisoformat(d["date"])
        v = [date , area]
        lCnt.append( v )
    
        
    
    dfc = pd.DataFrame(lCnt, columns=['date','HaCoupled'])
    dfc["HaCoupled"] = dfc["HaCoupled"]/10000
    
    dfc = dfc.set_index("date")
    dfc.sort_values(by=['date'], inplace=True)
    
    lCnt = []
    for contour in contours2:
        ffGSon= ffToGeoJson(contour,mode = "mnhPGD",baseDate=baseDate)
        d = ffGSon.parse()
        area = d["area"] 
        date = datetime.fromisoformat(d["date"])
        v = [date , area]
        lCnt.append( v )
            
    
    dnc = pd.DataFrame(lCnt, columns=['date','HaNoCoupled'])
    dnc["HaNoCoupled"] = dnc["HaNoCoupled"]/10000
    
    dnc = dnc.set_index("date")
    dnc.sort_values(by=['date'], inplace=True)
    
     
    
     
    totalBurnArea = np.count_nonzero(dBmap.arrival_time_of_front > 000000)
    mSize = (dBmap.domain.Lx/dBmap.arrival_time_of_front.shape[0]) * (dBmap.domain.Ly/dBmap.arrival_time_of_front.shape[1])
    lCnt = []
    for lDate in dnc.index:
        yetBurnt =( totalBurnArea - np.count_nonzero(dBmap.arrival_time_of_front > (lDate-baseDate).seconds))*mSize
        lCnt.append( [lDate, yetBurnt ])
     
    
    doc = pd.DataFrame(lCnt, columns=['date','HaObserved'])
    doc["HaObserved"] = doc["HaObserved"]/10000
    
    doc = doc.set_index("date")
    doc.sort_values(by=['date'], inplace=True)
    
    dac = pd.concat([dfc,dnc,doc], join='inner', axis=1)
    dac.plot()       


def geojson_to_placemark(geojson_obj, timestamp):
    # Convertir la géométrie GeoJSON en géométrie Shapely
    geom = shape(geojson_obj['geometry'])

    # Créer un Placemark KML
    pm = kml.Placemark()
    pm.geometry = geom
    pm.timestamp = timestamp

    # Créer et définir le style (contour orange de 2 points, pas de remplissage)
    style = styles.Style()
    style.linestyle.color = styles.Color.rgb(255, 165, 0)  # Orange
    style.linestyle.width = 2  # 2 points de large
    style.polystyle.fill = 0  # Pas de remplissage
    pm.style = style

    return pm

def geojsons_to_kml(geojson_objs):
    # Créer un KML
    k = kml.KML()

    # Créer un document pour contenir les Placemarks
    doc = kml.Document()
    k.append(doc)

    # Convertir chaque objet GeoJSON en Placemark et l'ajouter au document
    for gj in geojson_objs:
        timestamp_str = gj['properties']['timestamp']  # Assumer une propriété 'timestamp'
        timestamp = datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))  # Convertir la chaîne en datetime
        pm = geojson_to_placemark(gj, timestamp)
        doc.append(pm)

    return k.to_string(prettyprint=True)

def otherJsonToKml():
    # Utilisation :
    geojson_objs = [
        # Vos objets GeoJSON vont ici
    ]
    kml_str = geojsons_to_kml(geojson_objs)
    with open('output.kml', 'w') as f:
        f.write(kml_str)

def create_kml(west, south, east, north, Name, pngfile, outkml_path, pngcbarfile=None):
    kml_template = '''<?xml version="1.0" encoding="UTF-8"?>
    <kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">
        <Document id="1">'''
    
    if pngcbarfile:
        kml_template += '''
        <ScreenOverlay>
            <name>Légende</name>
            <Icon>
                <href>{pngcbarfile}</href>
            </Icon>
            <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>
            <screenXY x="0.05" y="0.05" xunits="fraction" yunits="fraction"/>
            <rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
            <size x="50" y="200" xunits="pixels" yunits="pixels"/>
        </ScreenOverlay>'''
    
    kml_template += '''
            <GroundOverlay id="2">
                <name>{Name}</name>
                <Icon id="3">
                    <href>{pngfile}</href>
                </Icon>
                <LatLonBox>
                    <north>{north}</north>
                    <south>{south}</south>
                    <east>{east}</east>
                    <west>{west}</west>
                </LatLonBox>
            </GroundOverlay>
        </Document>
    </kml>'''

    kml_content = kml_template.format(Name=Name, pngfile=pngfile, north=north, south=south, east=east, west=west, pngcbarfile=pngcbarfile)
    
    with open(outkml_path, "w") as f:
        f.write(kml_content)



def draw_arrows_on_image(image, vx, vy, scale=10):
    draw = ImageDraw.Draw(image)
    nj, ni = vx.shape
    dx_space = image.width / ni
    dy_space = image.height / nj
    
    for j in range(nj):
        for i in range(ni):
            position = (i * dx_space + dx_space/2, j * dy_space + dy_space/2)
            dx = vx[j, i] * scale
            dy = vy[j, i] * scale
            arrow_start = position
            arrow_end = (position[0] + dx, position[1] - dy)
            draw.line((arrow_start, arrow_end), fill='black', width=2)

            angle = math.atan2(dy, dx)
            arrow_head_length = 5 * (scale/10)
            arrow_point1 = (arrow_end[0] - arrow_head_length * math.cos(angle - math.pi/4),
                            arrow_end[1] + arrow_head_length * math.sin(angle - math.pi/4))
            arrow_point2 = (arrow_end[0] - arrow_head_length * math.cos(angle + math.pi/4),
                            arrow_end[1] + arrow_head_length * math.sin(angle + math.pi/4))

            draw.line((arrow_end, arrow_point1), fill='black', width=2)
            draw.line((arrow_end, arrow_point2), fill='black', width=2)




def arrayToPng(data, output_filename, output_cbar_png=None, cmap_str = "viridis",vx=None, vy=None, scale=10):


    # Normaliser les données entre 0 et 255
 
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
     
    
    imsize=(data.shape[0]*scale, data.shape[1]*scale)
    revdata = np.flipud(data)
    
    norm_data = ((revdata  - vmin) / (vmax - vmin))
    
    cmap = cm.get_cmap(cmap_str)
    rgba_data = (cmap(norm_data) * 255).astype(np.uint8)
    
    # Rendre les valeurs NaN transparentes
    rgba_data[np.isnan(revdata)] = (0, 0, 0, 0)
    
    # Créer et sauvegarder l'image
    image = Image.fromarray(rgba_data, "RGBA")
    image = image.resize(imsize)
    
    # Si les vecteurs vx et vy sont fournis, dessiner les flèches
    if vx is not None and vy is not None:
        draw_arrows_on_image(image, vx, vy, scale)
    # Enregistrer l'image
    image.save(output_filename)
    
    if output_cbar_png is not None:
        # Créer un gradient 1D de la colormap
        height, width = 200, 50
        gradient = np.linspace(1, 0, height).reshape(-1, 1)  # 1 at the top to 0 at the bottom
        gradient = np.hstack([gradient]*width)  # Repeat for each column
        
        # Apply the colormap to the gradient
        rgba_gradient = (cmap(gradient) * 255).astype(np.uint8)
        
        # Create an extended version of the gradient for ticks and labels
        #extra_width = 50  # 50 extra pixels for ticks and labels
        #extended_rgba_gradient = np.hstack([rgba_gradient, rgba_gradient[:, :extra_width]])
        
        # Create a PIL image from the extended array
        colorbar_image = Image.fromarray(rgba_gradient, "RGBA")
        
        # Prepare the object for drawing and the font
        draw = ImageDraw.Draw(colorbar_image)
        default_font = ImageFont.load_default()
        
        tick_positions = [0, height // 4, height // 2, 3 * height // 4, height - 1]
        tick_labels = [
            str(vmax), 
            str(vmin + 0.75 * (vmax - vmin)), 
            str((vmax + vmin) / 2), 
            str(vmin + 0.25 * (vmax - vmin)), 
            str(vmin)
        ]
        
        # Draw the ticks and labels
        for pos, label in zip(tick_positions, tick_labels):
            font = default_font
            fill = "black"
            vertical_adjustment = -5
            
            # Special formatting for top and bottom labels
            if pos == 0:
                fill = "white"
                vertical_adjustment = 2
            elif pos == height - 1:
                fill = "white"
                vertical_adjustment = -15
            
            # Draw a line for the tick
            draw.line([(5, pos), (5 + 10, pos)], fill="black", width=1)
            # Draw the text for the label
            draw.text((5 + 15, pos + vertical_adjustment), label, fill=fill, font=font)

        # Sauvegarder l'image modifiée
        colorbar_image.save(output_cbar_png)

 
# def bmap2kml(ds, wsen ,fnameKMLOUT,method="pil",maxSpeed=0.1):
  
#     import matplotlib.pyplot as plt
#     import matplotlib.colors as colors
#     import simplekml
#     import matplotlib.cm as cm
#     from PIL import Image, ImageOps
#     lon_min,  lat_min, lon_max, lat_max = wsen
#     # Créez une figure
#     fnameKMLPNGOUT="%s.png"%fnameKMLOUT
   
#     #minRealVal = float(ds['arrival_time_of_front'].where(ds['arrival_time_of_front'] > 0).min())
#     #max_val = float(ds['arrival_time_of_front'].max())
#     #norm_data = (ds['arrival_time_of_front'] -minRealVal)/  (max_val - minRealVal)
#     norm_data = ds.arrival_time_of_front.values
#     norm_data[norm_data < 0] =  np.nan
#     gradient_y, gradient_x = np.gradient(norm_data, 1)
#     BMapresolution = float(ds.domain.Lx / len(ds.DIMX))
#     speed = 100000 * (BMapresolution / np.sqrt(gradient_x**2 + gradient_y**2) )
#     speed[speed > 30000] = 30000
 
#     speed[np.isnan(speed)] = -1


#     array_with_nan = np.where(np.isnan(speed), -1, speed)
    
#     # Normaliser la matrice pour que les valeurs soient entre 0 et 1, à l'exception de -1 qui reste pour les NaN
#     min_val = np.nanmin(speed)
#     max_val = np.nanmax(speed)
#     normalized_array = np.where(array_with_nan == -1, -1, (array_with_nan - min_val) / (max_val - min_val))
    
#     # Appliquer la palette de couleurs "jet"
#     colored_array = (cm.jet(normalized_array) * 255).astype(np.uint8)
    
#     # Remettre les pixels correspondant aux valeurs NaN à transparent
#     colored_array[normalized_array == -1] = [0, 0, 0, 0]
    
#     # Créer une image PIL à partir de la matrice RGBA
#     pil_image = Image.fromarray(colored_array, 'RGBA')
     
 
#     pil_image = ImageOps.flip(pil_image)
    
#     # Sauvegardez l'image en PNG
#     pil_image.save(fnameKMLPNGOUT)

        
        
       
 
    # Supposons que 'lon_min', 'lat_min', 'lon_max', et 'lat_max' sont les coordonnées des coins de votre DataSet

    
 
    # kml = simplekml.Kml()

    # ground = kml.newgroundoverlay(name='GroundOverlay') 
    # ground.icon.href = fnameKMLPNGOUT  # Utilisez le fichier PNG enregistré
    # ground.latlonbox.north = lat_max
    # ground.latlonbox.south = lat_min
    # ground.latlonbox.east = lon_max
    # ground.latlonbox.west = lon_min
 
    
    # # Sauvegardez le kml
    # kml.save(fnameKMLOUT)



def genKMLFiles(PGDFILE, BMAPFILE, FFINPUTPATTERN, BMAPKMLOUT,frontsKMLOUT, everyNFronts=1,change_color_every=6):

    dBmap = xr.open_dataset(BMAPFILE)
    baseDate= datetime(int(dBmap.domain.refYear), 1, 1, 0, 0, 0, 0)
    baseDate= baseDate+timedelta(days=int(dBmap.domain.refDay))
    
    wsen, lbrt, ZS = get_WSEN_LBRT_ZS_From_Pgd(PGDFILE)
    west, south, east, north = wsen
     
    
    contours1  = glob.glob(FFINPUTPATTERN)
        
    lCnt=[]
    import os
    selectionSorted =  sorted(contours1)#[::6]    
    for contour in selectionSorted[::everyNFronts]:#[50:261]:
        f = ffData(contour,mode = "mnhPGD", baseDate=baseDate,wsen=wsen)
        lCnt.append((f.toGeoJson(),f.frontDate))
    with open(frontsKMLOUT, "w") as text_file:
        text_file.write(to_timed_kml(lCnt, croll=change_color_every))
    
    fdat = ffData(mode = "mnhPGD", wsen=wsen, lbrt = lbrt)
    
   
    fdat.setATimeMatrix(dBmap.arrival_time_of_front.values)
    ros = fdat.getROSMatrix()

    arrayToPng(ros, BMAPKMLOUT+"ROS.png",cmap_str="hsv", output_cbar_png=BMAPKMLOUT+"ROSCBAR.png")
    
    create_kml(west, south, east, north,"ROS", BMAPKMLOUT+"ROS.png",  BMAPKMLOUT, pngcbarfile=BMAPKMLOUT+"ROSCBAR.png")
    
    #bmap2kml(dBmap,wsen,BMAPKMLOUT+"2.kml")

def ffFromPgd(PGDFILE,domainDate=None,ignitions = None, FuelTest = -1):
    wsen, lbrt, ZS = get_WSEN_LBRT_ZS_From_Pgd(PGDFILE)
    fdat = ffData(mode = "mnhPGD", wsen=wsen, lbrt = lbrt)
    secsDomain = 0
    dom =""
    if domainDate is not None :
        secsDomain =  (domainDate - domainDate.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
        #dom+=f"setParameters[year={domainDate.year};month={domainDate.month};day={domainDate.day}]\n"
    dom+= fdat.getDomainInitString(secsDomain)
    
    if ignitions is not None :   
        for ig in ignitions:
            igT = datetime.fromisoformat(ig["when"])
            igSecs =  (igT - igT.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
            ix,iy,iz = fdat.lalo2xy((ig["latitude"],ig["longitude"],0))
            dom+=f"\nstartFire[loc=({ix},{iy},0.);t={igSecs}]"
            
    if FuelTest > 0:
        l,b,r,t = lbrt
        ix = np.linspace(l, r, FuelTest+2)
 
        iy = np.ones(ix.shape) * (b + (t - b) / 4)
        
        iz = np.zeros(ix.shape)
        
        # Si vous souhaitez les coordonnées comme des paires (ix[i], iy[i])
        points = np.array(list(zip(ix, iy, iz)))
        for ix, iy, iz in points[1:-1]:
            dom+=f"\nstartFire[loc=({ix},{iy},0.);t={secsDomain}]"
        
    return dom
   

    


# #pedrogao 
# BMAPFILE = "/Users/filippi_j/data/2022/pedrogao/REF_REPORT_BMAP.nc"
# PGDFILE = "/Users/filippi_j/data/2022/pedrogao/PGD_D80mA.nested.nc"
# FFINPUTPATTERN='/Users/filippi_j/data/2022/pedrogao/OutputsNoCoupled/output.0.*'

# #prunelli
# BMAPFILE = "//Users/filippi_j/Volumes/fcouto/KTEST_PEDROGAO/nested3FirePattern/006_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KTEST_PEDROGAO/nested3FirePattern/006_runff/PGD_D80mA.nested.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/fcouto/KTEST_PEDROGAO/nested3FirePattern/006_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/prunelli/4n/"
# BMAPKMLOUT = "%s/speedNC.kml"%outDir
# #BMAPKMLOUT = "%s/speedNC.kml"%outDir
# frontsKMLOUT = "%s/fronts.kml"%outDir
# #genKMLFiles(PGDFILE, BMAPFILE, FFINPUTPATTERN, BMAPKMLOUT,frontsKMLOUT, everyNFronts=30, change_color_every=30)
# #x=(182000,88496,0)
# #y = ffFromPgd(PGDFILE).xy2lalo(x)
# #y=(42.0037,9.3418,0)
# #z = ffFromPgd(PGDFILE).lalo2xy(y)

# #print(y,z,x)



# #prunelli local
# BMAPFILE = "/Users/filippi_j/data/2023/prunelli/frontData/coupled.nc"
# PGDFILE = "/Users/filippi_j/data/2023/prunelli/PGD_D80mA.nc"
# FFINPUTPATTERN= "/Users/filippi_j/data/2023/prunelli/frontData/cpl/output.0.*"



# #Villa1nest 
# BMAPFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_1nest/003_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_1nest/001_pgd/PGD_D120mA.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_1nest/003_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/baseFeux/villa1n/"

 

# #Villa2nest 
# BMAPFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_2nest/005_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_2nest/001_pgd/PGD_D160mA.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_2nest/005_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/baseFeux/villa2n/"

# #Villa3nest 
# BMAPFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_3nest/006_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_3nest/001_pgd/PGD_D80mA.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST/KTEST_3nest/006_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/baseFeux/villa3n/"

# #Oliveira2nest 
# BMAPFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_2nest/005_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_2nest/001_pgd/PGD_D160mA.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_2nest/005_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/baseFeux/KTEST_Oliveira/"

# #Oliveira3nest 
# BMAPFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_3nest/006_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_3nest/001_pgd/PGD_D80mA.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_3nest/006_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/baseFeux/KTEST_Oliveira3n/"



    
# #Pigna
# BMAPFILE = "//Users/filippi_j/Volumes/orsu/firecaster/2023/nest150Ref/006_runff/ForeFire/Outputs/ForeFire.0.nc"
# PGDFILE = "/Users/filippi_j/Volumes/orsu/firecaster/2023/nest150Ref/001_pgd/PGD_D80mA.nested.nc"
# FFINPUTPATTERN='/Users/filippi_j/Volumes/orsu/firecaster/2023/nest150Ref/006_runff/ForeFire/Outputs/output.0.*'
# outDir="/Users/filippi_j/data/2023/corbara20230727/"
# BMAPKMLOUT = "%s/speedNC.kml"%outDir

# frontsKMLOUT = "%s/OutputsNonCoupled.kml"%outDir


# #genKMLFiles(PGDFILE, BMAPFILE, FFINPUTPATTERN, BMAPKMLOUT,frontsKMLOUT,everyNFronts=10)

# PGDFILE = "/Users/filippi_j/data/2023/corbara20230727/mnhCase/001_pgd/PGD_D80mA.nested.nc"
# PGDFILE = "/Users/filippi_j/Volumes/fcouto/KDATABASE/KTEST_Oliveira/KTEST_2nest/001_pgd/PGD_D160mA.nested.nc"

# #x = ffFromPgd(PGDFILE).lalo2xy((42.0037,9.3418))
#x = ffFromPgd(PGDFILE).lalo2xy((40.5993,-8.1801))

#y = ffFromPgd(PGDFILE).xy2lalo(x)

#print(x,y)
 
