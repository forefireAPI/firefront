#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#dynamic of surface winds, a video .
#surface wind with clear thresholds#
#
#- report is more a text, the strory. 
#- a white column, .. expect downburst 15knots ahead of the front.
#- only for extended attack.. like if fire still burn next moning, data should be ready.


"""
Created on Mon Sep 11 11:34:06 2023

Sample tools usages for pre and post processing of coupled MesoNH-ForeFire

@author: Jean-Baptiste Filippi
@todo : timed kml with wind from preal -
@todo : pMNH2vtk - with pgd file, not with domains
@todo : pMNH2vtk : include bmaptovtk with the atimebmap in standard and ZS from PDG
@todo : mode output diagnowtics
@todo : 2Dmore frequent outputs
@todo : include the Mars request
@todo : a documentation that is also a presentation for MesoNH days
@todo : BMap reconstructions from contours... first anything to FireFront format, then make simulations, then compile results
@todo : intelligent analysis of reconstructed bmaps : find where there is the most coherent data and speed (filter by ROS, gradient like wind/slope...)
@How To simulate a wildfire with ForeFire-MesoNH :
@ 1) get the source : site github, compile : cmake, put lib in exe folder,compile MNH with "FF=1.0"
@ 2) have a python ready environment (see yaml file)) 
@ first step - generate the case - there is a sample use file in the folder - you need Lat-Lon-date and a pattern (in the template)
@ second step - get BC data - sample 
@ sautes de feu ! comme les 
"""

from datetime import datetime, timedelta


orig = "/Users/filippi_j/data/2024/barbaggio/nest150/"

Barbaggio_03012024 = {
    "run_info": {
        "start_time": "2024-01-03T03:00:00",
        "end_time": "2024-01-04T13:00:00",
        "latitude_center": 42.689245,
        "longitude_center": 9.3964812,
        "XOR1TO2": 46,
        "YOR1TO2": 46,
        "XOR2TO3": 60,
        "YOR2TO3": 60
    },
    "ignitions": [
        {"when": "2024-01-03T10:36:00", "latitude": 42.689245, "longitude": 9.3864812}
        ],
    "case_path": "%s/"%orig,
    "template_path": "/Users/filippi_j/soft/firefront/tools/3nestFFCASE/",
    "PGDFILES": [
        "/001_pgd/PGD_D2000mA.nc",
        "/001_pgd/PGD_D400mA.nc",
        "/001_pgd/PGD_D80mA.nc"
    ],
    "CoupledLandscape_path": "/006_runff/ForeFire/landcoupled.nc",
    "UnCoupledLandscape_path": "/006_runff/ForeFire/land.nc",
    "initff_path": "/006_runff/ForeFire/Init.ff",
    "BMAPFILE": "/006_runff/ForeFire/Outputs/ForeFire.0.nc",
    "FFINPUTPATTERN": "/006_runff/ForeFire/Outputs/output.0.*",
    "MODELOUTPATTERN": ["/006_runff/MODEL1/output","/006_runff/MODEL2/output","/006_runff/MODEL3/output"],
    
    "FFOUTVTKPATH": ["/RESULTS/vtkout1/","/RESULTS/vtkout2/","/RESULTS/vtkout3/"],
    "OUTKMLDOMAINFILE": "/RESULTS/domains.kml",
    "frontsKMLOUT": "/RESULTS/fronts.kml",
    "BMAPKMLOUT": "/RESULTS/ros.kml",
    "fuel_TIF_path": "/landscape/fuel.tif",
    "fuel_png_path": "/landscape/fuel.png",
    "fuel_kml_path": "/landscape/fuel.kml",
    "output_land_path": "/landscape/",
    "fuelff_path": "/006_runff/ForeFire/fuel.ff"
}




CAST = Barbaggio_03012024
 

gen_empty_FFMNH_case = False
gen_domain_kml = False
gen_fuel_map  = False
gen_fuel_map_from_landcover = False
gen_FAF_case = False
gen_AF_case = False
gen_2D_VTK_OUT = False
gen_3D_VTK_OUT = False
gen_KML_OUT = True 
gen_TESTFUEL_case = False

OUTKMLDOMAINFILE = CAST["case_path"]+CAST["OUTKMLDOMAINFILE"]
PGDFILES = [CAST['case_path'] + file for file in CAST['PGDFILES']]
FFOUTVTKPATH = [CAST['case_path'] + file for file in CAST['FFOUTVTKPATH']]
MODELOUTPATTERN = [CAST['case_path'] + file for file in CAST['MODELOUTPATTERN']]

fuel_TIF_path= CAST['case_path'] + CAST['fuel_TIF_path']
fuel_png_path= CAST['case_path'] + CAST['fuel_png_path']
fuel_kml_path= CAST['case_path'] + CAST['fuel_kml_path']
CoupledLandscape_path= CAST['case_path'] + CAST['CoupledLandscape_path']
UnCoupledLandscape_path= CAST['case_path'] + CAST['UnCoupledLandscape_path']
initff_path= CAST['case_path'] + CAST['initff_path']
BMAPFILE = CAST['case_path'] + CAST['BMAPFILE'] 
FFINPUTPATTERN= CAST['case_path'] + CAST['FFINPUTPATTERN']
BMAPKMLOUT= CAST['case_path'] + CAST['BMAPKMLOUT']
frontsKMLOUT= CAST['case_path'] + CAST['frontsKMLOUT']
output_land_path=  CAST['case_path'] + CAST['output_land_path']
fuelff_path=  CAST['case_path'] + CAST['fuelff_path']


# with date and location, define domains > generate template and copy files also from weather archive
if gen_empty_FFMNH_case:
    # Convertir les chaînes en objets datetime
    start_date = datetime.fromisoformat(CAST['run_info']['start_time'])
    end_date = datetime.fromisoformat(CAST['run_info']['end_time'])
    
    # Extraire la date de début au format YYYYMMDD
    start_date_str = start_date.strftime("%Y%m%d")
    
    # Calculer la différence en secondes entre la fin et le début
    delta_seconds = (end_date - start_date).total_seconds()
    
    # Trouver la valeur heure de départ en multiple de 3
    start_hour = start_date.hour
    start_hour = (start_hour // 3) * 3
    
    # Calculer la valeur heure de fin en multiple de 3
    delta_hours = delta_seconds // 3600  # Convertir la différence en heures
    end_hour = ((delta_hours + 2) // 3) * 3  # Arrondir à l'entier multiple de 3 supérieur
    end_hour += start_hour  # Ajouter la valeur heure de départ
    
    start_hour_str = "%02d"%start_hour
    end_hour_str = "%02d"%end_hour
      
    # Trouver l'ignition la plus tôt
    earliest_ignition_time = min(datetime.fromisoformat(ignition["when"]) for ignition in CAST["ignitions"])
   
    # Trouver l'heure arrondie au multiple de 3 immédiatement inférieure à cette différence, mais supérieure à l'heure de début
    MODEL_1_ALONE_max_hour =((earliest_ignition_time.hour // 3)+1) * 3
    MODEL_2_3_hour = earliest_ignition_time.hour - 1
    MNHFILE_2_3_hour = MODEL_2_3_hour - start_hour 
    
    M123_PREAL_START = ((earliest_ignition_time.hour // 3)) * 3
    M123_PREAL_END = int(end_hour)
    
    # Exécution de la fonction
    
    print(f"# Fire at date {start_date_str}, ignition at time {earliest_ignition_time}" )    
    print(f"# Run from {CAST['run_info']['start_time']} to {CAST['run_info']['end_time']}" )
    print(f"# First domain running alone from {start_hour_str} to {MODEL_1_ALONE_max_hour} " )
    print(f"# Domain 2 and 3 starting at time {MODEL_2_3_hour} that is hourly step {MNHFILE_2_3_hour} of run1 that started at {start_hour_str} " )
    print(f"# Domain 1-2-3 starting at time {MODEL_2_3_hour} (file init {MNHFILE_2_3_hour}) using PREAL BC from {M123_PREAL_START} to {M123_PREAL_END} inited at {start_hour_str} " )
    
    print(f"# Configuration files generation script :")
    print(f"cp -r {CAST['template_path']} {CAST['case_path']}")
    print(f"cd {CAST['case_path']}")
    print(f"cd 001_pgd/ ; bash MAKE_PGD {CAST['run_info']['latitude_center']} {CAST['run_info']['longitude_center']} {CAST['run_info']['XOR1TO2']} {CAST['run_info']['YOR1TO2']} {CAST['run_info']['XOR2TO3']} {CAST['run_info']['YOR2TO3']}") 
    print(f"cd ../002_real/ ; bash MAKE_PREAL {start_hour_str} {end_hour_str} {start_date_str} ") 
    print(f"cd ../003_run/ ; bash MAKE_RUN1 {start_hour} {MODEL_1_ALONE_max_hour} {start_date_str} ") 
    print(f"cd ../004_SpawnReal/ ; bash MAKE_SPAWNREAL {MNHFILE_2_3_hour} {start_date_str} {MODEL_2_3_hour}") 
    print(f"cd ../005_SpawnReal/ ; bash MAKE_SPAWNREAL23 {start_date_str} {MODEL_2_3_hour}")

    print(f"cd ../006_runff/ ; bash MAKE_RUNFF {start_date_str} {M123_PREAL_START} {M123_PREAL_END} {MNHFILE_2_3_hour} {MODEL_2_3_hour}")     
    
# now run PGD, and see if domains OK in google earth/ re-run and adjust if needed
if gen_domain_kml:
    from preprocessing.kmlDomain import pgds_to_KML
    pgds_to_KML(PGDFILES,OUTKMLDOMAINFILE)

# Now we can try to make a fuel distribution file using PNG's and databases.. approach here is from 
if gen_fuel_map:
    from preprocessing.learnFuelTifPng import makeFuelMapFromPgd
    print(PGDFILES[-1])
    makeFuelMapFromPgd(PGDFILES[-1],fuel_TIF_path,fuel_png_path,fuel_kml_path)   

if gen_fuel_map_from_landcover:
    from preprocessing.extract_from_eu_land_cover import landcover_roads_to_fuel
    from preprocessing.ffToGeoJson import get_WSEN_LBRT_ZS_From_Pgd

    land_in_path_files = "/Users/filippi_j/data/landEurope/"
    S2GLC_tif = f"{land_in_path_files}S2GLC_Europe_2017_v1.2.tif"
    legend_file_path = f"{land_in_path_files}S2GLC_Europe_2017_v1.2_legend_and_version_info.txt"
 
    WSEN, LBRT, ZS = get_WSEN_LBRT_ZS_From_Pgd(PGDFILES[-1])
    
    landcover_roads_to_fuel(S2GLC_tif,legend_file_path, WSEN, LBRT, output_land_path)
    
# IF OK, make all required files to ForeFire, also need tima and start
# FA : Fire-To-Atm one is AF : Atm-To-Fire last is FAF: Fire-To-Atm-To-Fire

if gen_FAF_case:
    from preprocessing.prealCF2Case import PGD2Case
    from preprocessing.ffToGeoJson import ffFromPgd
    domTime = datetime.fromisoformat(CAST['run_info']['start_time'])
    PGD2Case(PGDFILES[-1],fuel_png_path,UnCoupledLandscape_path, domTime,gen_wind=1)
    with open(initff_path, 'w') as file:
        file.write(ffFromPgd(PGDFILES[-1],domainDate=domTime,ignitions = CAST['ignitions']))

if gen_TESTFUEL_case:
    import csv
    from preprocessing.prealCF2Case import PGD2Case
    from preprocessing.ffToGeoJson import ffFromPgd
    # from preprocessing.extract_from_eu_land_cover import landcover_roads_to_fuel
    # from preprocessing.ffToGeoJson import get_WSEN_LBRT_ZS_From_Pgd

    # land_in_path_files = "/Users/filippi_j/data/landEurope/"
    # S2GLC_tif = f"{land_in_path_files}S2GLC_Europe_2017_v1.2.tif"
    # legend_file_path = f"{land_in_path_files}S2GLC_Europe_2017_v1.2_legend_and_version_info.txt"
 
    # WSEN, LBRT, ZS = get_WSEN_LBRT_ZS_From_Pgd(PGDFILES[-1])
    
    # landcover_roads_to_fuel(S2GLC_tif,legend_file_path, WSEN, LBRT, output_land_path)
    
    indices = []
    with open(fuelff_path, "r") as file:
        reader = csv.reader(file, delimiter=';')
        next(reader)  # Skip header
        for row in reader:
            indices.append(int(row[0]))
    print(indices)
    indices = [73, 75, 82, 83, 102, 103, 104]
    domTime = datetime.fromisoformat(CAST['run_info']['start_time'])
    PGD2Case(PGDFILES[-1],fuel_png_path,UnCoupledLandscape_path, domTime, fuel_test=indices)
    
    with open(initff_path, 'w') as file:
        file.write(ffFromPgd(PGDFILES[-1],domainDate=domTime,fuel_test=indices))

# for the atmosphere to fire only, Trying to use the PREAL files from Domain Nmax (smallest)
if gen_AF_case:
    PGD2Case(FIREDOMPGDFILE,IGNITIONS,OUTFFLANDSCAPE, OUTFFINIT)    

# now we assume everythng has run we can do the post-processing first the 2D VTKOUT
if gen_2D_VTK_OUT:
    print_vtkout_slurm_command_2D()
    
# now wa assume everythng has run we can do the post-processing the 3D VTKOUT
if gen_3D_VTK_OUT:
    from postprocessing.pMNHFF2VTK import ffmnhFileToVtk, ffFrontsToVtk
    ffmnhFileToVtk(inpattern = MODELOUTPATTERN[1],pgdFile = PGDFILES[1],outPath = FFOUTVTKPATH[1])

if gen_KML_OUT:
    from preprocessing.ffToGeoJson import genKMLFiles
    genKMLFiles(PGDFILES[-1], BMAPFILE, FFINPUTPATTERN, BMAPKMLOUT,frontsKMLOUT, everyNFronts=6, change_color_every=1,btslice=(40,-1))

from preprocessing.ffToGeoJson import plotRos
#plotRos(PGDFILES[-1], BMAPFILE,max_speed_filter=1000.0)

#wind params examples : 
#Southerly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(0.0,10.0,0.);t=0.]
#Westerly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(10.0,0.0,0.);t=0.]
#Northerly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(0.0,-10.0,0.);t=0.]
#Easterly 10m.s-1  trigger[wind;loc=(0.,0.,0.);vel=(-10.0,0.0,0.);t=0.]

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime



# define when to start
# define the domain
#

cost = [0.0,
0.0,
0.0,
210.17010116577148,
211.17010116577148,
215.17010116577148,
221.42010116577148,
229.67010116577148,
432.80348205566406,
1968.0598907470703,
2022.0598907470703,
2299.4276733398438,
2949.2632446289062,
4105.273971557617,
15482.962337493896,
20914.563625335693,
20958.063625335693,
22070.49808883667,
22102.74808883667,
22152.99808883667,
27016.63991165161,
31189.866333007812,
31412.628520965576,
33319.78584098816,
37538.17907142639,
42307.246057510376,
45969.47225379944,
47950.97472572327,
63450.30722236633,
98820.76568984985,
118991.89973831177,
165543.5011806488,
207644.5391407013,
258222.29362487793,
316281.03370666504,
359526.27380371094,
426399.09272644046,
471431.36740753177,
568360.2472064209,
627504.1055247497,
645412.9739634704,
659830.8771464538,
755253.0878855896,
804791.1836268615,
825648.5333583069,
865845.019369812,
888656.672583313,
955530.7864253235,
1026989.1899516296,
1120100.1706092071,
1237810.4050776672,
1356820.6299311828,
1514884.2269732666,
1594973.6059386444,
1651007.5541236114,
1683435.6126587105,
1731498.5213243675,
1739892.421066494,
1770438.2120802116,
1782285.448363514,
1818582.7460586738,
1861542.8454639625,
1894287.3541576576,
1905846.0894405555,
1929709.776379795,
1939955.0624038887,
1969903.5098488044,
1995118.1327755165,
2026537.5029518318,
2054908.0466148567,
2056578.453730793,
2093923.470576496,
2110355.661647053,
2131127.224284382,
2141263.9343216135,
2160388.001452656,
2167340.813575001,
2171339.1651680185,
2185996.7539436533,
2190421.6003219797,
2195660.7494308664,
2198974.6313583567,
2202753.9343979075,
2206021.3943664744,
2210215.363146992,
2221235.047381611,
2227777.806933613,
2235665.1993667795,
2237372.2564346506,
2247969.6464073374,
2249192.6416084482,
2250594.816216679,
2250638.816216679,
2265639.499686451,
2282383.091773243,
2282427.591773243,
2284826.036002369,
2286531.68581316,
2286576.93581316,
2290015.165564747,
2300200.7319557383,
2304126.8798134043]

def plotArea(larea,cf):
    areas = [area / 10000 for area, _ in larea]  # Convert from m² to hectares
    dates = [date for _, date in larea]
    modified_dates = [date + timedelta(hours=2) for date in dates]

    # Create the plot
    plt.figure(figsize=(20, 12))
    ax1 = plt.gca()
    ax1.plot(modified_dates, areas, marker='o', label='Surface')
    ax1.axhline(y=200, color='red', linestyle='-')
    ax1.axhline(y=35, color='green', linestyle='-')
    ax1.set_ylabel('Surface parcourue (en hectares)')

    # Adding the derivative on a secondary y-axis
    ax2 = ax1.twinx()
    ax2.plot(modified_dates, cf, marker='x', color='green', label='Coût du brulé, méthode LISA/A.Belgodere')
    ax2.set_ylabel('Coût du brulé, en million Euros, méthode LISA/A.Belgodere')

    # Format the x-axis
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    ax1.xaxis.set_major_locator(mdates.HourLocator())
    plt.xticks(rotation=45)

    # Activate the grid
    plt.grid(True)

    # Set the x-axis label
    plt.xlabel('Date (UTC)')

    # Set the title
    plt.title('Coût et Surface Parcourue (simulation) - en rouge 200ha / 4 Meur')

    # Show the plot
    plt.show()

def makeCost():
    from preprocessing.ffToGeoJson import genCostDict
    out_path = 'cost/results_badbaggio.pickle'
    import pickle
    #doto = genCostDict(PGDFILES[-1], BMAPFILE, FFINPUTPATTERN, everyNFronts=2,btslice=(35,-1))
    #with open(out_path, 'wb') as file:
    #    pickle.dump(doto,file)
    
    import numpy as np
    from datetime import datetime
    
    with open(out_path, 'rb') as file:
        data = pickle.load(file)
    
    
    costTimeFormat = '%Y%m%d-%H%M%S'
    
    areas = [np.sum(a) for a in data['area'][0]]
    cdates = [datetime.strptime(a, costTimeFormat) for a in data['contour_datetime']]
             
    burned = list(zip(areas,cdates))
    mycost =np.array(cost)/10000
    print(np.max(mycost))
    plotArea(burned,mycost)











