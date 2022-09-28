import sys
import os
import json
from pyproj import Transformer

def load_json(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
        return data

def get_inEpsg(data):
    return data['fronts'][0]['projection'].lower()

def get_coords(data):
    coordinates_3 = data['fronts'][0]['coordinates'].split()
    coords = []
    for c in coordinates_3:
        c2 = c.rpartition(',')[0].split(',')
        coords.append([float(c2[0]), float(c2[1])])
    return coords

def reproject(coords, inEpsg):
    [x1,y1] = coords
    outEpsg = 'epsg:4326'
    transformer = Transformer.from_crs(inEpsg, outEpsg, always_xy=True)
    x2, y2 = transformer.transform(x1, y1)
    return [x2,y2]

def get_coords_4326(coords, inEpsg):
    coords_4326 = []
    for c in coords:
        coords_4326.append(reproject(c, inEpsg))
    return coords_4326

def save_geojson(filepath, coords_4326):
    data = {}
    data['type'] = 'Feature'
    geometry = {}
    geometry['type'] = 'MultiPoint'
    data['geometry'] = geometry
    geometry['coordinates'] = coords_4326

    filenameJson = filepath.split('/')[-1]
    filename = filenameJson.split('.')[0]

    savePath = filepath.replace(filenameJson, '')
    completeName = os.path.join(savePath, f'{filename}.geojson')

    with open(completeName, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

def ffjson2geojson(filepath):
    data = load_json(filepath)
    inEpsg = get_inEpsg(data)
    coords = get_coords(data)
    coords_4326 = get_coords_4326(coords, inEpsg)
    save_geojson(filepath, coords_4326)

if __name__ == '__main__':
    try:
        filepath = sys.argv[1]
        ffjson2geojson(filepath)
        print(f'{filepath} converted to GeoJSON.')
    except Exception as e:
        if (str(e) == 'list index out of range'):
            print('Please provide the json file path')
        else:
            print(e)