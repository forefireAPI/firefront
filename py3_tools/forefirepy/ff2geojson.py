import sys
import json
from pyproj import Transformer

def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
        return data

def reproject(xy, inEpsg, outEpsg='epsg:4326'):
    [x1, y1] = xy
    transformer = Transformer.from_crs(inEpsg, outEpsg, always_xy=True)
    x2, y2 = transformer.transform(x1, y1)
    return [x2,y2]

def ffjson2geojson(filepath):
    ff_geojson = load_json(filepath)
    
    #reproject
    inEpsg = ff_geojson['projection'].lower()
    for feature in ff_geojson["features"]:
        reproj = [reproject(xy, inEpsg) for xy in feature['geometry']['coordinates'][0]]
        reproj.append(reproj[0])
        feature['geometry']['coordinates'][0] = reproj
    ff_geojson['projection'] = 'epsg:4326'

    #save
    fileextension = filepath.split('.')[-1]
    savePath = filepath.replace(f'.{fileextension}', '.geojson')
    with open(savePath, 'w', encoding='utf-8') as f:
        json.dump(ff_geojson, f, ensure_ascii=False, indent=4)
    

if __name__ == '__main__':
    try:
        filepath = sys.argv[1]
        ffjson2geojson(filepath)
    except Exception as e:
        print(e)