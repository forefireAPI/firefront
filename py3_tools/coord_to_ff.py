import argparse

from ff2geojson import *
from ForeFirepy3 import *

def main():
	ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	ap.add_argument(
        "--lat",
        help="latitude in aullene"
    )
	ap.add_argument(
        "--lon",
        help="longitude in aullene"
		)
	args = vars(ap.parse_args())
	lat = args.get('lat')
	lon = args.get('lon')

	[x, y] = reproject([lon, lat], inEpsg='epsg:4326', outEpsg='epsg:32632')

	output_path = create_ff(x,y)['output_path']
	filename = create_ff(x,y)['filename']

	os.system(f'cd /firefront/examples/aullene; ../../bin/CommandShell -i {filename}')
	
	ffjson2geojson(output_path + '0-2009-07-24T14-57-39Z.json')


if __name__ == '__main__':
    main()