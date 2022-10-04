import argparse
from turtle import st

from ff2geojson import *
from ForeFirepy3 import *

def main():
	ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	ap.add_argument(
        "--lat",
        help="latitude of fire start"
    )
	ap.add_argument(
        "--lon",
        help="longitude of fire start"
		)
	ap.add_argument(
        "--t",
        help="time step for simulation"
		)
	args = vars(ap.parse_args())
	lat = args.get('lat')
	lon = args.get('lon')
	t = args.get('t')

	[x, y] = reproject([lon, lat], inEpsg='epsg:4326', outEpsg='epsg:32632')

	dir_path = os.path.dirname(os.path.realpath(__file__))
	output_path = dir_path + '/../examples/aullene/'
	filename = f'{x}_{y}_{t}.ff'
	complete_path = output_path + filename

	ff = Forefire()
	ff.configBasicFf(lon=x, lat=y, t=t)
	ff.saveFf(complete_path)

	os.system(f'cd {output_path}; ../../bin/CommandShell -i {filename}')
	
	ffjson2geojson(output_path + '0-2009-07-24T14-57-39Z.json')


if __name__ == '__main__':
    main()