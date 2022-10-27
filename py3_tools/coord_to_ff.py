import argparse

from forefirepy.ForeFire import *

def main():
	ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	ap.add_argument(
        "--lat",
				default=41.6,
        help="latitude of fire start"
    )
	ap.add_argument(
        "--lon",
				default=9.2,
        help="longitude of fire start"
		)
	args = vars(ap.parse_args())
	lat = args.get('lat')
	lon = args.get('lon')

	dir_path = os.path.dirname(os.path.realpath(__file__))
	output_path = dir_path + '/../examples/aullene/'
	filename = f'aullene2.ff'
	complete_path = output_path + filename

	[x, y] = reproject([lon, lat], inEpsg='epsg:4326', outEpsg='epsg:32632')

	ff = Forefire()
	ff.configBasicFf(lon=x, lat=y)
	ff.saveFf(complete_path)

	os.system(f'cd {output_path}; forefire -i {filename}')
	
	ff.convert_to_geojson(output_path + '0-2009-07-24T14-57-39Z.ffgeojson')


if __name__ == '__main__':
    main()