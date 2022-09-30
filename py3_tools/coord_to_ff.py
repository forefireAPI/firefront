import argparse

from ff2geojson import *

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
	print('Coord', lat, lon)

	# reproj coords
	# save .ff file in examples DIR ../examples/aullene
	# os.system(forefire)
	# conver to geojson with specific name


if __name__ == '__main__':
    main()