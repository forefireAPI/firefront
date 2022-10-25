# ForeFire

![logo](./doc/forefire.jpg)

ForeFire is an [open-source code for wildland fire spread models](https://www.researchgate.net/publication/278769168_ForeFire_open-source_code_for_wildland_fire_spread_models), developed and maintained by Université de Corse Pascal Paoli.

Access the [demo simulator here](http://forefire.univ-corse.fr/sim/dev/)

![demo](./doc/sim-forefire.jpg)


It has been designed and run on Unix systems, three modules can be built with the source code.

The main binaries are
  
  - An interpreter (executable)
  - A dynamic library (shared, with C/C++/Java and Fortran bindings)

## 1. Requirements

The requirements and ForeFire can be installed by running `install-forefire.sh` (Ubuntu or Debian distributions)

```
cd forefire

sudo sh install-forefire.sh
```

The program will be built in: `./bin/forefire`

OR run the commands:

```
apt-get update

apt install build-essential -y

apt install libnetcdf-dev libnetcdf-cxx-legacy-dev -y

apt install cmake -y
```

To install
- The C++ compiler (that may come pre-built with your linux distribution)
- [NetCDF Library](https://www.unidata.ucar.edu/software/netcdf/) and [NetCDF-C++ legacy](https://www.unidata.ucar.edu/downloads/netcdf/netcdf-cxx/index.jsp), for compatibilities issues
- [cmake](https://cmake.org/), is used to build the executable and the python library

## 2. Build

### 2.1 Cmake

To build with cmake run the script
```
sh cmake-build.sh
```

To make the program [executable from eveywhere](https://unix.stackexchange.com/questions/3809/how-can-i-make-a-program-executable-from-everywhere) (during the session) Add the bin folder to path
```
export PATH=$PATH:`pwd`/bin
```
If you want to change it permanently add `export PATH=$PATH:</path/to/file>` to your ~/.bashrc file


### 2.2 Scons and Other build systems

Forefire can also be built with [scons](https://www.scons.org/). Install it with
```
apt install scons -y
```

A sample `SConstruct` file is included with the distribution.
Run it with
```
scons
```

To build with all warnings enabled
```
 scons -Q w=1 
```

Troubleshooting: If it does not work, try replacing the `Sconstruct` file with `./tools/Sconstruct`. Set the environment variables, and insert the path to the Netcdf (and Java headers for JNI bindings if required).

Make: A simple `makefile` is also available in the `tools` directory

More info on build systems can be found on [this issue](https://github.com/forefireAPI/firefront/issues/9)

## 3. Running an example

```
cd firefront/examples/aullene/

../../bin/forefire -i aullene.ff
```
The simulation result will be outputed in JSON format

### Converting output JSON to GeoJSON


Use the script `ff2geojson.py` with the .json file as argument.
```
python3 py3_tools/ff2geojson.py examples/aullene/1-2009-07-24T15-01-00Z.json
```
The JSON will be converted to GeoJSON (EPSG 4326) of geometry type MultiPoint and saved in the same directory.

### Running an example in a chosen location

Use the script `coord_to_ff.py` using `--lon` and `--lat` flags to pass coordinates. An example simulation will be outputted in GeoJSON in `/examples/aullene`.

## 4. Building python Lib
The `/swig` folder contains and `Sconstruct` file for python bindings.

Requires numpy (and numpy.i), swig, and matplotlib for testing. 

## 5. Building with Docker
A sample Dockerfile can allow to build a Docker image with
```
docker build . -t forefire
```

To run this image and interactively acces the continer use
```
docker run -it forefire bash
```