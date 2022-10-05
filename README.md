# ForeFire

![logo](./doc/forefire.jpg)

ForeFire is an [open-source code for wildland fire spread models](https://www.researchgate.net/publication/278769168_ForeFire_open-source_code_for_wildland_fire_spread_models), developed and maintained by Université de Corse Pascal Paoli.

Access the [demo simulator here](http://forefire.univ-corse.fr/sim/dev/)

![demo](./doc/sim-forefire.jpg)


It has been designed and run on Unix systems, three modules can be built with the source code.

The main binaries are
  
  - An interpreter (executable)
  - A dynamic library (shared, with C/C++/Java and Fortran bindings)
  - Pre-Post processing helping scripts

## 1. Installation

### 1.1 By script
The requirements and ForeFire can be installed by running `install-forefire.sh` (Ubuntu or Debian distributions)

```
cd forefire

sudo sh install-requirements.sh
```

The program will be built in: `/bin/forefire`

To be executable from anywhere add it to path
```
export PATH=$PATH:./bin
```

### 1.2 Manual Instalation - Requirements

Run

```
apt-get update

apt install build-essential -y

apt install libnetcdf-dev libnetcdf-cxx-legacy-dev -y

apt install scons -y
```
To install:

C++ compiler
- Compilation requires a c++ compiler, that usually comes pre-built with your linux distribution

NetCDF Library
- https://www.unidata.ucar.edu/software/netcdf/

- NetCDF-C++ `legacy` is required for compatibilities issues
https://www.unidata.ucar.edu/downloads/netcdf/netcdf-cxx/index.jsp

Scons
- The [SCons python tool](https://www.scons.org/) is used to make the executable and the python library

### 1.3 Manual installation - Building with scons

A sample `SConstruct` file is included with the distribution.
Run it with
```
cd firefront

scons
```

Troubleshooting: If it does not work, try using the `/tools/Sconstruct` file. Replace the `Sconstruct` file with `/tools/Sconstruct`. Set the environment variables, and insert the path to the Netcdf (and Java headers for JNI bindings if required).

## 2. Running an example

```
cd firefront/examples/aullene/

../../bin/CommandShell -i aullene.ff
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


## 3. Building python Lib
The "swig" repository contains python bindings requires numpy (and numpy.i), swig, and matplotlib for testing. 

## 4. Building with Docker
A sample Dockerfile can allow to build a Docker image with
```
docker build . -t forefire
```

To run this image and interactively acces the continer use
```
docker run -it forefire bash
```