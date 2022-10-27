# ForeFire

![logo](./doc/forefire.jpg)

ForeFire is an [open-source code for wildland fire spread models](https://www.researchgate.net/publication/278769168_ForeFire_open-source_code_for_wildland_fire_spread_models), developed and maintained by Université de Corse Pascal Paoli.

Access the [demo simulator here](http://forefire.univ-corse.fr/sim/dev/)

![demo](./doc/sim-forefire.jpg)


It has been designed and runs on Unix systems. Three modules can be built with the source code.

The main binaries are  
  - An interpreter (executable)
  - A shared library (with C/C++/Java and Fortran bindings)

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
- The C++ compiler
- [NetCDF Library](https://www.unidata.ucar.edu/software/netcdf/) and [NetCDF-C++ legacy](https://www.unidata.ucar.edu/downloads/netcdf/netcdf-cxx/index.jsp)
- [Cmake](https://cmake.org/) build tool

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
If you want to change it permanently, paste
```
export PATH=$PATH:</path/to/file>
```
at the end of your `~/.bashrc` file. The file can be edited with
```
nano ~/.bashrc
```


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

An example for the region of aullene in south France is provided. The example contains 3 files
- fuels.ff
- aullene.ff
- landscape.nc:

Run the example with

```
cd firefront/examples/aullene/

../../bin/forefire -i aullene.ff
```
The simulation result will be outputed in JSON format


### 4. Running with python

Installing requirements
```
cd py3_tools
pip install -r requirements.txt
```

You can use the script `coord_to_ff.py` to run the simulation in a default location

```
python coord_to_ff.py
```

For running in a chosen location, the script accepts latitude and longitude in epsg:4326 projection as inputs. It reprojects the coordinates into epsg:32632 projection, used in aullene's landscape.
```
python coord_to_ff.py --lat 41.6 --lon 9.1
```

The GeoJSON of geometry type Polygon will be saved in the `/examples/aullene` folder

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