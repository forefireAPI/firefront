# ForeFire

ForeFire is an open-source code for wildland fire spread models. The complete paper can be found [here](https://www.researchgate.net/publication/278769168_ForeFire_open-source_code_for_wildland_fire_spread_models)


It has been designed and run on Unix systems, three modules can be built with the source code.

The main binaries are
  
  - An interpreter (executable)
  - A dynamic library (shared, with C/C++/Java and Fortran bindings)
  - Pre-Post processing helping scripts

## 1. Requirements

The folowing comands are suited for Ubuntu or Debian distributions.

### Update apt repository
```
sudo apt update
```

### C++ compiler
Compilation requires a c++ compiler, that usually comes pre-built with your linux distribution, but can be installed with
```
sudo apt install build-essential
```

### NetCDF Library 

```
sudo apt install libnetcdf-dev libnetcdf-cxx-legacy-dev
```

- https://www.unidata.ucar.edu/software/netcdf/

- NetCDF-C++ >>LEGACY<< is required for compatibilities issues
https://www.unidata.ucar.edu/downloads/netcdf/netcdf-cxx/index.jsp


### Scons

The SCons python tool is used to make the executable and the python library
```
sudo apt get scons
```

More about at 
- https://www.scons.org/
- https://scons.org/documentation.html

## 2. Building the executable

Build  with
```
cd firefront

scons
```
A sample `SConstruct` file is included with the distribution, try it and if it does not work, set the environment variables, edit it and insert the path to the Netcdf (and Java headers for JNI bindings if required).

## 3. Running an example

```
cd firefront/Examples/aullene/

../../CommandShell -i aullene.ff
```
The simulation result will be outputed in Json format

## 4. Building python Lib
The "swig" repository contains python bindings requires numpy (and numpy.i), swig, and matplotlib for testing. 
