ForeFire has been designed and run on Unix systems, three modules can be built with the source code.

  - An interpreter (executable)
  - A dynamic library (shared, with C/C++/Java and Fortran bindings)

NetCDF  Library V3 or later must be installed on the system to build Forefire
Get it from http://www.unidata.ucar.edu/software/netcdf/

Compilation requires a c++ compiler, but it has only been tested on gcc/g++ compiler.
The SCons python tool is used to make the library and executable, get it from  http://www.scons.org
A sample SConstruct file is included with the distribution, try it and if it does not work, set the environment variables, edit it and insert the path to the Netcdf (and Java headers for JNI bindings if required).

to run it, type "./CommandShell -i examplescript" from the commandline

The "swig" repository contains python bindings requires numpy (and numpy.i), swig, and matplotlib for testing. 
