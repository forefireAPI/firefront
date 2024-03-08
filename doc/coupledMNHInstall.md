"""
Created on Fri Mar  8 16:00:33 2024

@author: filippi_j
How to install coupled cpde with Meso NH
"""

This should work on apple silicon with OSX Ventura


Requirements :

some packages installed with brew
gcc (you should have that already)
gfortran
openmpi
libxml2

get MesoNH preferably from depot at http://mesonh.aero.obs-mip.fr
but also from tarball if git troubles
http://mesonh.aero.obs-mip.fr/mesonh/dir_open/dir_MESONH/MNH-V5-7-0.tar.gz


!!!  Use bash, not ksh !!!! it may broke quite a bit of scripts

First export these variables :

export ARCH=LXgfortran
export VER_MPI=MPIAUTO
export OPTLEVEL=O2
export MNH_FOREFIRE=1.0

then enter the src directory in the MesoNH dir and do 
./configure

it should create -->  ../conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-7-0-FF-MPIAUTO-O2
load it by runing
. ../conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-7-0-FF-MPIAUTO-O2


SOME FIXUPS for OSX :
    libxml2 on configure
        if there is a problem while configure of Netcdf 4.9.2  (like requires libxml --disableDAP4-:
        a fixup : install autotools:
        brew install autoconf automake libtool
        go to the src/LIB/netcdf-c-4.9.2/ 
        edit configure.ac
        find line "	if test "x$ISOSX" = xyes && test "x$have_libxml2" = xno ; then " around l:1222
        remove "test"
        it should now look like 
        "	if test "x$ISOSX" = xyes && "x$have_libxml2" = xno ; then "
        run autoreconf to regenerate configure script from modified configure.ac
        test is now broken, configure should pass

    Netcdfcxx configure hangs at "-n 4.3.1" problem
        Meso-NH/src/LIB/netcdf-cxx4-4.3.1 there is a "VERSION" file that needs to be empty not "-n 4.4.3" 
        or whatever.... problem it is made during the "configure step" of netcdfcxx 
        as fixup.... you will need 2 terminals... one is in src/LIB/netcdf-cxx4-4.3.1 directory, the other mone is running the MesoNH make
        during the make, ith the netcdfcxx termnal, rerun frenetically : "cat VERSION; echo "" > VERSION " so it can be cleared in the time beween it is created and the one it is used..
        HINT : it is not comiming straight away after you launch make.. first is compiled HDF5 (long.. a few minutes),
        then after a more verbose netcdf9 (quick.. a few seconds.. start clearing), then eventually netcdfcxx

then run "make installmaster"
you should have a running set-up

you can then compile ForeFire, enven linking with the netcdfcxx from MesoNH if you need to
one forefire is installed and compiled with shared lib (it should be in firefront/lib/libforefireL.dylib)
go to the MesoNH  exe/ directory and link with 
ln -s PATHTOCOMPILEDFF/lib/libforefireL.dylib libForeFire.so
do in the exe dir you should have stuff like 
MESONH-LXgfortran-R8I4-MNH-V5-6-0-FF-MPIAUTO-O2 -> XXXX/src/dir_obj-LXgfortran-R8I4-MNH-V5-6-0-FF-MPIAUTO-O2/MASTER/MESONH
AND (for forefire) : 
libForeFire.so -> /Users/filippi_j/soft/firefront/lib/libforefireL.dylib

to run, get the case 016_FOREFIRE
go to the case directory
do :
export MPIRUN="mpirun -np 2" (for 2 cpus)

then 
RUN FIXUPS for OSX :
    Mpi run in sme kind of infinite loop
        go to the MesoNH bin dir and run :
        mv Mpirun NoMpirun   OSx is somehow not case sentitive in the path


