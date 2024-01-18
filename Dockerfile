FROM osgeo/gdal:ubuntu-small-3.3.1

# install requirements first to cache it
RUN apt-get update

RUN apt install build-essential -y

RUN apt install libnetcdf-dev libnetcdf-cxx-legacy-dev libnetcdf-c++4-dev -y

RUN apt install cmake -y

# copy files inside docker
WORKDIR /forefire

COPY . .

# install forefire
RUN sh cmake-build.sh

# add executable to PATH
RUN cp /forefire/bin/forefire /bin