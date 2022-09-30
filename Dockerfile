FROM osgeo/gdal:ubuntu-small-3.3.1

# install requirements first to cache it
RUN apt-get update

RUN apt install build-essential -y

RUN apt install libnetcdf-dev libnetcdf-cxx-legacy-dev -y

RUN apt install scons -y

# copy files inside docker
WORKDIR /forefire

COPY Sconstruct .

COPY examples ./examples

COPY src ./src

# install forefire
RUN scons

# add executable to PATH
RUN cp /forefire/bin/forefire /bin