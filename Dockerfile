FROM osgeo/gdal:ubuntu-small-3.3.1

RUN apt-get update

RUN apt install build-essential -y

RUN apt install libnetcdf-dev libnetcdf-cxx-legacy-dev -y

RUN apt install scons -y

WORKDIR /app

COPY Sconstruct /app/

COPY examples /app/examples

COPY src /app/src

RUN scons
