FROM osgeo/gdal:ubuntu-small-3.3.1

COPY install-requirements.sh .

RUN sh install-requirements.sh

COPY Sconstruct /app/

COPY examples /app/examples

COPY src /app/src

RUN scons
