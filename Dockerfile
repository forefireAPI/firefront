FROM osgeo/gdal:ubuntu-small-3.3.1

WORKDIR /forefire

# install requirements first to cache it
COPY install-requirements.sh .

RUN sh install-requirements.sh

# install forefire
COPY Sconstruct .

COPY examples ./examples

COPY src ./src

RUN scons

RUN mv /forefire/bin/forefire /bin