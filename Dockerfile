FROM osgeo/gdal:ubuntu-small-3.3.1

COPY install-requirements.sh .

RUN sh install-requirements.sh

WORKDIR /forefire

COPY Sconstruct .

COPY examples ./examples

COPY src ./src

RUN apt install build-essential -y

RUN scons

RUN mv /forefire/bin/forefire /bin