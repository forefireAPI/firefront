ARG UBUNTU_IMAGE=ubuntu:24.04

# Stage 1: builder
FROM ${UBUNTU_IMAGE} AS builder

ARG DEBIAN_FRONTEND=noninteractive

# Keep the BuildKit apt cache mount populated across builds: disable
# docker-clean's post-install hooks AND tell apt to keep downloaded .debs.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    rm -f /etc/apt/apt.conf.d/docker-clean \
 && echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' \
        > /etc/apt/apt.conf.d/keep-cache \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libnetcdf-c++4-dev \
        python3-dev \
        python3-pip

WORKDIR /build

# LICENSE is required by CPack at configure time (CMakeLists includes CPack).
COPY --link CMakeLists.txt LICENSE ./
COPY --link app/ ./app/
# tools/ is needed at configure time: CMakeLists.txt invokes
# tools/check_all_lfs.bash and references tools/runANN/ANNTest.cpp.
COPY --link tools/ ./tools/
COPY --link src/ ./src/

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
 && cmake --build build -j"$(nproc)"

COPY --link bindings/ ./bindings/

RUN --mount=type=cache,target=/root/.cache/pip \
    pip3 wheel --no-deps --wheel-dir /wheels ./bindings/python


# Stage 2: runtime
FROM ${UBUNTU_IMAGE}

ARG DEBIAN_FRONTEND=noninteractive

ENV FOREFIREHOME=/forefire \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    rm -f /etc/apt/apt.conf.d/docker-clean \
 && echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' \
        > /etc/apt/apt.conf.d/keep-cache \
 && apt-get update \
 && apt-get install -y --no-install-recommends \
        libnetcdf-c++4-1 \
        python3 \
        python3-pip \
        ca-certificates \
        curl

RUN --mount=type=cache,target=/root/.cache/pip \
    pip3 install --break-system-packages \
        lxml xarray netCDF4

COPY --from=builder --link /build/bin/forefire        /usr/local/bin/forefire
COPY --from=builder --link /build/lib/libforefireL.so /usr/local/lib/
COPY --from=builder --link /wheels/                   /tmp/wheels/

RUN --mount=type=cache,target=/root/.cache/pip \
    pip3 install --break-system-packages /tmp/wheels/*.whl \
 && rm -rf /tmp/wheels \
 && ldconfig

# Non-root user at UID 1000 so files written to mounted volumes get the
# host user's ownership. Ubuntu 24.04 ships a default `ubuntu` user at
# 1000, so we drop it first to free the slot. The /forefire/bin symlink
# keeps `../../bin/forefire` working in test scripts (tests/runff/ff-run.bash:16).
RUN userdel --remove ubuntu 2>/dev/null || true \
 && groupadd --gid 1000 forefire \
 && useradd  --gid 1000 --uid 1000 --create-home --shell /bin/bash forefire \
 && mkdir -p /forefire/bin \
 && ln -s /usr/local/bin/forefire /forefire/bin/forefire \
 && chown -R 1000:1000 /forefire

COPY --chown=1000:1000 app/      /forefire/app/
COPY --chown=1000:1000 bindings/ /forefire/bindings/
COPY --chown=1000:1000 tests/    /forefire/tests/
COPY --chown=1000:1000 tools/    /forefire/tools/

LABEL org.opencontainers.image.source="https://github.com/forefireAPI/forefire" \
      org.opencontainers.image.documentation="https://forefire.readthedocs.io/" \
      org.opencontainers.image.description="ForeFire is an open-source code for wildland fire spread models" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.title="forefire"

USER forefire
WORKDIR /forefire

EXPOSE 8000
CMD ["bash"]
