#!/usr/bin/env bash
#
# Install the NetCDF stack inside a manylinux_2_28 container, for cibuildwheel.
#
# EPEL provides the NetCDF C library and HDF5 for both x86_64 and aarch64, but
# no package for the legacy C++4 API that ForeFire's DataBroker includes as
# <netcdf>, so that one is built from source. Everything lands in /usr/local
# where auditwheel can find it and vendor it into the wheel.
set -euo pipefail

NETCDF_CXX4_VERSION="${NETCDF_CXX4_VERSION:-4.3.1}"
NETCDF_CXX4_SHA256="6a1189a181eed043b5859e15d5c080c30d0e107406fbb212c8fb9814e90f3445"
PREFIX="${PREFIX:-/usr/local}"

if [[ -f "${PREFIX}/include/netcdf" ]]; then
    echo "NetCDF C++4 already installed in ${PREFIX}, nothing to do."
    exit 0
fi

echo "::group::Installing NetCDF C library from EPEL"
dnf install -y epel-release
dnf install -y netcdf-devel hdf5-devel
echo "::endgroup::"

echo "::group::Building netcdf-cxx4 ${NETCDF_CXX4_VERSION}"
workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT
cd "${workdir}"

tarball="netcdf-cxx4-${NETCDF_CXX4_VERSION}.tar.gz"
curl -fsSL -o "${tarball}" \
    "https://downloads.unidata.ucar.edu/netcdf-cxx/${NETCDF_CXX4_VERSION}/${tarball}"
echo "${NETCDF_CXX4_SHA256}  ${tarball}" | sha256sum --check --strict

tar xzf "${tarball}"
cd "netcdf-cxx4-${NETCDF_CXX4_VERSION}"

# Autotools rather than the bundled CMake build: netcdf-cxx4 4.3.1 declares
# `cmake_minimum_required(VERSION 2.8.12)`, and CMake 4 refuses anything below
# 3.5. configure picks up the EPEL NetCDF through nc-config.
./configure \
    --prefix="${PREFIX}" \
    --enable-shared \
    --disable-static
make -j "$(nproc)"
make install
ldconfig
echo "::endgroup::"

echo "NetCDF C++4 installed:"
ls -l "${PREFIX}"/lib*/libnetcdf_c++4.so* "${PREFIX}/include/netcdf"
