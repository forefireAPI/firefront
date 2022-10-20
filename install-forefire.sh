#!/bin/bash
echo "====== UNIX REQUIREMENTS ========"

apt-get update

apt install build-essential -y

apt install libnetcdf-dev libnetcdf-cxx-legacy-dev -y

apt install scons cmake -y

echo "==========================="
echo "========= FOREFIRE ========"
echo "==========================="

scons