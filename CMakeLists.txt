cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -g -DNETCDF_NOT_LEGACY")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(netcdfHome /Users/filippi_j/soft/)
link_directories(${netcdfHome}/lib)
include_directories(${netcdfHome}/include)
project(forefire VERSION 1.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)

link_libraries("-lnetcdf_c++ -lnetcdf")

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(forefire ${SRC_FILES})

add_library(forefireL SHARED ${SRC_FILES})
