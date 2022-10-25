cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(forefire VERSION 1.0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/../bin)

link_libraries("-lnetcdf_c++ -lnetcdf")

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(forefire ${SRC_FILES})