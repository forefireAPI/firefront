
# NetCDF library paths to set here if required
#set(netcdfHome /Where/it/is/installed/)
#link_directories(${netcdfHome}/lib)
#include_directories(${netcdfHome}/include)


cmake_minimum_required(VERSION 3.10)
project(forefire VERSION 1.0)

# Set C++17 as the standard for the project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O0")
# Aggressive
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native -funroll-loops -fstrict-aliasing -ftree-vectorize")

# Output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)

# Link against NetCDF libraries globally (if all targets need it)
link_libraries("-lnetcdf_c++4 -lnetcdf")

# Main project executable
file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(forefire ${SRC_FILES})
target_link_libraries(forefire "-lnetcdf_c++4 -lnetcdf")

# Shared library
add_library(forefireL SHARED ${SRC_FILES})
target_link_libraries(forefireL "-lnetcdf_c++4 -lnetcdf")


# Specific main for the test if different
set(TEST_MAIN tests/ANNTest.cpp)
add_executable(ANN_test ${TEST_MAIN})