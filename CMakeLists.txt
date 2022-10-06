cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(forefire VERSION 1.0)

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(forefire ${SRC_FILES})

target_include_directories(forefire PUBLIC usr/include)

# NOT WORKING YET -> SHOULD include netcdf