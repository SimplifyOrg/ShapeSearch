cmake_minimum_required (VERSION 2.6)
project(ShapeSearch)

add_subdirectory(stlparser/src)
add_subdirectory(geomeng/src)

set(CMAKE_CXX_FLAGS "-O0 -g3 -Wall -std=c++11")
