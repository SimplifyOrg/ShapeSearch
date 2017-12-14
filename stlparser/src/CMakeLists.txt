include_directories(${PROJECT_SOURCE_DIR}/stlparser/src/)
set(CMAKE_CXX_FLAGS "-O0 -g3 -Wall -std=c++11")

add_library(stlparser stlparser.cpp stl.cpp)
