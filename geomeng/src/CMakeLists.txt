include_directories(${PROJECT_SOURCE_DIR}/geomeng/src/)
include_directories(${PROJECT_SOURCE_DIR}/stlparser/src/)
include_directories(/usr/include)
link_directories(/usr/lib)

set(CMAKE_CXX_FLAGS "-O0 -g3 -Wall -std=c++11")
add_executable(geomeng Source.cpp)
target_link_libraries(geomeng armadillo blas lapack m stlparser)
