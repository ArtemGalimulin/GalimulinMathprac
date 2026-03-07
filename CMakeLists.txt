cmake_minimum_required(VERSION 3.11)
project(gmshdemo)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

add_executable(example example.cpp)
add_executable(torus torus.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so.4.12"
        INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)

target_link_libraries(torus gmsh)