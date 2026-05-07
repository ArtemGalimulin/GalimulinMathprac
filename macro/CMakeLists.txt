cmake_minimum_required(VERSION 3.11)
project(mesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

add_executable(mesh mesh.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so.4.12"
        INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)

target_link_libraries(mesh gmsh)