cmake_minimum_required(VERSION 3.23)
project(MCF_sket)

include_directories("/home/boy/3D/vcglib")
include_directories("/home/boy/3D/vcglib/eigenlib")
include_directories("/home/boy/3D/qhull/lib/")
include_directories("/home/boy/3D/qhull/src/")
include_directories("/home/boy/3D/qhull/src/libqhullcpp")

link_directories("/home/boy/3D/qhull/build")

add_library(libqhullstatic_r STATIC IMPORTED)
set_property(TARGET libqhullstatic_r PROPERTY IMPORTED_LOCATION /home/boy/3D/qhull/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qhullstatic_r${CMAKE_STATIC_LIBRARY_SUFFIX})
add_dependencies(libqhullstatic_r qhull)

add_library(libqhullcpp STATIC IMPORTED)
set_property(TARGET libqhullcpp PROPERTY IMPORTED_LOCATION /home/boy/3D/qhull/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qhullcpp${CMAKE_STATIC_LIBRARY_SUFFIX})
set_property(TARGET libqhullcpp PROPERTY INTERFACE_LINK_LIBRARIES libqhullstatic_r)
add_dependencies(libqhullcpp qhull)

LINK_DIRECTORIES(/home/boy/3D/qhull/build)
set(CMAKE_CXX_STANDARD 14)

add_executable(MCF_sket main.cpp LSMAT.hpp sphereShrinking.hpp utils.hpp meanCurvatureFlow.hpp LaplaceHelper.hpp)
target_link_libraries(${PROJECT_NAME} libqhullstatic_r libqhullcpp)
