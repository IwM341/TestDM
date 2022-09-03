cmake_minimum_required(VERSION 3.5)

project(TestDm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories( "..\\..\\DMFramework\\utils" )

set(MAINNAME main0.cpp)
#add_definitions(-DTEST_INTERPOL)
#add_definitions(-DTEST_GRID_FUNCTION)
add_definitions(-DTEST_M_GRID_FUNCTION)
#add_definitions(-DTEST_STRUCT_INTERPOL)
add_executable(TestDm ${MAINNAME})
