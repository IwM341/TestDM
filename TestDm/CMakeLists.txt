cmake_minimum_required(VERSION 3.5)

project(TestDm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories( "D:\\Important\\articles\\DMFramework\\utils" )

set(MAINNAME main0.cpp)
add_executable(TestDm ${MAINNAME})
