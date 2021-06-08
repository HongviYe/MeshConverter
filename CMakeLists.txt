cmake_minimum_required (VERSION 2.8)
project(MeshConverter)
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCES
    src/meshIO.h
    src/meshIO.cpp
    src/MeshConverter.cpp)
include_directories(./extern/cli11)
include_directories(./extern/eigen)
add_executable(MeshConverter ${SOURCES})
