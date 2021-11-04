cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
project(MeshConverter)
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()
if(DOWNLOAD_WITH_CNPM)
    set (GITHUB_REPOSITE "github.com.cnpmjs.org")
elseif(DOWNLOAD_WITH_GITCLONE)
    set (GITHUB_REPOSITE "gitclone.com/github.com")
else()
    set (GITHUB_REPOSITE "github.com")
endif()
include(MeshConverterDependencies)
set(SOURCES
    src/meshIO.h
    src/meshIO.cpp
src/meshAlgorithm.cpp
src/meshAlgorithm.h
	src/MeshOrient.cpp
    src/MeshConverter.cpp)
include_directories(./extern/cli11)
include_directories(./extern/eigen)
add_executable(MeshConverter ${SOURCES})
