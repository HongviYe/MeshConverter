cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
project(MeshConverter)
set(CMAKE_CXX_FLAGS "-std=c++11" ${CMAKE_CXX_FLAGS})
# options
option(MESHORIENT_BUILD_EXECUTABLE "Build executable." OFF)
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
IF(MSVC)
 add_compile_options(-bigobj)
endif()
set(SOURCES
    src/meshIO.cpp
    src/meshAlgorithm.cpp
    src/MeshConverter.cpp
    src/remesh.cpp
    
    src/Halfedge/PolyMesh_Base.cpp
    src/Halfedge/PolyMesh.cpp
    src/Halfedge/PolyMeshIterators.cpp
    src/Halfedge/AABB_Tree.cpp
    
    )
include_directories(./src)
include_directories(./extern/cli11)
include_directories(${EIGEN_BASE})
include_directories(./extern/meshorient/include)
add_executable(MeshConverter ${SOURCES})

target_link_libraries(MeshConverter MeshOrient)
