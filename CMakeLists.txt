cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
project(MeshConverter)
set(CMAKE_CXX_FLAGS "-std=c++11" ${CMAKE_CXX_FLAGS})
# options
option(MESHORIENT_BUILD_EXECUTABLE "Build executable." OFF)
option(MESHORIENT_DOWNLOAD "Download extern dependent library." OFF)
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

if(NOT EXTERNAL_DIR)
    set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/extern)
endif()

if(MESHORIENT_DOWNLOAD)
    include(MeshConverterDependencies)
else()
endif()

set(EIGEN_BASE ${EXTERNAL_DIR}/eigen)
set(IGL_BASE ${EXTERNAL_DIR}/libigl)


IF(MSVC)
 add_compile_options(-bigobj)
endif()


get_directory_property(hasParent PARENT_DIRECTORY)
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
if(hasParent)
set(EXTERNAL_DIR ${hasParent}/third_party)
set(EIGEN_BASE ${EXTERNAL_DIR}/eigen)
set(IGL_BASE ${EXTERNAL_DIR}/libigl)
message("the include path =${EXTERNAL_DIR}/eigen")
endif()




include_directories(${EXTERNAL_DIR}/cli11)
include_directories(./src)
include_directories(${EIGEN_BASE})
include_directories(${IGL_BASE}/include)
include_directories(${EXTERNAL_DIR}/VTK)
include_directories(${EXTERNAL_DIR}/Common/Core)
include_directories(./src/MeshOrient)
add_executable(MeshConverter ${SOURCES})
add_subdirectory(src/MeshOrient)
target_link_libraries(MeshConverter MeshOrient)


# try to install to a target 
set(target_bin MeshConverter)
if(hasParent)
set(bin_dest ${hasParent}/bin)
add_subdirectory(src/MeshOrient)
install(TARGETS ${target_bin} DESTINATION "${bin_dest}")
endif() 
