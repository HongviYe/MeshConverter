set(CMAKE_CXX_STANDARD 11)
project(MeshOrient)

# options
option(MESHORIENT_BUILD_EXECUTABLE "Build executable." OFF)

set(SOURCES
./MeshOrient.h
./triMesh.h
./triMesh.cpp
./MeshOrient.cpp)
include_directories(${EIGEN_BASE})
include_directories(${IGL_BASE}/include)
message("the include eigen path =${EIGEN_BASE}")
# target_include_directories(MeshOrient PUBLIC ./)
add_library(MeshOrient STATIC ${SOURCES})
if(MSVC)
else()
target_link_libraries(MeshOrient PUBLIC pthread)
endif()

if(MESHORIENT_BUILD_EXECUTABLE)
    add_executable(MeshOrientExe MeshOrientExe.cpp)
    target_link_libraries(MeshOrientExe PUBLIC MeshOrient)
    set_target_properties(MeshOrientExe PROPERTIES OUTPUT_NAME MeshOrient)
endif()
