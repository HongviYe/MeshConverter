set(CMAKE_CXX_STANDARD 11)
project(MeshOrient)

# options
option(MESHORIENT_BUILD_EXECUTABLE "Build executable." OFF)

set(SOURCES
./MeshOrient.h
./triMesh.h
./triMesh.cpp
./MeshOrient.cpp)
add_library(MeshOrient STATIC ${SOURCES})
target_include_directories(MeshOrient PUBLIC ${EIGEN_BASE})
target_include_directories(MeshOrient PUBLIC ${IGL_BASE}/include)
# target_include_directories(MeshOrient PUBLIC ./)

if(MSVC)
else()
target_link_libraries(MeshOrient PUBLIC pthread)
endif()

if(MESHORIENT_BUILD_EXECUTABLE)
    add_executable(MeshOrientExe src/MeshOrientExe.cpp)
    target_link_libraries(MeshOrientExe PUBLIC MeshOrient)
    set_target_properties(MeshOrientExe PROPERTIES OUTPUT_NAME MeshOrient)
endif()
