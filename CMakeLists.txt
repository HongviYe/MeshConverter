cmake_minimum_required(VERSION 3.10)  # 升级最低版本要求，支持现代CMake特性

project(MeshConverter LANGUAGES CXX)

# 设置模块路径
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# 启用C++11标准
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # 禁用编译器特定扩展

# 配置选项
option(MESHORIENT_BUILD_EXECUTABLE "Build executable" OFF)
option(MESHORIENT_DOWNLOAD "Download external dependencies" OFF)
option(USE_GEOGRAM "Use Geogram library" OFF)
option(DOWNLOAD_WITH_CNPM "Use cnpmjs mirror for downloads" OFF)
option(DOWNLOAD_WITH_GITCLONE "Use gitclone mirror for downloads" OFF)

# 设置默认构建类型
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Defaulting to Release build type")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# 配置下载源
if(DOWNLOAD_WITH_CNPM)
  set(GITHUB_REPOSITE "github.com.cnpmjs.org")
elseif(DOWNLOAD_WITH_GITCLONE)
  set(GITHUB_REPOSITE "gitclone.com/github.com")
else()
  set(GITHUB_REPOSITE "github.com")
endif()

# 设置外部依赖目录
if(NOT EXTERNAL_DIR)
  set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/extern)
endif()

# 下载依赖
if(MESHORIENT_DOWNLOAD)
  include(MeshConverterDependencies)
endif()

# 设置依赖路径
set(EIGEN_BASE ${EXTERNAL_DIR}/eigen)
set(IGL_BASE ${EXTERNAL_DIR}/libigl)

# 平台特定配置
if(MSVC)
  add_compile_options(/bigobj)  # 修正MSVC编译选项语法
else()
  # 为非MSVC编译器添加常见优化选项
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 源文件列表
set(SOURCES
  src/meshIO.cpp
  src/meshAlgorithm.cpp
  src/MeshConverter.cpp
  src/remesh.cpp
  src/facet_classification.h
  src/facet_classification.cpp
  
  src/Halfedge/PolyMesh_Base.cpp
  src/Halfedge/PolyMesh.cpp
  src/Halfedge/PolyMeshIterators.cpp
  src/Halfedge/AABB_Tree.cpp
  src/SurfaceHoleFilling.cpp
  src/SurfaceRemesh.cpp
)

# 检查是否作为子项目
get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
  set(EXTERNAL_DIR ${hasParent}/third_party)
  set(EIGEN_BASE ${EXTERNAL_DIR}/eigen)
  set(IGL_BASE ${EXTERNAL_DIR}/libigl)
  message(STATUS "Include path: ${EXTERNAL_DIR}/eigen")
endif()

# 使用现代target_include_directories替代旧方法
add_executable(MeshConverter ${SOURCES})

target_include_directories(MeshConverter SYSTEM PRIVATE  # SYSTEM减少编译警告
  ${EXTERNAL_DIR}/cli11/include/CLI
  ${EIGEN_BASE}
  ${IGL_BASE}/include
  ${EXTERNAL_DIR}/VTK
  ${EXTERNAL_DIR}/Common/Core
)

target_include_directories(MeshConverter PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# 链接Geogram库(如果启用)
if(USE_GEOGRAM)
  target_link_libraries(MeshConverter PRIVATE geogram)
  target_compile_definitions(MeshConverter PUBLIC USE_GEOGRAM)
endif()



# 添加assert包含以修复编译错误
if(MSVC)
else()
target_compile_options(MeshConverter PRIVATE -include cassert)
endif()



# 安装配置
if(hasParent)
  set(bin_dest ${hasParent}/bin)
  install(TARGETS MeshConverter DESTINATION "${bin_dest}")
endif()
