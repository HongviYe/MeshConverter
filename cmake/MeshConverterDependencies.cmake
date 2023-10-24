###############################################################################
# CMake download helpers
###############################################################################

# download external dependencies
include(MeshConverterDownloadExternal)

###############################################################################
# Required dependencies
###############################################################################
if(USE_GEOGRAM)
    include(geogram)
endif()
# eigen
if(NOT EIGEN_BASE)
    meshconverter_download_eigen()
endif()

# libigl
if(NOT IGL_BASE)
    meshorient_download_libigl()
endif()
