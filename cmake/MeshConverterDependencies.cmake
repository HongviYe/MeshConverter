###############################################################################
# CMake download helpers
###############################################################################

# download external dependencies
include(MeshConverterDownloadExternal)

###############################################################################
# Required dependencies
###############################################################################

# eigen
if(NOT EIGEN_ROOT)
    meshconverter_download_eigen()
endif()
