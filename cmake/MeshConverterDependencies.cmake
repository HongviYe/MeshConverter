###############################################################################
# CMake download helpers
###############################################################################

# download external dependencies
include(MeshConverterDownloadExternal)

###############################################################################
# Required dependencies
###############################################################################

# eigen
if(NOT TARGET Eigen3::Eigen)
    meshconverter_download_eigen()
endif()
