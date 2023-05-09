###############################################################################
# CMake download helpers
###############################################################################

# download external dependencies
include(MeshConverterDownloadExternal)

###############################################################################
# Required dependencies
###############################################################################

# eigen
if(NOT EIGEN_BASE)
    meshconverter_download_eigen()
    set(EIGEN_BASE ${EXTERNAL_DIR}/eigen)
endif()

# libigl
if(NOT IGL_BASE)
    meshorient_download_libigl()
    set(IGL_BASE ${EXTERNAL_DIR}/libigl)
endif()
