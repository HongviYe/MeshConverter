################################################################################
include(DownloadProject)

if(NOT EXTERNAL_DIR)
    set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/extern)
endif()

if(NOT GITHUB_REPOSITE)
    set(GITHUB_REPOSITE "github.com")
endif()

# Shortcut function
function(meshconverter_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${EXTERNAL_DIR}/${name}
        DOWNLOAD_DIR ${EXTERNAL_DIR}/.cache/${name}
        ${ARGN}
    )
endfunction()

################################################################################

## eigen
function(meshconverter_download_eigen)
    meshconverter_download_project(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.3.7
    )
endfunction()

## MeshOrient
function(meshconverter_download_meshorient)
    meshconverter_download_project(meshorient
        GIT_REPOSITORY https://github.com/xq-meng/MeshOrient.git
    )
endfunction()
