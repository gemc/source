message(STATUS "Checking for XercesC")
if(XercesC_FOUND)
    message(STATUS "GEANT4 already seems to have found Xerces-C for with version ${XercesC_VERSION}" )
else(XercesC_FOUND)
    find_package(XercesC QUIET)
    if(NOT XercesC_FOUND)
        message(STATUS "XercesC was not found and will be installed")
        set(XercesC_VERSION 3.2.3 CACHE STRING "XercesC version" FORCE) # Is this necessary?
        set(XercesC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/XercesC)
        externalproject_add(
            XercesC
            URL                  "https://downloads.apache.org/xerces/c/3/sources/xerces-c-${XercesC_VERSION}.tar.gz"
            SOURCE_DIR           ${CMAKE_BINARY_DIR}/XercesC
            CONFIGURE_COMMAND    ./configure --enable-shared --disable-static --prefix=${XercesC_INSTALL_DIR}
            DOWNLOAD_NO_PROGRESS ON
            BUILD_IN_SOURCE      ON
        )
        set(XercesC_LIBRARY ${XercesC_INSTALL_DIR}/lib/libxerces-c${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE FILEPATH "XercesC libraries" FORCE)
        set(XercesC_INCLUDE_DIR ${XercesC_INSTALL_DIR}/include CACHE PATH "XercesC include dir" FORCE)
    else()
        add_custom_target(XercesC) # dummy target
        message(STATUS "XercesC library was found at: ${XercesC_LIBRARY}")
    endif()
endif(XercesC_FOUND)
add_dependencies(dependencies XercesC)
