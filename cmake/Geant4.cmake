#
# Check for GEANT4
#
# If GEANT4 is not found, add it to the list of packages that will be retreived and
# build using the ExternalProject mechanism.
#
message(STATUS "Checking for Geant4")
find_package(Geant4 QUIET COMPONENTS vis_all ui_all qt gdml)
if(XercesC_LIBRARY AND NOT XercesC_LIBRARY_RELEASE)            # An issue with some XercesC installs that do not set the _RELEAST and _DEBUG
    set(XercesC_LIBRARY_RELEASE ${XercesC_LIBRARY})
endif()
if(NOT Geant4_FOUND)

    set(XercesC_VERSION 3.2.3 CACHE STRING "XercesC version" FORCE)
    set(Geant4_VERSION 10.6.3 CACHE STRING "Geant4 version" FORCE)
    set(Geant4_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/Geant4)
    #
    # The GDML component of GEANT4 depends on XercesC. Look for it, and if not found, build it.
    #
    find_package(XercesC QUIET)
    if(NOT XercesC_FOUND)
        message(STATUS "XercesC was not found and will be installed before GEANT4")
        set(XercesC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
        externalproject_add(
                XercesC
                URL                  "https://downloads.apache.org/xerces/c/3/sources/xerces-c-${XercesC_VERSION}.tar.gz"
                SOURCE_DIR           ${CMAKE_BINARY_DIR}/XercesC
                CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} -j4
        )
    else()
        add_custom_target(XercesC) # dummy target
        message(STATUS "XercesC library was found at: ${XercesC_LIBRARY}")
    endif()

    message(STATUS "**********************************************************************")
    message(STATUS "** Geant4 was not found it will be installed (which is slow!)       **")
    message(STATUS "** OR set CMAKE_PREFIX_PATH to the location of Geant4Config.cmake   **")
    message(STATUS "** and rerun cmake. I.e. cmake -DCMAKE_PREFIX_PATH=<loc of G4/lib>  **")
    message(STATUS "**********************************************************************")

    add_dependencies(dependencies Geant4)

    externalproject_add(
        Geant4
        GIT_REPOSITORY   "https://github.com/Geant4/geant4"
        GIT_TAG          v${Geant4_VERSION}
        GIT_SHALLOW      ON
        SOURCE_DIR       ${CMAKE_BINARY_DIR}/Geant4
        CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${Geant4_INSTALL_DIR}
                         -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                         -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                         -DGEANT4_USE_GDML=ON
                         -DGEANT4_USE_QT=ON
                         -DGEANT4_INSTALL_DATA=ON
                         -DGEANT4_USE_SYSTEM_EXPAT=OFF
                         -DGEANT4_INSTALL_EXAMPLES=OFF
                         -DGEANT4_USE_QT=ON
        BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} -j4
        UPDATE_COMMAND   ""
    )

    add_dependencies(Geant4 XercesC)

    SET( DEPENDENCIES_RERUN_CMAKE "YES" FORCE )
else()
    # add_custom_target(Geant4) # dummy target
    message(STATUS "Geant4 was found at: ${Geant4_DIR}")
endif()
