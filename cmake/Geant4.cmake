#
# Check for GEANT4
#
# If GEANT4 is not found, add it to the list of packages that will be retreived and
# build using the ExternalProject mechanism.
#

# Check for the external CLHEP.
# The Build-in CLHEP with GEANT4 seems to not be sufficient for GEMC. (Missing RandomGaussT.h)
message(STATUS "Checking for CLHEP")
find_package(CLHEP QUIET)
if(NOT CLHEP_FOUND)
    set(CLHEP_VERSION 2.4.4.1 CACHE STRING "CLHEP version" FORCE)

    string(REPLACE "." "_" CLHEP_TAG "CLHEP_${CLHEP_VERSION}")
    message(STATUS "***********************************************************")
    message(STATUS "* CLHEP was not found and will be installed before GEANT4 *")
    message(STATUS "* CLHEP version ${CLHEP_VERSION}, tag: ${CLHEP_TAG}              *")
    message(STATUS "***********************************************************")

    externalproject_add(
            CLHEP
            GIT_REPOSITORY https://gitlab.cern.ch/CLHEP/CLHEP.git
            GIT_TAG   ${CLHEP_TAG}
            SOURCE_DIR ${CMAKE_BINARY_DIR}/CLHEP
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM}
    )
    set(CLHEP_ROOT_DIR ${CMAKE_INSTALL_PREFIX})

    SET( DEPENDENCIES_RERUN_CMAKE "YES" FORCE )
else()
    add_custom_target(CLHEP) # dummy target
    message(STATUS "CLHEP library was found at: ${CLHEP_INCLUDE_DIRS}, ${CLHEP_LIBRARIES}")
    get_filename_component(CLHEP_ROOT_DIR ${CLHEP_DIR} DIRECTORY)
endif()

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
        message(STATUS "*************************************************************")
        message(STATUS "* XercesC was not found and will be installed before GEANT4 *")
        message(STATUS "* XercesC version: ${XercesC_VERSION}                                    *")
        message(STATUS "*************************************************************")
        set(XercesC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
        externalproject_add(
                XercesC
                URL                  "https://downloads.apache.org/xerces/c/3/sources/xerces-c-${XercesC_VERSION}.tar.gz"
                SOURCE_DIR           ${CMAKE_BINARY_DIR}/XercesC
                CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM}
        )
        # Somewhat troublesome, XercesC insists on installing in lib64 on Linux.
        if(APPLE)
            set(XercesC_LIBRARY ${XercesC_INSTALL_DIR}/lib/libxerces-c${CMAKE_SHARED_LIBRARY_SUFFIX})
        else()
            set(XercesC_LIBRARY ${XercesC_INSTALL_DIR}/lib64/libxerces-c${CMAKE_SHARED_LIBRARY_SUFFIX})
        endif()
        set(XercesC_INCLUDE_DIR	${XercesC_INSTALL_DIR}/include)

    else()
        add_custom_target(XercesC) # dummy target
        message(STATUS "XercesC library was found at: ${XercesC_LIBRARY} and ${XercesC_INCLUDE_DIR}")
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
                         -DXercesC_LIBRARY=${XercesC_LIBRARY}
                         -DXercesC_INCLUDE_DIR=${XercesC_INCLUDE_DIR}
                         -DCLHEP_ROOT_DIR=${CLHEP_ROOT_DIR}
                         -DGEANT4_USE_QT=ON
                         -DGEANT4_INSTALL_DATA=ON
                         -DGEANT4_USE_SYSTEM_EXPAT=OFF
                         -DGEANT4_INSTALL_EXAMPLES=OFF
                         -DGEANT4_USE_QT=ON
        BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM}
        UPDATE_COMMAND   ""
    )
    add_dependencies(Geant4 XercesC)
    SET( DEPENDENCIES_RERUN_CMAKE "YES" FORCE )
else()
    add_custom_target(Geant4) # dummy target
    message(STATUS "Geant4 was found at: ${Geant4_DIR}")
endif()

add_dependencies(Geant4 CLHEP)

