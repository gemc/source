#
# Check for GEANT4
#
# If GEANT4 is not found, add it to the list of packages that will be retreived and
# build using the ExternalProject mechanism.
#
message(STATUS "Checking for Geant4")
find_package(Geant4 QUIET COMPONENTS vis_all ui_all qt gdml)
if(NOT Geant4_FOUND)
    message(STATUS "**********************************************************************")
    message(STATUS "** Geant4 was not found it will be installed (which is slow!)       **")
    message(STATUS "** OR set CMAKE_PREFIX_PATH to the location of Geant4Config.cmake   **")
    message(STATUS "** and rerun cmake. I.e. cmake -DCMAKE_PREFIX_PATH=<loc of G4/lib>  **")
    message(STATUS "**********************************************************************")
    set(Geant4_VERSION 10.6.3 CACHE STRING "Geant4 version" FORCE)
    set(Geant4_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/Geant4)
    add_dependencies(dependencies Geant4)
    externalproject_add(
        Geant4
        GIT_REPOSITORY   "https://github.com/Geant4/geant4"
        GIT_TAG          v${Geant4_VERSION}
        GIT_SHALLOW      ON
        SOURCE_DIR       ${CMAKE_BINARY_DIR}/Geant4
        CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${Geant4_INSTALL_DIR}
                         -DGEANT4_USE_GDML=ON
                         -DGEANT4_USE_QT=ON
                         -DGEANT4_INSTALL_DATA=ON
                         -DGEANT4_USE_SYSTEM_EXPAT=OFF
                         -DGEANT4_INSTALL_EXAMPLES=OFF
                         -DGEANT4_USE_QT=ON
        BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} -j4
        UPDATE_COMMAND   ""
    )
    SET( DEPENDENCIES_RERUN_CMAKE "YES" FORCE )
else()
    # add_custom_target(Geant4) # dummy target
    message(STATUS "Geant4 was found at: ${Geant4_DIR}")
endif()
