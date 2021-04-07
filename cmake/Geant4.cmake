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
        CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${Geant4_INSTALL_DIR} -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_SYSTEM_EXPAT=OFF -DGEANT4_INSTALL_EXAMPLES=OFF -DGEANT4_USE_QT=ON
        BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} -j4
        UPDATE_COMMAND   ""
    )
    # vis options left out for now
    # -DGEANT4_USE_OPENGL_X11=${GEANT4_USE_OPENGL_X11} -DGEANT4_USE_QT=${GEANT4_USE_QT}
    set(Geant4_DIR ${Geant4_INSTALL_DIR}/lib/Geant4-${Geant4_VERSION} CACHE PATH "Geant4 install dir" FORCE)
    set(Geant4_INCLUDE_DIRS ${Geant4_INSTALL_DIR}/include/Geant4 CACHE PATH "Geant4 include dirs" FORCE)
    set(Geant4_LIBRARIES "Geant4::G4Tree;Geant4::G4FR;Geant4::G4GMocren;Geant4::G4visHepRep;Geant4::G4RayTracer;Geant4::G4VRML;Geant4::G4OpenGL;Geant4::G4gl2ps;Geant4::G4vis_management;Geant4::G4modeling;Geant4::G4interfaces;Geant4::G4persistency;Geant4::G4analysis;Geant4::G4error_propagation;Geant4::G4readout;Geant4::G4physicslists;Geant4::G4run;Geant4::G4event;Geant4::G4tracking;Geant4::G4parmodels;Geant4::G4processes;Geant4::G4digits_hits;Geant4::G4track;Geant4::G4particles;Geant4::G4geometry;Geant4::G4materials;Geant4::G4graphics_reps;Geant4::G4intercoms;Geant4::G4global;Geant4::G4clhep;Geant4::G4expat;Geant4::G4zlib;Geant4::G4UIVisDefinitions")
else()
    # add_custom_target(Geant4) # dummy target
    message(STATUS "Geant4 was found at: ${Geant4_DIR}")
    foreach(_print_item Geant4_INCLUDE_DIRS Geant4_LIBRARIES XercesC_FOUND XercesC_INCLUDE_DIR XercesC_LIBRARY_DEBUG XercesC_LIBRARY_RELEASE)
        message(STATUS "${_print_item}  = ${${_print_item}}")
    endforeach()
endif()
