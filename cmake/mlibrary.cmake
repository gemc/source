#
message(STATUS "Checking for mlibrary")
find_package(mlibrary QUIET)
if(NOT mlibrary_FOUND)
    message(STATUS "**********************************************************************")
    message(STATUS "** mlibrary was not found, so it will be installed.                 **")
    message(STATUS "** OR set CMAKE_PREFIX_PATH to the location of mlibraryConfig.cmake **")
    message(STATUS "** and rerun cmake.                                                 **")
    message(STATUS "**********************************************************************")
    set(MLIBRARY_VERSION iss197 CACHE STRING "mlibrary version" FORCE)
    externalproject_add(
            mlibrary_external
            GIT_REPOSITORY   "https://github.com/mholtrop/mlibrary.git"
            GIT_TAG          ${MLIBRARY_VERSION}
            SOURCE_DIR       ${CMAKE_BINARY_DIR}/mlibrary
            INSTALL_DIR      ${CMAKE_INSTALL_PREFIX}
            CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} install
            UPDATE_COMMAND   ""
    )
    add_library(mlibrary SHARED IMPORTED)
    add_dependencies(dependencies mlibrary)
    add_dependencies(mlibrary mlibrary_external)
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/mlibrary)
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/assimp)
    set_target_properties(mlibrary PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/mlibrary;${CMAKE_INSTALL_PREFIX}/include"
                          INTERFACE_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib"
                          INTERFACE_LINK_LIBRARIES "Geant4::G4Tree;Geant4::G4FR;Geant4::G4GMocren;Geant4::G4visHepRep;Geant4::G4RayTracer;Geant4::G4VRML;Geant4::G4vis_management;Geant4::G4modeling;Geant4::G4interfaces;Geant4::G4persistency;Geant4::G4analysis;Geant4::G4error_propagation;Geant4::G4readout;Geant4::G4physicslists;Geant4::G4run;Geant4::G4event;Geant4::G4tracking;Geant4::G4parmodels;Geant4::G4processes;Geant4::G4digits_hits;Geant4::G4track;Geant4::G4particles;Geant4::G4geometry;Geant4::G4materials;Geant4::G4graphics_reps;Geant4::G4intercoms;Geant4::G4global;Geant4::G4clhep;Geant4::G4expat;Geant4::G4zlib;Geant4::G4UIVisDefinitions;Qt5::Xml;Qt5::Widgets;Qt5::Sql;cadmesh"
                          )
    set_property(TARGET mlibrary APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(mlibrary PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libmlibrary.dylib"
                          IMPORTED_SONAME_RELEASE "@rpath/libmlibrary.dylib"
                          )

    add_library(mlibrary_static STATIC IMPORTED)
    add_dependencies(dependencies mlibrary_static)
    add_dependencies(mlibrary_static mlibrary_external)


    set_target_properties(mlibrary_static PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/mlibrary;${CMAKE_INSTALL_PREFIX}/include"
                          INTERFACE_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib"
                          INTERFACE_LINK_LIBRARIES "Geant4::G4Tree;Geant4::G4FR;Geant4::G4GMocren;Geant4::G4visHepRep;Geant4::G4RayTracer;Geant4::G4VRML;Geant4::G4vis_management;Geant4::G4modeling;Geant4::G4interfaces;Geant4::G4persistency;Geant4::G4analysis;Geant4::G4error_propagation;Geant4::G4readout;Geant4::G4physicslists;Geant4::G4run;Geant4::G4event;Geant4::G4tracking;Geant4::G4parmodels;Geant4::G4processes;Geant4::G4digits_hits;Geant4::G4track;Geant4::G4particles;Geant4::G4geometry;Geant4::G4materials;Geant4::G4graphics_reps;Geant4::G4intercoms;Geant4::G4global;Geant4::G4clhep;Geant4::G4expat;Geant4::G4zlib;Geant4::G4UIVisDefinitions;Qt5::Xml;Qt5::Widgets;Qt5::Sql;cadmesh"
                          )
    set_property(TARGET mlibrary_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(mlibrary_static PROPERTIES
                          IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libmlibrary_static.a"
                          )


endif()
