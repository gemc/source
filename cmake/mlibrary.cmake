#
message(STATUS "Checking for mlibrary")
if(NOT Geant4_FOUND)
    message(STATUS "*********************************************************")
    message(STATUS "* mlibrary will not be build or configured until GEANT4 *")
    message(STATUS "* is found. Make sure GEANT4 builds, then rerun cmake   *")
    message(STATUS "*********************************************************")
else()
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
                                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
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
                              INTERFACE_LINK_LIBRARIES "${Geant4_LIBRARIES};Qt5::Xml;Qt5::Widgets;Qt5::Sql;cadmesh"
                              )
        set_property(TARGET mlibrary APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
        set_target_properties(mlibrary PROPERTIES
                              IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libmlibrary${CMAKE_SHARED_LIBRARY_SUFFIX}"
                              IMPORTED_SONAME_RELEASE "@rpath/libmlibrary${CMAKE_SHARED_LIBRARY_SUFFIX}"
                              )

        add_library(mlibrary_static STATIC IMPORTED)
        add_dependencies(dependencies mlibrary_static)
        add_dependencies(mlibrary_static mlibrary_external)


        set_target_properties(mlibrary_static PROPERTIES
                              INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/mlibrary;${CMAKE_INSTALL_PREFIX}/include"
                              INTERFACE_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib"
                              INTERFACE_LINK_LIBRARIES "${Geant4_LIBRARIES};Qt5::Xml;Qt5::Widgets;Qt5::Sql;cadmesh"
                              )
        set_property(TARGET mlibrary_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
        set_target_properties(mlibrary_static PROPERTIES
                              IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
                              IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libmlibrary_static.a"
                              )


    endif()
endif()
