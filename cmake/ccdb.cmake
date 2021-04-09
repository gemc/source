#
#  CCDB does not have a decent CMakeLists.txt -- this ought to be fixed.
#  None of the versions that are on GIT actually build correctly either if you use scons and Python 3
#  This means it needs a custom build script here
#
#
message(STATUS "Checking for ccdb")
find_package(ccdb QUIET)
if(NOT ccdb_FOUND)
    message(STATUS "********************************************************************")
    message(STATUS "* Warning: ccdb was not found.                                     *")
    message(STATUS "* We will get ccdb from GitHub and install it.                     *")
    message(STATUS "* If this is not what you want, make sure ccdb can be found        *")
    message(STATUS "* set CMAKE_PREFIX_PATH to include the ccdbConfig.cmake file       *")
    message(STATUS "********************************************************************")

    include(MYSQL) # We need to make sure MySQL is found.


    set(CCDB_VERSION cmake_v1.07 CACHE STRING "ccdb version" FORCE)
    externalproject_add(
            ccdb_external
            GIT_REPOSITORY   "https://github.com/mholtrop/ccdb"
            GIT_TAG          ${CCDB_VERSION}
            SOURCE_DIR       ${CMAKE_BINARY_DIR}/ccdb
            INSTALL_DIR      ${CMAKE_INSTALL_PREFIX}
            CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} -j8
            UPDATE_COMMAND   ""
    )
    add_library(ccdb_sqlite SHARED IMPORTED)
    add_library(ccdb SHARED IMPORTED)
    # Make (empty) include directories fir ccdb. These will be filled during "make"
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/ccdb)
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/ccdb/CCDB)
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/ccdb/SQLite)

    # Setup dependencies.
    add_dependencies(dependencies ccdb)
    add_dependencies(ccdb ccdb_external)

    set_target_properties(ccdb_sqlite PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/ccdb;${CMAKE_INSTALL_PREFIX}/include/ccdb/SQLite"
                          )

    set_target_properties(ccdb PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/ccdb;${CMAKE_INSTALL_PREFIX}/include/ccdb/CCDB"
                          INTERFACE_LINK_LIBRARIES "${MYSQL_LIBRARIES};ccdb_sqlite"
                          )

    set_property(TARGET ccdb_sqlite APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(ccdb_sqlite PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libccdb_sqlite${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          IMPORTED_SONAME_RELEASE "@rpath/libccdb_sqlite${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          )
    set_property(TARGET ccdb APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(ccdb PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libccdb${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          IMPORTED_SONAME_RELEASE "@rpath/libccdb${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          )

    add_library(ccdb_sqlite_static SHARED IMPORTED)
    add_library(ccdb_static SHARED IMPORTED)
    set_target_properties(ccdb_sqlite_static PROPERTIES
                          IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/ccdb"
                          IMPORTED_CONFIGURATIONS RELEASE
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/ccdb"
                          )
    set_target_properties(ccdb_static PROPERTIES
                          IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/ccdb"
                          IMPORTED_CONFIGURATIONS RELEASE
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/ccdb;${CMAKE_INSTALL_PREFIX}/include/ccdb/CCDB"
                          INTERFACE_LINK_LIBRARIES "${MYSQL_LIBRARIES};ccdb_sqlite"
                          )

    set_property(TARGET ccdb_sqlite_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(ccdb_sqlite_static PROPERTIES
                          IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libccdb_sqlite.a"
                          )
    set_property(TARGET ccdb_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(ccdb_static PROPERTIES
                          IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libccdb.a"
                          )


endif()
