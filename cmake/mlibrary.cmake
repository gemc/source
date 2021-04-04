#
message(STATUS "Checking for mlibrary")
find_package(mlibrary QUIET)
if(NOT mlibrary_FOUND)
    message(STATUS "**********************************************************************")
    message(STATUS "** mlibrary was not found, so it will be installed.                 **")
    message(STATUS "** OR set CMAKE_PREFIX_PATH to the location of mlibraryConfig.cmake **")
    message(STATUS "** and rerun cmake.                                                 **")
    message(STATUS "**********************************************************************")
    add_dependencies(dependencies mlibrary)
    externalproject_add(
            mlibrary
            GIT_REPOSITORY   "https://github.com/mholtrop/mlibrary.git"
            GIT_TAG          develop
            SOURCE_DIR       ${CMAKE_BINARY_DIR}/mlibrary
            CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} -j4 install
            UPDATE_COMMAND   ""
    )
endif()
