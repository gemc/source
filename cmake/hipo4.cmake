#
#  Hipo does not have a decent CMakeLists.txt that confirms to packages, or even builds the code properly.
#  This means it needs a custom build script here
#
#
message(STATUS "Checking for hipo4")
find_package(hipo4 QUIET)
if(NOT hipo4_FOUND)
    message(STATUS "**********************************************************************")
    message(STATUS "******* hipo4 was not found                                     *******")
    message(STATUS " Please install HIPO by hand and install in ${CMAKE_INSTALL_PREFIX} ")
    message(STATUS "***********************************************************************")
    #
    # 	https://github.com/gavalian/hipo4.git
    #
    find_package(LZ4 REQUIRED)
    if(NOT LZ4_FOUND)
        error("Could not find LZ4")
    endif()

    add_dependencies(dependencies HIPO)
    externalproject_add(
            HIPO
            GIT_REPOSITORY https://github.com/mholtrop/hipo.git
            GIT_TAG        master
            SOURCE_DIR     ${CMAKE_BINARY_DIR}/hipo
            INSTALL_DIR    ${CMAKE_INSTALL_PREFIX}
            CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            BUILD_COMMAND  ${CMAKE_MAKE_PROGRAM} install
    )
    set(HIPO_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libhipo4.a CACHE PATH "hipo library" FORCE)
    set(HIPO_LIBRARIES ${HIPO_LIBRARY} ${LZ4_LIBRARY})
    set(HIPO_INCLUDE_DIRS ${HIPO_INCLUDE_DIR})
endif()