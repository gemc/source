#
#  Hipo does not have a decent CMakeLists.txt that confirms to packages, or even builds the code properly.
#  This means it needs a custom build script here
#
#
message(STATUS "Checking for hipo4")
find_package(hipo4 QUIET)
if(NOT hipo4_FOUND)
    message(STATUS "**********************************************************************")
    message(STATUS "***    hipo4 was not found - we will install it automatically.     ***")
    message(STATUS "*** If you do not what that, please set CMAKE_INSTALL_PREFIX       ***")
    message(STATUS "***********************************************************************")
    #
    # 	https://github.com/gavalian/hipo4.git
    #
    set(HIPO_VERSION master CACHE STRING "hipo version" FORCE)
    externalproject_add(
            hipo4_external
            GIT_REPOSITORY https://github.com/mholtrop/hipo.git
            GIT_TAG        ${HIPO_VERSION}
            SOURCE_DIR     ${CMAKE_BINARY_DIR}/hipo
            INSTALL_DIR    ${CMAKE_INSTALL_PREFIX}
            CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                           -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                           -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            BUILD_COMMAND  ${CMAKE_MAKE_PROGRAM}
    )

    # We need so see where the LZ4 library will come from.
    find_package(LZ4 QUIET CONFIG PATHS ${PROJECT_SOURCE_DIR}/cmake )
    if(NOT LZ4_FOUND)
        message(STATUS "***********************************************************************")
        message(STATUS "* We did not find a system LZ4 package -- We will build this locally. *")
        message(STATUS "***********************************************************************")
        set(LZ4_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include)
        set(LZ4_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/liblz4.a)
        set(LZ4_LIBRARIES_STATIC ${CMAKE_INSTALL_PREFIX}/lib/liblz4.a)
        set(LZ4_LIBRARIES_SHARED ${CMAKE_INSTALL_PREFIX}/lib/liblz4${CMAKE_SHARED_LIBRARY_SUFFIX})
        file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
        file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)
    endif(NOT LZ4_FOUND)

    add_library(hipo4 SHARED IMPORTED)
    add_dependencies(dependencies hipo4)
    add_dependencies(hipo4 hipo4_external)
    file(MAKE_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/include/hipo4)
    set_target_properties(hipo4 PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include;${CMAKE_INSTALL_PREFIX}/include/hipo4"
                          INTERFACE_LINK_LIBRARIES "${LZ4_LIBRARIES_SHARED}"
                          )
    set_property(TARGET hipo4 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(hipo4 PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libhipo4${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          IMPORTED_SONAME_RELEASE "@rpath/libhipo4${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          )

    # Everything again for _static
    add_library(hipo4_static STATIC IMPORTED)
    #   add_dependencies(dependencies hipo4_static)    ## Putting these here again causes a LZ4 to be build twice?
    #   add_dependencies(hipo4_static hipo4_external)
    set_target_properties(hipo4_static PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/hipo4"
                          INTERFACE_LINK_LIBRARIES "${LZ4_LIBRARIES_STATIC}"
                          )
    set_property(TARGET hipo4_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(hipo4_static PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libhipo4.a"
                          )

else()
    get_target_property(_hipo4_found_dir hipo4 IMPORTED_LOCATION_RELEASE)
    message(STATUS "hipo4 found: ${_hipo4_found_dir}")

endif()
