#
# Use the https://github.com/mholtrop/evio-5.1.git version of EVIO so that evioConfig.cmake and eviocxxConfig.cmake
# are found correctly.  The evio target is a dependency of eviocxx, so we only need to find eviocxx
#
message(STATUS "Checking for eviocxx")
find_package(eviocxx QUIET)  # NO REQUIRED. If it is not found, we will install it.

if(NOT eviocxx_FOUND)   # We did not find it, so we will build it ourselves.
    message(STATUS "************************************************************************")
    message(STATUS "**  EVIO was not found. I will try to build it.                       **")
    message(STATUS "**  To avoid the build, place the EVIO library and include file in    **")
    message(STATUS "**  ${CMAKE_INSTALL_PREFIX}     **")
    message(STATUS "************************************************************************")
    #
    # Get: wget https://www.jlab.org/12gev_phys/packages/sources/evio/evio-5.1.tar.gz
    # unpack
    # Add in a CMakeLists.txt :-)
    # Run cmake, make, make install
    #
    add_dependencies(dependencies eviocxx)
    set(EVIO_VERSION main CACHE STRING "evio version" FORCE)
    externalproject_add(
            eviocxx_external
            GIT_REPOSITORY   https://github.com/mholtrop/evio-5.1.git
            GIT_TAG          ${EVIO_VERSION}
            SOURCE_DIR       ${CMAKE_BINARY_DIR}/evio
            INSTALL_DIR      ${CMAKE_INSTALL_PREFIX}
            CMAKE_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                             -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM}
            UPDATE_COMMAND   ""
    )
    add_library(eviocxx SHARED IMPORTED)
    add_library(evio    SHARED IMPORTED)
    add_dependencies(dependencies eviocxx)
    add_dependencies(eviocxx eviocxx_external)
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/eviocxx)
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include/evio)
    set_target_properties(evio PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/evio"
                          )
    set_target_properties(eviocxx PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/eviocxx"
                          INTERFACE_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib"
                          INTERFACE_LINK_LIBRARIES "evio"
                          )

    set_property(TARGET evio APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(evio PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libevio${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          IMPORTED_SONAME_RELEASE "@rpath/libevio${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          )

    set_property(TARGET eviocxx APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(eviocxx PROPERTIES
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libeviocxx${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          IMPORTED_SONAME_RELEASE "@rpath/libeviocxx${CMAKE_SHARED_LIBRARY_SUFFIX}"
                          )

    add_library(eviocxx_static STATIC IMPORTED)
    add_library(evio_static    STATIC IMPORTED)

    set_target_properties(evio_static PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/evio"
                          )
    set_target_properties(eviocxx_static PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/include/eviocxx"
                          INTERFACE_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib"
                          INTERFACE_LINK_LIBRARIES "evio"
                          )

    set_property(TARGET evio_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(evio_static PROPERTIES
                          IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libevio.a"
                          )

    set_property(TARGET eviocxx_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(eviocxx_static PROPERTIES
                          IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
                          IMPORTED_LOCATION_RELEASE "${CMAKE_INSTALL_PREFIX}/lib/libeviocxx.a"
                          )
endif()