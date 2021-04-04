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
    add_dependencies(dependencies evio)
    externalproject_add(
            evio
            GIT_REPOSITORY   https://github.com/mholtrop/evio-5.1.git
            GIT_TAG          main
            SOURCE_DIR       ${CMAKE_BINARY_DIR}/evio
            INSTALL_DIR    ${CMAKE_INSTALL_PREFIX}
            CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            BUILD_COMMAND    ${CMAKE_MAKE_PROGRAM} install
    )
endif()