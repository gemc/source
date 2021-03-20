#
#  Evio does not have a decent CMakeLists.txt that confirms to packages, or even builds the code properly.
#  This means it needs a custom build script here
#
# TODO: get Evio to use a proper CMakeLists.txt :-(
# TODO: get a build script here.
#
message(STATUS "Checking for evio")
find_package(EVIO QUIET)
if(NOT EVIO_FOUND)
    message(STATUS "******* EVIO was not found ******* BUT I don't know how to install this yet.. ")
    message(ERROR " Please install EVIO by hand and install in ${CMAKE_INSTALL_PREFIX} ")
    message(ERROR " or set the EVIO_INCLUDES and EVIO_LINRARIES variable from the command line. ")
    #
    # Get: wget https://www.jlab.org/12gev_phys/packages/sources/evio/evio-5.1.tar.gz
    # unpack
    # setup "scons" stuff and build with scons....
endif()