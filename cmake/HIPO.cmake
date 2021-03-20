#
#  Hipo does not have a decent CMakeLists.txt that confirms to packages, or even builds the code properly.
#  This means it needs a custom build script here
#
# TODO: get hipo to use a proper CMakeLists.txt :-(
# TODO: get a build script here.
#
message(STATUS "Checking for HIPO")
find_package(HIPO QUIET)
if(NOT HIPO_FOUND)
    message(STATUS "******* HIPO was not found ******* BUT I don't know how to install this yet.. ")
    message(ERROR " Please install HIPO by hand and install in ${CMAKE_INSTALL_PREFIX} ")
    message(ERROR " or set the HIPO_INCLUDES and HIPO_LIBRARIES variable from the command line. ")
    #
    # Get: wget https://www.jlab.org/12gev_phys/packages/sources/evio/evio-5.1.tar.gz
    # unpack
    # setup "scons" stuff and build with scons....
endif()