#
# Evio does not have a proper "Find" script, so we hand hack one here until that is fixed.
#
# TODO: get Evio to use a proper CMakeLists.txt that creates a FindEvio.cmake :-(
#
# Once done this will define
#  EVIO_FOUND - System has libevio.a
#  EVIO_INCLUDE_DIRS - The evio include directories
#  EVIO_LIBRARIES - The libraries needed to use libevio, including libevio itself.
#
find_path(EVIO_INCLUDE_DIR evio.h
          PATH_SUFFIXES include)

#find_path(EVIO_UTIL_INCLUDE_DIR evioUtil.hxx
#          PATH_SUFFIXES include)

#message("EVIO_INCLUDE_DIR = ${EVIO_INCLUDE_DIR}")

find_library(EVIO_LIBRARY NAMES evio
             PATH_SUFFIXES lib)
find_library(EVIO_LIBRARYXX NAMES evioxx
             PATH_SUFFIXES lib)

# message("EVIO_LIBRARYXX = ${EVIO_LIBRARYXX}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EVIO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(EVIO  DEFAULT_MSG
                                  EVIO_LIBRARY EVIO_INCLUDE_DIR)

mark_as_advanced(EVIO_INCLUDE_DIR EVIO_LIBRARY EVIO_LIBRARYXX)

set(EVIO_LIBRARIES ${EVIO_LIBRARY} ${EVIO_LIBRARYXX})
set(EVIO_INCLUDE_DIRS ${EVIO_INCLUDE_DIR})
