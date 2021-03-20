#
# Evio does not have a proper "Find" script, so we hand hack one here until that is fixed.
#
# TODO: get hipo to use a proper CMakeLists.txt that creates a FindHIPO.cmake :-(
#
# Once done this will define
#  HIPO_FOUND - System has libhipo4.a
#  HIPO_INCLUDE_DIRS - The hipo include directories
#  HIPO_LIBRARIES - The libraries needed to use libhipo4, including libhipo4 itself.
#
find_path(HIPO_INCLUDE_DIR hipoexceptions.h
          PATH_SUFFIXES include hipo4)

#message("HIPO_INCLUDE_DIR = ${HIPO_INCLUDE_DIR}")

find_library(HIPO_LIBRARY NAMES hipo4
             PATH_SUFFIXES lib)

# message("HIPO_LIBRARY = ${HIPO_LIBRARY}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set HIPO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(HIPO  DEFAULT_MSG
                                  HIPO_LIBRARY HIPO_INCLUDE_DIR)

find_package(LZ4 REQUIRED)
if(NOT LZ4_FOUND)
    error("Could not find LZ4")
endif()

mark_as_advanced(HIPO_INCLUDE_DIR HIPO_LIBRARY )

set(HIPO_LIBRARIES ${HIPO_LIBRARY} ${LZ4_LIBRARY})
set(HIPO_INCLUDE_DIRS ${HIPO_INCLUDE_DIR})
