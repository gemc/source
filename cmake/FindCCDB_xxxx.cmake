#
# CCDB does not have a proper "Find" script, so we hand hack one here until that is fixed.
#
# TODO: get CCDB to use a proper CMakeLists.txt that creates a FindCCDB.cmake :-(
#
# Once done this will define
#  CCDB_FOUND - System has libccdb.a
#  CCDB_INCLUDE_DIRS - The ccdb include directories
#  CCDB_LIBRARIES - The libraries needed to use libccdb, including libccdb itself.
#
if(DEFINED ENV{CCDB_HOME})
    set(CCDB_HOME $ENV{CCDB_HOME})
endif()

#message(STATUS "USING CCDB_HOME = ${CCDB_HOME}")
find_path(CCDB_INCLUDE_DIR CCDB/CCDBError.h
          HINTS ${CCDB_HOME}
          PATH_SUFFIXES include CCDB)

message("CCDB_INCLUDE_DIR = ${CCDB_INCLUDE_DIR}")

find_library(CCDB_LIBRARY NAMES ccdb CCDB
             HINTS ${CCDB_HOME}/lib )

message("CCDB_LIBRARY = ${CCDB_LIBRARY}")

include(FindPackageHandleStandardArgs)


# handle the QUIETLY and REQUIRED arguments and set CCDB_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CCDB  DEFAULT_MSG
                                  CCDB_LIBRARY CCDB_INCLUDE_DIR)

mark_as_advanced(CCDB_INCLUDE_DIR CCDB_LIBRARY )

set(CCDB_LIBRARIES ${CCDB_LIBRARY} )
set(CCDB_INCLUDE_DIRS ${CCDB_INCLUDE_DIR} )
