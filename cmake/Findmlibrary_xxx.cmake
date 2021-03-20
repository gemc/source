#
# Mlibrary does not have a proper "Find" script, so we hand hack one here until that is fixed.
#
# TODO: get mlibrary to use a proper CMakeLists.txt that creates a Findmlibrary_xxx.cmake :-(
#
# Once done this will define
#  mlibrary_FOUND - System has libmlibrary.a
#  mlibrary_INCLUDE_DIRS - The mlibrary include directories
#  mlibrary_LIBRARIES - The libraries needed to use libmlibrary, including libmlibrary itself.
#
if(DEFINED ENV{MLIBRARY})
    set(mlibrary_HOME $ENV{MLIBRARY})
endif()

#message(STATUS "USING mlibrary_HOME = ${mlibrary_HOME}")
find_path(mlibrary_INCLUDE_DIR goptions.h
          HINTS ${mlibrary_HOME}
          PATH_SUFFIXES include mlibrary)

# message("mlibrary_INCLUDE_DIR = ${mlibrary_INCLUDE_DIR}")

find_library(mlibrary_LIBRARY NAMES
              mlibrary
             HINTS ${mlibrary_HOME}/lib )

# message("mlibrary_LIBRARY = ${mlibrary_LIBRARY}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set mlibrary_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(mlibrary  DEFAULT_MSG
                                  mlibrary_LIBRARY mlibrary_INCLUDE_DIR)

mark_as_advanced(mlibrary_INCLUDE_DIR mlibrary_LIBRARY )

set(mlibrary_LIBRARIES ${mlibrary_LIBRARY} )
set(mlibrary_INCLUDE_DIRS ${mlibrary_INCLUDE_DIR} )
