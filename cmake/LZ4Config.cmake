#
#  Find the LZ4 library and include files and configure target for it.
#
get_filename_component(_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
find_path(LZ4_INCLUDE_DIR NAMES lz4.h)
find_library(LZ4_LIBRARY NAMES lz4)
if (LZ4_LIBRARY)
    include(CheckCSourceRuns)
    set(CMAKE_REQUIRED_INCLUDES ${LZ4_INCLUDE_DIR})
    set(CMAKE_REQUIRED_LIBRARIES ${LZ4_LIBRARY})
    check_c_source_runs("
#include <lz4.h>
int main() {
  int good = (LZ4_VERSION_MAJOR > 1) ||
    ((LZ4_VERSION_MAJOR == 1) && (LZ4_VERSION_MINOR >= 8));
return !good;
}" LZ4_GOOD_VERSION)
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
        LZ4 DEFAULT_MSG
        LZ4_LIBRARY LZ4_INCLUDE_DIR LZ4_GOOD_VERSION)

if (NOT LZ4_FOUND)
    message(STATUS "No system version of LZ4 found.")
else()
    if(NOT TARGET LZ4)
        add_library(LZ4 INTERFACE IMPORTED)
        set_target_properties(LZ4 PROPERTIES
                              INTERFACE_INCLUDE_DIRECTORIES "${LZ4_INCLUDE_DIR}"
                              )
    endif()
    message(STATUS "Found LZ4: ${LZ4_LIBRARY}")
    set(LZ4_INCLUDE_DIRS ${LZ4_INCLUDE_DIR})
    set(LZ4_LIBRARIES ${LZ4_LIBRARY})
    mark_as_advanced(LZ4_INCLUDE_DIRS LZ4_INCLUDE_DIR LZ4_LIBRARIES LZ4_LIBRARY)
endif (NOT LZ4_FOUND)
