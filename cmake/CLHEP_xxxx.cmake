#
#
#
message(STATUS "Checking for CLHEP")
find_package(CLHEP)
if(NOT CLHEP_FOUND)
    message(ERROR "CLHEP was not found and I do not know how to build it.")
endif()
