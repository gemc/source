#
message(STATUS "Checking for mlibrary")
find_package(mlibrary QUIET)
if(NOT mlibrary_FOUND)
    message(STATUS "******* mlibrary was not found ******* BUT I don't know how to install this yet.. ")
    message(ERROR " Please install mlibrary by hand. ")
    # git clone https://github.com/gemc/mlibrary.git
    # ...
endif()