message(STATUS "Checking for sqlite3")
find_package(SQLite3)
if(NOT Qt5_FOUND)
    message(FATAL_ERROR "********* NO System SQLite3 found, and no instructions to build it. ***********")

    #
    # TODO: Add 'get qt5' directions if qt5 is not found.
    #
endif()