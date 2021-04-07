# Check for Qt5
# This must be called AFTER Geant4, because Geant4's version must be the *same* as included here.
#
if(NOT Qt5_FOUND)  # Check if found already.
    if(Qt5Core_DIR)
        get_filename_component(Qt5_CMAKE_DIR ${Qt5Core_DIR} DIRECTORY)
        message(STATUS "Checking for qt5 in ${Qt5_CMAKE_DIR}")
        find_package(Qt5 COMPONENTS Core Gui Widgets OpenGL Xml Sql REQUIRED PATHS ${Qt5_CMAKE_DIR})
    else()
        message(STATUS "Checking for qt5")
        find_package(Qt5 COMPONENTS Core Gui Widgets OpenGL Xml Sql REQUIRED)
    endif()
    if(NOT Qt5_FOUND)
        message(FATAL_ERROR "********* NO System QT5 found, and no instructions to build it. ***********")
    #
    # TODO: Add 'get qt5' directions if qt5 is not found.
    #
    endif()
endif()
