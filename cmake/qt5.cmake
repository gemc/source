message(STATUS "Checking for qt5")
find_package(Qt5 COMPONENTS Core Gui Widgets OpenGL Xml Sql REQUIRED)
if(NOT Qt5_FOUND)
    message("********* NO System QT5 found, and no instructions to build it. ***********")
#
# TODO: Add 'get qt5' directions if qt5 is not found.
#
endif()