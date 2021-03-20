##########################################################
# CMake configuration for Gemc                           #
#                                                        #
# @author Maurik Holtrop, University of New Hampshire    #
#                                                        #
##########################################################

# minimum version of CMake, which is determined by Geant4's requirements
cmake_minimum_required(VERSION 3.8)

project(GEMC VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    set(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
        set(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
        message("CMAKE_INSTALL_RPATH set to ${CMAKE_INSTALL_RPATH} ")

    ENDIF("${isSystemDir}" STREQUAL "-1")
endif(APPLE)

# include CMake macros
include(ExternalProject)
include(FindPackageHandleStandardArgs)

# add path containing CMake macros for this project
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# external dependencies are added as dependencies of this target
add_custom_target(dependencies)

include(qt5)
include(XercesC)
include(Geant4)
include(CLHEP)
include(xercesc)
include(CCDB)
include(MYSQL)    ### This should eventually come from CCDB!!!
include(SQLite3)
include(EVIO)
include(mlibrary)
include(cadmesh)  ### Due to the way mlibrary currently gets cadmesh, we need to separately include it here.
include(HIPO)     ### Will be needed later.

add_definitions(${Geant4_DEFINITIONS})

message(STATUS "mlibrary_INCLUDE_DIRS = ${mlibrary_INCLUDE_DIRS}")
message(STATUS "mlibrary_LIBRARIES    = ${mlibrary_LIBRARIES}")

include_directories(
		${mlibrary_INCLUDE_DIRS}
		${EVIO_INCLUDE_DIRS}
		${CCDB_INCLUDE_DIRS}
		${Qt5Xml_INCLUDE_DIRS}
		${Qt5Widgets_INCLUDE_DIRS}
		${Qt5Sql_INCLUDE_DIRS}
		${Geant4_INCLUDE_DIRS}
		)

set(GEMC_ALL_LIBRARIES
	${mlibrary_LIBRARIES}
	${EVIO_LIBRARIES}
	${HIPO_LIBRARIES}
	${CCDB_LIBRARIES}
	${MYSQL_LIBRARIES}
	${SQLite3_LIBRARIES}
	${Qt5Xml_LIBRARIES}
	${Qt5Widgets_LIBRARIES}
	${Qt5Sql_LIBRARIES}
	${Geant4_LIBRARIES}
	${cadmesh_LIBRARIES}
	)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "INCLUDE_DIR ='${dir}'")
endforeach()

set(materials_sources
	materials/material_factory.cc
	materials/cpp_materials.cc
	materials/mysql_materials.cc
	materials/text_materials.cc)
include_directories(materials)
set(GEMC_ALL_SOURCES ${materials_sources})

set(mirrors_sources
	mirrors/mirrors_factory.cc
	mirrors/mysql_mirrors.cc
	mirrors/text_mirrors.cc)
include_directories(mirrors)
list(APPEND GEMC_ALL_SOURCES ${mirrors_sources})

set(parameters_sources
	parameters/parameter_factory.cc
	parameters/mysql_parameters.cc
	parameters/text_parameters.cc)
include_directories(parameters)
list(APPEND GEMC_ALL_SOURCES ${parameters_sources})

set(utilities_sources
	utilities/string_utilities.cc
	utilities/utils.cc
	utilities/lStdHep.cc
	utilities/lXDR.cc
	utilities/options.cc)
include_directories(utilities)
list(APPEND GEMC_ALL_SOURCES ${utilities_sources})


set(detector_sources
	detector/detector.cc
	detector/detector_factory.cc
	detector/mysql_det_factory.cc
	detector/gdml_det_factory.cc
	detector/cad_det_factory.cc
	detector/clara_det_factory.cc
	detector/text_det_factory.cc)
include_directories(detector)
list(APPEND GEMC_ALL_SOURCES ${detector_sources})


set(sensitivity_sources
	sensitivity/sensitiveDetector.cc
	sensitivity/identifier.cc
	sensitivity/Hit.cc
	sensitivity/backgroundHits.cc
	sensitivity/HitProcess.cc
	sensitivity/sensitiveID.cc)
include_directories(sensitivity)
list(APPEND GEMC_ALL_SOURCES ${sensitivity_sources})

set(physics_sources
	physics/PhysicsList.cc
	physics/GammaNuclearPhysics.cc
	physics/PhysicsListMessenger.cc)
include_directories(physics)
list(APPEND GEMC_ALL_SOURCES ${physics_sources})

set(fields_sources
	fields/field.cc
	fields/fieldFactory.cc
	fields/asciiField.cc
	fields/mappedField.cc
	fields/multipoleField.cc
	fields/symmetries/dipole.cc
	fields/symmetries/cylindrical.cc
	fields/symmetries/phi-segmented.cc
	fields/symmetries/cartesian_3d.cc)
include_directories(fields)
list(APPEND GEMC_ALL_SOURCES ${fields_sources})

set(hitprocess_sources
	hitprocess/HitProcess_MapRegister.cc
	hitprocess/flux_hitprocess.cc
	hitprocess/mirror_hitprocess.cc
	hitprocess/counter_hitprocess.cc
	hitprocess/clas12/micromegas/Lorentz.cc
	hitprocess/clas12/micromegas/FMT_hitprocess.cc
	hitprocess/clas12/micromegas/fmt_strip.cc
	hitprocess/clas12/micromegas/BMT_hitprocess.cc
	hitprocess/clas12/micromegas/bmt_strip.cc
	hitprocess/clas12/micromegas/ftm_hitprocess.cc
	hitprocess/clas12/micromegas/ftm_strip.cc
	hitprocess/clas12/rtpc_hitprocess.cc
	hitprocess/clas12/svt/bst_hitprocess.cc
	hitprocess/clas12/svt/bst_strip.cc
	hitprocess/clas12/alert/myatof_hitprocess.cc
	hitprocess/clas12/alert/ahdc_hitprocess.cc
	hitprocess/clas12/ctof_hitprocess.cc
	hitprocess/clas12/cnd_hitprocess.cc
	hitprocess/clas12/dc_hitprocess.cc
	hitprocess/clas12/ecal_hitprocess.cc
	hitprocess/clas12/ftof_hitprocess.cc
	hitprocess/clas12/ft_cal_hitprocess.cc
	hitprocess/clas12/ft_hodo_hitprocess.cc
	hitprocess/clas12/htcc_hitprocess.cc
	hitprocess/clas12/ltcc_hitprocess.cc
	hitprocess/clas12/rich_hitprocess.cc
	hitprocess/bdx/cormo_hitprocess.cc
	hitprocess/bdx/veto_hitprocess.cc
	hitprocess/bdx/crs_hitprocess.cc
	hitprocess/eic/eic_compton_hitprocess.cc
	hitprocess/eic/eic_dirc_hitprocess.cc
	hitprocess/eic/eic_ec_hitprocess.cc
	hitprocess/eic/eic_preshower_hitprocess.cc
	hitprocess/eic/eic_rich_hitprocess.cc
	hitprocess/injector/bubble_hitprocess.cc
	hitprocess/HPS/ECAL_hitprocess.cc
	hitprocess/HPS/SVT_hitprocess.cc
	hitprocess/HPS/muon_hodo_hitprocess.cc)
include_directories(hitprocess)
list(APPEND GEMC_ALL_SOURCES ${hitprocess_sources})

set(output_sources
	output/outputFactory.cc
	output/evio_output.cc
	output/hipo_output.cc
	output/hipoSchemas.cc
	output/txt_output.cc
	output/txt_simple_output.cc
	output/gbank.cc)
include_directories(output)
list(APPEND GEMC_ALL_SOURCES ${output_sources})

set(gui_sources
	gui/src/gemc_MainGui.cc
	gui/src/detector_editor.cc
	gui/src/runControl/run_control.cc
	gui/src/runControl/primaryTab.cc
	gui/src/runControl/momControls.cc
	gui/src/runControl/vtxControls.cc
	gui/src/camera_control.cc
	gui/src/detector_tree.cc
	gui/src/infos.cc
	gui/src/g4dialog.cc
	gui/src/gsignal.cc
	gui/src/physicsListGui.cc
	gui/src//gtrigger.cc
	utilities/graph.cc)
include_directories(gui/src)
list(APPEND GEMC_ALL_SOURCES ${gui_sources})

set(gemc_sources
	gemc.cc
	src/dmesg_init.cc
	src/run_conditions.cc
	src/gemc_options.cc
	src/MDetectorConstruction.cc
	src/MEventAction.cc
	src/MPrimaryGeneratorAction.cc
	src/ActionInitialization.cc
	src/MSteppingAction.cc)
include_directories(src)
list(APPEND GEMC_ALL_SOURCES ${gemc_sources})

add_executable(gemc ${GEMC_ALL_SOURCES})
target_link_libraries(gemc ${GEMC_ALL_LIBRARIES})