##########################################################
# CMake configuration for Gemc                           #
#                                                        #
# @author Maurik Holtrop, University of New Hampshire    #
#                                                        #
##########################################################

# minimum version of CMake, which is determined by Geant4's requirements
cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12)
if(${CMAKE_VERSION} VERSION_GREATER "3.17.0")
	cmake_policy(SET CMP0100 OLD) # Do not process *.hh files with moc, uic
endif()

project(GEMC VERSION 2.8)

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.3)
	message(STATUS "*************************************************************************************************************")
	message(STATUS "* ERROR:: Please use a later version of GCC.                                                                *")
	message(STATUS "* Your version of ${CMAKE_CXX_COMPILER_VERSION} will not compile GEANT4 code.                                                       *")
	message(STATUS "* Erase the current cache (just do rm -rf * :-) then rerun cmake adding the arguments:                      *")
	message(STATUS "*   -DCMAKE_C_COMPILER=/apps/gcc/9.2.0/bin/gcc -DCMAKE_CXX_COMPILER=/apps/gcc/9.2.0/bin/g++                 *")
	message(STATUS "*************************************************************************************************************")
	message(FATAL_ERROR "Compiler level too low.")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "COMPILERS Used:")
message(STATUS "CC = ${CMAKE_C_COMPILER}")
message(STATUS "CXX= ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_SHARED_LIBRARY_SUFFIX=${CMAKE_SHARED_LIBRARY_SUFFIX}")

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(STATUS  "GEMC requires an out-of-source build.")
	message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
	message(STATUS "CMakeCache.txt")
	message(STATUS "CMakeFiles")
	message(STATUS "Once these files are removed, create a separate directory, eg mkdir build")
	message(STATUS "and run CMake from there")
	message(FATAL_ERROR "in-source build detected")
endif()


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
	ENDIF("${isSystemDir}" STREQUAL "-1")
endif(APPLE)

# include CMake macros
include(ExternalProject)
include(FindPackageHandleStandardArgs)

# add path containing CMake macros for this project
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# external dependencies are added as dependencies of this target
add_custom_target(dependencies)

include(Geant4)
include(Qt5)
# include(XercesC)
include(ccdb)
include(eviocxx)
include(mlibrary)
include(hipo4)

if(DEPENDENCIES_RERUN_CMAKE)
	add_custom_target(build_ext_deps DEPENDS Geant4 ccdb evioxx hipo4)
	message(STATUS "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	message(STATUS "+ One or more of the dependency packages was not found.  +")
	message(STATUS "+ Run: make                                              +")
	message(STATUS "+ Then run cmake again                                   +")
	message(STATUS "+ Then finally run make once more to complete the build. +")
	message(STATUS "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
else()
	set(CMAKE_AUTOMOC ON)
	set(CMAKE_AUTORCC ON)
	set(CMAKE_AUTOUIC ON)

	add_definitions(${Geant4_DEFINITIONS})
	
	include_directories(
			${Qt5Xml_INCLUDE_DIRS}
			${Qt5Widgets_INCLUDE_DIRS}
			${Qt5Sql_INCLUDE_DIRS}
			${Geant4_INCLUDE_DIRS}
			)

	set(GEMC_ALL_LIBRARIES
		${Qt5Xml_LIBRARIES}
		${Qt5Widgets_LIBRARIES}
		${Qt5Sql_LIBRARIES}
		${Geant4_LIBRARIES}
		)

	set(GEMC_ALL_TARGETS
		ccdb
		eviocxx
		mlibrary
		hipo4
		)

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
		utilities/gemcUtils.cc
		utilities/lStdHep.cc
		utilities/lXDR.cc
		utilities/gemcOptions.cc)
	include_directories(utilities)
	list(APPEND GEMC_ALL_SOURCES ${utilities_sources})

	# Do not process these files with moc, setting this avoids the policy warning.
	set_property(SOURCE "utilities/lStdHep.hh" PROPERTY SKIP_AUTOMOC ON)
	set_property(SOURCE "utilities/lXDR.hh" PROPERTY SKIP_AUTOMOC ON)

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
	add_dependencies(gemc dependencies)
	target_link_libraries(gemc ${GEMC_ALL_TARGETS} ${GEMC_ALL_LIBRARIES})

	install(TARGETS gemc )
endif()