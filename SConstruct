from init_env import init_environment

# adding ccdb as temporary dependency
# will be removed once the hit process routines are plugins
env = init_environment("qt5 geant4 clhep evio xercesc ccdb mlibrary cadmesh")

# env.Replace(CXX = "/apps/gcc/4.7.2/bin/g++")
# env.Replace(CXX = "/usr/bin/clang++")
#env.Append(CXXFLAGS='-Wno-shorten-64-to-32')
#env.Append(CXXFLAGS='-Wno-sign-conversion')

# added because clhep is still behind clang (5/2015)
if env['PLATFORM'] == 'darwin':
	env.Append(CXXFLAGS='-Wno-absolute-value')

# Materials
env.Append(CPPPATH = 'materials')
materials_sources = Split("""
	materials/material_factory.cc
	materials/cpp_materials.cc
	materials/mysql_materials.cc
	materials/text_materials.cc""")
env.Library(source = materials_sources, target = "lib/materials")

# Mirrors
env.Append(CPPPATH = 'mirrors')
mirrors_sources = Split("""
	mirrors/mirrors_factory.cc
	mirrors/mysql_mirrors.cc
	mirrors/text_mirrors.cc""")
env.Library(source = mirrors_sources, target = "lib/mirrors")


# Parameters
env.Append(CPPPATH = 'parameters')
parameters_sources = Split("""
	parameters/parameter_factory.cc
	parameters/mysql_parameters.cc
	parameters/text_parameters.cc""")
env.Library(source = parameters_sources, target = "lib/parameters")


# Utilities
env.Append(CPPPATH = 'utilities')
util_sources = Split("""
	utilities/string_utilities.cc
	utilities/utils.cc
	utilities/lStdHep.cc
	utilities/lXDR.cc
	utilities/options.cc""")
env.Library(source = util_sources, target = "lib/utilities")


# Detector
env.Append(CPPPATH = 'detector')
det_sources = Split("""
	detector/detector.cc
	detector/detector_factory.cc
	detector/mysql_det_factory.cc
	detector/gdml_det_factory.cc
	detector/cad_det_factory.cc
	detector/clara_det_factory.cc
	detector/text_det_factory.cc""")
env.Library(source = det_sources, target = "lib/detector")


# Sensitivity
env.Append(CPPPATH = 'sensitivity')
sensi_sources = Split("""
	sensitivity/sensitiveDetector.cc
	sensitivity/identifier.cc
	sensitivity/Hit.cc
	sensitivity/HitProcess.cc
	sensitivity/sensitiveID.cc""")
env.Library(source = sensi_sources, target = "lib/sensitivity")


# Physics
env.Append(CPPPATH = 'physics')
phys_sources = Split("""
	physics/PhysicsList.cc
	physics/PhysicsListMessenger.cc""")
env.Library(source = phys_sources, target = "lib/physics")

# Fields
env.Append(CPPPATH = 'fields')
field_sources = Split("""
	fields/field.cc
	fields/fieldFactory.cc
	fields/asciiField.cc
	fields/mappedField.cc
	fields/multipoleField.cc 
	fields/symmetries/dipole.cc
	fields/symmetries/cylindrical.cc
	fields/symmetries/phi-segmented.cc""")
env.Library(source = field_sources, target = "lib/fields")


# Hit Processes
env.Append(CPPPATH = ['hitprocess', 'hitprocess/clas12', 'hitprocess/clas12/svt', 'hitprocess/clas12/micromegas'])
env.Append(CPPPATH = ['hitprocess/Aprime', 'hitprocess/GlueX', 'hitprocess/solid'])
hitp_sources = Split("""
	hitprocess/HitProcess_MapRegister.cc
	hitprocess/flux_hitprocess.cc
	hitprocess/clas12/micromegas/FMT_hitprocess.cc
	hitprocess/clas12/micromegas/fmt_strip.cc	
	hitprocess/clas12/micromegas/BMT_hitprocess.cc
	hitprocess/clas12/micromegas/bmt_strip.cc
	hitprocess/clas12/micromegas/ftm_hitprocess.cc
	hitprocess/clas12/micromegas/ftm_strip.cc
   hitprocess/clas12/bonus_hitprocess.cc
   hitprocess/clas12/svt/bst_hitprocess.cc
	hitprocess/clas12/svt/bst_strip.cc
	hitprocess/clas12/ctof_hitprocess.cc
	hitprocess/clas12/cnd_hitprocess.cc
	hitprocess/clas12/dc_hitprocess.cc
	hitprocess/clas12/ec_hitprocess.cc
	hitprocess/clas12/ecs_hitprocess.cc
	hitprocess/clas12/ftof_hitprocess.cc
	hitprocess/clas12/ft_cal_hitprocess.cc
	hitprocess/clas12/ft_hodo_hitprocess.cc
	hitprocess/clas12/htcc_hitprocess.cc
	hitprocess/clas12/ltcc_hitprocess.cc
	hitprocess/clas12/pcal_hitprocess.cc
	hitprocess/clas12/rich_hitprocess.cc
	hitprocess/bdx/cormo_hitprocess.cc
	hitprocess/bdx/veto_hitprocess.cc
	hitprocess/bdx/crs_hitprocess.cc
	hitprocess/injector/bubble_hitprocess.cc
	hitprocess/HPS/ECAL_hitprocess.cc
	hitprocess/HPS/SVT_hitprocess.cc
	hitprocess/HPS/muon_hodo_hitprocess.cc""")
env.Library(source = hitp_sources, target = "lib/hitprocess")


# Output
env.Append(CPPPATH = 'output')
output_sources = Split("""
	output/outputFactory.cc
	output/evio_output.cc
	output/txt_output.cc
	output/gbank.cc""")
env.Library(source = output_sources, target = "lib/output")

# GUI
env.Append(CPPPATH = 'gui/src')
gui_sources = Split("""
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
	utilities/graph.cc""")
env.Library(source = gui_sources, target = "lib/gui")


# gemc
env.Append(CPPPATH = 'src')
gemc_sources = Split("""
	gemc.cc
	src/dmesg_init.cc
	src/run_conditions.cc
	src/gemc_options.cc
	src/MDetectorConstruction.cc
	src/MDetectorMessenger.cc
	src/MEventAction.cc
	src/MPrimaryGeneratorAction.cc
	src/ActionInitialization.cc
	src/MSteppingAction.cc""")

env.Append(LIBPATH = ['lib'])
env.Prepend(LIBS =  ['materials', 'mirrors', 'parameters', 'utilities', 'detector', 'sensitivity', 'physics', 'fields', 'hitprocess', 'output', 'gui'])
env.Program(source = gemc_sources, target = "gemc")


if env['LIBRARY'] == "static":
	env.Library(source = gemc_sources, target = "gemc")


if env['LIBRARY'] == "shared":
	env.SharedLibrary(source = gemc_sources, target = "gemc")







