#=======================================
#	gemc mirrors definition
#
#	This file defines a MyMirror class that holds values to define an optical boundary.  In gemc, any type of optical boundary
#	is described as a "mirror", regardless of its use or reflective quality.
#
#	Each mirror definition is an instance of the class. All of the instance variables are strings.
#
#	Class members (all members are text strings):
#	name			- The name of the mirror definition
#	description			-  A description of the mirror definition
#	type				- The surface type (see below).  Example:  "dielectric_dielectric" or "dielectric_metal"
#	finish			- The type of finish of the optical surface (see below for options). Example:  "polishedfrontpainted"
#	model			- The optical model to use for the surface (see below for options)
#	border			- "SkinSurface" if the optical boundary represents the entire outside surface of a volume
#					- For a border surface defined as the contact area between two neighboring volumes, use the
#					- name of the bordering volume.  Optical surfaces are ordered in GEANT4 and gemc.
#	matOptProps		- use the name of a material with optical properties to use as the boundary properties
#					- The material does not have to be the same as either bordering volume, i.e., a thin paint
#					- if "notDefined" then use the properties defined in the following variables to define the boundary
#	photonEnergy		- a list of photon energies at which to evaluate the following optical properties
#	indexOfRefraction	- a list of boundary material refractive indices evaluated at the energies in photonEnergy
#	reflectivity			- a list of boundary material reflectivities evaluated at the energies in photonEnergy
#	efficiency			- a list of photoelectric absorption efficiencies evaluated at the energies in photonEnergy
#					- efficiency is used in "dielectric_metal" boundaries where the photon is either reflected or
#					- absorbed by the metal with this efficiency.  Can be used as the quantum efficiency in a PMT.
#	specularlobe		- defines scattering properties of a rough optical surface
#	specularspike		- defines scattering properties of a rough optical surface
#	backscatter			- defines scattering properties of a rough optical surface

# =================================================================
# Available finish in materials/include/G4OpticalSurface.hh:
#
# polished,                    // smooth perfectly polished surface
# polishedfrontpainted,        // smooth top-layer (front) paint
# polishedbackpainted,         // same is 'polished' but with a back-paint
#
# ground,                      // rough surface
# groundfrontpainted,          // rough top-layer (front) paint
# groundbackpainted,           // same as 'ground' but with a back-paint
#
# polishedlumirrorair,         // mechanically polished surface, with lumirror
# polishedlumirrorglue,        // mechanically polished surface, with lumirror & meltmount
# polishedair,                 // mechanically polished surface
# polishedteflonair,           // mechanically polished surface, with teflon
# polishedtioair,              // mechanically polished surface, with tio paint
# polishedtyvekair,            // mechanically polished surface, with tyvek
# polishedvm2000air,           // mechanically polished surface, with esr film
# polishedvm2000glue,          // mechanically polished surface, with esr film & meltmount
#
# etchedlumirrorair,           // chemically etched surface, with lumirror
# etchedlumirrorglue,          // chemically etched surface, with lumirror & meltmount
# etchedair,                   // chemically etched surface
# etchedteflonair,             // chemically etched surface, with teflon
# etchedtioair,                // chemically etched surface, with tio paint
# etchedtyvekair,              // chemically etched surface, with tyvek
# etchedvm2000air,             // chemically etched surface, with esr film
# etchedvm2000glue,            // chemically etched surface, with esr film & meltmount
#
# groundlumirrorair,           // rough-cut surface, with lumirror
# groundlumirrorglue,          // rough-cut surface, with lumirror & meltmount
# groundair,                   // rough-cut surface
# groundteflonair,             // rough-cut surface, with teflon
# groundtioair,                // rough-cut surface, with tio paint
# groundtyvekair,              // rough-cut surface, with tyvek
# groundvm2000air,             // rough-cut surface, with esr film
# groundvm2000glue             // rough-cut surface, with esr film & meltmount

# Available models in materials/include/G4OpticalSurface.hh:
#
# glisur,                      // original GEANT3 model
# unified,                     // UNIFIED model
# LUT                          // Look-Up-Table model

# Available surface types in materials/include/G4SurfaceProperty.hh
#
# dielectric_metal,            // dielectric-metal interface
# dielectric_dielectric,       // dielectric-dielectric interface
# dielectric_LUT,              // dielectric-Look-Up-Table interface
# firsov,                      // for Firsov Process
# x_ray                        // for x-ray mirror process

# Border Volume Types:
#
# SkinSurface: surface of a volume
# Border Surface: surface between two volumes (second volume must exist)

# =================================================================

# Mirror class definition
class MyMirror():
	
	def __init__(self, name="none", description="none", \
	type="notDefined", finish="notDefined", model="notDefined", border="notDefined", \
	matOptProps="none", photonEnergy="notDefined", indexOfRefraction="none", reflectivity="none", \
	efficiency="none", specularlobe="none", specularspike="none", backscatter="none"):
		
		self.name = name		# Mandatory
		self.description = description
		self.type = type			# Mandatory
		self.finish = finish		# Mandatory
		self.model = model		# Mandatory
		self.border = border		# Mandatory
		self.matOptProps = matOptProps
		self.photonEnergy = photonEnergy
		self.indexOfRefraction = indexOfRefraction
		self.reflectivity = reflectivity
		self.efficiency = efficiency
		self.specularlobe = specularlobe
		self.specularspike = specularspike
		self.backscatter = backscatter


# Function to initialize (overwrite) any existing mirror file so that any new mirrors can simply be appended in this project	
def init_mirrors_file(configuration):
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__mirrors_"+str(configuration.variation)+".txt"
		open(fileName, "w+")


# Function to write out a material definition to the material file for use as input by gemc
def print_mir(configuration, mirror):	
	if mirror.matOptProps == "notDefined":
		if mirror.photonEnergy == "none":
			print(" !! Error: there is no material with optical properties associated with this mirror.\n")
			print(" !! Optical properties must be defined.\n")

	# TEXT Factory
	if configuration.factory == "TEXT":
		fName = configuration.detector_name+"__mirrors_"+configuration.variation+".txt"
		with open(fName, "a+") as fn:
			if mirror.type == "notDefined":
				print("Error: type undefined.\n")
			if mirror.finish == "notDefined":
				print("Error: finish undefined.\n")
			if mirror.model == "notDefined":
				print("Error: model undefined.\n")
			if mirror.border == "notDefined":
				print("Error: border undefined.\n")

			lstr = ""
			
			lstr += "%20s  |" % str(mirror.name)
			lstr += "%30s  |" % str(mirror.description)
			lstr += "%24s  |" % str(mirror.type)
			lstr += "%20s  |" % str(mirror.finish)
			lstr += "%10s  |" % str(mirror.model)
			lstr += "%25s  |" % str(mirror.border)
			lstr += "%25s  |" % str(mirror.matOptProps)
			lstr += "%5s |" % str(mirror.photonEnergy)
			lstr += "%5s |" % str(mirror.indexOfRefraction)
			lstr += "%5s |" % str(mirror.reflectivity)
			lstr += "%5s |" % str(mirror.efficiency)
			lstr += "%5s |" % str(mirror.specularlobe)
			lstr += "%5s |" % str(mirror.specularspike)
			lstr += "%5s " % str(mirror.backscatter)
			
			lstr += "\n"
			fn.write(lstr)


	if int(configuration.verbosity) > 0:
		print("  + Mirror %s uploaded successfully for variation %s \n" %(configuration.detector_name, configuration.variation))
