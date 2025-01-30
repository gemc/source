#=======================================
#	gemc materials definition
#
#	This file defines a MyMaterials class that holds all the parameters needed to define a material in gemc or GEANT4.
#	Any material in the project is an instance of this class.  The "print_mat" function writes out the material parameters to the
#	materials definition text file that gemc takes as input
#
#	Class members (all members are text strings):
#	name			- The name of the material
#	description			-  A description of the material
#	density			- The material density in g/cm3
#	ncomponents		- The number of of components of the material, e.g., water (H2O) has 2 components: H and O
#	components		- A list of the components and their relative amounts in the material, e.g. "H 2 O 1"
#
#	*****  The next values set optical properties for materials.  They can be ignored (left to default values) if not using optical physics.
#
#	photonEnergy		- A list of photon energies at which any other optical parameters will be provided
#					- Not required (leave as default "none" if not using optical physics)
#					- if any optical parameter (indexOfRefraction, reflectivity, etc.) is defined, photonEnergy MUST also be defined
#					- provide as a list of energies with units, for example:  "1.0*eV 2.0*eV 3.0*eV 4.0*eV 5.0*eV 6.0*eV"
#	indexOfRefraction	- A list of the refractive index evaluated at the energies named in photonEnergy "1.40 1.39 1.38 1.37 1.36"
#					- must have same number of elements in list as in photonEnergy - same for all optical parameters
#	absorptionLength	- A list of the material absorption length evaluated at the energies in photonEnergy
#					- includes units, for example:  "72.8*m 53.2*cm 39.1*cm"
#	reflectivity			- A list of reflectivity values evaluated at the energies in photonEnergy
#	efficiency			- A list of absorption efficiency evaluated at the energies in photonEnergy
#					- efficiency is only used for a dielectric-metal optical boundary where there is no refraction
#					- At this boundary the photon is either reflected or absorbed by the metal with this efficiency
#					- This parameter can be used to define a quantum efficiency for a PMT, for example
#       mie                     - Sets Mie scattering length evaluated at the energies in photonEnergy
#       mieforward              - Constant, describing angular distribution of forward Mie scattering
#       miebackward             - Constant, describing angular distribution of backward Mie scattering
#       mieratio                - Constant, describing ratio of forward to backward Mie scattering (1 - all forward, based on g4 user guide)
#	*****  The next values are about defining scintillators.  They can be ignored (left to default values) if not using a scintillator
#	***** Scintillators are assumed to have a fast and slow response component, defined by relative spectra
#
#	fastcomponent		- A list of the fast component relative spectra values evaluated at the energies in photonEnergy
#	slowcomponent		- A list of the fast component relative spectra values evaluated at the energies in photonEnergy
#	scintillationyield		- Characteristic light yield in photons/MeV e-, given as a single number not a list:  "8400."
#	resolutionscale		- Resolution scale broadens the statistical distribution of generated photons
#					- due to impurities typical of doped crystals like NaI(Tl) and CsI(Tl).  Can be narrower
#					- when the Fano factor plays a role.  Actual number of emitted photons in a step fluctuates
#					- around the mean number with width (ResolutionScale*sqrt(MeanNumberOfPhotons)
#					- Resolution scale is given as a single number, not a list:  "2.0"
#	fasttimeconstant		- (??) believe this is related to the scintillator pulse rise time.  Given as number with units: "1.6*ns"
#	slowtimeconstant	- (??) believe this is related to scintillator slow decay time. Given as number with units: "3.2*ns"
#	yieldratio			- relative strength of the fast component as a fraction of total scintillation yield:  "0.8"
#	rayleigh			- A list of the Rayleigh scattering attenuation coefficient evaluated at the energies in photonEnergy
#
#
#	******	Note that photon energies can be obtained from the wavelength:
#			lambda * nu = c	where lambda is wavelength, c is the speed of light, and nu is frequency
#			E = h * nu		where h is Plank's constant
#			A handy relation for estimating is that h*c ~ 197 eV*nm

# Material class definition
class MyMaterial():
	def __init__(self, name="none", description="none", density="none", ncomponents="none", components="none",
	photonEnergy="none", indexOfRefraction="none", absorptionLength="none", reflectivity="none", efficiency="none",
	fastcomponent="none", slowcomponent="none", scintillationyield="-1", resolutionscale="-1", 
        fasttimeconstant="-1", slowtimeconstant="-1", yieldratio="-1", rayleigh="none", 
        mie="none", mieforward="-1", miebackward="-1", mieratio="-1"):

		self.name = name
		self.description = description
		self.density = density
		self.ncomponents = ncomponents
		self.components = components
		self.photonEnergy = photonEnergy
		self.indexOfRefraction = indexOfRefraction
		self.absorptionLength = absorptionLength
		self.reflectivity = reflectivity
		self.efficiency = efficiency
		self.fastcomponent = fastcomponent
		self.slowcomponent = slowcomponent
		self.scintillationyield = scintillationyield
		self.resolutionscale = resolutionscale
		self.fasttimeconstant = fasttimeconstant
		self.slowtimeconstant = slowtimeconstant
		self.yieldratio = yieldratio
		self.rayleigh = rayleigh
                self.mie = mie
                self.mieforward = mieforward
                self.miebackward = miebackward
                self.mieratio = mieratio

# Function to initialize (overwrite) any existing material file so that any new materials can simply be appended in this project	
def init_materials_file(configuration):
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__materials_"+str(configuration.variation)+".txt"
		open(fileName, "w+")

# Function to write out a material definition to the material file for use as input by gemc
def print_mat(configuration, material):
	# TEXT Factory
	if configuration.factory == "TEXT":
		fName = configuration.detector_name+"__materials_"+configuration.variation+".txt"
		with open(fName, "a+") as fn:
			lstr = ""
			
			lstr += "%20s  |" % str(material.name)
			lstr += "%30s  |" % str(material.description)
			lstr += "%10s  |" % str(material.density)
			lstr += "%10s  |" % str(material.ncomponents)
			lstr += "%50s  |" % str(material.components)
			# optical parameters
			lstr += "%5s  |" % str(material.photonEnergy)
			lstr += "%5s  |" % str(material.indexOfRefraction)
			lstr += "%5s  |" % str(material.absorptionLength)
			lstr += "%5s  |" % str(material.reflectivity)
			lstr += "%5s  |" % str(material.efficiency)
			# scintillation parameters
			lstr += "%5s  |" % str(material.fastcomponent)
			lstr += "%5s  |" % str(material.slowcomponent)
			lstr += "%5s  |" % str(material.scintillationyield)
			lstr += "%5s  |" % str(material.resolutionscale)
			lstr += "%5s  |" % str(material.fasttimeconstant)
			lstr += "%5s  |" % str(material.slowtimeconstant)
			lstr += "%5s  |" % str(material.yieldratio)
			# other optical processes
			lstr += "%5s  |" % str(material.rayleigh)
			lstr += "%5s  |" % str(material.mie)
			lstr += "%5s  |" % str(material.mieforward)
			lstr += "%5s  |" % str(material.miebackward)
                        lstr += "%5s  \n" % str(material.mieratio)

			fn.write(lstr)
