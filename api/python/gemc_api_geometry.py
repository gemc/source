#=======================================
#	gemc geometry definition
#
#	This file defines a MyDetector class that holds all the parameters needed to define a physical volume in gemc or GEANT4.
#	Any folume in the project is an instance of this class.  The "print_det" function writes out the volume parameters to the
#	geometry definition text file that gemc takes as input.
#
#	The instance variables are given default values when a volume is created.  Any variables can be defined with named  arguments
#	in the call to MyDetector() or updated later.  Any variable that will keep its default value can be ignored in the detector definition
#
#	Class members (all members are text strings):
#	name			- The name of the volume
#	mother			- The name of the parent volume, defaults to the world volume "root"
#	description			-  A description of the volume
#	type				- The volume type, e.g., shape.  For example, "Box" or "Tube" or "Polycone", from the GEANT4 list of solids
#	dimensions		- A string giving the dimensions (with units) required to define the shape, from the GEANT4 list of solids
#					- This will vary by the type.  A "Box" type needs 3 dimensions:  length, width, height "2*cm 5*m 24*cm"
#					- Note gemc uses the GEANT4 shape definitions that use half-lengths from the center of the shape
#					- A cylinder ("Tube") shape has inner radius, outer radius, height, starting angle, total angle:
#						"10*cm 10.4*cm 30.*cm 0.*deg 360.*deg" for a long thin-walled, hollow tube
#	pos				- A string giving the location of this volume relative to the center of the parent volume (including "root")
#	rotation			- A string giving the rotation of this volume relative to the coordinate system of the parent
#					- The rotation string has 3 values, defining rotations around the x-, y-, and z- axes
#					- For example "0*deg 0*deg 45*deg" defines a 45-degree rotation about the positive x-axis of the mother volume
#	color				- A hexidecimal color value string. For example, "0000ff" is blue.  An optional 7th digit from 0-5
#					- sets the transparency value where 0 is completely opaque and 5 is fully transparent (invisible)
#					- "0000ff4" gives the volume a mostly transparent blue color
#	material			- text string defining the volume's material.  This can be a GEANT4 element definition given by "G4_Fe"
#					- (replace Fe with the element symbol you want).  It can be any material in the GEANT4 NIST materials list.
#					- For example, "G4_MYLAR" or "G4_PLEXIGLASS".  It can also be the name of a material defined
#					- using the MyMaterial class in this API as long as it appears in the materials text file that gemc takes as input
#	mfield			- The name of a magnetic field attached to the volume ("no" is default). The field is defined in a separate
#					- file in a folder pointed to by the "FIELD_DIR" option in the gcard file.
#	ncopy
#	pMany
#	exist				- This is an integer field:  1 if the volume exists, 0 if not.  It is a way to turn off volumes in a geometry
#	visible			- This is an integer field: 1 if the volume should be visible when the geometry is displayed, 0 if not
#	style				- This is an integer field: 1 means display the volume as a solid, 0 means display as wireframe
#	sensitivity			- A string giving the name of the hit definition to associate with this volume if it is to be a sensitive detector
#					- Set to "flux" to record the true information from the track without digitization
#					- Set to the name of the digitization factory to get signals after the influence of the readout electronics
#					- At present, digitizations are hard-coded into gemc so users must either adapt one of the available
#					- factories ("ctof", "ec", etc) or recompile gemc to add their own.  The goal is to make digitization a
#					- plugin like the materials or the geometry to make it easy to specify your own. See the gemc roadmap
#					- A special "mirror: " sensitivity is used for optical boundaries on volumes.
#	hit_type			- A string specifying the type of hit to record.  At present this is a similar value to sensitivity.
#					- For a FLUX detector this will be "flux", for an optical boundary this will be "mirror"
#	identifiers			- A string with identifiers used to ID the volume for digitization.  For the "ctof" factory, for example,
#					- this might include "layer strip paddle" or others.  KEYWORD "manual" also means... ???
#					- For a FLUX detector, the identifier is just an integer number in the string to show where the hit came from

# Detector class definition
class MyDetector():
	def __init__(self, name="none", mother="root", description="no description", type="none", dimensions="0", 
	pos="0 0 0", rotation="0 0 0", color="999999", material="none", mfield="no", ncopy=1, pMany=1, 
	exist=1, visible=1, style=0, sensitivity="no", hit_type="no", identifiers="no"):
		
		self.name = name
		self.mother = mother
		self.description = description
		self.type = type
		self.dimensions = dimensions
		self.pos = pos
		self.rotation = rotation
		self.color = color
		self.material = material
		self.mfield = mfield
		self.ncopy = ncopy
		self.pMany = pMany
		self.exist = exist
		self.visible = visible			# 0 is invisible, 1 is visible
		self.style = style			# 0 is wireframe, 1 is solid
		self.sensitivity = sensitivity
		self.hit_type = hit_type
		self.identifiers = identifiers

# Function to initialize (overwrite) any existing detector file so that any new detectors can simply be appended in this project	
def init_geom_file(configuration):
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__geometry_"+str(configuration.variation)+".txt"
		open(fileName, "w+")

# Function to write out a detector definition to the detector file for use as input by gemc
def print_det(configuration, detector):
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__geometry_"+str(configuration.variation)+".txt"
		with open(fileName, "a+") as dn:
			lstr = ""
			lstr += "%20s  |" % detector.name
			lstr += "%20s  |" % detector.mother
			lstr += "%30s  |" % detector.description
			lstr += "%50s  |" % detector.pos
			lstr += "%40s  |" % detector.rotation
			lstr += "%7s   |" % detector.color
			lstr += "%20s  |" % detector.type
			lstr += "%60s  |" % detector.dimensions
			lstr += "%20s  |" % detector.material
			lstr += "%20s  |" % detector.mfield
			lstr += "%6s   |" % detector.ncopy
			lstr += "%6s   |" % detector.pMany
			lstr += "%6s   |" % detector.exist
			lstr += "%4s   |" % detector.visible
			lstr += "%4s   |" % detector.style
			lstr += "%20s  |" % detector.sensitivity
			lstr += "%20s  |" % detector.hit_type
			lstr += "%40s \n" % detector.identifiers
			
			dn.write(lstr)

