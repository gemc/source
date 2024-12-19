#=======================================
#	gemc utils definition
#
#	This file defines a MyConfiguration class that holds the configuration options read in from the config.dat configuration file.
#	Any folume in the project is an instance of this class.  The "print_det" function writes out the volume parameters to the
#	geometry definition text file that gemc takes as input.
#
#	The instance variables are given default values when a volume is created.  Any variables can be defined with named  arguments
#	in the call to MyDetector() or updated later.  Any variable that will keep its default value can be ignored in the detector definition
#
#	Class members (all members are text strings):
#	name			- The name of the detector (in the broad sense).  Think project name here.
#	variation			- The name of the project variation.  For example, one could have variations of the project where
#					- a volume has a different size or material.  The variation defaults to 'default'
#	factory			- The configuration factory defines how the generated files that gemc uses are stored.  Currently this
#					- API only supports "TEXT".  Other choices might include "MYSQL" or a java format
#	dbhost			- The hostname of the mysql database server where gemc detectors, materials, etc. may be stored
#					- Currently unsupported in this API
#	rmin				- The first run number for gemc to generate
#	rmax				- The last run number for gemc.  Used when running gemc to make a batch of runs for later analysis
#	comment			- Any comment regarding the project that the author would like to include
#	verbosity			- The default verbosity level for gemc at runtime.  Individual category verbosities may be set with
#					- command line options to gemc.command or in the .gcard file
#	

# Configuration class definition
class MyConfiguration():
	def __init__(self, detector_name="none", variation="default", factory="TEXT", dbhost="",
	rmin=0, rmax=0, comment="none", verbosity=1):
		self.detector_name = detector_name
		self.variation = variation
		self.factory = factory
		self.dbhost = dbhost
		self.rmin = rmin
		self.rmax = rmax
		self.comment = comment
		self.verbosity = verbosity

# Function to load the configuration data from the configuration file.
def load_configuration(cFile):
	configuration = MyConfiguration()
	with open(cFile,"r") as cn:
		for line in cn:
			if line.startswith("#") or line.isspace():
				continue
			key, value = line.split(":")
			setattr(configuration, key.strip(), value.strip() )
	
	if float(configuration.verbosity) > 0:
		print("\n  * Loading configuration from " + cFile)
		print("   > Detector Name: " + str(configuration.detector_name))
		print("   > Variation:     " + str(configuration.variation))
		print("   > Factory:       " + str(configuration.factory))
		print("   > DB host:       " + str(configuration.dbhost))
		print("   > Run Min:       " + str(configuration.rmin))
		print("   > Run Max:       " + str(configuration.rmin))
		print("   > Comment:       " + str(configuration.comment))
		print("   > Verbosity:     " + str(configuration.verbosity))

	return configuration

# The following code allows this module to be executed as a main python script for the purpose of testing the functions
#  To test, type:  'gemc_api_utils.py config.dat' on the command line, where config.dat is the name of the configuration file
if __name__ == "__main__":	
	import argparse, sys

	parser = argparse.ArgumentParser()
	parser.add_argument("config_filename", help="The name of the experiment configuration file")
	if len(sys.argv)==1:
	    parser.print_help()
	    sys.exit(1)
	args = parser.parse_args()
	print(args.config_filename)
	
	cfg = load_configuration(args.config_filename)
