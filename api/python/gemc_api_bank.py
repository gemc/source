#=======================================
#	gemc hit definition
#
#	This file writes a bank variable to the bank text file.  This is used along with a digitization factory compiled into gemc.
#	This part of the API was written to match the perl api function and does work with the ctof example.  It has not been
#	well tested or documented at this point.
#
#	This file does not currently define a class for the bank.  It simply write the bank variables to the gemc text file.
#
#	Bank variables (all members are text strings):
#	name			- The name of the bank
#	varname			- The name of the variable being inserted into the bank
#	description			- A description of the bank variable
#	num
#	type


# Function to initialize (overwrite) any existing bank file so that any new bank variables can simply be appended in this project	
def init_bank_file(configuration):
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__bank.txt"
		open(fileName, "w+")


# Print bank variable to TEXT file
def insert_bank_variable(configuration, bankname=None, varname=None, num=None, type=None, description=None):
	for item in [bankname, varname, num, type, description]:
		if item == None:
			print(" ERROR: To define a bank variable 4 arguments should be passed to <insert_bank_variable> \n")
			sys.exit(1)
	
	# TEXT Factory
	if configuration.factory == "TEXT":
		fName = configuration.detector_name+"__bank.txt"
		with open(fName, "a+") as fn:
			lstr = ""
			lstr += "%20s  |" % bankname
			lstr += "%20s  |" % varname
			lstr += "%50s  |" % description
			lstr += "%5s   |" % num
			lstr += "%20s  \n" % type
			fn.write (lstr)

	
	
	if (float(configuration.verbosity) > 0):
		print("  + variable %s uploaded successfully for variation %s"%(varname,str(configuration.variation)) )

		
