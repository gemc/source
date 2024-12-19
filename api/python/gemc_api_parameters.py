#=======================================
#	gemc parameters definition
#
#	This file defines a MyParameter class that holds values that can be input to gemc for use by the simulation.
#	An example of a parameter might be the dimension of a detector element or the number of copies.  There are places
#	in the structure to define where the parameter came from, say a CAD drawing or who the author is.
#
#	Each parameter is the instance of the class. All of the instance variables are strings that default to "none".
#	Any variable that will keep its default value can be ignored in the detector definition.  Only the name and value are really required
#
#	Class members (all members are text strings):
#	name			- The name of the parameter
#	value				- The value of the parameter
#	units				- The units for the value
#	description			-  A description of the parameter
#	author			- The name of the author of the parameter
#	author_email		- The email address for contacting the author of the parameter
#	link				- A link to a drawing containing the parameter, for example, a CAD drawing
#	drawing_varname	- The name of the parameter as it is used in the linked drawing
#	drawing_author		- The name of the author of the linked drawing

# Parameter class definition
class MyParameter():
	def __init__(self, name="none", value="none", units="none", description="none", 
	author="none", author_email="none", link="none", drawing_varname="none", drawing_author="none"):
		self.name = name
		self.value = value
		self.units = units
		self.description = description
		self.author = author
		self.author_email = author_email
		self.link = link
		self.drawing_varname = drawing_varname
		self.drawing_author = drawing_author


# Function to retrieve the list of parameters from the parameter text file.  The function returns a python dictionary (parameterDict) of
# (key, value) pairs where 'key' is the parameter name and 'value' is the parameter value.  The full class structure of each parameter
# is also stored in a list (parameterList) that is currently not returned.
def get_parameters(configuration):
	parameterDict = {}
	parameterList = []
	
	# Text Factory. The parameter file is assumed to be present
	# and named "parameters.txt"
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__parameters_"+str(configuration.variation)+".txt"
		with open(fileName, "r") as f:
			for line in f:
				param = MyParameter()
				sline = line.split("|")
				aline = list(map(str.strip, sline))
				#~ print(aline)
				param.name = aline[0]
				param.value = float(aline[1])
				param.units = aline[2]
				param.description = aline[3]
				param.author = aline[4]
				param.author_email = aline[5]
				param.link = aline[6]
				param.drawing_varname = aline[7]
				param.drawing_author = aline[8]
				parameterDict[param.name] = param.value
				parameterList.append(param)
	
	if float(configuration.verbosity) > 0:
		for key in parameterDict.keys():
			print(" * Parameter \"%s\" loaded with value: %s \n" % (key, parameterDict[key]) )
	return parameterDict

# Function to print out the full set of fields for each stored parameter.  This is currently used to test whether the parameters
# are being read in from the file correctly.
def print_parameters(parameterList):
	for param in parameterList:
		print("-----------------------------------------")
		print("  - Parameter name: %s" % param.name)
		print("  - Value: %s %s" %(param.value, param.units) )
		print("  - Description: %s" % param.description)
		print("  - Author: %s  email: %s" % (param.author, param.author_email) )
		print("  - Link to Drawing or Document: %s" % param.link)
		print("  - Variable name on the drawing: %s" % param.drawing_varname)
		print("  - Drawing Author: %s" % param.drawing_author)
		print("-----------------------------------------")

