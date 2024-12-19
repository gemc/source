#=======================================
#	gemc hit definition
#
#	This file defines a MyHit class that holds values that describe a detector hit in gemc. The hit definition is used to
#	include the effects of integration over a time window aw well as thresholds and signal propagation in readout electronics.
#	The hit will be used with digitization systems that are currently compiled into gemc.  Eventually these will become a
#	plugin that can easily be defined by the user.
#
#	Each parameter is the instance of the class. All of the instance variables are strings.
#
#	Class members (all members are text strings):
#	name			- The name of the hit
#	description			- A description of the hit
#	identifiers			- A string identifying the detector element recording the hit (see the detector api script)
#	SignalThreshold		- A string containing the minimum energy that must be deposited in the detector to register a hit
#	TimeWindow		- The length of time over which to integrate the detector deposit as part of a single hit
#	ProdThreshold		- (unsure of use in gemc)  May be confusing it with SignalThreshold
#	MaxStep			-
#	riseTime			- rise time of the voltage signal (likely used for the V(t) signal from flash ADC readout)
#	fallTime			- fall time of the voltage signal (likely used to describe the V(t) signal from flash ADC readout)
#	mvToMeV			- essentially the gain of the front end electronics (?)
#	pedestal			- adds a number of ADC counts to account for the base level in the readout electronics
#	delay				- adds a time delay to the signal to account for signal propagation in the readout electronics

# Hit class definition
class MyHit():
	def __init__(self, name="none", description="none", identifiers="id", SignalThreshold="none", TimeWindow="none", ProdThreshold=None,
	MaxStep="none", riseTime="none", fallTime="none", mvToMeV="none", pedestal="none", delay="none"):
		self.name = name
		self.description = description
		self.identifiers = identifiers
		self.SignalThreshold = SignalThreshold
		self.TimeWindow = TimeWindow
		self.ProdThreshold = ProdThreshold
		self.MaxStep = MaxStep
		self.riseTime = riseTime
		self.fallTime = fallTime
		self.mvToMeV = mvToMeV
		self.pedestal = pedestal
		self.delay = delay

# Function to initialize (overwrite) any existing hits file so that any new hits can simply be appended in this project	
def init_hits_file(configuration):
	if configuration.factory == "TEXT":
		fileName = configuration.detector_name + "__hit_"+str(configuration.variation)+".txt"
		open(fileName, "w+")
		

# Function to write out a hit definition to the hit file for use as input by gemc
def print_hit(configuration, hit):
	# TEXT Factory
	if configuration.factory == "TEXT":
		fName = configuration.detector_name+"__hit_"+configuration.variation+".txt"
		with open(fName, "a+") as fn:
			lstr = ""
			lstr += "%20s  |" % hit.name
			lstr += "%30s  |" % hit.description
			lstr += "%40s  |" % hit.identifiers
			lstr += "%8s   |" % hit.SignalThreshold
			lstr += "%8s   |" % hit.TimeWindow
			lstr += "%8s   |" % hit.ProdThreshold
			lstr += "%8s   |" % hit.MaxStep
			lstr += "%8s   |" % hit.riseTime
			lstr += "%8s   |" % hit.fallTime
			lstr += "%8s   |" % hit.mvToMeV
			lstr += "%8s   |" % hit.pedestal
			lstr += "%8s  \n" % hit.delay
			fn.write (lstr)
	
	if float(configuration.verbosity) > 0:
		print("  + Hit %s uploaded successfully for variation %s \n" % (hit.name, str(configuration.variation)))
