#ifndef detector_factory_H
#define detector_factory_H 1

// gemc headers
#include "options.h"
#include "detector.h"
#include "run_conditions.h"
#include "utils.h"

// c++ headers
#include <map>
#include <iostream>

class detectorFactory
{
	public:
		// Pure Virtual Method to create the map<string, detectors>
		virtual map<string, detector> loadDetectors() = 0;         
		virtual ~detectorFactory(){}

		// initialize factorytype, option and runcondition classes
		void initFactory(goptions, runConditions, string) ;    
		
		string factoryType;
		goptions gemcOpt;
		runConditions RC;
};

// Define detectorFactoryInMap as a pointer to a function that returns a pointer 
typedef detectorFactory *(*detectorFactoryInMap)();                                                 

// returns detectorFactory from Factory Map
detectorFactory *getDetectorFactory(map<string, detectorFactoryInMap> *detectorFactoryMap, string);  

// Registers detectorFactories in detectorFactoryMap
map<string, detectorFactoryInMap> registerDetectorFactory();                                        

// build detectors according to their factory
map<string, detector> buildDetector(map<string, detectorFactoryInMap> detectorFactoryMap, goptions go, runConditions rc);

string check_factory_existance(map<string, detectorFactoryInMap> detectorFactoryMap, runConditions rc);

// load detector from gtable
detector get_detector(gtable, goptions go, runConditions rc);

#endif
