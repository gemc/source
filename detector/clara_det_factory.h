#ifndef CLARA_DET_FACTORY_H
#define CLARA_DET_FACTORY_H 1

// gemc headers
#include "detector_factory.h"

class clara_det_factory : public detectorFactory
{
	public:
		clara_det_factory(){;}
	
		// load all detectors that matches factorytype
		map<string, detector> loadDetectors();      
		
		// initialize factorytype, option and runcondition classes
		void initFactory(goptions, runConditions, string);		
		
		static detectorFactory *createFactory() 
		{
			return new clara_det_factory; 
		}
};


#endif
