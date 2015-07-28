#ifndef MYSQL_DET_FACTORY_H
#define MYSQL_DET_FACTORY_H 1

// gemc headers
#include "detector_factory.h"

class mysql_det_factory : public detectorFactory
{
	public:
		mysql_det_factory(){;}
	
		// load all detectors that matches factorytype
		map<string, detector> loadDetectors();      
		
		// initialize factorytype, option and runcondition classes
		void initFactory(goptions, runConditions, string);		
		
		static detectorFactory *createFactory() 
		{
			return new mysql_det_factory; 
		}
};


#endif
