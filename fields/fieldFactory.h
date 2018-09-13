#ifndef FIELD_FACTORY_H
#define FIELD_FACTORY_H

// C++ headers
#include <map>
#include <iostream>
using namespace std;

// gemc headers
#include "field.h"


class fieldFactory
{
	public:
		virtual ~fieldFactory(){}
		
		virtual bool isEligible(string)                              = 0; // check if field object contain a valid gfield XML header
		virtual gfield loadField(string, goptions)                   = 0; // load field definitions 
		virtual void loadFieldMap(gMappedField*, double)             = 0; // load field map. This is a dispatcher function for the various types of fields below
		
		virtual void loadFieldMap_Dipole(gMappedField*, double)      = 0; // load 1D dipole field depending on 2 coordinates (transverse and longitudinal)
		virtual void loadFieldMap_Cylindrical(gMappedField*, double) = 0; // load cylindrical field depending on 2 coordinates (transverse and longitudinal)
		virtual void loadFieldMap_phiSegmented(gMappedField*, double) = 0; // load phi-symmetric field depending on 3 coordinates (azimuthal,transverse and longitudinal)
		virtual void loadFieldMap_cartesian3d(gMappedField*, double) = 0; // load cartesian_3D field depending on 3 coordinates (XX,YY,ZZ)				
		
		string factoryType;

		void initFactory(string ft)
		{
			cout << "   >> gemc Init: " << ft << " Field Factory is Initialized "  << endl;
			factoryType = ft;
		}
};

typedef fieldFactory *(*fieldFactoryInMap)();                                        // Define fieldFactoryInMap as a pointer to a function that returns a pointer 

fieldFactory *getFieldFactory(map<string, fieldFactoryInMap> *, string);             // returns fieldFactory Function from Factory Map

map<string, fieldFactoryInMap> registerFieldFactories();                             // Registers FieldFactory in Factory Map

map<string, gfield> loadAllFields(map<string, fieldFactoryInMap>, goptions opts );  // returns a map of all available gfields

#endif
