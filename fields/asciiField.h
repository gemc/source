#ifndef ASCII_FIELDS_H
#define ASCII_FIELDS_H

#include "fieldFactory.h"

class asciiField : public fieldFactory
{
public:
	~asciiField(){}
	
	// check if field object contains a gfield XML header
	bool isEligible(string);
	
	bool isSymmetric;
	
	// load field definitions
	gfield loadField(string, goptions);
	
	// load field map. This is a dispatcher function for the various types of fields below
	void loadFieldMap(gMappedField*, double);
	
	void loadFieldMap_Dipole(gMappedField*, double);       // load dipole field map
	void loadFieldMap_Cylindrical(gMappedField*, double);  // load cylindrical field map
	void loadFieldMap_phiSegmented(gMappedField*, double); // load phiSegmented field map
	void loadFieldMap_cartesian3d(gMappedField*, double);  // load cartesian3d field map

	// clas12 specific. Notice: this should be in a dedicated binary factory
	// however:
	// - clas12 would be the only map using it. Not really worth it, also considering that:
	// - gemc is moving to gemc3 with a better (plugin) mechanism to load fields
	virtual void loadFieldMap(gclas12BinaryMappedField*, double) {} // empty implementation

	static fieldFactory *createFieldFactory() {
		return new asciiField;
	}
	
};


#endif
