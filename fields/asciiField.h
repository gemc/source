#ifndef ASCII_FIELDS_H
#define ASCII_FIELDS_H

#include "fieldFactory.h"

class asciiField : public fieldFactory
{
	public:
	~asciiField(){}
	
	// check if field object contains a gfield XML header
	bool isEligible(string);
	
	// load field definitions
	gfield loadField(string, goptions);
	
	// load field map. This is a dispatcher function for the various types of fields below
	void loadFieldMap(gMappedField*, double);
	
	void loadFieldMap_Dipole(gMappedField*, double);       // load dipole field map
	void loadFieldMap_Cylindrical(gMappedField*, double);  // load cylindrical field map
	void loadFieldMap_phiSegmented(gMappedField*, double); // load cylindrical field map
	
	static fieldFactory *createFieldFactory()
	{
		return new asciiField;
	}
	
};


#endif
