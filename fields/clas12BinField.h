#ifndef C12BIN_FIELDS_H
#define C12BIN_FIELDS_H

#include "fieldFactory.h"

// uses David's Heddle cmag: https://github.com/JeffersonLab/clas12-cmag
// example of fields used in real data:
// https://github.com/JeffersonLab/clas12-offline-software/blob/development/etc/services/data.yaml

class clas12BinField : public fieldFactory
{
public:
	~clas12BinField(){}

	// save map list for loadFieldMap
	map<string, vector<string> > validC12MapNames;


	// implementing virtual method
	bool isEligible(string);

	// load field definitions
	gfield loadField(string, goptions);
	
	// load field map. This is a dispatcher function for the various types of fields below
	void loadFieldMap(gMappedField*, double)              {} // empty implementation
	void loadFieldMap_Dipole(gMappedField*, double)       {} // empty implementation
	void loadFieldMap_Cylindrical(gMappedField*, double)  {} // empty implementation
	void loadFieldMap_phiSegmented(gMappedField*, double) {} // empty implementation
	void loadFieldMap_cartesian3d(gMappedField*, double)  {} // empty implementation

	// clas12 specific. Notice: this should be in a dedicated binary factory
	// however:
	// - clas12 would be the only map using it. Not really worth it, also considering that:
	// - gemc is moving to gemc3 with a better (plugin) mechanism to load fields
	virtual void loadFieldMap(gclas12BinaryMappedField*, double); // load clas12 binary field map.
		
	static fieldFactory *createFieldFactory() {
		return new clas12BinField;
	}
	
	double solenoidScale = 1;
	double torusScale    = 1;

	double solenoidOrigin[3] = {0, 0, 0};
	double torusOrigin[3]    = {0, 0, 0};


};


#endif
