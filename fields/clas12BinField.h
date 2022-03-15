#ifndef C12BIN_FIELDS_H
#define C12BIN_FIELDS_H

#include "fieldFactory.h"

extern "C" {
#include "magfield.h"	
#include "magfieldio.h"
#include "munittest.h"
#include "magfieldutil.h"
}

class clas12BinField : public fieldFactory
{
public:
	~clas12BinField(){}

	// save map list for loadFieldMap
	map<string, vector<string> > validC12MapNames;

	// constructor:  create map list with hardcoded names
	clas12BinField() {
		validC12MapNames["c12BinaryTorusSymmSolenoid2018"]  = {"Symm_Solenoid_r601...", "Symm_torus_r2501..."};
		validC12MapNames["c12BinaryTorusASymmSolenoid2018"] = {"Symm_Solenoid_r601...", "Full_torus_r251..."};
	}

	// check if the binary map filename is a match for a pre-defined list
	bool isEligible(string);


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
	virtual void loadFieldMap(gclas12BinaryMappedField*, double); // load clas12 binary field map.
	
	MagneticFieldPtr symmetricTorus;
	MagneticFieldPtr solenoid;
	MagneticFieldPtr fullTorus;
	
	static fieldFactory *createFieldFactory() {
		return new asciiField;
	}
	
};


#endif
