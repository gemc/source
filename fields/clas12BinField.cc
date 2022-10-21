// gemc headers
#include "fieldFactory.h"
#include "clas12BinField.h"
#include "string_utilities.h"
#include "gemcUtils.h"

// c++
#include <string>
using namespace std;

extern "C" {
#include "magfield.h"
#include "magfieldio.h"
#include "munittest.h"
#include "magfieldutil.h"
}

// mlibrary
#include "gstring.h"
using namespace gstring;


bool clas12BinField::isEligible(string compositeFieldsName)
{
	if (compositeFieldsName ==  TorusSymmSolenoid2018 ) {
		return 1;
	} else if ( compositeFieldsName == TorusASymmSolenoid2018 ) {
		return 1;
	} else {
		return 0;
	}
}


// load field definitions
gfield clas12BinField::loadField(string file, goptions opts)
{
	gfield gf(opts);

	gf.name        = file;
	gf.description = "Field from David Heddle cMag library: " + file;
	gf.format      = "bc12map";
	gf.factory     = "CLAS12BIN";
	gf.integration = "G4ClassicalRK4";
	gf.minStep     = 0.01;
	gf.unit        = "kilogauss";
	gf.symmetry    = "cMag";

	if(gf.bc12map == nullptr) {
		gf.bc12map = new gclas12BinaryMappedField(file);
		if(getenv("FIELD_DIR") != nullptr) {
			string fieldDir=getenv("FIELD_DIR");
			gf.bc12map->symmetricTorusFileName = fieldDir + "/" + validC12MapNames[TorusSymmSolenoid2018][1];
			gf.bc12map->solenoidFileName       = fieldDir + "/" + validC12MapNames[TorusSymmSolenoid2018][0];
			gf.bc12map->fullTorusFileName      = fieldDir + "/" + validC12MapNames[TorusASymmSolenoid2018][1];
		}
	}
	
	// initialize field and bc12map field map
	gf.initialize(opts);

	if(gf.bc12map->interpolation == "none") {
		setAlgorithm(NEAREST_NEIGHBOR);
	} else {
		setAlgorithm(INTERPOLATION);
	}

	
	
	return gf;
}


// called by create_MFM
void clas12BinField::loadFieldMap(gclas12BinaryMappedField* b12map, double v) {

	
	if (v > 0 ) {
		cout << endl << "  #### Loading Binary Field Maps for " << b12map->identifier << endl;
	}

	// initialize map pointers
	b12map->combinedValuePtr = (FieldValuePtr) malloc(sizeof (FieldValue));
	if (b12map->identifier == TorusSymmSolenoid2018 ) {
	
		b12map->solenoid       = initializeSolenoid(b12map->solenoidFileName.c_str());
		b12map->symmetricTorus = initializeTorus(b12map->symmetricTorusFileName.c_str());
		
	} else if (b12map->identifier == TorusASymmSolenoid2018 ) {
	
		b12map->solenoid  = initializeSolenoid(b12map->solenoidFileName.c_str());
		b12map->fullTorus = initializeTorus(b12map->fullTorusFileName.c_str());
	}

	cout << endl << "  ####  Binary Field Maps for " << b12map->identifier << " loading complete." << endl << endl;

}





