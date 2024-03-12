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
	// we assume the files specified in the gcard exist
	return 1;
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
			gf.bc12map->defineNamesAndType(fieldDir);
		}
	}
	
	// initialize field and bc12map field map
	gf.initialize(opts);

	if(gf.bc12map->interpolation == "none") {
		setAlgorithm(NEAREST_NEIGHBOR);
	} else {
		setAlgorithm(INTERPOLATION);
	}

	
	string hardcodedBinaryTorusOptionName = "binary_torus";
	string hardcodedBinarySolenOptionName = "binary_solenoid";

	vector<aopt> FIELD_SCALES_OPTION = opts.getArgs("SCALE_FIELD");
	for (unsigned int f = 0; f < FIELD_SCALES_OPTION.size(); f++) {
		vector < string > scales = getStringVectorFromStringWithDelimiter(FIELD_SCALES_OPTION[f].args, ",");
		if(scales.size() == 2) {
			if (scales[0].find(hardcodedBinaryTorusOptionName) != string::npos) {
				torusScale = get_number(scales[1]);
			} else if (scales[0].find(hardcodedBinarySolenOptionName) != string::npos) {
				solenoidScale = get_number(scales[1]);
			}
		}
	}
	
	vector<aopt> FIELD_DISPLACEMENT_OPTION = opts.getArgs("DISPLACE_FIELDMAP");
	for (unsigned int f = 0; f < FIELD_DISPLACEMENT_OPTION.size(); f++) {
		vector < string > displacement = getStringVectorFromStringWithDelimiter(FIELD_DISPLACEMENT_OPTION[f].args, ",");
		if(displacement.size() == 4) {
			if (displacement[0].find(hardcodedBinaryTorusOptionName) != string::npos) {
				torusOrigin[0] = get_number(displacement[1]);
				torusOrigin[1] = get_number(displacement[2]);
				torusOrigin[2] = get_number(displacement[3]);
			} else if (displacement[0].find(hardcodedBinarySolenOptionName) != string::npos) {
				solenoidOrigin[0] = get_number(displacement[1]);
				solenoidOrigin[1] = get_number(displacement[2]);
				solenoidOrigin[2] = get_number(displacement[3]);
			}
		}
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

	b12map->solenoidPtr = initializeSolenoid(b12map->solenoidMapFileName.c_str());
	b12map->torusPtr    = initializeTorus(b12map->torusMapFileName.c_str());

	b12map->solenoidPtr->scale = solenoidScale;
	b12map->torusPtr->scale    = torusScale;

	b12map->solenoidPtr->shiftX = solenoidOrigin[0]/cm;
	b12map->solenoidPtr->shiftY = solenoidOrigin[1]/cm;
	b12map->solenoidPtr->shiftZ = solenoidOrigin[2]/cm;

	b12map->torusPtr->shiftX  = torusOrigin[0]/cm;
	b12map->torusPtr->shiftY  = torusOrigin[1]/cm;
	b12map->torusPtr->shiftZ  = torusOrigin[2]/cm;
		

	cout << endl << "  ####  Binary Field Maps for " << b12map->identifier << " loading complete." << endl << endl;

}





