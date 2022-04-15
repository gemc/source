// gemc headers
#include "fieldFactory.h"
#include "clas12BinField.h"
#include "string_utilities.h"
#include "gemcUtils.h"
#include "magfieldio.h"

// c++
#include <string>
using namespace std;


// cMag
#include "magfieldio.h"

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

	return gf;
}


// called by create_MFM
void clas12BinField::loadFieldMap(gclas12BinaryMappedField* b12map, double v) {

	
	if (v > 0 ) {
		cout << "  > Loading Binary Field Maps for " << TorusSymmSolenoid2018 << endl;
	}

	if (b12map->identifier == TorusSymmSolenoid2018 ) {
	
		
		
		
	} else if (b12map->identifier == TorusSymmSolenoid2018 ) {
	
	
	}

	// use validC12MapNames instead of solenoidPath and torusSymmetricPath


//	char *solenoidPath = (char*) malloc(255);
//	char *torusSymmetricPath = (char*) malloc(255);
//	char *torusFullPath = (char*) malloc(255);
//
//	const char *dataDir;
//
//	dataDir = "/w/hallb_scshelf2102/clas12/jnewton/binary/data/fieldmaps";
//
//	sprintf(solenoidPath, "%s/Symm_solenoid_r601_phi1_z1201_13June2018.dat", dataDir);

//	if(isSymmetric==true)  {
//		sprintf(torusSymmetricPath, "%s/Symm_torus_r2501_phi16_z251_24Apr2018.dat",dataDir);//Absolute Path To Symmetric Torus
//		symmetricTorus = initializeTorus("/w/hallb_scshelf2102/clas12/jnewton/binary/data/fieldmaps/Symm_torus_r2501_phi16_z251_24Apr2018.dat");
//	}
//
//	else  {
//		sprintf(torusFullPath, "%s/Full_torus_r251_phi181_z251_03March2020.dat",dataDir);//Absolute Path To Full Torus
//		fullTorus = initializeTorus("/w/hallb_scshelf2102/clas12/jnewton/binary/data/fieldmaps/Full_torus_r251_phi181_z251_03March2020.dat");
//	}


//	vector<string> c12symm = validC12MapNames["c12BinaryTorusSymmSolenoid2018"];
//	vector<string> c12asymm = validC12MapNames["c12BinaryTorusASymmSolenoid2018"];
//
//	string a = c12symm[1];
//	string b = c12asymm[2];
//	string c = c12symm[2];
//
//	const char *dir1 = a.c_str();
//	const char *dir2 = b.c_str();
//	const char *dir3 = c.c_str();
//
//	solenoid = initializeSolenoid(dir1);
//	symmetricTorus = initializeTorus(dir2);
//	fullTorus = initializeTorus(dir3);

}





