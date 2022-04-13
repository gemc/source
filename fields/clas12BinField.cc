// gemc headers
#include "fieldFactory.h"
#include "clas12BinField.h"
#include "string_utilities.h"
#include "gemcUtils.h"
#include "magfieldio.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "magfieldio.h"

// mlibrary
#include "gstring.h"
using namespace gstring;


bool clas12BinField::isEligible(string compositeFieldsName)
{

	if (compositeFieldsName ==  TorusSymmSolenoid2018 || compositeFieldsName == TorusASymmSolenoid2018) {
		return 1;
	}

	return 0;
}


// load field definitions
gfield clas12BinField::loadField(string file, goptions opts)
{
	gfield gf(opts);

	gf.name        = file;
	gf.description = "Field from David Heddle cMag library";
	gf.format      = "bc12map";
	gf.factory     = "CLAS12BIN";
	gf.integration = "G4ClassicalRK4";
	gf.minStep     = 0.01;
	gf.unit        = "kilogauss";
	gf.symmetry    = "cMag";

	// initialize field and field map
	gf.initialize(opts);

	return gf;
}



// load field map
void clas12BinField::loadFieldMap(gMappedField* map, double v)
{
	cout << "  > Loading field map from " << map->identifier << " with symmetry: " << map->symmetry << endl;

	// actual load of map: call function below
//	loadFieldMap(map, v);
}


void clas12BinField::loadFieldMap(gclas12BinaryMappedField* map, double v) {


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


	vector<string> c12symm = validC12MapNames["c12BinaryTorusSymmSolenoid2018"];
	vector<string> c12asymm = validC12MapNames["c12BinaryTorusASymmSolenoid2018"];

	string a = c12symm[1];
	string b = c12asymm[2];
	string c = c12symm[2];

	const char *dir1 = a.c_str();
	const char *dir2 = b.c_str();
	const char *dir3 = c.c_str();

	solenoid = initializeSolenoid(dir1);
	symmetricTorus = initializeTorus(dir2);
	fullTorus = initializeTorus(dir3);

}





