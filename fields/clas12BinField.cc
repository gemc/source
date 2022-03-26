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


bool clas12BinField::isEligible(string file)
{

	// check that file name is one of the items in the list
	// loop over map check that file is one of the keys of the map.

	
	  //Check if file name is consistent with any of the symmetric Torus files in FIELD_DIR                                                                                                                                                     
  	  for(int i = 0; i < validC12MapNames["c12BinaryTorusSymmSolenoid2018"].size();i++)  {
		  vector<string> c12vec = validC12MapNames["c12BinaryTorusSymmSolenoid2018];
		  string c12string = c12vec[i];
                  if(file == c12sgring) return 0;
          }

          //Check if file name is consistent with any of the asymmetric Torus files in FIELD_DIR                                                                                                                                                    
          for(int i = 0; i < validC12MapNames["c12BinaryTorusASymSolenoid2018"].size(); i++) {
		  vector<string> c12vec = validC12MapNames["c12BinaryTorusASymmSolenoid2018"];
		  string c12string = c12vec[i];
                  if(file == c12string) return 0;
          }

	return 1;
}


// load field definitions
gfield clas12BinField::loadField(string file, goptions opts)
{
	gfield gf(opts);

	gf.name        = assignAttribute(e, "name", "na");
	gf.factory     = assignAttribute(e, "factory", "na");
	gf.description = assignAttribute(e, "comment", "no comment");
	gf.symmetry    = assignAttribute(e, "type",   "na");
	gf.format      = assignAttribute(e, "format", "bc12map");
	gf.factory     = assignAttribute(e, "factory", "CLAS12BIN");
	gf.unit        = assignAttribute(e, "unit", "kG");

	// fill all these properties according to David's definitions

	//	string integration;     ///< Integration Method
	//	double verbosity;       ///< Log verbosity
	//	double minStep;         ///< Minimum Step for the G4ChordFinder
	//	string unit;            ///< Field Unit

	// initialize field and field map
	gf.initialize(opts);

	return gf;
}



// load field map
void clas12BinField::loadFieldMap(gMappedField* map, double v)
{
	cout << "  > Loading field map from " << map->identifier << " with symmetry: " << map->symmetry << endl;

	// actual load of map: call function below

}


void clas12BinField::loadFieldMap(gclas12BinaryMappedField* map, double v) {


	// use validC12MapNames instead of solenoidPath and torusSymmetricPath

	const char *dataDir;
	
	vector<string> c12symm = validC12MapNames["c12BinaryTorusSymmSolenoid2018];
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



// Examples

/*
 
 clas12 solenoid:
 <mfield>
 <description name="clas12-solenoid" factory="ASCII" comment="clas12 superconducting solenoid"/>
 <symmetry type="cylindrical-z" format="map"/>
 <map>
 <coordinate>
 <first  name="transverse"    npoints="601"   min="0"  max="3" units="m"/>
 <second name="longitudinal"  npoints="1201"  min="-3" max="3" units="m"/>
 </coordinate>
 <field unit="T"/>
 </map>
 </mfield>
 clas12 torus:
 <mfield>
 <description name="clas12-torus" factory="ASCII" comment="clas12 superconducting torus"/>
 <symmetry type="phi-segmented" format="map""/>
 <map>
 <coordinate>
 <first  name="azimuthal"     npoints="61"   min="0"   max="30"  units="deg"/>
 <second name="transverse"    npoints="126"  min="0"   max="500" units="cm"/>
 <third  name="longitudinal"  npoints="126"  min="100" max="600" units="cm"/>
 </coordinate>
 <field unit="kilogauss"/>
 </map>
 </mfield>
 Example of uniform field:
 <mfield>
 <description name="uniform" factory="ASCII" comment="Uniform 10 T Magnetic Field along x-axis"/>
 <symmetry type="uniform" format="simple"/>
 <dimension bx="10" by="0" bz="0" units="T"/>
 </mfield>
 */



