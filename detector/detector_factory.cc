// gemc headers
#include "detector_factory.h"
#include "mysql_det_factory.h"
#include "text_det_factory.h"
#include "gdml_det_factory.h"
#include "clara_det_factory.h"

// c++ headers
#include <set>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

detectorFactory *getDetectorFactory(map<string, detectorFactoryInMap> *detectorFactoryMap, string fname)
{

	if(detectorFactoryMap->find(fname) == detectorFactoryMap->end())
	{
		cout << " *** WARNING: " << fname << " NOT FOUND IN  Detector Factory Map." << endl;
		return NULL;
	}
	
	return (*detectorFactoryMap)[fname]();
}


// factory registration and initialization
map<string, detectorFactoryInMap> registerDetectorFactory()
{
	// the key of the STL map also set factoryType in initFactory
	map<string, detectorFactoryInMap> dFactoryMap;

	// mysql factory
	dFactoryMap["MYSQL"] = &mysql_det_factory::createFactory;

	// text factory
	dFactoryMap["TEXT"]  = &text_det_factory::createFactory;

	// gdml factory
	dFactoryMap["GDML"]  = &gdml_det_factory::createFactory;
	
	// clara factory
	dFactoryMap["CLARA"]  = clara_det_factory::createFactory;
	
	return dFactoryMap;
}



map<string, detector> buildDetector(map<string, detectorFactoryInMap> detectorFactoryMap, goptions go, runConditions rc)
{
	map<string, detector> hallMap;
	
	// getting detector factories one by one
	for(map<string, detectorFactoryInMap>::iterator it = detectorFactoryMap.begin(); it != detectorFactoryMap.end(); it++)
	{
		// building detectors from this factory
		detectorFactory *thisFactory = getDetectorFactory(&detectorFactoryMap, it->first);
		
		// initialize factory with the key of the STL map
		thisFactory->initFactory(go, rc, it->first);
		
		// loading detectors
		map<string, detector> thisDMap = thisFactory->loadDetectors();
			
		// merging these detectors to hallMap
		for(map<string, detector>::iterator idet = thisDMap.begin(); idet != thisDMap.end(); idet++)
		{
			// big warning if detector already exist
			// detector is NOT loaded if already existing
			if(hallMap.find(idet->first) != hallMap.end())
			{
				cout << " *** WARNING! A detector >" << idet->first << "< in factory " << it->first << " exists already " << endl;
			}
			// loading detector if not present yet
			// assigning member factory to it
			else
			{
				hallMap[idet->first] = idet->second;
			}
		}
		
		// done with the factory, deleting factory pointer
		delete thisFactory;
	}

	// adding root mother volume
	string hall_mat   = go.optMap["HALL_MATERIAL"].args;
	string hall_field = go.optMap["HALL_FIELD"].args;
	vector<string> hall_dims  = get_strings(go.optMap["HALL_DIMENSIONS"].args, ",");
	if(hall_dims.size() != 3)
		cout << "   !!! Error: Hall dimensions is not a vector of 3 numbers. Example of dimensions: \"20*m, 20*m, 20*m\"" << endl;
	
	detector queen;
	queen.name        = "root";
	queen.mother      = "akasha";
	queen.description = "mother of us all";
	queen.pos         =  G4ThreeVector(0, 0, 0);
	queen.rot         =  G4RotationMatrix(G4ThreeVector(1, 0, 0),
                                        G4ThreeVector(0, 1, 0),
                                        G4ThreeVector(0, 0, 1));
	queen.type        =  "Box";
	queen.dimensions.push_back(get_number(hall_dims[0],1));
	queen.dimensions.push_back(get_number(hall_dims[1],1));
	queen.dimensions.push_back(get_number(hall_dims[2],1));
	queen.material    = hall_mat;
	queen.magfield    = hall_field;
	queen.exist       = 1;
	queen.visible     = 0;
	queen.ncopy       = 0;
	queen.scanned     = 1;
	hallMap[queen.name] = queen;

	// Transmitting the magnetic field properties along the genealogy
	// if they have mfield set to "no" and if the option argument NO_FIELD doesn't apply
	string nofield    = go.optMap["NO_FIELD"].args;
	
	// Transmitting the magnetic field properties along the genealogy if they are set to "no"
	for( map<string, detector>::iterator it =  hallMap.begin(); it!=hallMap.end() && nofield != "all" ; it++)
	{
		// if this is tagged for no field, continue
		if(it->first  == nofield)
			continue;
		
		if(it->second.magfield == "no")
		{
			// looking up the whole genealogy, until the first field is found
			string mother = it->second.mother;
			string firstAncestorFieldFound = "no";
			while(mother != "akasha" && firstAncestorFieldFound == "no")
			{
				if(hallMap[mother].magfield != "no")
				{
					firstAncestorFieldFound = hallMap[mother].magfield;
					it->second.magfield = hallMap[mother].magfield;
				}
				// moving up in genealogy
				mother = hallMap[mother].mother;
			}
		}
	}
	
	return hallMap;
}

// returns a string that log if the factory requested are present or not
string check_factory_existance(map<string, detectorFactoryInMap> detectorFactoryMap, runConditions rc)
{
		
	set<string> requested;
	set<string> present;
	
	string frequested;
	string fnotfound;

 	// logging requested factories
	for(map<string, detectorCondition>::iterator it = rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++)
	{
		requested.insert(it->second.get_factory());
	}

	// logging present factories
	for(map<string, detectorFactoryInMap>::iterator it = detectorFactoryMap.begin(); it != detectorFactoryMap.end(); it++)
		present.insert(it->first);

	int found_all = 1;
	for(set<string>::iterator it = requested.begin(); it != requested.end(); it++ )
	{
		frequested.append(*it + " ");
		if(present.find(*it) == present.end())
		{
			found_all = 0;
			fnotfound.append(*it + " ");
		}
	}


	if(!found_all)
		return string(" *** WARNING: These detector factories were requested but not found: " + fnotfound);
		
	return string(" >> All detector factories requested are found: " + frequested);

}

void detectorFactory::initFactory(goptions go, runConditions rc, string ft)
{
	if(gemcOpt.optMap["LOG_VERBOSITY"].arg > 0)
		cout << "  > " << ft << " Detector Factory is Initialized "  << endl;
	
	factoryType = ft;
	gemcOpt     = go;
	RC          = rc;
}




// returns detector from a gtable
detector get_detector(gtable gt, goptions go, runConditions RC)
{
	if(gt.data.size() < 18)
	{
		cout << " !!! ERROR: Detector data size should be at least 18. There are " << gt.data.size() << " items on the line for " << gt.data[0] << endl;
		exit(0);
	}	
	string hd_msg     = " >> TEXT Factory: ";
	double verbosity  = go.optMap["GEO_VERBOSITY"].arg;
	
	string catch_v    = go.optMap["CATCH"].args;

	detector det;
	
	// 0,1,2: Id, Mother, Description
	det.name        = gt.data[0];
	det.mother      = gt.data[1];
	det.description = gt.data[2];

	// 3: Position Vector
	det.pos = calc_position(gt.data[3]);

	// 4: Rotation Vector
	det.rot = calc_rotation(gt.data[4], det.name);
				
	// Checking for displacements and rotation from nominal position
	if(RC.detectorConditionsMap.find(det.name) != RC.detectorConditionsMap.end())
	{
		// Adding gcard displacement for this detector if non zero
		G4ThreeVector shiftp = RC.detectorConditionsMap[det.name].get_position();
					
		if(shiftp.mag2() != 0)
		{
			if(verbosity > 3 || det.name.find(catch_v))
				cout << hd_msg << " Detector " << det.name << " is displaced by: " << shiftp/cm << " cm" << endl;
						
			det.pos += shiftp;
		}
						
    	// Adding gcard rotation for this detector if if non zero
		G4ThreeVector more_rot = RC.detectorConditionsMap[det.name].get_vrotation();
		if(more_rot.mag2() != 0)
		{
			if(verbosity > 3 || det.name.find(catch_v))
				cout << hd_msg << " Detector " << det.name << " is rotated by: " << more_rot/deg << " deg" << endl;
						
			det.rot.rotateX(more_rot.x());
			det.rot.rotateY(more_rot.y());
			det.rot.rotateZ(more_rot.z());
		}
	}  // end of checking displacement
				
				
	// 5: Color, opacity
	if(gt.data[5].size() != 6 && gt.data[5].size() != 7)
	{
		cout << hd_msg << " Color Attributes for " << det.name << "<" << gt.data[5] << ">  have wrong size: " << gt.data[5].size()
									 << ". It should be 6 or 7 digits  rrggbb[t]  (red, green, blue hexadecimals + optional transparency)." << endl;
		exit(9);
	}

	G4Colour thisCol = gcol(gt.data[5]);

	// 6: Solid Type
	det.type = gt.data[6];

	// 7: Dimensions
	stringstream vars(gt.data[7]);
	string var;
	while(!vars.eof())
	{
		vars >> var;
		det.dimensions.push_back(get_number(var,1));
	}

	// 8: Material
	det.material =  gt.data[8];
	// resetting Material if asked
	vector<aopt> changeMatOptions = go.getArgs("CHANGEVOLUMEMATERIALTO");
	for (unsigned int f = 0; f < changeMatOptions.size(); f++)
	{
		vector < string > VolumeNewMats = get_strings(changeMatOptions[f].args, ",");
		if(VolumeNewMats.size() == 2)
		{
			// VolumeNewMats[0] = volume name
			// VolumeNewMats[1] = new material
			if(det.name == TrimSpaces(VolumeNewMats[0]))
				det.material = TrimSpaces(VolumeNewMats[1]);
		}
	}



	// 9: Magnetic Field
	det.magfield = gt.data[9];

	// 10: copy number
	det.ncopy   = atoi(gt.data[10].c_str());

	// 11: pMany
	det.pMany   = atoi(gt.data[11].c_str());

	// 12: Activation flag
	det.exist   = atoi(gt.data[12].c_str());
			
	// Overwriting existance is set in the gcard
	if(RC.detectorConditionsMap.find(det.name) != RC.detectorConditionsMap.end())
	{
		det.exist = RC.detectorConditionsMap[det.name].get_existance();
		if(verbosity > 3 || det.name.find(catch_v))
			cout << hd_msg << " Detector " << det.name << " has existance set to: " << det.exist << endl;
	}
				
	// 13:  Visibility
	det.visible = atoi(gt.data[13].c_str());

	// 14: Style
	det.style   = atoi(gt.data[14].c_str());

	// Setting the Visual Atributes Color, Visibility, Style
	det.VAtts = G4VisAttributes(thisCol);
	det.visible ? det.VAtts.SetVisibility(true) : det.VAtts.SetVisibility(false);
	if(det.visible)
		det.style   ? det.VAtts.SetForceSolid(true) : det.VAtts.SetForceWireframe(true);


	// 15: sensitivity
	det.sensitivity = gt.data[15];

	if(det.sensitivity != "no")
	{
		// 16: hitType
		det.hitType = gt.data[16];

		// 17: identity
		det.identity = get_identifiers(gt.data[17]);
	}
	// removing Sensitivity if asked
	vector < string > vnames = get_strings(go.optMap["REMOVESENSITIVITY"].args, ",");
	for(unsigned vtr = 0; vtr < vnames.size(); vtr++)
	{
		if(det.name == TrimSpaces(vnames[vtr]))
		{
			det.sensitivity = "no";
			det.hitType = "";
			det.identity.clear();
		}
	}



	// 18 detector system
	det.system  = gt.data[18];
	det.factory = gt.data[19];


	if(gt.data.size() > 20)
		det.variation = gt.data[20];

	if(gt.data.size() > 21)
		det.run     = atoi(gt.data[21].c_str());

	if(verbosity>2)
		cout << det;	
	return det;
}





