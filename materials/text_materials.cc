// gemc headers
#include "material_factory.h"
#include "text_materials.h"
#include "string_utilities.h"
#include "utils.h"

// G4 headers
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"


map<string, G4Material*> text_materials::initMaterials(runConditions rc, goptions opts)
{
	
	string hd_msg    = opts.optMap["LOG_MSG"].args + " TEXT Materials Factory: >> ";
	double verbosity = opts.optMap["MATERIAL_VERBOSITY"].arg;
	
	map<string, material> mymats;  // material map

	// first check if there's at least one detector with TEXT factory
	if(!check_if_factory_is_needed(rc.detectorConditionsMap, "TEXT"))
		return materialsFromMap(mymats);
	
	
	// Looping over detectorConditionsMap for detector names
	// To each detector is associated a material and (optional) opt properties
	for(map<string, detectorCondition>::iterator it=rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++)
	{
		// building materials belonging to detectors that are tagged with MYSQL factory
		if(it->second.get_factory() != "TEXT")
			continue;

		if(verbosity)
			cout << hd_msg << " Initializing " << it->second.get_factory() << " for detector " << it->first << endl;
		
		// only add "main" if it's the main variation
		string dname     = it->first ;
		string variation = get_variation(it->second.get_variation());
		string filename  =  dname + "__materials_" + variation + ".txt";
		
		ifstream IN(filename.c_str());
		if(!IN)
		{
			// if file is not found, maybe it's in the GEMC_DATA_DIR directory
			if(getenv("GEMC_DATA_DIR")  != NULL)
			{
				
				string maybeHere = (string) getenv("GEMC_DATA_DIR") + "/" + filename;
				
				IN.open(maybeHere.c_str());
				if(!IN)
				{
					if(verbosity>1)
						cout << hd_msg << "Warning: The system >" << dname
						<< "< does not have a material file associated with it. "
						<< "Probably it's using default Geant4 parameters." << endl;
					continue;
				}
			}

			
			if(!IN)
			{
				if(verbosity>1)
					cout << hd_msg << "Warning: The system >" << dname
			    		           << "< does not have a material file associated with it. "
								   << "Probably it's using default Geant4 parameters." << endl;
				continue;
			}
		}
		
		// else loading parameters from file
		while(!IN.eof())
		{
			string dbline;
			getline(IN, dbline);
			
			if(!dbline.size())
			continue;
			
			gtable gt(get_strings(dbline, "|"));
			
			material thisMat(TrimSpaces(             gt.data[0])); // name
			thisMat.desc                =            gt.data[1];   // description
			thisMat.density             = get_number(gt.data[2]);  // density
			thisMat.ncomponents         = get_number(gt.data[3]);  // number of components
			thisMat.componentsFromString(            gt.data[4]);  // component + quantity list
			thisMat.opticalsFromString(              gt.data[5], "photonEnergy");
			thisMat.opticalsFromString(              gt.data[6], "indexOfRefraction");
			thisMat.opticalsFromString(              gt.data[7], "absorptionLength");
			thisMat.opticalsFromString(              gt.data[8], "reflectivity");
			thisMat.opticalsFromString(              gt.data[9], "efficiency");
			
			// scintillation
			// this condition is for backward compatibility,
			// scintillation was added with gemc 2.3
			if( gt.data.size() == 18)
			{
				thisMat.opticalsFromString(             gt.data[10], "fastcomponent");
				thisMat.opticalsFromString(             gt.data[11], "slowcomponent");
				thisMat.scintillationyield = get_number(gt.data[12]);
				thisMat.resolutionscale    = get_number(gt.data[13]);
				thisMat.fasttimeconstant   = get_number(gt.data[14]);
				thisMat.slowtimeconstant   = get_number(gt.data[15]);
				thisMat.yieldratio         = get_number(gt.data[16]);
				thisMat.opticalsFromString(             gt.data[17], "rayleigh");
			}
			mymats[thisMat.name] = thisMat;
			
		}
	}
	cout << endl;
	
	map<string, G4Material*> returnMap = materialsFromMap(mymats);
 	if(verbosity>0) printMaterials(returnMap);

	return returnMap;
}




