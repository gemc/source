// gemc headers
#include "mirrors_factory.h"
#include "text_mirrors.h"
#include "string_utilities.h"
#include "gemcUtils.h"

// mlibrary
#include "gstring.h"
using namespace gstring;


map<string, mirror*> text_mirrors::initMirrors(runConditions rc, goptions opts)
{
	string hd_msg    = opts.optMap["LOG_MSG"].args + " TEXT Mirrors Factory: >> ";
	double verbosity = opts.optMap["MIRROR_VERBOSITY"].arg;
	
	map<string, mirror*> mymirs;  // mirror map

	// first check if there's at least one detector with TEXT factory
	if(!check_if_factory_is_needed(rc.detectorConditionsMap, "TEXT"))
		return mymirs;
		
	// Looping over detectorConditionsMap for detector names
	// To each detector maybe associated a mirror
	for(map<string, detectorCondition>::iterator it=rc.detectorConditionsMap.begin(); it != rc.detectorConditionsMap.end(); it++)
	{
		// building mirror belonging to detectors that are tagged with MYSQL factory
		if(it->second.get_factory() != "TEXT")
			continue;

		if(verbosity)
			cout << hd_msg << " Initializing " << it->second.get_factory() << " for detector " << it->first << endl;
		
		// only add "main" if it's the main variation
		string dname     = it->first ;
		string variation = get_variation(it->second.get_variation());
		string filename  =  dname + "__mirrors_" + variation + ".txt";
		
		ifstream IN(filename.c_str());
		if(!IN) {
			// if file is not found, maybe it's in the GEMC_DATA_DIR directory
			if(getenv("GEMC_DATA_DIR")  != nullptr) {
				
				string maybeHere = (string) getenv("GEMC_DATA_DIR") + "/" + filename;
				
				IN.open(maybeHere.c_str());
				if(!IN) {
					if(verbosity>1)
						cout << hd_msg << "Warning: The system >" << dname
						     << "< does not have a mirror file associated with it. " << endl;
					continue;
				}
			}
			
			if(!IN) {
				if(verbosity>1)
					cout << hd_msg << "Warning: The system >" << dname
			    		           << "< does not have a mirror file associated with it. " << endl;
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
			
			gtable gt(getStringVectorFromStringWithDelimiter(dbline, "|"));
			if( gt.data.size() >13 ) {
				mirror *thisMir = new mirror(trimSpacesFromString( gt.data[0])); // name
				thisMir->desc         =                  gt.data[1];   // description
				thisMir->type         =       trimSpacesFromString(gt.data[2]);  // type
				thisMir->finish       =       trimSpacesFromString(gt.data[3]);  // finish
				thisMir->model        =       trimSpacesFromString(gt.data[4]);  // model
				thisMir->border       =       trimSpacesFromString(gt.data[5]);  // border
				thisMir->maptOptProps =       trimSpacesFromString(gt.data[6]);  // material with optical properties
				thisMir->opticalsFromString(             gt.data[7],  "photonEnergy");
				thisMir->opticalsFromString(             gt.data[8],  "indexOfRefraction");
				thisMir->opticalsFromString(             gt.data[9],  "reflectivity");
				thisMir->opticalsFromString(             gt.data[10], "efficiency");
				thisMir->opticalsFromString(             gt.data[11], "specularlobe");
				thisMir->opticalsFromString(             gt.data[12], "specularspike");
				thisMir->opticalsFromString(             gt.data[13], "backscatter");

				// this condition is for backward compatibility,
				// sigma Alpha were added with gemc 2.6
				// 15 exact quantities
				if( gt.data.size() > 14) {
					thisMir->sigmaAlpha = get_number(gt.data[14]);
				}

				mymirs[thisMir->name] = thisMir;
			}





		}
	}
	cout << endl;
	
 	if(verbosity>0) printMirrors(mymirs);

	return mymirs;
}
