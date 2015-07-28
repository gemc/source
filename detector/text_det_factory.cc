// Qt headers
#include <QtSql>

// gemc headers
#include "text_det_factory.h"
#include "utils.h"

map<string, detector> text_det_factory::loadDetectors()
{
	string hd_msg     = gemcOpt.optMap["LOG_MSG"].args + " TEXT Factory: >> ";
	double verbosity  = gemcOpt.optMap["GEO_VERBOSITY"].arg;
	
	map<string, detector> dets;
	// first check if there's at least one detector with TEXT factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
	return dets;
	
	// building detectors that are tagged with TEXT factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++)
	{
		if(it->second.get_factory() != factoryType )
			continue;
		
		string dname     = it->first;
		string fname     = dname + "__geometry";
		string variation = get_variation(it->second.get_variation());
		
		if(verbosity)
			cout <<  hd_msg << " Importing Detector: " <<  dname << " with " << factoryType << " factory, variation " << variation << endl;
		
		fname += "_" + variation + ".txt";
		ifstream IN(fname.c_str());
		if(!IN)
		{
			
			// if file is not found, maybe it's in the GEMC_DATA_DIR directory
			if(getenv("GEMC_DATA_DIR")  != NULL)
			{
				string maybeHere = (string) getenv("GEMC_DATA_DIR") + "/" + fname;

				IN.open(maybeHere.c_str());
				if(!IN)
				{
					
					cout << hd_msg << "  Failed to open geometry file " << maybeHere << " for system: " << dname
					     << ". Maybe the filename doesn't exist? Exiting." << endl;
					exit(0);
				}
			}
			
			if(!IN)
			{
				
				cout << hd_msg << "  Failed to open geometry file " << fname << " for system: " << dname
					 << ". Maybe the filename doesn't exist? Exiting." << endl;
				exit(0);
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
			if( gt.data.size() != 18)
				cout << "ERROR: Incorrect number of geometry items (" << gt.data.size() << ") for " << gt.data[0] << endl;
			
			gt.add_data(dname);
			gt.add_data((string)"TEXT");
			gt.add_data(variation);
			
			// big warning if detector already exist
			// detector is NOT loaded if already existing
			if(dets.find(gt.data[0]) != dets.end())
			{
				cout << endl <<  " *** WARNING! A detector >" << gt.data[0]
				     << " exists already. Keeping original, not loading this instance. " << endl << endl;
			}
			else
			{
				dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
			}
			
		}
		IN.close();
	}
	
	return dets;
}





