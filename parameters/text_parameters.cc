// gemc headers
#include "parameter_factory.h"
#include "text_parameters.h"
#include "string_utilities.h"
#include "utils.h"


map<string, double> text_parameters::loadParameters(goptions opts, runConditions RC)
{
	string hd_msg    = opts.optMap["LOG_MSG"].args + " TEXT Parameters: >> ";
	double verbosity = opts.optMap["PARAMETER_VERBOSITY"].arg;
	
	map<string, double> GParameters;   // parameters maps
	// first check if there's at least one detector with MYSQL factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return GParameters;
	
	// Looping over detectorConditionsMap for detector names
	// To each detector is associated a material table and an opt properties table
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++)
	{
		if(it->second.get_factory() != factoryType)
			continue;
		
		string dname     = it->first;
		string fname     = dname + "__parameters";
		string variation = get_variation(it->second.get_variation());
		
		if(verbosity)
			cout <<  hd_msg << " Importing Parameters for Detector: " <<  dname << " with " << factoryType << " factory, variation " << variation << endl;
		
		fname += "_" + variation + ".txt";
		ifstream IN(fname.c_str());
		if(!IN)
		{
			// file doesn't exist for this system.
			if(verbosity > 1)
				cout << hd_msg << "  Warning: failed to open " << fname << " for system: " << dname << ". Maybe the filename doesn't exist?" << endl;
			continue;
		}
		
		while(!IN.eof())
		{
			string dbline;
			getline(IN, dbline);
			
			if(!dbline.size())
				continue;
			
			gtable gt(get_strings(dbline, "|"));
			
			GParameters[gt.data[0]] = get_par_value(gt);
			
			if(verbosity > 1)
				log_value(gt, factoryType);
		}
		
		IN.close();
	}
	
	return GParameters;
}
