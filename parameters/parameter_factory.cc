// GEMC headers
#include "parameter_factory.h"
#include "mysql_parameters.h"
#include "text_parameters.h"
#include "string_utilities.h"
#include "utils.h"


parametersFactory *getParameterFactory(map<string, parameterFactoryInMap> *parametersFactoryMap, string parametersMethod)
{
	
	if(parametersFactoryMap->find(parametersMethod) == parametersFactoryMap->end())
	{
		cout << "  ** WARNING: " << parametersMethod << " NOT FOUND IN parameter Factory Map." << endl;
		return NULL;
	}
	
	return (*parametersFactoryMap)[parametersMethod]();
}

map<string, parameterFactoryInMap> registerParameterFactories()
{
	map<string, parameterFactoryInMap> parameterFactoryMap;
	
	// mysql factory
	parameterFactoryMap["MYSQL"] = &mysql_parameters::createParametersFactory;
	
	// text factory
	parameterFactoryMap["TEXT"]  = &text_parameters::createParametersFactory;
	
	return parameterFactoryMap;
}


map<string, string> getParametersMap(map<string, double> dataMap)
{
	map<string, string> parmap;
	
 	for(map<string, double>::iterator it = dataMap.begin(); it != dataMap.end(); it++)
	{
		string key = "parameter " + it->first;
		parmap[key] = stringify(it->second);
	}
	
	return parmap;
}


map<string, double> loadAllParameters(map<string, parameterFactoryInMap> parameterFactoryMap, goptions go, runConditions rc)
{
	double verbosity = go.optMap["PARAMETER_VERBOSITY"].arg;
	map<string, double> gParameters;
	
	// getting parameters factories one by one
	for(map<string, parameterFactoryInMap>::iterator it = parameterFactoryMap.begin(); it != parameterFactoryMap.end(); it++)
	{
		// building parameters from this factory
		parametersFactory *thisFactory = getParameterFactory(&parameterFactoryMap, it->first);
		
		// initialize factory
		thisFactory->initFactory(it->first);
		
		// loading parameters
		map<string, double> thisPMap = thisFactory->loadParameters(go, rc);
		
		// merging these detectors to hallMap
		for(map<string, double>::iterator ipar = thisPMap.begin(); ipar != thisPMap.end(); ipar++)
		{
			// big warning if parameter already exist
			// detector is NOT loaded if already existing
			if(gParameters.find(ipar->first) != gParameters.end() && verbosity)
			{
				cout << "  ** WARNING! Parameter " << ipar->first << " in factory " << it->first << " exists already! " << endl;
			}
			// loading parameter if not present yet
			else
				gParameters[ipar->first] = ipar->second;
		}
		
		// done with the factory, deleting factory pointer
		delete thisFactory;
	}
	
	return gParameters;
}

double get_par_value(gtable gt)
{
	string pvalue = gt.data[1];
	string punits = gt.data[2];
	return get_number(pvalue + "*" + punits);
}

void log_value(gtable gt, string factory)
{
	cout <<  "  > gemc " << factory << " Parameters: \"" + gt.data[0] +  "\" loaded with value " + gt.data[1] + gt.data[2] + ". Description: " + gt.data[3] << endl;
}





