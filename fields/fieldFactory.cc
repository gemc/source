// gemc headers
#include "fieldFactory.h"
#include "asciiField.h"
#include "utils.h"


fieldFactory *getFieldFactory(map<string, fieldFactoryInMap> *fieldsFactoryMap, string fieldsMethod)
{

	if(fieldsFactoryMap->find(fieldsMethod) == fieldsFactoryMap->end())
	{
		cout << endl << endl << "  >>> WARNING: " << fieldsMethod << " NOT FOUND IN Field Factory Map." << endl;
		return NULL;
	}
	
	return (*fieldsFactoryMap)[fieldsMethod]();
}

map<string, fieldFactoryInMap> registerFieldFactories()
{
	map<string, fieldFactoryInMap> fieldFactoryMap;
	
	// ASCII factory
	fieldFactoryMap["ASCII"] = &asciiField::createFieldFactory;
		
	return fieldFactoryMap;
}


map<string, gfield> loadAllFields(map<string, fieldFactoryInMap> fieldFactoryMap, goptions opts)
{
	double verbosity = opts.optMap["FIELD_VERBOSITY"].arg ;
	// get list of files in directories in:
	// 
	// - JLAB_ROOT/noarch/data
	// - GEMC_DATA_DIR if exists
	// - FIELD_DIR gemc option if set
	map<string, string> filesMap;

	if(getenv("JLAB_ROOT")     != NULL)        mergeMaps(filesMap, getFilesInDirectory((string) getenv("JLAB_ROOT") + "/noarch/data/" ));
	if(getenv("GEMC_DATA_DIR") != NULL)        mergeMaps(filesMap, getFilesInDirectory((string) getenv("GEMC_DATA_DIR") ));
	if(opts.optMap["FIELD_DIR"].args != "env") mergeMaps(filesMap, getFilesInDirectory(opts.optMap["FIELD_DIR"].args));


	// checking eligibility of each file
	// if eligible, load field definitions
	map<string, gfield> gfields;
	
	for(map<string, string>::iterator	it = filesMap.begin(); it != filesMap.end(); it++)
	{
		// if factory exist, calling isEligible
		fieldFactory *thisFactory = getFieldFactory(&fieldFactoryMap, it->second);

		if(thisFactory != NULL)
		{
			if(thisFactory->isEligible(it->first))
			{
				gfield gf = thisFactory->loadField(it->first, opts);
				gf.fFactory = thisFactory;
								
				// if symmetry is set, it's probably a good field
				if(gf.symmetry != "na")
				{
					gfields[gf.name] = gf;			
					if(verbosity > 0) cout << gfields[gf.name] << endl;
				}
			}
		}	
		// not done with the factory, cannot delete factory pointer
		// it's needed later for loading field maps
	}
	
	return gfields;
}

