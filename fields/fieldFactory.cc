// gemc headers
#include "fieldFactory.h"
#include "asciiField.h"
#include "clas12BinField.h"
#include "gemcUtils.h"


fieldFactory *getFieldFactory(map<string, fieldFactoryInMap> *fieldsFactoryMap, string fieldsMethod)
{
	if(fieldsFactoryMap->find(fieldsMethod) == fieldsFactoryMap->end()) {
		cout << endl << endl << "  >>> WARNING: " << fieldsMethod << " NOT FOUND IN Field Factory Map." << endl;
		return nullptr;
	}
	
	return (*fieldsFactoryMap)[fieldsMethod]();
}

map<string, fieldFactoryInMap> registerFieldFactories()
{
	map<string, fieldFactoryInMap> fieldFactoryMap;
	
	// ASCII factory
	fieldFactoryMap["ASCII"] = &asciiField::createFieldFactory;

	// CLAS12BinaryMap factory
	fieldFactoryMap["CLAS12BIN"] = &clas12BinField::createFieldFactory;

	return fieldFactoryMap;
}


map<string, gfield> loadAllFields(map<string, fieldFactoryInMap> fieldFactoryMap, goptions opts)
{
	double verbosity = opts.optMap["FIELD_VERBOSITY"].arg ;
	// get list of files in directories in:
	//
	// - GEMC_DATA_DIR environment variable if exists
	// - FIELD_DIR gemc option if set other than "env", otherwise FIELD_DIR environment variable if exists
	map<string, string> filesMap;

	if(getenv("GEMC_DATA_DIR") != nullptr) mergeMaps(filesMap, getFilesInDirectory((string) getenv("GEMC_DATA_DIR") ));
	if(opts.optMap["FIELD_DIR"].args != "env") {
		mergeMaps(filesMap, getFilesInDirectory(opts.optMap["FIELD_DIR"].args));
	} else {
		if(getenv("FIELD_DIR") != nullptr) mergeMaps(filesMap, getFilesInDirectory((string) getenv("FIELD_DIR") ));
	}

	// checking eligibility of each file
	// if eligible, load field definitions
	map<string, gfield> gfields;

	for(map<string, string>::iterator it = filesMap.begin(); it != filesMap.end(); it++) {

		// if factory exist, calling isEligible
		fieldFactory *thisFactory = getFieldFactory(&fieldFactoryMap, it->second);

		if ( verbosity > 0 ) {
			if(thisFactory != nullptr) {
				if(thisFactory->isEligible(it->first)) {
					cout << "  > Field Factory: Loading file " << it->first << " using factory: " << it->second << endl;
				} else {
					cout << "  > Field Factory: Candidate file " << it->first << " cannot be loaded by factory: " << it->second << endl;
				}
			}
		}

		if(thisFactory != nullptr) {
			if(thisFactory->isEligible(it->first)) {
				gfield gf = thisFactory->loadField(it->first, opts);
				gf.fFactory = thisFactory;

				// if symmetry is set, it's probably a good field
				if(gf.symmetry != "na") {
					gfields[gf.name] = gf;
					if(verbosity > 0) cout << gfields[gf.name] << endl;
				}
			}
		}
		
	}
	
	return gfields;
}

