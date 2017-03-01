// c++ headers
#include <set>

// gemc headers
#include "mirrors_factory.h"
#include "mysql_mirrors.h"
#include "text_mirrors.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;


mirrors *getMirrorFactory(map<string, mirrorFactory> *factory, string mirrorsMethod)
{
	
	if(factory->find(mirrorsMethod) == factory->end())
	{
		cout << endl << endl << "  >>> WARNING: " << mirrorsMethod << " NOT FOUND IN Mirror Factory Map." << endl;
		return NULL;
	}
	
	return (*factory)[mirrorsMethod]();
}

map<string, mirrorFactory> registerMirrorFactories()
{
	
	map<string, mirrorFactory> mirrorMethodMap;
	
	// MYSQL initialization
	mirrorMethodMap["MYSQL"] = &mysql_mirrors::createMirrors;
	
	// TEXT initialization
	mirrorMethodMap["TEXT"] = &text_mirrors::createMirrors;
		
	return mirrorMethodMap;
}

void printMirrors(map<string, mirror*> mirMap)
{
	for(map<string, mirror*>::iterator it = mirMap.begin(); it != mirMap.end(); it++)
		cout << "    - mirror: >" << it->first << "< >" << it->second->name <<"<" << endl;
}


void mirror::opticalsFromString(string s, string what)
{
	stringstream comps(s);
	
	while(!comps.eof())
	{
		string c;
		comps >> c ;
		string trimmedC = trimSpacesFromString(c);
		if(what != "none")
		{
			if(what == "photonEnergy")
				photonEnergy.push_back(get_number(trimmedC));
			
			if(what == "indexOfRefraction")
				indexOfRefraction.push_back(get_number(trimmedC));
			
			if(what == "reflectivity")
				reflectivity.push_back(get_number(trimmedC));
			
			if(what == "efficiency")
				efficiency.push_back(get_number(trimmedC));
			
			if(what == "specularlobe")
				specularlobe.push_back(get_number(trimmedC));
			
			if(what == "specularspike")
				specularspike.push_back(get_number(trimmedC));
			
			if(what == "backscatter")
				backscatter.push_back(get_number(trimmedC));

			if(what == "sigmaAlpha")
				sigmaAlpha = get_number(trimmedC);
			else
				sigmaAlpha = 0;
		}
	}
	
	if(what == "backscatter")
	{
		// backscatter is the last vector to be loaded
		// now we can check the vector sizes for comparison
		// if no match, resetting quantities
		if(indexOfRefraction.size() != photonEnergy.size()) indexOfRefraction.clear();
		if(reflectivity.size()      != photonEnergy.size())	reflectivity.clear();
		if(efficiency.size()        != photonEnergy.size())	efficiency.clear();
		if(specularlobe.size()      != photonEnergy.size())	specularlobe.clear();
		if(specularspike.size()     != photonEnergy.size())	specularspike.clear();
		if(backscatter.size()       != photonEnergy.size())	backscatter.clear();
	}
}




// Load all mirrors coming from MYSQL, TEXT factories
map<string, mirror*> buildMirrors(map<string, mirrorFactory> mirrorFactoryMap, goptions go, runConditions rc)
{
	// Loading MYSQL def
	mirrors *mirrorSelectedFactory = getMirrorFactory(&mirrorFactoryMap, "MYSQL");
	map<string, mirror*> mirs = mirrorSelectedFactory->initMirrors(rc, go);

	// adding TEXT
	mirrors *textFactory = getMirrorFactory(&mirrorFactoryMap, "TEXT");
	map<string, mirror*> textMirs = textFactory->initMirrors(rc, go);
	for(map<string, mirror*>::iterator it = textMirs.begin(); it != textMirs.end(); it++)
		mirs[it->first] = it->second;

	return mirs;
}


ostream &operator<<(ostream &stream, mirror Mirror)
{
	cout  << endl;
	cout << "   Mirror name:  "  << Mirror.name        << "  -  " <<  Mirror.desc << endl;
	cout << "   type:         "  << Mirror.type        << endl;
	cout << "   finish:       "  << Mirror.finish      << endl;
	cout << "   model:        "  << Mirror.model       << endl;
	cout << "   border:       "  << Mirror.border      << endl;
	cout << "   sigmaAlpha:   "  << Mirror.sigmaAlpha  << endl;
	cout << endl;
	
	return stream;
}


