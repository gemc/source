#ifndef MIRRORS_FACTORY_H
#define MIRRORS_FACTORY_H

// C++ headers
#include <map>
#include <iostream>
using namespace std;

// gemc headers
#include "run_conditions.h"
#include "options.h"

class mirror
{
	
	public:
		mirror(){;}
		mirror(string n){name = n;}
		~mirror(){;}
		
		string name;
		string desc;
		string type;
		string finish;
		string model;
		string border;
		double sigmaAlpha;

		string maptOptProps;
	
		vector<double> photonEnergy;
		vector<double> indexOfRefraction;
		vector<double> reflectivity;
		vector<double> efficiency;
		vector<double> specularlobe;
		vector<double> specularspike;
		vector<double> backscatter;
	
		// load optical properties from DB entry
		void opticalsFromString(string, string);

		///< Overloaded "<<" for mirror class. Dumps infos on screen.
		friend ostream &operator<<(ostream &stream, mirror);
	
};


class mirrors
{
	public:
		// Pure Virtual Method to initialize gemc mirrors
		virtual map<string, mirror*> initMirrors(runConditions, goptions) = 0;
		virtual ~mirrors(){}
};

// Define mirrorFactory as a pointer to a function that returns a pointer
typedef mirrors *(*mirrorFactory)();

// returns mirrorFactory Function from Factory Map
mirrors *getMirrorFactory(map<string, mirrorFactory> *, string);

// Registers mirrorFactory in Factory Map
map<string, mirrorFactory> registerMirrorFactories();

// build all mirrors from all factories
map<string, mirror*> buildMirrors(map<string, mirrorFactory> mirrorFactoryMap, goptions go, runConditions rc);

void printMirrors(map<string, mirror*> mirMap);

#endif

