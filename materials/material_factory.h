#ifndef MATERIAL_FACTORY_H
#define MATERIAL_FACTORY_H

// G4 headers
#include "G4Material.hh"

// C++ headers
#include <map>
#include <iostream>
using namespace std;

// gemc headers
#include "run_conditions.h"
#include "options.h"

class material
{
	
	public:
		material(){;}
		material(string n){name = n;}
		~material(){;}
		
		string name;
		string desc;
		double density;
		int ncomponents;
		
		vector<string> components;
		vector<double> fracs;
	
		vector<double> photonEnergy;
		vector<double> indexOfRefraction;
		vector<double> absorptionLength;
		vector<double> reflectivity;
		vector<double> efficiency;

		// scintillation
		vector<double> fastcomponent;
		vector<double> slowcomponent;
		double scintillationyield;
		double resolutionscale;
		double fasttimeconstant;
		double slowtimeconstant;
		double yieldratio;
		vector<double> rayleigh;
		double birkConstant;
	
		// load material components from DB entry
		void componentsFromString(string);

		// load optical properties from DB entry
		void opticalsFromString(string, string);

};


class materials
{
	public:
		// Pure Virtual Method to initialize G4 Materials
		virtual map<string, G4Material*> initMaterials(runConditions, goptions) = 0;
		map<string, G4Material*> materialsFromMap(map<string, material>);
		virtual ~materials(){}
};

typedef materials *(*materialFactory)();                                // Define materialFactory as a pointer to a function that returns a pointer 

materials *getMaterialFactory(map<string, materialFactory> *, string);  // returns materialFactory Function from Factory Map

map<string, materialFactory> registerMaterialFactories();               // Registers materialFactory in Factory Map

map<string, G4Material*> buildMaterials(map<string, materialFactory> materialFactoryMap, goptions go, runConditions rc);

// build material with standard isotopes.
map<string, G4Material*> materialsWithIsotopes();

void printMaterials(map<string, G4Material*> matMap);

#endif

