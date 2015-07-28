#ifndef PHYSICS_LIST_H
#define PHYSICS_LIST_H 1

// gemc headers
#include "options.h"

// geant4 headers
#include "G4VModularPhysicsList.hh"

// c++ headers
#include <string>
using namespace std;


// the ingredient list is a combination of names
// such as "FTFP_BERT + em_opt1 + optical"
//
// valid physics ingredients
// - any G4 phys list
// - em_std, em_opt1, em_opt2, em_opt3, em_opt4
// - optical
class PhysicsList: public G4VModularPhysicsList
{
public:
	
	PhysicsList(goptions);
	~PhysicsList();

	void list();
	
	vector<string> physIngredients;
	bool validateIngredients();
	
	// SetCuts is required by G4VModularPhysicsList
	// will actually call the individual particles sets
	void SetCuts();
	void SetCutForGamma(double);
	void SetCutForElectron(double);
	void SetCutForPositron(double);
	void SetCutForProton(double);

	
private:

	double cutForGamma;
	double cutForElectron;
	double cutForPositron;
	double cutForProton;

	goptions gemcOpt;   ///< gemc options map

	double verbosity;
	string ingredientsList;
	string hadronicPhys;
	string EMPhys;
	string opticalPhys;
	string HPPhys;
	
	vector<G4String> g4HadronicList;
	vector<G4String> g4EMList;
	
	
	G4VPhysicsConstructor*  g4EMPhysics;
	G4VPhysicsConstructor*  g4ParticleList;
	vector<G4VPhysicsConstructor*>  g4HadronicPhysics;

	// build the geant4 physics
	void cookPhysics();
	void ConstructParticle();
	void ConstructProcess();

};



#endif
