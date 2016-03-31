#ifndef ActionInitialization_h
#define ActionInitialization_h 1

// gemc
#include "MPrimaryGeneratorAction.h"
#include "MEventAction.h"
#include "MSteppingAction.h"
#include "options.h"


// geant4
#include "G4VUserActionInitialization.hh"


/// Action initialization class.
class ActionInitialization : public G4VUserActionInitialization
{
public:
	ActionInitialization(goptions *go, map<string, double> *gPars);
	virtual ~ActionInitialization();
	
	virtual void BuildForMaster() const;
	virtual void Build() const;
	
	
	MPrimaryGeneratorAction *genAction;
	MEventAction            *evtAction;
	MSteppingAction         *stpAction;
};


#endif


