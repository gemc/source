/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

// gemc
#include "ActionInitialization.h"

ActionInitialization::ActionInitialization(goptions* go, map<string, double> *gPars) : G4VUserActionInitialization()
{
	genAction = new MPrimaryGeneratorAction(go);
	evtAction = new MEventAction(*go, *gPars);
	stpAction = new MSteppingAction(*go);
}


ActionInitialization::~ActionInitialization()
{}


void ActionInitialization::BuildForMaster() const
{
	// SetUserAction(new RunAction);
}


void ActionInitialization::Build() const
{
	SetUserAction(genAction);
	SetUserAction(evtAction);
	SetUserAction(stpAction);
}

