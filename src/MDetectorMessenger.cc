// gemc headers
#include "MDetectorMessenger.h"


MDetectorMessenger::MDetectorMessenger(MDetectorConstruction* Det) : MDetector(Det)
{
	gemcDir = new G4UIdirectory("/gemc/");
	gemcDir->SetGuidance("UI commands of gemc");
	
	UpdateGeoCmd = new G4UIcmdWithoutParameter("/gemc/updateGeo",this);
	UpdateGeoCmd->SetGuidance("Update detector geometry.");
	UpdateGeoCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
	UpdateGeoCmd->SetGuidance("if you changed geometrical value(s).");
	UpdateGeoCmd->AvailableForStates(G4State_Idle);
}

MDetectorMessenger::~MDetectorMessenger()
{
	delete gemcDir;
	delete UpdateGeoCmd;
}


void MDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == UpdateGeoCmd )
	{
		MDetector->updateGeometry();
	}
}

