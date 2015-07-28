// gemc headers
#include "PhysicsListMessenger.h"


// geant4 headers
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UImanager.hh"


PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys) :fPhysicsList(pPhys)
{
	fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutGamma",this);
	fGammaCutCmd->SetGuidance("Set gamma cut.");
	fGammaCutCmd->SetParameterName("Gcut",false);
	fGammaCutCmd->SetUnitCategory("Length");
	fGammaCutCmd->SetRange("Gcut>=0.0");
	fGammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	fElectCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutEl",this);
	fElectCutCmd->SetGuidance("Set electron cut.");
	fElectCutCmd->SetParameterName("Ecut",false);
	fElectCutCmd->SetUnitCategory("Length");
	fElectCutCmd->SetRange("Ecut>=0.0");
	fElectCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	fPosCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutPos",this);
	fPosCutCmd->SetGuidance("Set positron cut.");
	fPosCutCmd->SetParameterName("Pcut",false);
	fPosCutCmd->SetUnitCategory("Length");
	fPosCutCmd->SetRange("Pcut>=0.0");
	fPosCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	fCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutProt",this);
	fCutCmd->SetGuidance("Set proton cut.");
	fCutCmd->SetParameterName("ProtCut",false);
	fCutCmd->SetUnitCategory("Length");
	fCutCmd->SetRange("ProtCut>=0.0");
	fCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	fAllCutCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/CutsAll",this);
	fAllCutCmd->SetGuidance("Set cut for all.");
	fAllCutCmd->SetParameterName("cut",false);
	fAllCutCmd->SetUnitCategory("Length");
	fAllCutCmd->SetRange("cut>=0.0");
	fAllCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	fPListCmd = new G4UIcmdWithAString("/testhadr/Physics",this);
	fPListCmd->SetGuidance("Add modula physics list.");
	fPListCmd->SetParameterName("PList",false);
	fPListCmd->AvailableForStates(G4State_PreInit);
	
	fListCmd = new G4UIcmdWithoutParameter("/testhadr/ListPhysics",this);
	fListCmd->SetGuidance("Available Physics Lists");
	fListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

PhysicsListMessenger::~PhysicsListMessenger()
{
	delete fGammaCutCmd;
	delete fElectCutCmd;
	delete fPosCutCmd;
	delete fCutCmd;
	delete fAllCutCmd;
	delete fPListCmd;
	delete fListCmd;
}

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
}

