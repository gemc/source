#include "G4BertiniElectroNuclearBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4GammaGeneralProcess.hh"
#include "G4LossTableManager.hh"

#include "G4HadronicParameters.hh"

#include "GammaNuclearPhysics.h"


GammaNuclearPhysics::GammaNuclearPhysics(const G4String& name) :  G4VPhysicsConstructor(name) { }
GammaNuclearPhysics::~GammaNuclearPhysics() { }

void GammaNuclearPhysics::ConstructProcess()
{
   G4PhotoNuclearProcess* thePhotoNuclearProcess = new G4PhotoNuclearProcess();

	G4ElectronNuclearProcess *theElectronNuclearProcess = new G4ElectronNuclearProcess;
	G4ElectroVDNuclearModel *theElectroReaction = new G4ElectroVDNuclearModel;
	theElectronNuclearProcess->RegisterMe(theElectroReaction);

	G4PositronNuclearProcess* thePositronNuclearProcess = new G4PositronNuclearProcess;
	thePositronNuclearProcess->RegisterMe(theElectroReaction);


	G4CascadeInterface *theGammaReaction = new G4CascadeInterface;
	theGammaReaction->SetMaxEnergy(3.5*GeV);
	thePhotoNuclearProcess->RegisterMe(theGammaReaction);

	G4QGSMFragmentation *theFragmentation = new G4QGSMFragmentation;

	G4ExcitedStringDecay *theStringDecay = new G4ExcitedStringDecay(theFragmentation);
	G4QGSModel< G4GammaParticipants >  *theStringModel = new G4QGSModel< G4GammaParticipants >;
	theStringModel->SetFragmentationModel(theStringDecay);

	G4GeneratorPrecompoundInterface *theCascade = new G4GeneratorPrecompoundInterface;

	G4TheoFSGenerator *theModel = new G4TheoFSGenerator;
	theModel->SetTransport(theCascade);
	theModel->SetHighEnergyGenerator(theStringModel);

	G4ProcessManager * aProcMan = nullptr;

	theModel->SetMinEnergy(1.5*GeV);
	theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
	thePhotoNuclearProcess->RegisterMe(theModel);

	G4GammaGeneralProcess* sp = (G4GammaGeneralProcess*)G4LossTableManager::Instance()->GetGammaGeneralProcess();
	if(sp) {
		sp->AddHadProcess(thePhotoNuclearProcess);
	} else {
		aProcMan = G4Gamma::Gamma()->GetProcessManager();
		aProcMan->AddDiscreteProcess(thePhotoNuclearProcess);
	}

	aProcMan = G4Electron::Electron()->GetProcessManager();
	aProcMan->AddDiscreteProcess(theElectronNuclearProcess);

	aProcMan = G4Positron::Positron()->GetProcessManager();
	aProcMan->AddDiscreteProcess(thePositronNuclearProcess);


}
