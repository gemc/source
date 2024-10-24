// c++ headers
#include <iostream>
using namespace std;

// gemc headers
#include "PhysicsList.h"
#include "PhysicsListMessenger.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// geant4 headers
#include "G4LossTableManager.hh"
#include "G4PhysListFactory.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4DecayTable.hh"
#include "G4ProcessTable.hh"

// geant4 physics headers
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4MuonDecayChannelWithSpin.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

PhysicsList::PhysicsList(goptions opts) : G4VModularPhysicsList()
{

	gemcOpt = opts;
	verbosity       = gemcOpt.optMap["PHYS_VERBOSITY"].arg ;
	ingredientsList = gemcOpt.optMap["PHYSICS"].args ;
	int fastmcMode  = gemcOpt.optMap["FASTMCMODE"].arg;
	defaultCutValue = gemcOpt.optMap["PRODUCTIONCUT"].arg;
	if(fastmcMode > 0) {
		defaultCutValue = 5000;
	}

	// copy of G4String ss[23] in G4PhysListFactory.cc
	availableHadronicPhysic = {
		"FTFP_BERT","FTFP_BERT_TRV","FTFP_BERT_ATL","FTFP_BERT_HP","FTFQGSP_BERT",
		"FTFP_INCLXX","FTFP_INCLXX_HP","FTF_BIC", "LBE","QBBC",
		"QGSP_BERT","QGSP_BERT_HP","QGSP_BIC","QGSP_BIC_HP","QGSP_BIC_AllHP",
		"QGSP_FTFP_BERT","QGSP_INCLXX","QGSP_INCLXX_HP","QGS_BIC",
		"Shielding","ShieldingLEND","ShieldingM","NuBeam"};


	// default physics lists
	hadronicPhys  = "none";
	EMPhys        = "none";
	opticalPhys   = "none";
	HPPhys        = "none";
	photoNuclear  = "none";

	// making sure vectors are empty (although, by construction they should be)
	g4HadronicPhysics.clear();
	g4EMList.clear();

	// loading hadronic and em physics lists
	G4PhysListFactory factory;
	g4HadronicList = factory.AvailablePhysLists();

	// default EM list
	g4EMList.push_back("STD");

	// em comes with "_", stripping it
	vector<G4String> allEM = factory.AvailablePhysListsEM();
	for(unsigned i=0; i<allEM.size(); i++) {
		vector<string> emstripped = getStringVectorFromStringWithDelimiter(allEM[i], "_");
		string stripped;

		if(emstripped.size() == 2)
		stripped = trimSpacesFromString(emstripped[1]);
		else if (emstripped.size() == 1)
		stripped = emstripped[0];
		else
		continue;

		g4EMList.push_back(stripped);
	}

	G4LossTableManager::Instance();

	cutForGamma     = defaultCutValue;
	cutForElectron  = defaultCutValue;
	cutForPositron  = defaultCutValue;
	cutForProton    = defaultCutValue;

	physIngredients = getStringVectorFromStringWithDelimiter(ingredientsList, "+");

	// G4VPhysicsConstructor
	g4EMPhysics = nullptr;
	g4DecayPhys = nullptr;

	// validateIngredients will also set hadronicPhys, EMPhys, opticalPhys
	if(!validateIngredients()) {
		cout << "  !!! Error: physics ingredients list not valid: >" << ingredientsList << "<" << endl;
		list();

		cout << " Exiting." << endl;
		exit(1);
	} else {
		cookPhysics();
	}
}



PhysicsList::~PhysicsList() { ; }


void PhysicsList::list()
{
	cout << "     > Available hadronic physics list: " << endl;
	for(unsigned i=0; i<g4HadronicList.size(); i++)
		cout << "       - " << g4HadronicList[i] << endl;
	cout << endl;

	cout << "     > Available EM physics list: " << endl;
	for(unsigned i=0; i<g4EMList.size(); i++)
		cout << "       - " << g4EMList[i] << endl;

	cout << "   > Optical physics: " << opticalPhys << endl;
	cout << "   > Photo-Nuclear: "   << photoNuclear << endl;
}



// check the ingredients consistency
bool PhysicsList::validateIngredients()
{
	unsigned isHadronicLegit     = 0;
	unsigned isEMLegit           = 0;
	unsigned isOpticalLegit      = 0;
	unsigned isHPLegit           = 0;
	unsigned isPhotoNuclearLegit = 0;

	G4PhysListFactory factory;
	for(unsigned i=0; i<physIngredients.size(); i++) {

		// checking if the hadronic list is registered in the factory
		string ingredient = trimSpacesFromString(physIngredients[i]);

		if(factory.IsReferencePhysList(ingredient)) {
			isHadronicLegit = 1;
			hadronicPhys = ingredient;
		}

		// checking that the EM is is registered in the factory
		// through our stripped EMList array
		for(unsigned i=0; i<g4EMList.size(); i++) {
			if(ingredient == g4EMList[i]) {
				isEMLegit = 1;
				EMPhys = ingredient;
			}
		}

		// if optical is present, activate it
		if(ingredient == "Optical") {
			isOpticalLegit = 1;
			opticalPhys = "yes";
		}

		// if HP is present, activate it
		if(ingredient == "HP") {
			HPPhys    = "yes";
			isHPLegit = 1;
		}

		// if photoNuclear is present, activate it
		if(ingredient == "PhotoNuclear") {
			isPhotoNuclearLegit = 1;
			photoNuclear = "yes";
		}

	}

	if(verbosity > 0) {
		cout << "  >> Physics: " << ingredientsList << endl;
		cout << "   > Hadronic: " << hadronicPhys  << endl;
		cout << "   > EM: " << EMPhys  << endl;
		cout << "   > HP: " << HPPhys  << endl;
		cout << "   > Optical: " << opticalPhys << endl;
		cout << "   > Photo-Nuclear: " << photoNuclear << endl << endl;
	}

	if(physIngredients.size() == (isHadronicLegit + isEMLegit + isOpticalLegit + isHPLegit + isPhotoNuclearLegit)) {
		return TRUE;
	}

	return FALSE;
}


void PhysicsList::SetCuts()
{

	if (verbosity >0) {
		cout << "PhysicsList::SetCuts:";
		cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << endl;
	}

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma,    "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");
	SetCutValue(cutForProton,   "proton");

	if (verbosity>0) DumpCutValuesTable();
}


// required in G4VModularPhysicsList
void PhysicsList::SetCutForGamma(double cut)
{
	cutForGamma = cut;
	SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}
void PhysicsList::SetCutForElectron(double cut)
{
	cutForElectron = cut;
	SetParticleCuts(cutForElectron, G4Electron::Electron());
}
void PhysicsList::SetCutForPositron(double cut)
{
	cutForPositron = cut;
	SetParticleCuts(cutForPositron, G4Positron::Positron());
}
void PhysicsList::SetCutForProton(double cut)
{
	cutForProton = cut;
	SetParticleCuts(cutForProton, G4Proton::Proton());
}


// the list of physics list can be found in G4PhysListFactory.cc
// the includes should be loaded from the entries in each entry in open /source/physics_lists/constructors/hadron_inelastic/src/G4HadronPhysics*

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "G4HadronPhysicsFTFP_BERT_ATL.hh"
#include "G4HadronPhysicsFTFQGSP_BERT.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4HadronPhysicsQGS_BIC.hh"
#include "G4HadronPhysicsShielding.hh"
#include "G4HadronPhysicsShieldingLEND.hh"
#include "G4HadronPhysicsNuBeam.hh"

// Stopping physics
#include "G4StoppingPhysics.hh"


// Photonuclear processes
#include "GammaNuclearPhysics.h"

// Optical Physics
#include "G4OpticalPhysics.hh"
#include "G4SynchrotronRadiation.hh"
#include "G4SynchrotronRadiationInMat.hh"

#include "G4StepLimiter.hh"



void PhysicsList::cookPhysics()
{
	// Particles
	g4DecayPhys = new G4DecayPhysics("decays");

	// EM Physics
	// see also https://geant4-userdoc.web.cern.ch/UsersGuides/PhysicsListGuide/html/electromagnetic/index.html
	if(g4EMPhysics) delete  g4EMPhysics;
	if(EMPhys == "STD")  g4EMPhysics = new G4EmStandardPhysics();
	else if(EMPhys == "EMV")  g4EMPhysics = new G4EmStandardPhysics_option1();
	else if(EMPhys == "EMX")  g4EMPhysics = new G4EmStandardPhysics_option2();
	else if(EMPhys == "EMY")  g4EMPhysics = new G4EmStandardPhysics_option3();
	else if(EMPhys == "EMZ")  g4EMPhysics = new G4EmStandardPhysics_option4();
	else if(EMPhys == "LIV")  g4EMPhysics = new G4EmLivermorePhysics();
	else if(EMPhys == "LVP")  g4EMPhysics = new G4EmLivermorePolarizedPhysics();
	else if(EMPhys == "PEN")  g4EMPhysics = new G4EmPenelopePhysics();
	else if(EMPhys == "LEM")  g4EMPhysics = new G4EmLowEPPhysics();
	else if(EMPhys != "none") {
		cout << " !! Wrong EMPhys " << EMPhys << endl << "Exiting." << endl;
		exit(1);
	}

	// Hadronic Physics
	// The lists in FTFP_BERT() functions also include the EM, hadron elastic and so on. Here we separate them.
	// This is a general version of Hadr01 example

	// em extra physics always there
	if(hadronicPhys != "none") {
		g4HadronicPhysics.push_back( new G4EmExtraPhysics(verbosity));

		// Hadron Elastic Physics
		if(HPPhys == "yes") {
			cout << "   >  Loading High Precision Cross Sections... this may take a while..." << endl;
			g4HadronicPhysics.push_back( new G4HadronElasticPhysicsHP(verbosity) );
		} else {
 			g4HadronicPhysics.push_back( new G4HadronElasticPhysics(verbosity) );
		}

		// binary cascade, bertini models, or standard
		// ion physics
		if(hadronicPhys.find("BIC") != string::npos) {
			g4HadronicPhysics.push_back( new G4IonBinaryCascadePhysics(verbosity));
		} else if(hadronicPhys.find("BERT") != string::npos) {
			g4HadronicPhysics.push_back( new G4IonPhysics(verbosity));
		}

		g4HadronicPhysics.push_back( new G4StoppingPhysics(verbosity));
		g4HadronicPhysics.push_back( new G4NeutronTrackingCut(verbosity));
	}

	// adding the hadronic physics list
	// the complete list is in G4PhysListFactory.cc, however those are G4VModularPhysicsList
	// here we use the the G4VPhysicsConstructor directly. The includes are in "include/Geant4"
	// To get the constructor: open /source/physics_lists/constructors/hadron_inelastic/src/G4HadronPhysics*
	// PRAGMA TODO:
	// Notice that the (int) constructor sets the QuasiElastic physics to false


	if(     hadronicPhys == "FTFP_BERT")      {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
	else if(hadronicPhys == "FTFP_BERT_HP")   {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT_HP(verbosity));}
	else if(hadronicPhys == "FTFP_BERT_TRV")  {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT_TRV(verbosity));}
	else if(hadronicPhys == "FTFP_BERT_ATL")  {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT_ATL(verbosity));}
	else if(hadronicPhys == "FTFQGSP_BERT")   {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFQGSP_BERT(verbosity));}
//	else if(hadronicPhys == "FTFP_INCLXX")    {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
//	else if(hadronicPhys == "FTFP_INCLXX_HP") {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
	else if(hadronicPhys == "FTF_BIC")        {g4HadronicPhysics.push_back( new G4HadronPhysicsFTF_BIC(verbosity));}
//	else if(hadronicPhys == "LBE")            {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
//	else if(hadronicPhys == "QBBC")           {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
	else if(hadronicPhys == "QGSP_BERT")      {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BERT(verbosity));}
	else if(hadronicPhys == "QGSP_BERT_HP")   {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BERT_HP(verbosity));}
	else if(hadronicPhys == "QGSP_BIC")       {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BIC(verbosity));}
	else if(hadronicPhys == "QGSP_BIC_HP")    {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BIC_HP(verbosity));}
	else if(hadronicPhys == "QGSP_BIC_AllHP") {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BIC_AllHP(verbosity));}
	else if(hadronicPhys == "QGSP_FTFP_BERT") {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_FTFP_BERT(verbosity));}
//	else if(hadronicPhys == "QGSP_INCLXX")    {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
//	else if(hadronicPhys == "QGSP_INCLXX_HP") {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
	else if(hadronicPhys == "QGS_BIC")        {g4HadronicPhysics.push_back( new G4HadronPhysicsQGS_BIC(verbosity));}
	else if(hadronicPhys == "Shielding")      {g4HadronicPhysics.push_back( new G4HadronPhysicsShielding(verbosity));}
	else if(hadronicPhys == "ShieldingLEND")  {g4HadronicPhysics.push_back( new G4HadronPhysicsShieldingLEND(verbosity));}
//	else if(hadronicPhys == "ShieldingM")     {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT(verbosity));}
	else if(hadronicPhys == "NuBeam")         {g4HadronicPhysics.push_back( new G4HadronPhysicsNuBeam(verbosity));}
	else if(hadronicPhys == "none")           {;}
	else {
		cout << " > " << hadronicPhys << " is not supported in this version of GEMC yet. Exiting." << endl;
		exit(1);
	}


	// optical physics
	// taken from example: extended/optical/LXe
	if(opticalPhys == "yes") {
		// verbosity is set to zero at the constructor level by default
		// see G4OpticalPhysics.hh

		G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
		// enable Mie scattering (do not think it is enabled by default)
		// make this optional?
		opticalPhysics->Configure(kMieHG, true);
		
		// this was deprecated?
		// Method G4OpticalPhysics::SetWLSTimeProfile is deprecated.
		// Use G4OpticalParameters::SetWLSTimeProfile(G4String) instead.
		// but we did try it? Maybe a bug?
		//opticalPhysics->SetWLSTimeProfile("delta");

		g4HadronicPhysics.push_back(opticalPhysics);
	}

	if(photoNuclear == "yes") {
		g4HadronicPhysics.push_back(new GammaNuclearPhysics("GammaNuclearPhysics"));
	}
}


void PhysicsList::ConstructParticle()
{
	muonRadDecay = 0;
	muonRadDecay = gemcOpt.optMap["FORCE_MUON_RADIATIVE_DECAY"].arg;

	g4DecayPhys->ConstructParticle();

	if(muonRadDecay){
		G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
		G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
		MuonPlusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+",1.00));
		MuonMinusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-",1.00));
		G4MuonPlus::MuonPlusDefinition()->SetDecayTable(MuonPlusDecayTable);
		G4MuonMinus::MuonMinusDefinition()->SetDecayTable(MuonMinusDecayTable);
	}

}


void PhysicsList::ConstructProcess()
{
	AddTransportation();
	int fastmcMode = gemcOpt.optMap["FASTMCMODE"].arg;
	int synrad     = gemcOpt.optMap["SYNRAD"].arg;

	if(fastmcMode%10 < 2) {
		G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
		G4VProcess* decay;

		if(g4EMPhysics) {
			g4EMPhysics->ConstructProcess();
		}

		g4DecayPhys->ConstructProcess();

		for(size_t i=0; i<g4HadronicPhysics.size(); i++)
			g4HadronicPhysics[i]->ConstructProcess();

		// sync radiation
		G4SynchrotronRadiation*      fSync    = nullptr;
		G4SynchrotronRadiationInMat* fSyncMat = nullptr;

		if (synrad == 1) fSync    = new G4SynchrotronRadiation();
		else if (synrad == 2) fSyncMat = new G4SynchrotronRadiationInMat();
		//G4AutoDelete::Register(fSync);

		auto theParticleIterator = GetParticleIterator();

		// PhysicsList contains theParticleIterator
		theParticleIterator->reset();

		while( (*theParticleIterator)() ) {
			G4ParticleDefinition* particle = theParticleIterator->value();
			G4ProcessManager*     pmanager = particle->GetProcessManager();
			string                pname    = particle->GetParticleName();

			// Adding Step Limiter
			if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0)) {
				if(verbosity > 2) {
					cout << "   >  Adding Step Limiter for " << pname << endl;
				}

				pmanager->AddProcess(new G4StepLimiter,       -1,-1, 3);
			}

			// G4SynchrotronRadiation if requested
			if (synrad == 1) {

				if (pname == "e-") {
					//electron
					pmanager->AddProcess(fSync,               -1,-1, 4);
					pmanager->AddProcess(new G4StepLimiter,   -1,-1, 5);

				} else if (pname == "e+") {
					//positron
					pmanager->AddProcess(fSync,              -1,-1, 5);
					pmanager->AddProcess(new G4StepLimiter,  -1,-1, 6);
				}
			} else if (synrad == 2) {

				if (pname == "e-") {
					//electron
					pmanager->AddProcess(fSyncMat,            -1,-1, 4);
					pmanager->AddProcess(new G4StepLimiter,   -1,-1, 5);

				} else if (pname == "e+") {
					//positron
					pmanager->AddProcess(fSyncMat,           -1,-1, 5);
					pmanager->AddProcess(new G4StepLimiter,  -1,-1, 6);
				}
			}

			if(muonRadDecay){
				G4DecayWithSpin* theDecayProcess = new G4DecayWithSpin();
				decay = processTable->FindProcess("Decay",particle);
				if (theDecayProcess->IsApplicable(*particle)) {
					if(decay) pmanager->RemoveProcess(decay);
					pmanager->AddProcess(theDecayProcess);
					pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
					pmanager->SetProcessOrderingToLast(theDecayProcess, idxAtRest);
				}
			}
		}
	}
}
