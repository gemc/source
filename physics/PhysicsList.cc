// GEMC PHYSICS LIST
//
// -----------------
//
// The



// c++ headers
#include <iostream>
using namespace std;

// gemc headers
#include "PhysicsList.h"
#include "PhysicsListMessenger.h"
#include "string_utilities.h"

// geant4 headers
#include "G4LossTableManager.hh"
#include "G4PhysListFactory.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
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
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4DecayWithSpin.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

PhysicsList::PhysicsList(goptions opts) : G4VModularPhysicsList()
{

	gemcOpt = opts;
	verbosity       = gemcOpt.optMap["PHYS_VERBOSITY"].arg ;
	ingredientsList = gemcOpt.optMap["PHYSICS"].args ;
	
	// default physics lists
	hadronicPhys = "none";
	EMPhys       = "none";
	opticalPhys  = "none";
	HPPhys       = "none";
	
	// loading hadronic and em physics lists
	G4PhysListFactory factory;
	g4HadronicList = factory.AvailablePhysLists();
	// em comes with "_", stripping it
	vector<G4String> allEM       = factory.AvailablePhysListsEM();
	g4EMList.push_back("STD");
	for(unsigned i=0; i<allEM.size(); i++)
	{
		vector<string> emstripped = get_strings(allEM[i], "_");
		string stripped;
		
		if(emstripped.size() == 2) stripped = TrimSpaces(emstripped[1]);
		if(stripped != "")
			g4EMList.push_back(stripped);
		
	}
	
	G4LossTableManager::Instance();
	defaultCutValue = gemcOpt.optMap["PRODUCTIONCUT"].arg;
	cutForGamma     = defaultCutValue;
	cutForElectron  = defaultCutValue;
	cutForPositron  = defaultCutValue;
	cutForProton    = defaultCutValue;
	
	physIngredients = get_strings(ingredientsList, "+");
	
	g4EMPhysics    = NULL;
	g4ParticleList = NULL;
	g4HadronicPhysics.clear();
	
	// validateIngredients will also set hadronicPhys, EMPhys, opticalPhys
	if(!validateIngredients())
	{
		cout << "  !!! Error: physics ingredients list not valid: >" << ingredientsList << "<" << endl;
		list();
		
		cout << " Exiting." << endl;
		exit(0);
	}
	else
		cookPhysics();

}




PhysicsList::~PhysicsList()
{
	;
}


void PhysicsList::list()
{	
	cout << "     > Available hadronic physics list: " << endl;
	for(unsigned i=0; i<g4HadronicList.size(); i++)
		cout << "       - " << g4HadronicList[i] << endl;
	cout << endl;
	
	cout << "     > Available EM physics list: " << endl;
	for(unsigned i=0; i<g4EMList.size(); i++)
		cout << "       - " << g4EMList[i] << endl;
	
	cout << "   > Optica: optical" << endl;
}



// check the ingredients consistency
bool PhysicsList::validateIngredients()
{
	unsigned isHadronicLegit = 0;
	unsigned isEMLegit       = 0;
	unsigned isOpticalLegit  = 0;
	unsigned isHPLegit       = 0;
	
	G4PhysListFactory factory;
	for(unsigned i=0; i<physIngredients.size(); i++)
	{
		string ingredient = TrimSpaces(physIngredients[i]);
		
		if(factory.IsReferencePhysList(ingredient))
		{
			isHadronicLegit = 1;
			hadronicPhys = ingredient;
		}

		for(unsigned i=0; i<g4EMList.size(); i++)
			if(ingredient == g4EMList[i])
			{
				isEMLegit = 1;
				EMPhys = ingredient;
			}
		
		if(ingredient == "Optical")
		{
			isOpticalLegit = 1;
			opticalPhys = "yes";
		}
		if(ingredient == "HP")
		{
			HPPhys    = "yes";
			isHPLegit = 1;
		}
	}
	
	if(verbosity > 0)
	{
		cout << "  >> Physics: " << ingredientsList << endl;
		cout << "   > Hadronic: " << hadronicPhys  << endl;
		cout << "   > EM: " << EMPhys  << endl;
		cout << "   > HP: " << HPPhys  << endl;
		cout << "   > Optical: " << opticalPhys << endl << endl;
	}
	
	if(physIngredients.size() == (isHadronicLegit + isEMLegit + isOpticalLegit + isHPLegit))
		return TRUE;
	
	return FALSE;
}




void PhysicsList::SetCuts()
{
	
	if (verbosity >0)
	{
		cout << "PhysicsList::SetCuts:";
		cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << endl;
	}
	
	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");
	SetCutValue(cutForProton, "proton");
	
	if (verbosity>0) DumpCutValuesTable();
}


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


#include "FTFP_BERT.hh"
#include "FTFP_BERT_TRV.hh"
#include "FTFP_BERT_HP.hh"
#include "FTF_BIC.hh"
#include "LBE.hh"
#include "QBBC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_FTFP_BERT.hh"
#include "QGS_BIC.hh"
#include "QGSP_INCLXX.hh"
#include "Shielding.hh"
#include "G4OpticalPhysics.hh"

#include "G4StepLimiter.hh"

void PhysicsList::cookPhysics()
{
	// Particles
	g4ParticleList = new G4DecayPhysics("decays");
	
	// EM Physics
	if(g4EMPhysics) delete  g4EMPhysics;
	     if(EMPhys == "STD")  g4EMPhysics = new G4EmStandardPhysics();
	else if(EMPhys == "EMV")  g4EMPhysics = new G4EmStandardPhysics_option1();
	else if(EMPhys == "EMX")  g4EMPhysics = new G4EmStandardPhysics_option2();
	else if(EMPhys == "EMY")  g4EMPhysics = new G4EmStandardPhysics_option3();
	else if(EMPhys == "EMZ")  g4EMPhysics = new G4EmStandardPhysics_option4();
	else if(EMPhys == "LIV")  g4EMPhysics = new G4EmLivermorePhysics();
	else if(EMPhys == "PEN")  g4EMPhysics = new G4EmPenelopePhysics();				
	else if(EMPhys != "none")
	{
		cout << " !! Wrong EMPhys " << EMPhys << endl << "Exiting." << endl;
		exit(0);
	}
	  
	// Hadronic Physics
	// always adding these
	// this is a general version of
	// Hadr01 example
	// LBE and QBBC were removed
	
	// em extra physics always there
	if(hadronicPhys != "none")
	{
		g4HadronicPhysics.push_back( new G4EmExtraPhysics(verbosity));
		
		if(HPPhys == "yes")
		{
			cout << "   >  Loading High Precision Cross Sections... this may take a while..." << endl;
			g4HadronicPhysics.push_back( new G4HadronElasticPhysicsHP(verbosity) );
		}
		else
			g4HadronicPhysics.push_back( new G4HadronElasticPhysics(verbosity) );
		
		// binary cascade, bertini models, or standard
		// ion physics
		if(hadronicPhys.find("BIC") != string::npos)
			g4HadronicPhysics.push_back( new G4IonBinaryCascadePhysics(verbosity));
		else if(hadronicPhys.find("BERT") != string::npos)
			g4HadronicPhysics.push_back( new G4IonPhysics(verbosity));
		
		
		g4HadronicPhysics.push_back( new G4NeutronTrackingCut(verbosity));
	}
	// adding the hadronic physics list
 	// I don't understand why there isn't a factory method for this
	// had to hardcode the cases
	// in any case they are in G4PhysListFactory
	
	     if(hadronicPhys == "FTFP_BERT")      {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT());}
	else if(hadronicPhys == "FTFP_BERT_TRV")  {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT_TRV());}
	else if(hadronicPhys == "FTFP_BERT_HP")   {g4HadronicPhysics.push_back( new G4HadronPhysicsFTFP_BERT_HP());}
	else if(hadronicPhys == "FTF_BIC")        {g4HadronicPhysics.push_back( new G4HadronPhysicsFTF_BIC());}
	else if(hadronicPhys == "QGSP_BERT")      {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BERT());}
	else if(hadronicPhys == "QGSP_BERT_HP")   {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BERT_HP());}
	else if(hadronicPhys == "QGSP_BIC")       {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BIC());}
	else if(hadronicPhys == "QGSP_BIC_HP")    {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_BIC_HP());}
	else if(hadronicPhys == "QGSP_FTFP_BERT") {g4HadronicPhysics.push_back( new G4HadronPhysicsQGSP_FTFP_BERT());}
	else if(hadronicPhys == "QGS_BIC")        {g4HadronicPhysics.push_back( new G4HadronPhysicsQGS_BIC());}
	else if(hadronicPhys == "none")            {;}
	//else if(hadronicPhys == "Shielding")      {g4HadronicPhysics.push_back( new Shielding());}
	else {cout << "Wrong hadronicPhys " << hadronicPhys << endl << "Exiting." << endl; exit(0);}
	
	// optical physics
	// taken from example: optical/LXe
	if(opticalPhys == "yes")
	{
		// verbosity is set to zero at the constructor level by default
		// see G4OpticalPhysics.hh
		G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
		opticalPhysics->SetWLSTimeProfile("delta");
		
		g4HadronicPhysics.push_back(opticalPhysics);
	}
	
}


void PhysicsList::ConstructParticle()
{
	string cosmics =  gemcOpt.optMap["COSMICRAYS"].args;
	vector<string> csettings = get_info(cosmics, string(",\""));
	string decayType = "default";
	string particleType;
	int len = csettings.size();
	if(csettings[0] == "default"){
	  if(len>4){
	    particleType = csettings[3];
	    decayType = csettings[4];
	  }
	}else{
	  if(len>6){
	    particleType = csettings[5];
	    decayType = csettings[6];
	  }
	}
	// warn if muon radiative decay is selected but the simulated
        // cosmic rays is not a muon
	if(decayType=="radiative" && particleType!="muon") 
	  cout << "!!! Check COSMICRAYS data card, muon radiative decay required but no muon being simulated " << endl;
	
	g4ParticleList->ConstructParticle();
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
	
	G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
	G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
	
	if(decayType == "radiative"){
	  MuonPlusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+",1.00));
	  MuonMinusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-",1.00));
	}else{ // default
	  MuonPlusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu+",0.986));
	  MuonPlusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
	  MuonMinusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu-",0.986));
	  MuonMinusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-",1.00));
  }
	G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);
	G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);
}


void PhysicsList::ConstructProcess()
{
	AddTransportation();
	theDecayProcess = new G4DecayWithSpin();
	G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
	G4VProcess* decay;

	if(g4EMPhysics)
		g4EMPhysics->ConstructProcess();
	
	g4ParticleList->ConstructProcess();
	
	for(size_t i=0; i<g4HadronicPhysics.size(); i++)
		g4HadronicPhysics[i]->ConstructProcess();



	// PhysicsList contains theParticleIterator
	theParticleIterator->reset();

	while( (*theParticleIterator)() )
	{

		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager*     pmanager = particle->GetProcessManager();
		string                pname    = particle->GetParticleName();
		decay = processTable->FindProcess("Decay",particle);      
	
		// Adding Step Limiter
		if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0) && (pname != "chargedgeantino"))
		{
			if(verbosity > 2)
				cout << "   >  Adding Step Limiter for " << pname << endl;
			
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
		}

		if (theDecayProcess->IsApplicable(*particle)) {
		  if(decay) pmanager->RemoveProcess(decay);
		  pmanager->AddProcess(theDecayProcess);
		  pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
		  pmanager->SetProcessOrderingToLast(theDecayProcess, idxAtRest);
		}
	}
}










