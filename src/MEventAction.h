/// \file MEventAction.h
/// Derived from G4UserEventAction.\n
/// The two functions:
/// - BeginOfEventAction
/// - EndOfEventAction
/// define the routines at the begin
/// and at the end of each event.\n
/// In EndOfEventAction the output is processed and
/// written out (if the output option is selected).\n
/// The event number is initialized and increased.
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef MEventAction_H
#define MEventAction_H 1

// G4 headers
#include "G4UserEventAction.hh"

// gemc headers
#include "outputFactory.h"
#include "gbank.h"
#include "HitProcess.h"
#include "sensitiveDetector.h"
#include "options.h"
#include "MPrimaryGeneratorAction.h"


/// \class TInfos
/// <b> TInfos: Track information of particles in sensitive detectors1 </b>\n\n
/// Contains:
/// - particle ID (int from PDG encoding).
/// - mother particle ID.
/// - vertex
/// - mother particle vertex.\n
/// For ions:
/// Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.\n
/// For a nucleus consisting of np protons and nn neutrons \n
/// A = np + nn and Z = np.\n
/// I gives the isomer level, with I = 0 corresponding \n
/// to the ground state and I >0 to excitations \n
class TInfos
{
	public:
		TInfos(){;}
		TInfos(int MTID)
		{
			mtid = MTID;
			mpid = 0;
			mv   = G4ThreeVector(0.,0.,0.);
		}
		~TInfos(){;}

	public:
		int mtid;  // mother track id
		int mpid;  // mother PID
	
		G4ThreeVector mv;
};

vector<int> vector_mtids(  map<int, TInfos> tinfos, vector<int> tids);
vector<int> vector_mpids(  map<int, TInfos> tinfos, vector<int> tids);
vector<G4ThreeVector> vector_mvert(  map<int, TInfos> tinfos, vector<int> tids);
vector<int>           vector_zint(  int size);  ///< provides a vector of 0
vector<G4ThreeVector> vector_zthre( int size);  ///< provides a vector of (0,0,0)


/// \class MEventAction
/// <b> MEventAction </b>\n\n
/// Derived from G4UserEventAction.\n
/// The two functions:
/// - BeginOfEventAction
/// - EndOfEventAction
/// define the routines at the begin
/// and at the end of each event.\n
/// In EndOfEventAction the output is written out
/// (if the output option is selected)
class MEventAction : public G4UserEventAction
{
	public:
		MEventAction(goptions, map<string, double>);       ///< Constructor copies gemc options
		~MEventAction();                                   ///< Destructor
		
		goptions gemcOpt;                                  ///< gemc options
		
		outputContainer                  *outContainer;     ///< outputContainer class - contains the output format.
		map<string, outputFactoryInMap>  *outputFactoryMap; ///< outputFactory map
		map<string, sensitiveDetector*>  SeDe_Map;          ///< Sensitive detector Map
		map<string, HitProcess_Factory>  *hitProcessMap;    ///< Hit Process Routine Factory Map
		map<string, gBank>               *banksMap;         ///< Bank Map
		map<string, double>               gPars;            ///< Parameters Map
		MPrimaryGeneratorAction          *gen_action;       ///< Generator Action
    
		map<int, int> hierarchy;                     ///< Hierarchy map
		map<int, int> momDaugther;                   ///< mom - daughter relationship
		vector<int> vector_otids(vector<int> tids);  ///< return original track id of a vector of tid

	
		int    evtN;            ///< Event Number
		string hd_msg;          ///< Event Action Message
		int    Modulo;          ///< Print Log Event every Modulo
		double VERB;            ///< Event Verbosity
		string catch_v;         ///< Print Log for volume
		int   SAVE_ALL_MOTHERS; ///< Loops over the stored trajectories to store mother vertex and pid in the output
		int   MAXP;             ///< Max number of generated particles to save on output stream
		string WRITE_ALLRAW;    ///< List of detectors for which geant4 all raw info need to be saved
		string WRITE_INTRAW;    ///< List of detectors for which geant4 raw integrated info need to be saved
		string WRITE_INTDGT;    ///< List of detectors for which digitized integrated info need to be NOT saved
		string SIGNALVT;        ///< List of detectors for which voltage versus time need to be saved
		
	public:
		void BeginOfEventAction(const G4Event*);            ///< Routine at the start of each event
		void EndOfEventAction(const G4Event*);              ///< Routine at the end of each event
		void SetEvtNumber(int N){evtN = N;}                 ///< Sets Event Number
	
		runWeights rw;
	
};

#endif




