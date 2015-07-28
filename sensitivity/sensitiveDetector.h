/// \file sensitiveDetector.h
/// Defines the gemc Sensitive Detector class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef sensitiveDetector_H
#define sensitiveDetector_H 1

// G4 headers
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"

// gemc headers
#include "sensitiveID.h"
#include "detector.h"
#include "Hit.h"
#include "HitProcess.h"

// C++ headers
#include <iostream>
#include <string>
#include <set>
using namespace std;


/// \class sensitiveDetector
/// <b> sensitiveDetector </b>\n\n
/// This is the gemc Sensitive Detector.\n
/// When a track enters a volume associated with this SD, the
/// ProcessHits routine is called. ProcessHits builds the MHit and
/// adds it in the MHitCollection.\n
/// At the end of event, each MHitCollection instantiate the
/// HitProcess_Factory relative to the Hit Process and calls its method ProcessHit.\n
class sensitiveDetector : public G4VSensitiveDetector
{
	public:
		sensitiveDetector(G4String, goptions, string factory, int run, string variation, string system);       ///< Constructor
		virtual ~sensitiveDetector();
		
		virtual void Initialize(G4HCofThisEvent*);                   ///< Virtual Method called at the beginning of each hit event
		virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);    ///< Virtual Method called for each step of each hit
		virtual void EndOfEvent(G4HCofThisEvent*);                   ///< Virtual Method called at the end of each hit event
		
		G4String HCname;                                             ///< Sensitive Detector/Hit Collection Name
		map<string, detector>           *hallMap;                    ///< detector map
		map<string, HitProcess_Factory> *hitProcessMap;              ///< Hit Process Routine Factory Map
		set<vector<identifier> >        Id_Set;                      ///< Identifier Set. Used to determine if a step is inside a new/existing element.
		
		goptions    gemcOpt;   ///< gemc option class
		sensitiveID SDID;      ///< sensitiveID used for identification, hit properties and digitization
		
	private:
		MHitCollection *hitCollection;                               ///< G4THitsCollection<MHit>
		HitProcess     *ProcessHitRoutine;                           ///< To call PID
		int             HCID;                                        ///< HCID increases every new hit collection.
		
		string hd_msg1;        ///< New Hit message
		string hd_msg2;        ///< Normal Message
		string hd_msg3;        ///< End of hit Collection message
		string catch_v;        ///< Volume Name for Verbosity
		double verbosity;      ///< Hit Verbosity
		double RECORD_PASSBY;  ///< If set to one, records particles even if they do not leave any energy
		double RECORD_MIRROR;  ///< If set to one, records particles in the mirror detectors
		
		
	public:
		vector<identifier> GetDetectorIdentifier(string name) {return (*hallMap)[name].identity;} ///< returns detector identity
		string GetDetectorHitType(string name)                {return (*hallMap)[name].hitType;}  ///< returns detector hitType
		MHitCollection* GetMHitCollection()                   {if(hitCollection) return hitCollection; else return NULL;}              ///< returns hit collection
		MHit* find_existing_hit(vector<identifier>);                                               ///< returns hit collection hit inside identifer
};


#endif










