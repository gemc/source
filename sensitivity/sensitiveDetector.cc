// G4 headers
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"

// gemc headers
#include "identifier.h"
#include "sensitiveDetector.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

sensitiveDetector::sensitiveDetector(G4String name, goptions opt, string factory, int run, string variation, string system):G4VSensitiveDetector(name), gemcOpt(opt), HCID(-1)
{
	HCname = name;
	collectionName.insert(HCname);
	hitCollection = NULL;
	
	hd_msg1 = gemcOpt.optMap["LOG_MSG"].args + " New Hit: <<< ";
	hd_msg2 = gemcOpt.optMap["LOG_MSG"].args + " > ";
	hd_msg3 = gemcOpt.optMap["LOG_MSG"].args + " End of Hit Collections >> ";
	catch_v = gemcOpt.optMap["CATCH"].args;
	verbosity = gemcOpt.optMap["HIT_VERBOSITY"].arg;
	RECORD_PASSBY = gemcOpt.optMap["RECORD_PASSBY"].arg;
	RECORD_MIRROR = gemcOpt.optMap["RECORD_MIRRORS"].arg;
	ELECTRONICNOISE  = replaceCharInStringWithChars(gemcOpt.optMap["ELECTRONICNOISE"].args, ",", "  ");
	fastMCMode       = gemcOpt.optMap["FASTMCMODE"].arg;  // fast mc = 2 will increase prodThreshold and maxStep to 5m

	// when background is being saved, all tracks passing by detectors
	// are saved even if they do not deposit energy
	// for FASTMC mode = 2 this is necessary so it does record the hit
	if(gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg == 3 || fastMCMode == 2)
		RECORD_PASSBY = 1;
		
	SDID = sensitiveID(HCname, gemcOpt, factory, variation, system);

}

sensitiveDetector::~sensitiveDetector(){}


void sensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
	Id_Set.clear();
	hitCollection = new MHitCollection(HCname, collectionName[0]);
	if(HCID < 0)  HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	HCE->AddHitsCollection( HCID, hitCollection );
	ProcessHitRoutine = NULL;
	if(verbosity > 1)
		cout << "   > " << collectionName[0] << " initialized." << endl;


}


G4bool sensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	// First check on energy deposited
	double depe = aStep->GetTotalEnergyDeposit();
	// don't enter if RECORD_PASSBY is not set and it's not an optical photon
	// Notice: a gamma will not directly release energy on a scintillator
	// but will convert, and the pair will release energy
	// so by default right now gammas are not recorded and the hit belongs to the pair
	
	G4VTouchable* TH =  (G4VTouchable*) aStep->GetPreStepPoint()->GetTouchable();

	if(depe == 0 && RECORD_PASSBY == 0 && aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
		return false;

	// do not record Mirrors unless specified
	if(aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() && RECORD_MIRROR == 0) return false;
	
	G4Track *trk = aStep->GetTrack();
	// this is expensive should we really check?
	if(trk->GetDefinition()->GetParticleName().find("unknown") != string::npos) return false;
	
	G4StepPoint   *prestep     = aStep->GetPreStepPoint();
	G4StepPoint   *poststep    = aStep->GetPostStepPoint();
	string         processName = "na";
	

	///< Hit informations
	///< The hit position is taken from PostStepPoint (inside the sensitive volume)
	///< Transformation to local coordinates has to be done with prestep
	double         Dx      = aStep->GetStepLength();
	string         name    = TH->GetVolume(0)->GetName();                                     ///< Volume name
	G4ThreeVector   xyz    = poststep->GetPosition();                                         ///< Global Coordinates of interaction
	G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                      ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);
	G4ThreeVector  vert    = trk->GetVertexPosition();
	double         ctime   = poststep->GetGlobalTime();                                       ///< Time of step
	G4ThreeVector  pxyz    = prestep->GetMomentum();                                          ///< Track Momentum (before entering the volume)
	double         ene     = prestep->GetTotalEnergy();                                       ///< Track Energy (before entering the volume)
	int            tid     = trk->GetTrackID();                                               ///< Track ID
	int            pid     = trk->GetDefinition()->GetPDGEncoding();                          ///< Track PID
	int            q       = (int) trk->GetDefinition()->GetPDGCharge();                      ///< Track Charge
	if(trk->GetCreatorProcess()) {
		processName     = trk->GetCreatorProcess()->GetProcessName();                      ///< Process that originated the track
	}
	string materialName    = poststep->GetMaterial()->GetName();                              ///< Material name in this step
	vector<identifier> VID = SetId(GetDetectorIdentifier(name), TH, ctime, SDID.timeWindow, tid);                              ///< Identifier at the geant4 level, using the G4 hierarchy to set the copies
	
	// Get the ProcessHitRoutine to calculate the new vector<identifier>
	if(ProcessHitRoutine == NULL)
	{
		ProcessHitRoutine = getHitProcess(hitProcessMap, GetDetectorHitType(name));
	}
	
	// if not existing, exit
	// this should never happen though
	if(ProcessHitRoutine == NULL)
	{
		cout << endl << "  !!! Error: >" <<  GetDetectorHitType(name) << "< NOT FOUND IN  ProcessHit Map for volume: " << name << " - exiting." << endl;
		return false;
	}

    // getting magnetic field
    const double point[4] = {xyz.x(), xyz.y(), xyz.z(), 10};
    double fieldValue[3] = {0, 0, 0};
    double hitFieldValue = 0;

    G4FieldManager *fmanager = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetFieldManager();

    // if no field manager, the field is zero
    if(fmanager) {
        fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);

        hitFieldValue = sqrt(fieldValue[0]*fieldValue[0] + fieldValue[1]*fieldValue[1] + fieldValue[2]*fieldValue[2]);

    }

	///< Process VID: getting Identifier at the ProcessHitRoutine level
	///< A process routine can generate hit sharing
	vector<identifier> PID = ProcessHitRoutine->processID(VID, aStep, (*hallMap)[name]);
	int singl_hit_size = VID.size();
	int multi_hit_size = PID.size()/singl_hit_size;
	
	// splitting PIDs into an array
	for(int mh = 0; mh<multi_hit_size; mh++)
	{
		vector<identifier> mhPID;
		
		for(int this_id = 0; this_id<singl_hit_size; this_id++)
		{
			identifier this_shit; // adding this single hit
			identifier thisPID = PID[this_id + mh*singl_hit_size];
			this_shit.name       = thisPID.name;
			this_shit.rule       = thisPID.rule;
			this_shit.id         = thisPID.id;
			this_shit.time       = thisPID.time;
			this_shit.TimeWindow = thisPID.TimeWindow;
			this_shit.TrackId    = thisPID.TrackId;
			this_shit.id_sharing = thisPID.id_sharing;
			mhPID.push_back(this_shit);
		}
		
		if(verbosity > 9 || name.find(catch_v) != string::npos)
			cout << endl << hd_msg2 << " Before hit Process Identification:"  << endl << VID
			     << hd_msg2 << " After:  hit Process Identification:" << endl << mhPID << endl;
		
		///< Checking if it's new hit or existing hit. Use the overloaded "=="
		if(verbosity > 9) cout << endl << endl << " BEGIN SEARCH for same hit in Identifier Set..." << endl;
		
		set<vector<identifier> > :: iterator itid;
		int hit_found = 0;
		
		for(itid = Id_Set.begin(); itid!= Id_Set.end() && !hit_found; itid++)
		{
			if(*itid == mhPID)  hit_found=1;
			if(verbosity > 9 )
				cout << "   >> Current Step:  " << mhPID
				<< "   >> Set Hit Index: " << *itid
				<< (hit_found ? "   >> FOUND at this Set Entry. " : "   >> Not found yet. ") << endl;
		}
		if(verbosity > 10) cout << " SEARCH ENDED." << (hit_found ? " 1 " : " No ") << "hit found in the Set." << endl << endl;
		
		
		
		// New Hit
		if(!hit_found)
		{
			MHit *thisHit = new MHit();
			thisHit->SetPos(xyz);
			thisHit->SetLPos(Lxyz);
			thisHit->SetVert(vert);
			thisHit->SetTime(ctime);
			thisHit->SetEdep(depe*mhPID[singl_hit_size-1].id_sharing);
			thisHit->SetDx(Dx);
			thisHit->SetMom(pxyz);
			thisHit->SetE(ene);
			thisHit->SetTrackId(tid);
			thisHit->SetDetector((*hallMap)[name]);
			thisHit->SetId(mhPID);
			thisHit->SetPID(pid);
			thisHit->SetCharge(q);
			thisHit->SetMatName(materialName);
			thisHit->SetProcID(processID(processName));
            thisHit->SetSDID(SDID);
            thisHit->SetMgnf(hitFieldValue);
			hitCollection->insert(thisHit);
			Id_Set.insert(mhPID);
			
			if(verbosity > 6 || name.find(catch_v) != string::npos)
			{
				string pid    = aStep->GetTrack()->GetDefinition()->GetParticleName();
				cout << endl << hd_msg1 << endl
				<< "  > This element was not hit yet in this event. Identity:" << endl << thisHit->GetId()
				<< "  > Creating new hit by a E=" <<  ene/MeV << ", p=" <<  pxyz.mag()/MeV << " MeV "  << pid
				<< ", track ID = " << tid << ", inside Hit Collection <" << HCname << ">." << endl
				<< "  > Energy Deposited this step: " << depe/MeV << " MeV" << endl
				<< "  > Time of this step: " << G4BestUnit(ctime, "Time") << endl
				<< "  > Position of this step:   " << xyz/cm  << " cm"  << endl
				<< "  > Local Position in volume: " << Lxyz/cm  << " cm" << endl;
			}
		}
		else
		{
			// Adding hit info only if the poststeppint remains in the volume?
			// if( aStep->GetPreStepPoint()->GetTouchable()->GetVolume(0) == aStep->GetPostStepPoint()->GetTouchable()->GetVolume(0))
			{
				MHit *thisHit = find_existing_hit(mhPID);
				if(!thisHit)
				{
					cout << " Hit not found in collection but found in PID. This should never happen. Exiting." << endl;
					exit(0);
				}
				else
				{
					thisHit->SetPos(xyz);
					thisHit->SetLPos(Lxyz);
					thisHit->SetVert(vert);
					thisHit->SetTime(ctime);
					thisHit->SetEdep(depe*mhPID[singl_hit_size-1].id_sharing);
					thisHit->SetDx(Dx);
					thisHit->SetMom(pxyz);
					thisHit->SetE(ene);
					thisHit->SetTrackId(tid);
					thisHit->SetPID(pid);
					thisHit->SetCharge(q);
					thisHit->SetMatName(materialName);
					thisHit->SetProcID(processID(processName));
					thisHit->SetDetector((*hallMap)[name]);
                    thisHit->SetMgnf(hitFieldValue);

					if(verbosity > 6 || name.find(catch_v) != string::npos)
					{
						string pid    = aStep->GetTrack()->GetDefinition()->GetParticleName();
						cout << hd_msg2 << " Step Number " << thisHit->GetPos().size()
						<< " inside Identity: "  << endl << thisHit->GetId()
						<< "            >  Adding hit inside Hit Collection <" << HCname << ">."
						<< " by a E=" <<  ene/MeV << ", p=" <<  pxyz.mag()/MeV << " MeV "
						<< pid << ", track ID = " << tid << endl
						<< "            >  Energy Deposited this step: " << depe/MeV << " MeV" << endl
						<< "            >  Time of this step: " << G4BestUnit(ctime, "Time")
						<< " is within this element Time Window of " << SDID.timeWindow/ns << " ns. " << endl
						<< "            >  Position of this step:   " << xyz/cm << " cm" << endl
						<< "            >  Local Position in volume: " << Lxyz/cm  << " cm" << endl;
					}
				}
			}
		}
	}
	
	
	return true;
}


void sensitiveDetector::EndOfEvent(G4HCofThisEvent *HCE)
{
	int nhitC = hitCollection->GetSize();
	if(!nhitC) return;

	// adding electronic noise to hits
	// only if requested by user
	// notice: there should be routine to decide if hit is in the same TW
	// (in that case it is not a new hit)
	if(ELECTRONICNOISE.find(HCname) != string::npos)
	{
		vector<MHit*> noiseHits = ProcessHitRoutine->electronicNoise();
		for(unsigned int h=0; h<noiseHits.size(); h++)
			hitCollection->insert(noiseHits[h]);
	}

//	cout << HCname << " Before: HITC " << hitCollection << " size: " << nhitC << endl;

//	// adding background noise to hits
//	for(auto bgh: currentBackground) {
//		if(bgh != nullptr) {
//			hitCollection->insert(new MHit(bgh->energy, bgh->timeFromEventStart, bgh->nphe, bgh->identity));
//		}
//	}
//
//	nhitC = hitCollection->GetSize();
//	cout << HCname << " After: HITC " << hitCollection << " size: " << nhitC << endl;



	MHit *aHit;
	double Etot;
	if(verbosity > 2 && nhitC)
	{
		cout << endl;
		cout << hd_msg3 << " Hit Collections <" << HCname << ">: " << nhitC << " hits." << endl;
		
		for (int i=0; i<nhitC; i++)
		{
			aHit = (*hitCollection)[i];
			string vname = aHit->GetId()[aHit->GetId().size()-1].name;
			if(verbosity > 5 || vname.find(catch_v) != string::npos)
			{
				cout << hd_msg3 << " Hit " << i + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
				cout << aHit->GetId();
				Etot = 0;
				for(unsigned int e=0; e<aHit->GetPos().size(); e++) Etot = Etot + aHit->GetEdep()[e];
				cout << "   Total energy deposited: " << Etot/MeV << " MeV" << endl;
			}
		}
	}




	if(ProcessHitRoutine) delete ProcessHitRoutine; // not needed anymore
}


MHit*  sensitiveDetector::find_existing_hit(vector<identifier> PID)  ///< returns hit collection hit inside identifer
{
	for(unsigned int i=0; i<hitCollection->GetSize(); i++)
		if((*hitCollection)[i]->GetId() == PID) return (*hitCollection)[i];
	
	return NULL;
}

// to check process name go to $G4ROOT/$GEANT4_VERSION/source/geant$GEANT4_VERSION/source/processes/
// mgrep "const G4String&" | grep process
// notice: this map should be written out in the output stream
int sensitiveDetector::processID(string procName)
{
	if(procName == "eIoni")                 return 1;
	if(procName == "compt")                 return 2;
	if(procName == "eBrem")                 return 3;
	if(procName == "phot")                  return 4;
	if(procName == "conv")                  return 5;
	if(procName == "annihil")               return 6;
	if(procName == "photonNuclear")         return 7;
	if(procName == "electronNuclear")       return 8;
	if(procName == "positronNuclear")       return 9;
	if(procName == "CoulombScat")           return 10;
	if(procName == "Cerenkov")              return 11;
	if(procName == "Scintillation")			return 12;
	if(procName == "SynRad")           		return 13;
	if(procName == "SynchrotronRadiation")	return 14;
	if(procName == "hadElastic")            return 20;
	if(procName == "hBrems")                return 21;
	if(procName == "hIoni")                 return 22;
	if(procName == "hPairProd")             return 23;
	if(procName == "protonInelastic")       return 30;
	if(procName == "neutronInelastic")      return 31;
	if(procName == "nCapture")              return 32;
	if(procName == "anti_neutronInelastic") return 33;
	if(procName == "pi-Inelastic")          return 40;
	if(procName == "pi+Inelastic")          return 41;
	if(procName == "Decay")                 return 50;
	if(procName == "DecayWithSpin")         return 51;
	if(procName == "muIoni")                return 60;
	if(procName == "muPairProd")            return 61;
	if(procName == "muBrems")               return 62;
	if(procName == "muonNuclear")           return 63;
	if(procName == "kaon-Inelastic")        return 70;
	if(procName == "kaon+Inelastic")        return 71;
	if(procName == "kaon0Inelastic")        return 72;
	if(procName == "kaon0LInelastic")       return 73;
	if(procName == "kaon0SInelastic")       return 74;
	if(procName == "alphaInelastic")        return 80;
	if(procName == "lambdaInelastic")       return 90;
	if(procName == "sigma-Inelastic")       return 100;
	if(procName == "dInelastic")            return 110;
	if(procName == "ionIoni")               return 120;
	if(procName == "ionInelastic")          return 121;
	if(procName == "tInelastic")            return 130;
	if(procName == "xi0Inelastic")          return 140;
	if(procName == "omega-Inelastic")       return 150;

	if(procName == "na")                    return 999;

	if(verbosity > 0)
		cout << " Warning: process name " << procName << " not catalogued." << endl;

	return 999;
}




