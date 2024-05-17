// G4 headers
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"
#include "G4UImanager.hh"

// gemc headers
#include "MEventAction.h"
#include "Hit.h"

// mlibrary
#include "frequencySyncSignal.h"

// c++
#include <iostream>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// return original track id of a vector of tid
vector<int> MEventAction::vector_otids(vector<int> tids)
{
	vector<int> otids;
	for(unsigned int t=0; t<tids.size(); t++)
	{
		if(hierarchy.find(tids[t]) != hierarchy.end())
			otids.push_back(hierarchy[tids[t]]);
		else
			otids.push_back(0);
	}
	return otids;
}


// the following functions return the pid, mpid and charges vector<int> from a map<int, TInfos>
vector<int> vector_zint(int size)
{
	vector<int> zints;
	for(int t=0; t<size; t++)
		zints.push_back(0);
	
	return zints;
}

vector<G4ThreeVector> vector_zthre(int size)
{
	vector<G4ThreeVector> zthre;
	for(int t=0; t<size; t++)
		zthre.push_back(G4ThreeVector(0,0,0));
	
	return zthre;
}

vector<int> vector_mtids(map<int, TInfos> tinfos, vector<int> tids)
{
	vector<int> mtids;
	for(unsigned int t=0; t<tids.size(); t++)
		mtids.push_back(tinfos[tids[t]].mtid);
	
	return mtids;
}

vector<int> vector_mpids(map<int, TInfos> tinfos, vector<int> tids)
{
	vector<int> mpids;
	for(unsigned int t=0; t<tids.size(); t++)
		mpids.push_back(tinfos[tids[t]].mpid);
	
	return mpids;
}

vector<G4ThreeVector> vector_mvert(map<int, TInfos> tinfos, vector<int> tids)
{
	vector<G4ThreeVector> mvert;
	for(unsigned int t=0; t<tids.size(); t++)
		mvert.push_back(tinfos[tids[t]].mv);
	
	return mvert;
}

MEventAction::MEventAction(goptions opts, map<string, double> gpars)
{
	gemcOpt          = opts;
	hd_msg           = gemcOpt.optMap["LOG_MSG"].args + " Event Action: >> ";
	Modulo           = (int) gemcOpt.optMap["PRINT_EVENT"].arg ;
	VERB             = gemcOpt.optMap["EVENT_VERBOSITY"].arg ;
	catch_v          = gemcOpt.optMap["CATCH"].args;
	SAVE_ALL_MOTHERS = (int) gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg ;
	SAVE_ALL_ANCESTORS = (int) gemcOpt.optMap["SAVE_ALL_ANCESTORS"].arg ;
	gPars            = gpars;
	MAXP             = (int) gemcOpt.optMap["NGENP"].arg;
	FILTER_HITS      = (int) gemcOpt.optMap["FILTER_HITS"].arg;
	FILTER_HADRONS   = (int) gemcOpt.optMap["FILTER_HADRONS"].arg;
	FILTER_HIGHMOM   = (int) gemcOpt.optMap["FILTER_HIGHMOM"].arg;
	SKIPREJECTEDHITS = (int) gemcOpt.optMap["SKIPREJECTEDHITS"].arg;
	rw               = runWeights(opts);
	
	WRITE_ALLRAW     = replaceCharInStringWithChars(gemcOpt.optMap["ALLRAWS"].args, ",", "  ");
	WRITE_INTRAW     = replaceCharInStringWithChars(gemcOpt.optMap["INTEGRATEDRAW"].args, ",", "  ");
	WRITE_INTDGT     = replaceCharInStringWithChars(gemcOpt.optMap["INTEGRATEDDGT"].args, ",", "  ");
	SIGNALVT         = replaceCharInStringWithChars(gemcOpt.optMap["SIGNALVT"].args, ",", "  ");
	RFSETUP          = replaceCharInStringWithChars(gemcOpt.optMap["RFSETUP"].args, ",", "  ");
	RFSTART          = replaceCharInStringWithChars(gemcOpt.optMap["RFSTART"].args, ",", "  ");
	fastMCMode       = gemcOpt.optMap["FASTMCMODE"].arg;  // fast mc = 2 will increase prodThreshold and maxStep to 5m
	
	requestedNevents = (long int) gemcOpt.optMap["N"].arg ;
	ntoskip        = gemcOpt.optMap["SKIPNGEN"].arg;
	last_runno     = -99;

	
	// fastMC mode will set SAVE_ALL_MOTHERS to 1
	// a bit cluncky for now
	if(fastMCMode>0) {
		SAVE_ALL_MOTHERS = 1;
	}
	
	// SAVE_ALL_ANCESTORS will set SAVE_ALL_MOTHERS to nonzero
	if (SAVE_ALL_ANCESTORS && (SAVE_ALL_MOTHERS == 0)) {
		SAVE_ALL_MOTHERS = 1;
	}
	tsampling  = get_number(get_info(gemcOpt.optMap["TSAMPLING"].args).front());
	nsamplings = get_number(get_info(gemcOpt.optMap["TSAMPLING"].args).back());
	
	if(SAVE_ALL_MOTHERS>1) {
		lundOutput = new ofstream("background.dat");
		cout << " > Opening background.dat file to save background particles in LUND format." << endl;
	}
	
	evtN = gemcOpt.optMap["EVTN"].arg;
	
	// background hits
	backgroundHits = nullptr;
	BGFILE = gemcOpt.optMap["MERGE_BGHITS"].args;
	
	// there's no check that the map is built correctly
	if(BGFILE != "no") {
		backgroundHits = new GBackgroundHits(BGFILE, requestedNevents, VERB);
	}
	
	backgroundEventNumber.clear();
	
	// SAVE_SELECTED parameters
	string arg = gemcOpt.optMap["SAVE_SELECTED"].args;
	if (arg == "" || arg == "no") {
		ssp.enabled = false;
	} else {
		vector<string> values;
		string units;
		values       = get_info(gemcOpt.optMap["SAVE_SELECTED"].args, string(",\""));
		if (values.size() == 5 || values.size() == 6) {
			ssp.enabled = true;
			ssp.targetId  = values[0];
			ssp.tIdsize   = ssp.targetId.size();
			ssp.targetPid = get_number(values[1]);
			ssp.lowLim    = get_number(values[2]);
			ssp.hiLim     = get_number(values[3]);
			ssp.variable  = values[4];
			if (values.size() == 5)
				ssp.dir = "./";
			else
			{
				ssp.dir = values[5];
				if (ssp.dir[ssp.dir.size()-1] != '/' ) ssp.dir += "/";
			}
#ifdef WIN32
			std::replace(ssp.dir.begin(), ssp.dir.end(),'/','\\');
#endif
			G4RunManager::GetRunManager()->SetRandomNumberStore(true);
			G4RunManager::GetRunManager()->SetRandomNumberStoreDir(ssp.dir);
		}
	}
	
	
	if(RFSETUP == "clas12_ccdb") {
		setup_clas12_RF(rw.getRunNumber(evtN));

	} else if(RFSETUP != "no")  {
		rfvalue_strings = getStringVectorFromString(RFSETUP);
		set_and_show_rf_setup();
	}
	
}

MEventAction::~MEventAction()
{
	if(SAVE_ALL_MOTHERS>1)
		lundOutput->close();
}

void MEventAction::BeginOfEventAction(const G4Event* evt)
{
	G4RunManager *runManager = G4RunManager::GetRunManager();;
	if(gen_action->isFileOpen() == false) {
		runManager->AbortRun();
		cout << " No more events in the input file." << endl;
		return;
	}
	
	MPrimaryGeneratorAction* pga = (MPrimaryGeneratorAction*)(runManager->GetUserPrimaryGeneratorAction());
	if (pga->doneRerun())
		return;
	if (pga->isRerun())
		evtN = pga->rerunEvent();
	
	rw.getRunNumber(evtN);
	bgMap.clear();
		
	static int lastEvtN = -1;
	if(evtN > lastEvtN && evtN%Modulo == 0 ) {
		cout << hd_msg << " Begin of event " << evtN << "  Run Number: " << rw.runNo;
		if(rw.isNewRun) {
			cout << " (new) ";
			G4Random::getTheGenerator ()->showEngineStatus() ;
		}
		if ( VERB > 1 ) {
			G4Random::getTheGenerator ()->showEngineStatus() ;
		}
		cout << endl;
		lastEvtN = evtN;
	}
	
	
	// background hits:
	// checking the whole hit map
	if(VERB > 4) {
		if(backgroundHits != nullptr) {
			for(auto sDet: SeDe_Map) {
				map<int, vector<BackgroundHit*> > *backgroundHitsEventMap = backgroundHits->getBackgroundForSystem(sDet.first);
				if(backgroundHitsEventMap != nullptr) {
					for(auto bgHits: (*backgroundHitsEventMap)) {
						if(bgHits.first == evtN) {
							cout << " >>> Background hits for detector " << sDet.first << ", event number: " << bgHits.first <<  endl;
							for(auto bgh: bgHits.second) {
								cout << *bgh << endl;
							}
							
						}
					}
				}
			}
		}
	}
	ssp.decision = false; // default is don't save RNG
}

void MEventAction::EndOfEventAction(const G4Event* evt)
{

	if ((gen_action->isFileOpen() == false) ||
		 (gen_action->doneRerun() == true))
		return;

	// completely skip the event
	// (still need to increase event number)
	if(evtN <= ntoskip) {
		evtN++;
		return;
	}
	
	MHitCollection* MHC;
	int nhits;
	
	
	// if FILTER_HITS is set, checking if there are any hits
	if(FILTER_HITS) {
		int anyHit = 0;
		for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++) {
			MHC = it->second->GetMHitCollection();
			if (MHC) anyHit += MHC->GetSize();
		}
		
		// stop here if there are no hits and FILTER_HITS is set
		if(anyHit==0) return;
	}
	
	
	// if FILTER_HADRONS is set, checking if there are any (matching) hadrons
	if (FILTER_HADRONS == 1 || abs(FILTER_HADRONS) > 99) {
		int foundHad = 0;
		for (map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++) {
			MHC = it->second->GetMHitCollection();
			if (MHC) {
				nhits = MHC->GetSize();
			} else {
				nhits = 0;
			}
			for (int h=0; h<nhits; h++) {
				vector<int>           pids = (*MHC)[h]->GetPIDs();
				for (vector<int>::const_iterator pit = pids.begin(); pit != pids.end(); pit++) {
					if ((FILTER_HADRONS == 1 && abs(*pit) > 99) || *pit == FILTER_HADRONS) {
						foundHad = 1;
						break;
					}
				}
			}
		}
		
		// stop here if there are no hadrons and FILTER_HADRONS is set
		if (not foundHad) return;
	}
	
	// if FILTER_HIGHMOM is set, checking if there is  any high mom hits
	if (FILTER_HIGHMOM != 0) {
		int foundHighmom = 0;
		for (map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++) {
			MHC = it->second->GetMHitCollection();
			if (MHC) nhits = MHC->GetSize();
			else nhits = 0;
			for (int h=0; h<nhits; h++) {
				vector<G4ThreeVector> mmts = (*MHC)[h]->GetMoms();
				for(unsigned int t=0; t<mmts.size(); t++) {
					// 		      	cout << "mom " << mmts[t].mag() << endl;
					if (mmts[t].mag()>FILTER_HIGHMOM) {
						// 		      	cout << "mom " << mmts[t].mag() << endl;
						foundHighmom = 1;
						break;
					}
				}

			}
		}

		// stop here if there are no  any high mom hits
		if (not foundHighmom) return;
	}

	if(evtN%Modulo == 0 && VERB > 0) {
		cout << hd_msg << " EndOfEventAction for event number " << evtN << ",  Run Number: " << rw.runNo << endl;
	}
	
	// building the tracks set database with all the tracks in all the hits
	// if SAVE_ALL_MOTHERS is set
	set<int> track_db;
	if(SAVE_ALL_MOTHERS) {
		for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++) {
			MHC = it->second->GetMHitCollection();
			if (MHC) nhits = MHC->GetSize();
			else nhits = 0;
			
			for(int h=0; h<nhits; h++) {
				vector<int> tids = (*MHC)[h]->GetTIds();
				
				for(unsigned int t=0; t<tids.size(); t++) {
					if((*MHC)[h]->isBackgroundHit == 0)
						track_db.insert(tids[t]);
				}
				
				if(SAVE_ALL_MOTHERS>1) {
					vector<int>           pids = (*MHC)[h]->GetPIDs();
					// getting the position of the hit, not vertex of track
					vector<G4ThreeVector> vtxs = (*MHC)[h]->GetPos();
					vector<G4ThreeVector> mmts = (*MHC)[h]->GetMoms();
					vector<double>        tims = (*MHC)[h]->GetTime();
					// only put the first step of a particular track
					// (don't fill if track exist already)
					for(unsigned int t=0; t<tids.size(); t++)
						if(bgMap.find(tids[t]) == bgMap.end())
							bgMap[tids[t]] = BGParts(pids[t], tims[t], vtxs[t], mmts[t]);
				}
			}
		}
	}
	
	// now filling the map of tinfos with tracks infos from the track_db database
	// this won't get the mother particle infos except for their track ID
	map<int, TInfos> tinfos;
	set<int>::iterator it;
	
	// the container is full only if /tracking/storeTrajectory 2
	G4TrajectoryContainer *trajectoryContainer = nullptr;
	
	set<int> track_db2 = track_db;
	if(SAVE_ALL_MOTHERS) {
		trajectoryContainer = evt->GetTrajectoryContainer();
		momDaughter.clear();
		
		if(VERB>3) {
			cout << " >> Total number of tracks " << trajectoryContainer->size() << endl;
		}
		
		while(trajectoryContainer && track_db.size()) {
			// looping over all tracks
			for(unsigned int i=0; i< trajectoryContainer->size(); i++)
			{
				// cout << " index track " << i << endl;
				G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
				int tid = trj->GetTrackID();
				
				// adding track in mom daughter relationship
				// for all tracks
				if(momDaughter.find(tid) == momDaughter.end())
					momDaughter[tid] = trj->GetParentID();
				
				// is this track involved in a hit?
				// if yes, add it to the track map db
				it = track_db.find(tid);
				if(it != track_db.end())
				{
					int           mtid = trj->GetParentID();
					tinfos[tid]        = TInfos(mtid);
					// cout << " At map level: " << tid << " " <<mtid << " " << pid << "   charge: " << q << endl;
					// remove the track from the database so we don't have to loop over it next time
					track_db.erase(it);
				}
			}
		}
		
		// building the hierarchy map
		// for all secondary tracks
		for(map<int, int>::iterator itm = momDaughter.begin(); itm != momDaughter.end(); itm++) {
			int ancestor = itm->first;
			if(momDaughter[ancestor] == 0)
				hierarchy[itm->first] = itm->first;
			
			while (momDaughter[ancestor] != 0) {
				hierarchy[itm->first] = momDaughter[ancestor];
				ancestor = momDaughter[ancestor];
			}
		}
		
		// now accessing the mother particle infos
		for(map<int, TInfos>::iterator itm = tinfos.begin(); itm != tinfos.end(); itm++) {
			int mtid = (*itm).second.mtid;
			// looking for mtid infos in the trajectoryContainer
			for(unsigned int i=0; i< trajectoryContainer->size() && mtid != 0; i++) {
				G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
				int tid = trj->GetTrackID();
				if(tid == mtid) {
					tinfos[(*itm).first].mpid   = trj->GetPDGEncoding();
					tinfos[(*itm).first].mv     = trj->GetPoint(0)->GetPosition();
				}
			}
		}
		
		// removing daughter tracks if mother also gave hits
		if(SAVE_ALL_MOTHERS>2) {
			vector<int> bgtIDs;
			for(map<int, BGParts>::iterator it = bgMap.begin(); it != bgMap.end(); it++)
				bgtIDs.push_back(it->first);
			
			for(unsigned i=0; i<bgtIDs.size(); i++) {
				int daughter = bgtIDs[i];
				int momCheck = 0;
				if(momDaughter.find(daughter) != momDaughter.end()) {
					momCheck = momDaughter[daughter];
				}

				while(momCheck != 0) {
					// mom has a hit
					if(bgMap.find(momCheck) != bgMap.end()) {
						// so if daughter is found, delete it
						// daughter may already be deleted before
						if(bgMap.find(daughter) != bgMap.end())
							bgMap.erase(bgMap.find(daughter));
					}
					// going up one generation
					if(momDaughter.find(momCheck) != momDaughter.end())
						momCheck = momDaughter[momCheck];
					else
						momCheck = 0;
				}
			}
		}
	}
	
	
	// Making sure the output routine exists in the ProcessOutput map
	// If no output selected, or HitProcess not found, end current event
	map<string, outputFactoryInMap>::iterator ito = outputFactoryMap->find(outContainer->outType);
	if(ito == outputFactoryMap->end()) {
		if(outContainer->outType != "no") {
			cout << hd_msg << " Warning: output type <" << outContainer->outType
			<< "> is not registered in the outContainerput Factory. " << endl
			<< "      This event will not be written out." << endl;
		}
		evtN++;
		return;
	}
	outputFactory *processOutputFactory = getOutputFactory(outputFactoryMap, outContainer->outType);

	// configuration contains:
	// number of hits in the hit collection for each sensitive detector
	map<string, double> configuration;
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++) {
		MHC = it->second->GetMHitCollection();

		if (MHC) {
			if (MHC->GetSize() > 0) {
				string hitType = it->first;

				if(WRITE_INTRAW.find(hitType) != string::npos || WRITE_INTRAW == "*") {
					configuration[hitType] = MHC->GetSize();
				}
			}
		}
	}


	processOutputFactory->prepareEvent(outContainer, &configuration);
	
	// Header Bank contains event number
	// Need to change this to read DB header bank
	map<string, double> header;
	header["runNo"]    = rw.runNo;
	header["evn"]      = evtN;
	header["evn_type"] = -1;  // physics event. Negative is MonteCarlo event
	header["beamPol"]  = gen_action->getBeamPol();
	
	// write event header bank
	processOutputFactory->writeHeader(outContainer, header, getBankFromMap("header", banksMap));
	
	// user header should be in a different tag than the normal header
	// for now, we're ok
	// assuming 100 user vars max
	map<string, double> userHeader;
	for(unsigned i=0; i<gen_action->headerUserDefined.size(); i++) {
		string tmp = "userVar" ;
		if(i<9)        tmp +="00";
		else if (i<99) tmp +="0";
		
		tmp += to_string(i+1);
		
		userHeader[tmp] = gen_action->headerUserDefined[i];
	}
	
	// write event header bank
	processOutputFactory->writeUserInfoseHeader(outContainer, userHeader);
	
	// write RF bank if present
	// do not write in FASTMC mode
	if(RFSETUP!= "no" && fastMCMode == 0) {
		
		double additionalTime = 0;
		
		vector<string> rfstartoption = getStringVectorFromString(RFSTART);
		// 4 options means we're in eventVertex case
		if(rfstartoption.size() == 4) {
			// verifying that we're in the eventVertex case
			if(rfstartoption[0] == "eventVertex") {
				G4ThreeVector referenceRFPosition(get_number(rfstartoption[1]), get_number(rfstartoption[2]), get_number(rfstartoption[3]));
				G4ThreeVector firstParticleVertex = evt->GetPrimaryVertex(0)->GetPosition();
				additionalTime = (firstParticleVertex - referenceRFPosition).mag()/mm / 299.792; // speed of light in mm / ns
				if(firstParticleVertex.z() > referenceRFPosition.z()) additionalTime = -additionalTime;
			}
		}

        CLHEP::HepRandomEngine* currentEngine = CLHEP::HepRandom::getTheEngine();

        // this should be an array of longs describing the engine status
        // the first two elements are
        const long *engineStatus = currentEngine->getSeeds();

        int g4rseed = engineStatus[2];
        if (g4rseed == 0) {
            cout << "Error: engineStatus third element is 0" << endl;
            exit(1);
        }

		// getting time window
		string rfsetup_string = to_string(g4rseed) + " " + to_string(gen_action->getTimeWindow()) + " " ;
		
		// getting start time of the event
		rfsetup_string += to_string(gen_action->getStartTime() + additionalTime) + " " ;
		
		if(RFSETUP == "clas12_ccdb"){
			setup_clas12_RF(rw.runNo);
		}

		for(unsigned i=0; i<rfvalue_strings.size(); i++) {
			rfsetup_string += rfvalue_strings[i] + " " ;
		}
		
		FrequencySyncSignal rfs(rfsetup_string);
		processOutputFactory->writeRFSignal(outContainer, rfs, getBankFromMap("rf", banksMap));
		
		if(VERB > 1) {
			cout << rfs << endl;
		}
	}
	
	// Getting Generated Particles info
	// Are these loops necessary, revisit later1
	vector<generatedParticle> MPrimaries;
	for(int pv=0; pv<evt->GetNumberOfPrimaryVertex() && pv<MAXP; pv++) {
		generatedParticle Mparticle;
		G4PrimaryVertex* MPV = evt->GetPrimaryVertex(pv);
		Mparticle.vertex = MPV->GetPosition();
		double thisTime = MPV->GetT0();
		int thisMult    = MPV->GetNumberOfParticle();
		
		for(int pp = 0; pp<MPV->GetNumberOfParticle() && pv<MAXP; pp++) {
			G4PrimaryParticle *PP  = MPV->GetPrimary(pp);
			Mparticle.momentum     = PP->GetMomentum();
			Mparticle.PID          = PP->GetPDGcode();
			Mparticle.time         = thisTime;
			Mparticle.multiplicity = thisMult;
		}
		MPrimaries.push_back(Mparticle)  ;
	}
	
	if(SAVE_ALL_MOTHERS>1) {
		saveBGPartsToLund();
	}

	if(VERB > 4) {
		for(unsigned pi = 0; pi<MPrimaries.size(); pi++) {
			cout << " Particle " << pi + 1 << " has " << MPrimaries[pi].pSum.size() << " particle summaries:" << endl;
			for(unsigned ss =0; ss<MPrimaries[pi].pSum.size(); ss++)
			{
				cout << " \t det: " << MPrimaries[pi].pSum[ss].dname <<
				"  Etot: "  <<  MPrimaries[pi].pSum[ss].etot <<
				"  time: "  <<  MPrimaries[pi].pSum[ss].t << endl;
			}
			cout << endl;
		}
	}

	
	map<int, vector<hitOutput> > hit_outputs_from_AllSD;
	
	// loop over sensitive detectors
	// if there are hits, process them and/or write true infos out
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++) {

		MHC = it->second->GetMHitCollection();
		
		// adding background if existing
		// adding background noise to hits
		vector<BackgroundHit*> currentBackground = getNextBackgroundEvent(it->first);
		for(auto bgh: currentBackground) {
			if(bgh != nullptr) {
				if(MHC)
					MHC->insert(new MHit(bgh->energy, bgh->timeFromEventStart, bgh->nphe, bgh->identity));
			}
		}
		
		if (MHC) nhits = MHC->GetSize();
		else nhits = 0;
		
		// The same ProcessHit Routine must apply to all the hits  in this HitCollection.
		// Instantiating the ProcessHitRoutine only once for the first hit.
		// this conditions applies to digitization and true information processing
		if(nhits) {

			string hitType = it->first;
			
			HitProcess *hitProcessRoutine = getHitProcess(hitProcessMap, hitType);
			if(!hitProcessRoutine)
				return;
			
			if(fastMCMode == 0 || fastMCMode > 9) {
				hitProcessRoutine->init(hitType, gemcOpt, gPars);
			}

			// geant4 integrated digitized information
			// by default they are all ENABLED
			// user can disable them one by one
			// using the INTEGRATEDDGT option
			// for FASTMC mode, do not digitize the info
			vector<hitOutput> allDgtOutput;

			// the digitization routine may decide to skip writing events.
			// keeping the hit number in a vector so we can skip the event writing for the true information as well
			vector<int> hitsToSkip;

			// creating summary information for each generated particle
			for(unsigned pi = 0; pi<MPrimaries.size(); pi++) {
				MPrimaries[pi].pSum.push_back(summaryForParticle("na"));
				if(fastMCMode > 0) {
					MPrimaries[pi].fastMC.push_back(fastMCForParticle("na"));
				}
			}

			
			for(int h=0; h<nhits; h++) {
				MHit* aHit = (*MHC)[h];

				// mother particle infos
				if(SAVE_ALL_MOTHERS) {


					// setting track infos before processing the hit
					vector<int> tids = aHit->GetTIds();
					vector<int> otids = vector_otids(tids);
					aHit->SetmTrackIds(vector_mtids(tinfos, tids));
					aHit->SetoTrackIds(otids);
					aHit->SetmPIDs(    vector_mpids(tinfos, tids));
					aHit->SetmVerts(   vector_mvert(tinfos, tids));

					// for every particle initializing a vector
					// the index is the primary particle index
					// the int will increase for each step
					// if the int > 0
					// then we count it as ONE hit by the track
					vector<int> hitByPrimary;
					for(unsigned pi = 0; pi<MPrimaries.size(); pi++) {
						hitByPrimary.push_back(0);
					}

					// all these vector have the same length.
					for(unsigned pi = 0; pi<MPrimaries.size(); pi++) {
						vector<double> edeps = aHit->GetEdep();
						vector<double> times = aHit->GetTime();
						MPrimaries[pi].pSum.back().nphe = aHit->GetTIds().size();
						if(fastMCMode > 0) {
							MPrimaries[pi].fastMC.back().pOrig  = aHit->GetMom();
							MPrimaries[pi].fastMC.back().pSmear = hitProcessRoutine->psmear(aHit->GetMom());
						}
						for(unsigned ss =0; ss<edeps.size(); ss++) {
							if(otids[ss] == (int) pi+1) {
								MPrimaries[pi].pSum.back().etot += edeps[ss];
								hitByPrimary[pi]++;
								// getting fastest time - should we put threshold here?
								if(MPrimaries[pi].pSum.back().t < 0 || MPrimaries[pi].pSum.back().t > times[ss])
									MPrimaries[pi].pSum.back().t = times[ss];

							}
						}

						if(hitByPrimary[pi]) MPrimaries[pi].pSum.back().stat++;

						if(MPrimaries[pi].pSum.back().etot > 0 || MPrimaries[pi].pSum.back().nphe > 0)
							MPrimaries[pi].pSum.back().dname = hitType;
					}
				} else {
					// filling mother infos with zeros
					int thisHitSize = aHit->GetId().size();
					vector<int>           zint  = vector_zint(thisHitSize);
					vector<G4ThreeVector> zthre = vector_zthre(thisHitSize);
					aHit->SetoTrackIds(zint);
					aHit->SetmTrackIds(zint);
					aHit->SetmPIDs(    zint);
					aHit->SetmVerts(   zthre);
				}
			}

			if(WRITE_INTDGT.find(hitType) == string::npos && (fastMCMode == 0 ||fastMCMode > 9)) {

				hitProcessRoutine->initWithRunNumber(rw.runNo);
				
				for(int h=0; h<nhits; h++) {

					hitOutput thisHitOutput;
					MHit* aHit = (*MHC)[h];
					
					// calling integrateDgt will also set writeHit
					thisHitOutput.setDgtz(hitProcessRoutine->integrateDgt(aHit, h+1));
					
					// include this hit. Users can set writeHit to false to avoid writing the hit
					// the hitProcessRoutine variable detectorThreshold could be used in integrateDgt
					if(hitProcessRoutine->writeHit) {
						allDgtOutput.push_back(thisHitOutput);
					} else {
						if(VERB > 4 ) {
							cout << " Event Action: hit " << h + 1 << " was rejected in " << hitType << " digitization routine." << endl;
						}
						hitsToSkip.push_back(h);
					}
					
					string vname = aHit->GetId()[aHit->GetId().size()-1].name;
					if(VERB > 6 || vname.find(catch_v) != string::npos)
					{
						cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
						cout << aHit->GetId();
						double Etot = 0;
						for(unsigned int e=0; e<aHit->GetPos().size(); e++) Etot = Etot + aHit->GetEdep()[e];
						cout << "   Total energy deposited: " << Etot/MeV << " MeV" << endl;
					}
				}
				processOutputFactory->writeG4DgtIntegrated(outContainer, allDgtOutput, hitType, banksMap);
				
			} // end of geant4 integrated digitized information
			
			// geant4 integrated raw information
			// by default they are all DISABLED
			// user can enable them one by one
			// using the INTEGRATEDRAW option

			bool WRITE_TRUE_INTEGRATED = 0;
			bool WRITE_TRUE_ALL = 0;

			if(WRITE_INTRAW.find(hitType) != string::npos || WRITE_INTRAW == "*") WRITE_TRUE_INTEGRATED = 1;
			if(WRITE_ALLRAW.find(hitType) != string::npos || WRITE_ALLRAW == "*") WRITE_TRUE_ALL = 1;

			vector<hitOutput> allRawOutput;


			for(int h=0; h<nhits; h++) {
				MHit* aHit = (*MHC)[h];

				// electronic noise hits disable? Why? TODO
				if(aHit->isElectronicNoise) {
					continue;
				}

				hitOutput thisHitOutput;


				if(fastMCMode == 0 || fastMCMode > 9) {
					thisHitOutput.setRaws(hitProcessRoutine->integrateRaw(aHit, h+1, WRITE_TRUE_INTEGRATED));
				}

				if(WRITE_TRUE_ALL && (fastMCMode == 0 || fastMCMode > 9)) {
					thisHitOutput.setAllRaws(hitProcessRoutine->allRaws(aHit, h+1));
				}

				if (SKIPREJECTEDHITS == 0) {
					allRawOutput.push_back(thisHitOutput);
				} else {
					if ( find(hitsToSkip.begin(), hitsToSkip.end(), h) == hitsToSkip.end() ) {
						allRawOutput.push_back(thisHitOutput);
					} else {
						if(VERB > 4) {
							cout << " hit number " << h + 1 << " is rejected in " << hitType << endl;
						}
					}
				}


				string vname = aHit->GetId()[aHit->GetId().size()-1].name;
				if(VERB > 6 || vname.find(catch_v) != string::npos) {
					cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
					cout << aHit->GetId();
					double Etot = 0;
					for(unsigned int e=0; e<aHit->GetPos().size(); e++) Etot = Etot + aHit->GetEdep()[e];
					cout << "   Total energy deposited: " << Etot/MeV << " MeV" << endl;
				}
			}




			if(WRITE_TRUE_INTEGRATED) {
				processOutputFactory->writeG4RawIntegrated(outContainer, allRawOutput, hitType, banksMap);
			}

			// geant4 all raw information
			// by default they are all DISABLED
			// user can enable them one by one
			// using the ALLRAWS option
			if(WRITE_TRUE_ALL) {
				processOutputFactory->writeG4RawAll(outContainer, allRawOutput, hitType, banksMap);
			}

			// geant4 voltage versus time
			// by default they are all DISABLED
			// user can enable them one by one
			// using the SIGNALVT option
			if(SIGNALVT.find(hitType) != string::npos) {
				vector<hitOutput> allVTOutput;
				
				
				for(int h=0; h<nhits; h++) {
					
					hitOutput thisHitOutput;
					
					MHit* aHit = (*MHC)[h];
					
					// process each step to produce a charge/time digitized information / step
					thisHitOutput.setChargeTime(hitProcessRoutine->chargeTime(aHit, h));
					
					vector<double> stepTimes   = thisHitOutput.getChargeTime()[3]; // time at electronics
					vector<double> stepCharges = thisHitOutput.getChargeTime()[2]; // charge at electronics
					vector<double> hardware    = thisHitOutput.getChargeTime()[5]; // crate/slot/channel
					
					map<int, int> vSignal;
					
					// crate, slot, channels as from translation table
					vSignal[0] = hardware[0];           // crate
					vSignal[1] = hardware[1];           // slot
					vSignal[2] = hardware[2];           // channel
					
					// Add comments what are these
					double pedestal_mean = hardware[3];
					double pedestal_sigm = hardware[4];
					
					for(unsigned ts = 0; ts<nsamplings; ts++) {
						double forTime = ts*tsampling;
						double voltage = 0;
						
						// create the voltage output based on the hit process
						// routine voltage(double charge, double time, double forTime)
						for(unsigned s=0; s<stepTimes.size(); s++) {
							
							double stepTime   = stepTimes[s];
							double stepCharge = stepCharges[s];
							
							// cout<<"StepTime = "<<stepTime<<endl;
							// cout<<"StepCharge =  = "<<stepCharge<<endl;
							
							voltage += hitProcessRoutine->voltage(stepCharge, stepTime, forTime);
							
							//cout<<"setpCharge = "<<stepCharge<<endl;
							
							// cout << " hit " << h <<    "step " << s << "  time: " << stepTime
							// << "   charge " << stepCharge << "  voltage " << voltage << "  for time bunch " << ts << endl;
						}
						
						
						// Now pedestal should be calculated, Assume it is a Gaussian
						double pedestal = G4RandGauss::shoot(pedestal_mean, pedestal_sigm);
						
						// need conversion factor from double to int
						// the first 3 entries are crate/slot/channels above
						// the total signal is the pedestal + voltage (from actuall hit), here voltage is actually represents
						// FADC counts
						vSignal[ts+3] = int(pedestal) + (int) voltage;
					}

					thisHitOutput.createQuantumS(vSignal);
					
					hit_outputs_from_AllSD[vSignal[0]].push_back(thisHitOutput);
					allVTOutput.push_back(thisHitOutput);
					
					// this is not written out yet
					
					string vname = aHit->GetId()[aHit->GetId().size()-1].name;
					
					if(VERB > 6 || vname.find(catch_v) != string::npos)
					{
						cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
						cout << aHit->GetId();
						double Etot = 0;
						for(unsigned int e=0; e<aHit->GetPos().size(); e++) Etot = Etot + aHit->GetEdep()[e];
						cout << "   Total energy deposited: " << Etot/MeV << " MeV" << endl;
					}
				}
				processOutputFactory->writeChargeTime(outContainer, allVTOutput, hitType, banksMap);
				
				// Event number (evtN) is needed in FADCMode1, therefore this is also passed as an argument
				//processOutputFactory->writeFADCMode1(outContainer, allVTOutput, evtN);
				
			}
			
			
			// Check whether to save RNG
			if (ssp.enabled && ssp.decision == false)
				for (int h = 0; h < nhits; ++h)
				{
					// Check if masked ID matches targetId
					int id = allDgtOutput[h].getIntDgtVar ("id");
					int id2 = id;
					int j = ssp.tIdsize-1;

					for (; j >= 0; --j)
					{
						if (ssp.targetId[j] != 'x' && id2 % 10 != atoi (ssp.targetId.substr(j,1).c_str()))
							break;
						id2 /= 10;
					}
					if (j >= 0)
						continue;

					// Check pid
					int pid = allRawOutput[h].getIntRawVar ("pid");
					if (pid != ssp.targetPid)
						continue;

					// Check given variable
					double varval = allRawOutput[h].getIntRawVar (ssp.variable);
					if (varval == -99)
						varval = allDgtOutput[h].getIntDgtVar (ssp.variable);
					if (varval == -99)
					{
						cout << "Unknown variable " << ssp.variable << " for SAVE_SELECTED, exiting" << endl;
						exit (0);
					}

					if (varval >= ssp.lowLim && varval <= ssp.hiLim)
					{
						ssp.decision = true;
						break;
					}
				}
			
			delete hitProcessRoutine;
		}
	}
	
	processOutputFactory->writeFADCMode1( hit_outputs_from_AllSD, evtN);
	
	// writing out generated particle infos
	processOutputFactory->writeGenerated(outContainer, MPrimaries, banksMap, gen_action->userInfo);
	
	// For hits, store all ancestors
	
	if (SAVE_ALL_ANCESTORS)
	{
		vector<ancestorInfo> ainfo;
		set<int> storedTraj;
		for (unsigned int i = 0; i < trajectoryContainer->size(); i++)
		{
			G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
			int tid = trj->GetTrackID();
			it = track_db2.find (tid);
			if (it != track_db2.end())
			{
				// This trajectory tid is involved in a hit, store it and its ancestors
				while (tid > 0 && storedTraj.find (tid) == storedTraj.end())
				{
					int mtid = trj->GetParentID();
					ancestorInfo ai;
					ai.pid = trj->GetPDGEncoding();
					ai.tid = tid;
					ai.mtid = mtid;
					ai.trackE = trj->GetInitialKineticEnergy();
					ai.p = trj->GetInitialMomentum();
					ai.vtx = trj->GetPoint(0)->GetPosition();
					ainfo.push_back (ai);
					storedTraj.insert (tid);
					
					// Find trajectory of the mother
					tid = 0;
					if (mtid > 0)
					{
						for (unsigned int ii = 0; ii < trajectoryContainer->size(); ii++)
						{
							trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[ii];
							if (trj->GetTrackID() == mtid)
							{
								tid = mtid;
								break;
							}
						}
					}
				}
			}
		}
		// write out ancestral trajectories
		processOutputFactory->writeAncestors (outContainer, ainfo, getBankFromMap("ancestors", banksMap));
	}
	
	
	processOutputFactory->writeEvent(outContainer);
	delete processOutputFactory;
	
	// Save RNG; can't use G4RunManager::GetRunManager()->rndmSaveThisEvent()
	// because GEANT doesn't know about GEMC run/event numbers
	if (ssp.decision) {
		G4String fileIn  = ssp.dir + "currentEvent.rndm";
		
		std::ostringstream os;
		os << "run" << rw.runNo << "evt" << evtN
		<< ".rndm" << '\0';
		G4String fileOut = ssp.dir + os.str();
		
		G4String copCmd = "/control/shell cp "+fileIn+" "+fileOut;
		G4UImanager::GetUIpointer()->ApplyCommand(copCmd);
	}
	
	if(evtN%Modulo == 0 && VERB > 2 ) {
		cout << hd_msg << " End of Event " << evtN << " Routine..." << endl << endl;
	}
	
	// Increase event number. Notice: this is different than evt->GetEventID()
	evtN++;

	return;
}





void MEventAction::saveBGPartsToLund()
{
	// for lund format see documentation at gemc.jlab.org:
	// https://gemc.jlab.org/gemc/html/documentation/generator/lund.html
	
	// this should work also for no hits, map size will be zero.
	
	*lundOutput << (int) bgMap.size() << "\t" << evtN <<  "\t 0 0 0 0 0 0 0 0 " << endl;
	
	int i = 1;
	for(map<int, BGParts>::iterator it = bgMap.begin(); it != bgMap.end(); it++)
		*lundOutput << i++ << "\t0\t1\t" << it->second.pid << "\t0\t" << it->first << "\t"
		<< it->second.p.x()/GeV << "\t" << it->second.p.y()/GeV << "\t" << it->second.p.z()/GeV << "\t" << it->second.time << "\t0\t"
		<< it->second.v.x()/cm  << "\t" << it->second.v.y()/cm  << "\t" << it->second.v.z()/cm << endl;
}




vector<BackgroundHit*> MEventAction::getNextBackgroundEvent(string forSystem)
{
	if(backgroundHits != nullptr) {
		
		if(backgroundEventNumber.find(forSystem) == backgroundEventNumber.end()) {
			backgroundEventNumber[forSystem] = 1;
		} else {
			backgroundEventNumber[forSystem]++;
		}
		
		map<int, vector<BackgroundHit*> > *backgroundHitsEventMap = backgroundHits->getBackgroundForSystem(forSystem);
		
		if(backgroundHitsEventMap != nullptr) {
			
			// resetting numbering to first one if we reach end of map
			if(backgroundEventNumber[forSystem] == (int) backgroundHitsEventMap->size()) {
				backgroundEventNumber[forSystem] = 1;
			}
			
			if(backgroundHitsEventMap->find(backgroundEventNumber[forSystem]) != backgroundHitsEventMap->end()) {
				return (*backgroundHitsEventMap)[backgroundEventNumber[forSystem]];
			}
		}
	}
	return {nullptr};
}











void MEventAction::setup_clas12_RF(int runno) {
	
	if(last_runno != runno) {
		string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
		string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
		string timestamp = "";
		if(digiSnapshotTime != "no") {
			timestamp = ":"+digiSnapshotTime;
		}
		vector<vector<double> > data;
		// int isec, ilay, istr;

		cout << " RF Setup: Run Number change detected:  " << runno << endl;
		last_runno = runno;
		
		string connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
		if(getenv ("CCDB_CONNECTION") != nullptr) {
			connection = (string) getenv("CCDB_CONNECTION");
		}
						
		unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(connection));
		char   database[80];

		snprintf(database, sizeof(database), "%s:%d:%s%s", "/calibration/eb/rf/config", runno, digiVariation.c_str(), timestamp.c_str());
		cout << " Connecting to " << database << endl;

		data.clear(); calib->GetCalib(data, database);
		
		double clock, prescale;
		
		// taking first entry only
		for(unsigned row = 0; row < 1; row++) {
//		for(unsigned row = 0; row < data.size(); row++) {
//			isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
			clock = data[row][4];
			prescale =  data[row][5];
		}
		
		
		rfvalue_strings = {to_string(clock), to_string(prescale)};
		set_and_show_rf_setup();
        for(auto& rfv: rfvalue_strings) {
            cout << "    - " << rfv << endl;
        }
	}
	

}


void MEventAction::set_and_show_rf_setup() {
	double rf_frquency_from_period = 1.0 / get_number(rfvalue_strings[0]);

	cout << " RF Setup: Period, Frequency [GHz], [prescales]" << endl;
    cout << "    - " << get_number(rfvalue_strings[0]) << endl;
    rfvalue_strings[0] = to_string(rf_frquency_from_period) + " ";
//	for(auto& rfv: rfvalue_strings) {
//		cout << "    - " << rfv << endl;
//	}
}
