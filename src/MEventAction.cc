// G4 headers
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"

// gemc headers
#include "MEventAction.h"
#include "Hit.h"

#include <iostream>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


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
	VERB             = gemcOpt.optMap["BANK_VERBOSITY"].arg ;
	catch_v          = gemcOpt.optMap["CATCH"].args;
	SAVE_ALL_MOTHERS = (int) gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg ;
	gPars            = gpars;
	MAXP             = (int) gemcOpt.optMap["NGENP"].arg;
	rw               = runWeights(opts);
	
	WRITE_ALLRAW     = replaceCharWithChars(gemcOpt.optMap["ALLRAWS"].args, ",", "  ");
	WRITE_INTRAW     = replaceCharWithChars(gemcOpt.optMap["INTEGRATEDRAW"].args, ",", "  ");
	WRITE_INTDGT     = replaceCharWithChars(gemcOpt.optMap["INTEGRATEDDGT"].args, ",", "  ");
	SIGNALVT         = replaceCharWithChars(gemcOpt.optMap["SIGNALVT"].args, ",", "  ");
	
	if(SAVE_ALL_MOTHERS>1)
	{
		lundOutput = new ofstream("background.dat");
		cout << " > Opening background.dat file to save background particles in LUND format." << endl;
	}
	
	evtN = gemcOpt.optMap["EVTN"].arg;
}

MEventAction::~MEventAction()
{
	if(SAVE_ALL_MOTHERS>1)
		lundOutput->close();
}

void MEventAction::BeginOfEventAction(const G4Event* evt)
{
	rw.getRunNumber(evtN);
	bgMap.clear();
	
	if(evtN%Modulo == 0 )
	{
		cout << hd_msg << " Begin of event " << evtN << "  Run Number: " << rw.runNo;
		if(rw.isNewRun) cout << " (new) ";
		cout << endl;
		cout << hd_msg << " Random Number: " << G4UniformRand() << endl;
		// CLHEP::HepRandom::showEngineStatus();
	}
}

void MEventAction::EndOfEventAction(const G4Event* evt)
{
	
	MHitCollection* MHC;
	int nhits;
	
	if(evtN%Modulo == 0 )
		cout << hd_msg << " Starting Event Action Routine " << evtN << "  Run Number: " << rw.runNo << endl;
	
	
	// building the tracks set database with all the tracks in all the hits
	// if SAVE_ALL_MOTHERS is set
	set<int> track_db;
	if(SAVE_ALL_MOTHERS)
		for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
		{
			string hitType;
			
			MHC = it->second->GetMHitCollection();
			nhits = MHC->GetSize();
			
			for(int h=0; h<nhits; h++)
			{
				vector<int>           tids = (*MHC)[h]->GetTIds();
				
				for(unsigned int t=0; t<tids.size(); t++)
					track_db.insert(tids[t]);
				
				
				if(SAVE_ALL_MOTHERS>1)
				{
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
	
	
	// now filling the map of tinfos with tracks infos from the track_db database
	// this won't get the mother particle infos except for their track ID
	map<int, TInfos> tinfos;
	set<int>::iterator it;
	
	// the container is full only if /tracking/storeTrajectory 2
	G4TrajectoryContainer *trajectoryContainer;
	
	if(SAVE_ALL_MOTHERS)
	{
		trajectoryContainer = evt->GetTrajectoryContainer();
		momDaughter.clear();
		
		if(VERB>3)
			cout << " >> Total number of tracks " << trajectoryContainer->size() << endl;
		
		while(trajectoryContainer && track_db.size())
		{
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
		for(map<int, int>::iterator itm = momDaughter.begin(); itm != momDaughter.end(); itm++)
		{
			int ancestor = itm->first;
			if(momDaughter[ancestor] == 0)
				hierarchy[itm->first] = itm->first;
			
			while (momDaughter[ancestor] != 0)
			{
				hierarchy[itm->first] = momDaughter[ancestor];
				ancestor = momDaughter[ancestor];
			}
		}
		
		// now accessing the mother particle infos
		for(map<int, TInfos>::iterator itm = tinfos.begin(); itm != tinfos.end(); itm++)
		{
			int mtid = (*itm).second.mtid;
			// looking for mtid infos in the trajectoryContainer
			for(unsigned int i=0; i< trajectoryContainer->size() && mtid != 0; i++)
			{
				G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
				int tid = trj->GetTrackID();
				if(tid == mtid)
				{
					tinfos[(*itm).first].mpid   = trj->GetPDGEncoding();
					tinfos[(*itm).first].mv     = trj->GetPoint(0)->GetPosition();
				}
			}
		}
		// removing daughter tracks if mother also gave hits
		if(SAVE_ALL_MOTHERS>2)
		{
			vector<int> bgtIDs;
			for(map<int, BGParts>::iterator it = bgMap.begin(); it != bgMap.end(); it++)
				bgtIDs.push_back(it->first);
		
			for(unsigned i=0; i<bgtIDs.size(); i++)
			{
				int daughter = bgtIDs[i];
				int momCheck = 0;
				if(momDaughter.find(daughter) != momDaughter.end())
					momCheck = momDaughter[daughter];
				
				while(momCheck != 0)
				{
					// mom has a hit
					if(bgMap.find(momCheck) != bgMap.end())
					{
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
	if(ito == outputFactoryMap->end())
	{
		if(outContainer->outType != "no")
			cout << hd_msg << " Warning: output type <" << outContainer->outType
			<< "> is not registered in the outContainerput Factory. " << endl
			<< "      This event will not be written out." << endl;
		evtN++;
		return;
	}
	outputFactory *processOutputFactory = getOutputFactory(outputFactoryMap, outContainer->outType);
	
	
	
	
	// Header Bank contains event number
	// Need to change this to read DB header bank
	map<string, double> header;
	header["runNo"]    = rw.runNo;
	header["evn"]      = evtN;
	header["evn_type"] = -1;  // physics event. Negative is MonteCarlo event
	header["beamPol"]  = gen_action->getBeamPol();
	
	for(unsigned i=0; i<gen_action->lundUserDefined.size(); i++)
	{
		string tmp = "var" + stringify((int) i+1);
		header[tmp] = gen_action->lundUserDefined[i];
	}
	
	// write event header bank
	processOutputFactory->writeHeader(outContainer, header, getBankFromMap("header", banksMap));
	
	// Getting Generated Particles info
	// Are these loops necessary, revisit later1
	vector<generatedParticle> MPrimaries;
	for(int pv=0; pv<evt->GetNumberOfPrimaryVertex() && pv<MAXP; pv++)
	{
		generatedParticle Mparticle;
		G4PrimaryVertex* MPV = evt->GetPrimaryVertex(pv);
		Mparticle.vertex = MPV->GetPosition();
		double thisTime = MPV->GetT0();
		int thisMult    = MPV->GetNumberOfParticle();

		for(int pp = 0; pp<MPV->GetNumberOfParticle() && pv<MAXP; pp++)
		{
			G4PrimaryParticle *PP  = MPV->GetPrimary(pp);
			Mparticle.momentum     = PP->GetMomentum();
			Mparticle.PID          = PP->GetPDGcode();
			Mparticle.time         = thisTime;
			Mparticle.multiplicity = thisMult;
		}
		MPrimaries.push_back(Mparticle)  ;
	}
	
	if(SAVE_ALL_MOTHERS>1)
		saveBGPartsToLund();
	
	for(map<string, sensitiveDetector*>::iterator it = SeDe_Map.begin(); it!= SeDe_Map.end(); it++)
	{
		MHC = it->second->GetMHitCollection();
		nhits = MHC->GetSize();
		
		// The same ProcessHit Routine must apply to all the hits  in this HitCollection.
		// Instantiating the ProcessHitRoutine only once for the first hit.
		if(nhits)
		{
			// the bank idtag is the one that corresponds to the hitType
			//MHit* aHit = (*MHC)[0];
			string vname = (*MHC)[0]->GetDetector().name;
			string hitType = it->second->GetDetectorHitType(vname);
			
			HitProcess *hitProcessRoutine = getHitProcess(hitProcessMap, hitType);
			if(!hitProcessRoutine)
				return;
			
			hitProcessRoutine->init(hitType, gemcOpt, gPars);
			
			bool WRITE_TRUE_INTEGRATED = 0;
			bool WRITE_TRUE_ALL = 0;
			if(WRITE_INTRAW.find(hitType) != string::npos) WRITE_TRUE_INTEGRATED = 1;
			if(WRITE_ALLRAW.find(hitType) != string::npos) WRITE_TRUE_ALL = 1;
			
			vector<hitOutput> allRawOutput;
			
			// creating summary information for each generated particle
			for(unsigned pi = 0; pi<MPrimaries.size(); pi++)
				MPrimaries[pi].pSum.push_back(summaryForParticle("na"));
			
			
			for(int h=0; h<nhits; h++)
			{
				MHit* aHit = (*MHC)[h];
				if(aHit->isElectronicNoise)
					continue;

				hitOutput thisHitOutput;
				
				// mother particle infos
				if(SAVE_ALL_MOTHERS)
				{
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
					for(unsigned pi = 0; pi<MPrimaries.size(); pi++)
						hitByPrimary.push_back(0);
					
					// all these vector have the same length.
					for(unsigned pi = 0; pi<MPrimaries.size(); pi++)
					{
						vector<double> edeps = aHit->GetEdep();
						vector<double> times = aHit->GetTime();
						MPrimaries[pi].pSum.back().nphe = aHit->GetTIds().size();
						for(unsigned ss =0; ss<edeps.size(); ss++)
						{
							if(otids[ss] == (int) pi+1)
							{
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
				}
				else
				{
					// filling mother infos with zeros
					int thisHitSize = aHit->GetId().size();
					vector<int>           zint  = vector_zint(thisHitSize);
					vector<G4ThreeVector> zthre = vector_zthre(thisHitSize);
					aHit->SetoTrackIds(zint);
					aHit->SetmTrackIds(zint);
					aHit->SetmPIDs(    zint);
					aHit->SetmVerts(   zthre);
				}
				
				thisHitOutput.setRaws(hitProcessRoutine->integrateRaw(aHit, h+1, WRITE_TRUE_INTEGRATED));
				
				if(WRITE_TRUE_ALL)
					thisHitOutput.setAllRaws(hitProcessRoutine->allRaws(aHit, h+1));
				
				
				allRawOutput.push_back(thisHitOutput);
				
				string vname = aHit->GetId()[aHit->GetId().size()-1].name;
				if(VERB > 4 || vname.find(catch_v) != string::npos)
				{
					cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
					cout << aHit->GetId();
					double Etot = 0;
					for(unsigned int e=0; e<aHit->GetPos().size(); e++) Etot = Etot + aHit->GetEdep()[e];
					cout << "   Total energy deposited: " << Etot/MeV << " MeV" << endl;
				}
			}
			
			// geant4 integrated raw information
			// by default they are all DISABLED
			// user can enable them one by one
			// using the INTEGRATEDRAW option
			if(WRITE_TRUE_INTEGRATED)
				processOutputFactory->writeG4RawIntegrated(outContainer, allRawOutput, hitType, banksMap);
			
			// geant4 all raw information
			// by default they are all DISABLED
			// user can enable them one by one
			// using the ALLRAWS option
			if(WRITE_TRUE_ALL)
				processOutputFactory->writeG4RawAll(outContainer, allRawOutput, hitType, banksMap);
			
			
			
			if(VERB > 4)
				for(unsigned pi = 0; pi<MPrimaries.size(); pi++)
				{
					cout << " Particle " << pi + 1 << " has " << MPrimaries[pi].pSum.size() << " particle summaries:" << endl;
					for(unsigned ss =0; ss<MPrimaries[pi].pSum.size(); ss++)
					{
						cout << " \t det: " << MPrimaries[pi].pSum[ss].dname <<
						"  Etot: "  <<  MPrimaries[pi].pSum[ss].etot <<
						"  time: "  <<  MPrimaries[pi].pSum[ss].t << endl;
					}
					cout << endl;
					
				}
			
			// geant4 integrated digitized information
			// by default they are all ENABLED
			// user can disable them one by one
			// using the INTEGRATEDDGT option
			if(WRITE_INTDGT.find(hitType) == string::npos)
			{
				hitProcessRoutine->initWithRunNumber(rw.runNo);
				
				vector<hitOutput> allDgtOutput;
				for(int h=0; h<nhits; h++)
				{
					
					hitOutput thisHitOutput;
					MHit* aHit = (*MHC)[h];
					
					thisHitOutput.setDgtz(hitProcessRoutine->integrateDgt(aHit, h+1));
					allDgtOutput.push_back(thisHitOutput);
					
					string vname = aHit->GetId()[aHit->GetId().size()-1].name;
					if(VERB > 4 || vname.find(catch_v) != string::npos)
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
			
			
			// geant4 voltage versus time
			// by default they are all DISABLED
			// user can enable them one by one
			// using the SIGNALVT option
			if(SIGNALVT.find(hitType) != string::npos)
			{
				vector<hitOutput> allVTOutput;
				
				for(int h=0; h<nhits; h++)
				{
					hitOutput thisHitOutput;
					
					MHit* aHit = (*MHC)[h];
					
					thisHitOutput.setSignal(hitProcessRoutine->signalVT(aHit));
					thisHitOutput.setQuantumS(hitProcessRoutine->quantumS(thisHitOutput.getSignalVT(), aHit));
					
					allVTOutput.push_back(thisHitOutput);
					
					// this is not written out yet
					
					string vname = aHit->GetId()[aHit->GetId().size()-1].name;
					if(VERB > 4 || vname.find(catch_v) != string::npos)
					{
						cout << hd_msg << " Hit " << h + 1 << " --  total number of steps this hit: " << aHit->GetPos().size() << endl;
						cout << aHit->GetId();
						double Etot = 0;
						for(unsigned int e=0; e<aHit->GetPos().size(); e++) Etot = Etot + aHit->GetEdep()[e];
						cout << "   Total energy deposited: " << Etot/MeV << " MeV" << endl;
					}
				}
			}
			
			delete hitProcessRoutine;
		}
	}
	// writing out generated particle infos
	processOutputFactory->writeGenerated(outContainer, MPrimaries, banksMap);
	
	processOutputFactory->writeEvent(outContainer);
	delete processOutputFactory;
	
	
	if(evtN%Modulo == 0 )
		cout << hd_msg << " End of Event " << evtN << " Routine..." << endl << endl;
	
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















