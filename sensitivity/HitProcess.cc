// gemc headers
#include "HitProcess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

HitProcess *getHitProcess(map<string, HitProcess_Factory> *hitProcessMap, string HCname)
{
	if(HCname == "no")
		return NULL;

	if(hitProcessMap->find(HCname) == hitProcessMap->end())
	{
		cout << endl << "  !!! Error: >" << HCname << "< NOT FOUND IN  ProcessHit Map - exiting." << endl;
		exit(0);
		return NULL;
	}
	return (*hitProcessMap)[HCname]();
}

// returns the list of Hit Factories registered
set<string> getListOfHitProcessHit(map<string, HitProcess_Factory> hitProcessMap)
{
	set<string> listF;
	
	for(map<string, HitProcess_Factory>::iterator it = hitProcessMap.begin(); it != hitProcessMap.end(); it++)
		listF.insert(it->first);
	
	return listF;
}


// - integrateRaw: returns geant4 raw information integrated over the hit
map<string, double> HitProcess::integrateRaw(MHit* aHit, int hitn, bool WRITEBANK)
{
	map<string, double> raws;

	if(WRITEBANK) {

		if(aHit->isBackgroundHit == 1) {
			raws["hitn"]    = hitn;
			raws["totEdep"] = aHit->GetEdep().front();
			raws["avg_t"]   = aHit->GetTime().front();
			raws["procID"]  = -1;
			raws["nsteps"]  = 1;

		} else {

			trueInfos tInfos(aHit);

			raws["hitn"]    = hitn;
			raws["pid"]     = (double) aHit->GetPID();
			raws["mpid"]    = (double) aHit->GetmPID();
			raws["tid"]     = (double) aHit->GetTId();
			raws["mtid"]    = (double) aHit->GetmTrackId();
			raws["otid"]    = (double) aHit->GetoTrackId();
			raws["trackE"]  = aHit->GetE();
			raws["totEdep"] = tInfos.eTot;
			raws["avg_x"]   = tInfos.x;
			raws["avg_y"]   = tInfos.y;
			raws["avg_z"]   = tInfos.z;
			raws["avg_lx"]  = tInfos.lx;
			raws["avg_ly"]  = tInfos.ly;
			raws["avg_lz"]  = tInfos.lz;
			raws["px"]      = aHit->GetMom().getX();
			raws["py"]      = aHit->GetMom().getY();
			raws["pz"]      = aHit->GetMom().getZ();
			raws["vx"]      = aHit->GetVert().getX();
			raws["vy"]      = aHit->GetVert().getY();
			raws["vz"]      = aHit->GetVert().getZ();
			raws["mvx"]     = aHit->GetmVert().getX();
			raws["mvy"]     = aHit->GetmVert().getY();
			raws["mvz"]     = aHit->GetmVert().getZ();
			raws["avg_t"]   = tInfos.time;
			raws["procID"]  = aHit->GetProcID();
			raws["nsteps"]  = aHit->GetPIDs().size();
		}
	}
	return raws;
}



map< string, vector <double> > HitProcess::allRaws(MHit* aHit, int hitn)
{
	map< string, vector <double> > allRaws;
	
	vector<int>    pids = aHit->GetPIDs();
	vector<double> pids_d(pids.begin(), pids.end());

	vector<int>    mpids = aHit->GetmPIDs();
	vector<double> mpids_d(mpids.begin(), mpids.end());

	vector<int>    tids = aHit->GetTIds();
	vector<double> tids_d(tids.begin(), tids.end());

	vector<int>    mtids = aHit->GetmTrackIds();
	vector<double> mtids_d(mtids.begin(), mtids.end());

	vector<int>    otids = aHit->GetoTrackIds();
	vector<double> otids_d(otids.begin(), otids.end());


	vector<double> hitnd(pids.size(), (double) hitn);
	
	vector<G4ThreeVector>  gpos = aHit->GetPos();
	vector<G4ThreeVector>  lpos = aHit->GetLPos();
	vector<G4ThreeVector>  mome = aHit->GetMoms();
	vector<G4ThreeVector>  vert = aHit->GetVerts();
	vector<G4ThreeVector>  mver = aHit->GetmVerts();


	vector<double> x, y, z;
	vector<double> lx, ly, lz;
	vector<double> px, py, pz;
	vector<double> vx, vy, vz;
	vector<double> mvx, mvy, mvz;
	vector<double> stepi;
	
	for(unsigned s=0; s<gpos.size(); s++)
	{
	
		stepi.push_back(s+1);

		x.push_back(gpos[s].getX());
		y.push_back(gpos[s].getY());
		z.push_back(gpos[s].getZ());
		
		lx.push_back(lpos[s].getX());
		ly.push_back(lpos[s].getY());
		lz.push_back(lpos[s].getZ());
		
		px.push_back(mome[s].getX());
		py.push_back(mome[s].getY());
		pz.push_back(mome[s].getZ());
		
		vx.push_back(vert[s].getX());
		vy.push_back(vert[s].getY());
		vz.push_back(vert[s].getZ());
		
		mvx.push_back(mver[s].getX());
		mvy.push_back(mver[s].getY());
		mvz.push_back(mver[s].getZ());
	}

	
	allRaws["stepn"]   = stepi;
	allRaws["hitn"]    = hitnd;
	allRaws["pid"]     = pids_d;
	allRaws["mpid"]    = mpids_d;
	allRaws["tid"]     = tids_d;
	allRaws["mtid"]    = mtids_d;
	allRaws["otid"]    = otids_d;
	allRaws["trackE"]  = aHit->GetEs();
	allRaws["edep"]    = aHit->GetEdep();
	allRaws["t"]       = aHit->GetTime();
	allRaws["x"]       = x;
	allRaws["y"]       = y;
	allRaws["z"]       = z;
	allRaws["lx"]      = lx;
	allRaws["ly"]      = ly;
	allRaws["lz"]      = lz;
	allRaws["px"]      = px;
	allRaws["py"]      = py;
	allRaws["pz"]      = pz;
	allRaws["vx"]      = vx;
	allRaws["vy"]      = vy;
	allRaws["vz"]      = vz;
	allRaws["mvx"]     = mvx;
	allRaws["mvy"]     = mvy;
	allRaws["mvz"]     = mvz;
	
	return allRaws;
}




trueInfos::trueInfos(MHit* aHit)
{
	eTot = 0;
	time = 0;
	x = y = z = lx = ly = lz = 0;

	// getting vectors of energy deposited, positions and times
	// nsteps is the size

	vector<G4double> Edep       = aHit->GetEdep();
	vector<G4ThreeVector> pos   = aHit->GetPos();
	vector<G4ThreeVector> Lpos  = aHit->GetLPos();
	vector<G4double>      times = aHit->GetTime();

	nsteps = Edep.size();
	for(unsigned int s=0; s<nsteps; s++)
		eTot += Edep[s];

	// if energy deposited, averaging quantities over the hit
	// weighted by energy deposited
	if(eTot)
	{
		for(unsigned int s=0; s<nsteps; s++)
		{
			x    +=  pos[s].x()*Edep[s];
			y    +=  pos[s].y()*Edep[s];
			z    +=  pos[s].z()*Edep[s];
			lx   += Lpos[s].x()*Edep[s];
			ly   += Lpos[s].y()*Edep[s];
			lz   += Lpos[s].z()*Edep[s];
			time += times[s]*Edep[s];
		}
		x    = x/eTot;
		y    = y/eTot;
		z    = z/eTot;
		lx   = lx/eTot;
		ly   = ly/eTot;
		lz   = lz/eTot;
		time = time/eTot;
	}
	else
	{
		for(unsigned int s=0; s<nsteps; s++)
		{
			x    +=  pos[s].x();
			y    +=  pos[s].y();
			z    +=  pos[s].z();
			lx   += Lpos[s].x();
			ly   += Lpos[s].y();
			lz   += Lpos[s].z();
			time += times[s];
		}
		x    = x/nsteps;
		y    = y/nsteps;
		z    = z/nsteps;
		lx   = lx/nsteps;
		ly   = ly/nsteps;
		lz   = lz/nsteps;
		time = time/nsteps;
	}
}



