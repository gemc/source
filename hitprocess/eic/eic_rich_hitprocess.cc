// G4 headers
#include "G4UnitsTable.hh"
#include "G4Poisson.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "eic_rich_hitprocess.h"

map<string, double> eic_rich_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	trueInfos tInfos(aHit);
//	predefined variable Etot, x, y, z, lx, ly, lz, time
	
	int nsteps = aHit->GetPos().size();
		
	vector<G4ThreeVector> pos  = aHit->GetPos();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double> times = aHit->GetTime();
	vector<G4ThreeVector> p = aHit->GetMoms();	

	dgtz["nsteps"] = nsteps;	
	dgtz["in_px"] = p[0].x();
	dgtz["in_py"] = p[0].y();
	dgtz["in_pz"] = p[0].z();
	dgtz["in_x"] = pos[0].x();
	dgtz["in_y"] = pos[0].y();
	dgtz["in_z"] = pos[0].z();
	dgtz["in_t"] = times[0];
	dgtz["out_px"] = p[nsteps-1].x();
	dgtz["out_py"] = p[nsteps-1].y();
	dgtz["out_pz"] = p[nsteps-1].z();	
	dgtz["out_x"] = pos[nsteps-1].x();
	dgtz["out_y"] = pos[nsteps-1].y();
	dgtz["out_z"] = pos[nsteps-1].z();
	dgtz["out_t"] = times[nsteps-1];
	
	dgtz["hitn"]    = hitn;
	dgtz["pid"]     = (double) aHit->GetPID();
	dgtz["mpid"]    = (double) aHit->GetmPID();
	dgtz["tid"]     = (double) aHit->GetTId();
	dgtz["mtid"]    = (double) aHit->GetmTrackId();
	dgtz["otid"]    = (double) aHit->GetoTrackId();
	dgtz["trackE"]  = aHit->GetE();
	dgtz["totEdep"] = tInfos.eTot;
	dgtz["avg_x"]   = tInfos.x;
	dgtz["avg_y"]   = tInfos.y;
	dgtz["avg_z"]   = tInfos.z;
	dgtz["avg_lx"]  = tInfos.lx;
	dgtz["avg_ly"]  = tInfos.ly;
	dgtz["avg_lz"]  = tInfos.lz;
	dgtz["avg_t"]   = tInfos.time;
	dgtz["px"]      = aHit->GetMom().getX();
	dgtz["py"]      = aHit->GetMom().getY();
	dgtz["pz"]      = aHit->GetMom().getZ();
	dgtz["vx"]      = aHit->GetVert().getX();
	dgtz["vy"]      = aHit->GetVert().getY();
	dgtz["vz"]      = aHit->GetVert().getZ();
	dgtz["mvx"]     = aHit->GetmVert().getX();
	dgtz["mvy"]     = aHit->GetmVert().getY();
	dgtz["mvz"]     = aHit->GetmVert().getZ();	
	
	dgtz["id"]  =  identity[0].id;	
	dgtz["hitn"] = hitn;
	
	return dgtz;  
}

vector<identifier>  eic_rich_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> eic_rich_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;

	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}


map< string, vector <int> >  eic_rich_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > eic_rich_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double eic_rich_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

