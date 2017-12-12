// gemc headers
#include "muon_hodo_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> muon_hodo_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;

	vector<identifier> identity = aHit->GetId();

	int idx = identity[0].id;
	int idy = identity[1].id;
	int idz = identity[2].id;
  	trueInfos tInfos(aHit);

	// Get the paddle length: in HS  paddles are boxes, it's the x
	double length = aHit->GetDetector().dimensions[0];
	
	// Distances from left, right
	double dLeft  = length + tInfos.lx;
	double dRight = length - tInfos.lx;
	
	// dummy formulas for now, parameters could come from DB
	int adcl = (int) (100*tInfos.eTot*exp(-dLeft/length/2  ));
	int adcr = (int) (100*tInfos.eTot*exp(-dRight/length/2 ));
	
	// speed of light is 30 cm/s
	int tdcl = (int) (100*(tInfos.time/ns + dLeft/cm/30.0));
	int tdcr = (int) (100*(tInfos.time/ns + dRight/cm/30.0));
	
	dgtz["hitn"]   = hitn;
	dgtz["idx"]    = idx;
	dgtz["idy"]    = idy;
	dgtz["idz"]    = idz;
	dgtz["adcl"]   = adcl;
	dgtz["adcr"]   = adcr;
	dgtz["tdcl"]   = tdcl;
	dgtz["tdcr"]   = tdcr;
	
	return dgtz;
}

vector<identifier>  muon_hodo_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> muon_hodo_HitProcess :: electronicNoise()
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



map< string, vector <int> >  muon_hodo_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}



// - charge: returns charge/time digitized information / step
map< int, vector <double> > muon_hodo_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double muon_hodo_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}








