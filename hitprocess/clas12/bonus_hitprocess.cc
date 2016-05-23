// gemc headers
#include "bonus_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> bonus_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// true information
	// for example tInfos.eTot is total energy deposited
	trueInfos tInfos(aHit);
	
	// local variable for each step
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	
	// energy at each step
	// so tInfos.eTot is the sum of all steps s of Edep[s] 
	vector<double>      Edep = aHit->GetEdep();

	
	// Distances from left, right
	double dLeft  =  tInfos.ly;
	double dRight =  tInfos.ly;
	
	// dummy formulas for now, parameters could come from DB
	int ADCL = (int) (100*tInfos.eTot*exp(-dLeft/2  ));
	int ADCR = (int) (100*tInfos.eTot*exp(-dRight/2 ));
	
	// speed of light is 30 cm/s
	int TDCL = (int) (100*(tInfos.time/ns + dLeft/cm/30.0));
	int TDCR = (int) (100*(tInfos.time/ns + dRight/cm/30.0));
	
	dgtz["hitn"]   = hitn;
	dgtz["ADCL"]   = ADCL;
	dgtz["ADCR"]   = ADCR;
	dgtz["TDCL"]   = TDCL;
	dgtz["TDCR"]   = TDCR;
	
	return dgtz;
}

vector<identifier>  bonus_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> bonus_HitProcess :: electronicNoise()
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


map< string, vector <int> >  bonus_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

