// gemc headers
#include "CVRT_hitprocess.h"

map<string, double> CVRT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	return dgtz;
}

vector<identifier>  CVRT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	return id;
}




map< string, vector <int> >  CVRT_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> CVRT_HitProcess :: electronicNoise()
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









