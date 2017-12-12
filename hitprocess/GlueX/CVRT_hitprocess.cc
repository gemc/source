// gemc headers
#include "CVRT_hitprocess.h"

map<string, double> CVRT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;

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




// - charge: returns charge/time digitized information / step
map< int, vector <double> > CVRT_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double CVRT_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}







