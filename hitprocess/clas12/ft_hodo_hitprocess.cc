// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "ft_hodo_hitprocess.h"

map<string, double> ft_hodo_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);
	
	// use Crystal ID to define IDX and IDY
	int ID    = identity[0].id;
	int Layer = identity[1].id;
	
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = 8191;

	if(tInfos.eTot>0)
	{
	  ADC = (int) (tInfos.eTot*100);
	  TDC = (int) (tInfos.time *100);
	}
	  
	
	dgtz["hitn"]  = hitn;
	dgtz["id"]    = ID;
	dgtz["layer"] = Layer;
	dgtz["adc"]   = ADC;
	dgtz["tdc"]   = TDC;
		
	return dgtz;
}

vector<identifier>  ft_hodo_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ft_hodo_HitProcess :: electronicNoise()
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



map< string, vector <int> >  ft_hodo_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}












