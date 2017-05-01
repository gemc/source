// gemc headers
#include "counter_hitprocess.h"

map<string, double> counter_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int id  = identity[0].id;
	
	if(verbosity>4)
		cout << log_msg << " counter detector id: " << id << endl;
	
	dgtz["hitn"] = hitn;
	dgtz["id"]   = id;
	
	
	// now counting the particles
	vector<int> pids = aHit->GetPIDs();
	vector<int> tids = aHit->GetTIds();

	int ngamma, nep, nem, npip, npim, npi0, nkp, nkm, nk0, nproton, nneutron, noptphoton;
	ngamma     = 0;
	nep        = 0;
	nem        = 0;
	npip       = 0;
	npim       = 0;
	npi0       = 0;
	nkp        = 0;
	nkm        = 0;
	nk0        = 0;
	nproton    = 0;
	nneutron   = 0;
	noptphoton = 0;
	
	set<int> trackIds;
	
	for(unsigned int p=0; p<pids.size(); p++) {
		int pid = pids[p];
		int tid = tids[p];
		
		if(trackIds.find(tid) == trackIds.end() ) {
			
				 if(pid == 22)   { ngamma++; }
			else if(pid == 11)   { nem++; }
			else if(pid == -11)  { nep++; }
			else if(pid == 211)  { npip++; }
			else if(pid == -211) { npim++; }
			else if(pid == 111)  { npi0++; }
			else if(pid == 321)  { ngamma++; }
			else if(pid == -321) { ngamma++; }
			else if(pid == 311)  { ngamma++; }
			else if(pid == 2212) { nproton++; }
			else if(pid == 2112) { nneutron++; }
			else if(pid == 0)    { noptphoton++; }
			
			trackIds.insert(tid);
		}
		
	}

	dgtz["ngamma"]     = ngamma;
	dgtz["nep"]        = nep;
	dgtz["nem"]        = nem;
	dgtz["npip"]       = npip;
	dgtz["npim"]       = npim;
	dgtz["npi0"]       = npi0;
	dgtz["nkp"]        = nkp;
	dgtz["nkm"]        = nkm;
	dgtz["nk0"]        = nk0;
	dgtz["nproton"]    = nproton;
	dgtz["nneutron"]   = nneutron;
	dgtz["noptphoton"] = noptphoton;

	return dgtz;
}

vector<identifier>  counter_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> counter_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;
	return noiseHits;
}



map< string, vector <int> >  counter_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}



// - charge: returns charge/time digitized information / step
map< int, vector <double> > counter_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double counter_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}












