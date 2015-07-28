// gemc headers
#include "LTCC_hitprocess.h"

// C++ headers
#include <set>

map<string, double> LTCC_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// if the particle is not an opticalphoton
	// assuming it will create one photo-electron anyway
	// if(aHit->GetPID() != 0) return out;
	
	vector<int> tids = aHit->GetTIds();
	vector<int> pids = aHit->GetPIDs();
	set<int> TIDS;
	
	
	for(unsigned int s=0; s<nsteps; s++)
	{
		// no conditions to count it as one photo-electron
		// a particle hitting on of the anodes directly
		// may produce in average one photo-electron?
		// if(pids[s] == 0)
		TIDS.insert(tids[s]);
	}
	
	int sector = identity[0].id;
	int side   = identity[1].id;
	int pmt    = identity[2].id;
	
	if(verbosity>4)
		cout <<  log_msg << " pmt: " << pmt << " x=" << x/cm << " y=" << y/cm << " z=" << z/cm << endl;
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["side"]   = side;
	dgtz["pmt"]    = pmt;
	dgtz["nphe"]   = TIDS.size();
	
	return dgtz;
}


vector<identifier>  LTCC_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}





map< string, vector <int> >  LTCC_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}






