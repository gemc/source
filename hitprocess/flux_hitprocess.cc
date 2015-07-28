// gemc headers
#include "flux_hitprocess.h"

map<string, double> flux_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int id  = identity[0].id;
	
	if(verbosity>4)
		cout << log_msg << " flux detector id: " << id << endl;
	
	dgtz["hitn"] = hitn;
	dgtz["id"]   = id;
	
	return dgtz;
}

vector<identifier>  flux_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}





map< string, vector <int> >  flux_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}













