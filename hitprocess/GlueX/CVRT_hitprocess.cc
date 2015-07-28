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











