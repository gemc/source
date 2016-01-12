// gemc headers
#include "bubble_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> bubble_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	int thisPid = aHit->GetPID();
	double totEnergy = aHit->GetE();
	
	if(thisPid == 11 || thisPid == -11) totEnergy -= electron_mass_c2;

	dgtz["detId"] = identity[0].id;
	dgtz["kinE"]  = totEnergy;
	dgtz["pid"]   = thisPid;
	dgtz["hitn"]  = hitn;
	
	return dgtz;
}

vector<identifier>  bubble_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


map< string, vector <int> >  bubble_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

