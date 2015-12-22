// gemc headers
#include "bubble_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> bubble_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// true information
	// for example tInfos.eTot is total energy deposited
	trueInfos tInfos(aHit);
	
	dgtz["detId"]   = identity[0].id;;
	dgtz["totEDep"] = tInfos.eTot;
	dgtz["pid"]     = aHit->GetPID();
	dgtz["hitn"]    = hitn;
	
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

