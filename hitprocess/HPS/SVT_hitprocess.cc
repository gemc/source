// %%%%%%%%%%%%
// gemc headers
// %%%%%%%%%%%%
#include "SVT_hitprocess.h"

map<string, double> SVT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// STR ID:
	// layer, type, sector, module, strip
	// class str_strip strs;
	// strs.fill_infos();
	
	// double checking dimensions
	//  double CardLength = 2.0*aHit->GetDetector().dimensions[2]/mm;  // length of 1 card
	// double CardWidth  = 2.0*aHit->GetDetector().dimensions[0]/mm;  // width 1 card
	//if(CardLength != strs.CardLength || CardWidth != strs.CardWidth)
	//cout << hd_msg << "  Warning: dimensions mismatch between card reconstruction dimensions and gemc card dimensions." << endl << endl;
	
	int slayer  = identity[0].id;
	int stype   = identity[1].id;
	int segment = identity[2].id;
	int module  = identity[3].id;
	int strip   = identity[4].id;
	//
	if(verbosity>4)
	{
		trueInfos tInfos(aHit);

		cout <<  log_msg << " layer: " << slayer << "  type: " << stype << "  segment: " << segment << "  module: "<< module
			 << "  Strip: " << strip << " x=" << tInfos.x << " y=" << tInfos.y << " z=" << tInfos.z << endl;
	}
 	dgtz["hitn"]    = hitn;
	dgtz["slayer"]  = slayer;
	dgtz["stype"]   = stype;
	dgtz["segment"] = segment;
	dgtz["module"]  = module;
	dgtz["strip"]   = strip;
	
	return dgtz;
}

#define ABS_(x) (x < 0 ? -x : x)

vector<identifier> SVT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	
	enum STR_identifiers {layer=0,type=1,segment=2,module=3,strip=4};
	
	//int slayer  = yid[0].id;
	//int stype   = yid[1].id;
	//int segment = yid[2].id;
	//int module  = yid[3].id;
	//int strip   = yid[4].id;
	
	
	
	double active_width  = 38.34; // mm  X-direction
	double active_height = 98.33; // mm  Y-direction
	double readout_pitch = 0.06; // mm, = 60 micron
	
	G4ThreeVector   xyz    = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector  Lxyz    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz); // Local coordinates
	
	double x = Lxyz.x();
	double y = Lxyz.y();
	
	if( x < -active_width/2)      yid[strip].id = -10;
	else if ( x > active_width/2) yid[strip].id = -9;
	else if( y < -active_height/2)yid[strip].id = -8;
	else if( y > active_height/2) yid[strip].id = -7;
	else
	{
		int nstrip = (int)((x + active_width/2 )/readout_pitch);
		yid[strip].id = nstrip;
	}
	yid[id.size()-1].id_sharing = 1;
	return yid;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> SVT_HitProcess :: electronicNoise()
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



map< string, vector <int> >  SVT_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > SVT_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double SVT_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}












