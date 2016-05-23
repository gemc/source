// gemc headers
#include "bst_hitprocess.h"
#include "bst_strip.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// geant4
#include "Randomize.hh"

map<string, double> bst_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	class bst_strip bsts;
	bsts.fill_infos();
	
	// double checking dimensions
	double SensorLength = 2.0*aHit->GetDetector().dimensions[2]/mm;  // length of 1 card
	double SensorWidth  = 2.0*aHit->GetDetector().dimensions[0]/mm;  // width 1 card
	
	if(SensorLength != bsts.SensorLength || SensorWidth != bsts.SensorWidth)
		cout << "  Warning: dimensions mismatch between sensor reconstruction dimensions and gemc card dimensions." << endl << endl;
	
	int layer  = 2*identity[0].id + identity[1].id - 2 ;
	int sector = identity[2].id;
	int card   = identity[3].id;
	int strip  = identity[4].id;
	
	trueInfos tInfos(aHit);
	
   // 1 DAC is 0.87 keV
   // If deposited energy is below 26.1 keV there is no hit.
   // If deposited energy is above 117.47 keV the hit will be in the last ADC bin
   // (from 117.47 keV to infinity), which means overflow.
   // I.e. there are 7 bins plus an overflow bin.
   
   double minHit = 0.0261*MeV;
   double maxHit = 0.11747*MeV;
   double deltaADC = maxHit - minHit;
	int adc     = floor(   7*(tInfos.eTot - minHit)/deltaADC);
   int adchd   = floor(8196*(tInfos.eTot - minHit)/deltaADC);
	
   if(tInfos.eTot>maxHit)
   {
      adc   = 7;
      adchd = 8196;
   }
   if(tInfos.eTot<minHit)
   {
      adc   = -5;
      adchd = -5000;
   }
		
	if(verbosity>4)
	{
		cout <<  log_msg << " layer: " << layer << "  sector: " << sector << "  Card: " << card <<  "  Strip: " << strip
			 << " x=" << tInfos.x << " y=" << tInfos.y << " z=" << tInfos.z << endl;
	}
	
	dgtz["hitn"]   = hitn;
	dgtz["layer"]  = layer;
	dgtz["sector"] = sector;
	dgtz["strip"]  = strip;
	dgtz["ADC"]    = adc;
   dgtz["ADCHD"]  = adchd;
	dgtz["time"]   = tInfos.time;
	dgtz["bco"]    = (int) 255*G4UniformRand();
	
	return dgtz;
}



vector<identifier> bst_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	
	G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();
	
	G4ThreeVector  Lxyz    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()  ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);
		
	class bst_strip bsts;
	bsts.fill_infos();
	
	int layer   = 2*yid[0].id + yid[1].id - 2 ;
	int sector  = yid[2].id;
	int isensor = yid[3].id;
	
	vector<double> multi_hit = bsts.FindStrip(layer-1, sector-1, isensor, Lxyz);
	
	int n_multi_hits = multi_hit.size()/2;
	
	// closest strip
	yid[4].id = (int) multi_hit[0];
	
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
	yid[4].id_sharing = multi_hit[1];
	
	// additional strip
	for(int h=1; h<n_multi_hits; h++)
	{
		for(int j=0; j<4; j++)
		{
			identifier this_id;
			this_id.name       = yid[j].name;
			this_id.rule       = yid[j].rule;
			this_id.id         = yid[j].id;
			this_id.time       = yid[j].time;
			this_id.TimeWindow = yid[j].TimeWindow;
			this_id.TrackId    = yid[j].TrackId;
			this_id.id_sharing = multi_hit[3];
			yid.push_back(this_id);
		}
		// last id is strip
		identifier this_id;
		this_id.name       = yid[4].name;
		this_id.rule       = yid[4].rule;
		this_id.id         = (int) multi_hit[2];
		this_id.time       = yid[4].time;
		this_id.TimeWindow = yid[4].TimeWindow;
		this_id.TrackId    = yid[4].TrackId;
		this_id.id_sharing = multi_hit[3];
		yid.push_back(this_id);
		
	}
		
	return yid;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> bst_HitProcess :: electronicNoise()
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


map< string, vector <int> >  bst_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}













