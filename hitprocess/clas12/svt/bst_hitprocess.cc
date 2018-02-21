// gemc headers
#include "bst_hitprocess.h"
#include "bst_strip.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// geant4
#include "Randomize.hh"

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static bstConstants initializeBSTConstants(int runno)
{
	// all these constants should be read from CCDB
	bstConstants bstc;


	// database
	bstc.runNo = runno;
	bstc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
	bstc.connection = (string) getenv("CCDB_CONNECTION");
	else
	bstc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	bstc.variation  = "main";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(bstc.connection));

	return bstc;
}

map<string, double> bst_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	double minHit = 0.0261*MeV;
	double maxHit = 0.11747*MeV;
	double deltaADC = maxHit - minHit;


	if(aHit->isBackgroundHit == 1) {

		vector<double> eDep = aHit->GetEdep();
		vector<double> stepTime = aHit->GetTime();

		// background hit has all the energy in the first step
		double totEdep = eDep[0];

		int adc     = floor(   7*(totEdep - minHit)/deltaADC);
		int adchd   = floor(8196*(totEdep - minHit)/deltaADC);

		dgtz["hitn"]   = hitn;
		dgtz["sector"] = identity[0].id;
		dgtz["layer"]  = identity[1].id;
		dgtz["strip"]  = identity[2].id;
		dgtz["ADC"]    = adc;
		dgtz["ADCHD"]  = adchd;
		dgtz["time"]   = stepTime[0];
		dgtz["bco"]    = (int) 255*G4UniformRand();

		return dgtz;
	}

	class bst_strip bsts;
	bsts.fill_infos();

	if(!aHit->isElectronicNoise) {
	// double checking dimensions
	double SensorLength = 2.0*aHit->GetDetector().dimensions[2]/mm;  // length of 1 card
	double SensorWidth  = 2.0*aHit->GetDetector().dimensions[0]/mm;  // width 1 card

	double diffLen = SensorLength - bsts.SensorLength;
	double diffWid = SensorWidth  - bsts.SensorWidth;


	// there may be precision issues that's why it is not compared to zero
	if( diffLen > 1E-3 || diffWid > 1E-3 )
		cout << "  Warning: dimensions mismatch between sensor reconstruction dimensions and gemc card dimensions." << endl << endl;
	}

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

    // decide if write an hit or not
    writeHit = true;
    // define conditions to reject hit
    bool rejectHitConditions = false;
    if(rejectHitConditions) {
        writeHit = false;
    }

	return dgtz;
}



vector<identifier> bst_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	// yid is the current strip identifier.
	// it has 5 dimensions:
	// 0. superlayer
	// 1. region
	// 2. sector
	// 3. sensor
	// 4. strip id
    // the first 4 are set in the identifer:
	// superlayer manual 1 type manual 1 segment manual 3 module manual 1 strip manual 1
	// the strip id is set by this routine
	//
	// this routine ALSO identifies contiguos strips and their energy sharing


	// this strip original identifier - the dimension is ONE identifier
	// later we'll push_back the contiguos strips
	vector<identifier> yid = id;

	G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();

	G4ThreeVector  Lxyz    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()  ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);

	class bst_strip bsts;
	bsts.fill_infos();

	int layer   = 2*yid[0].id + yid[1].id - 2 ;
	int sector  = yid[2].id;
	int isensor = yid[3].id;

	// FindStrip returns the strips index that are hit
	// this vector odd values are the strip numbers, the second number is the energy sharing
	// for bst there can be:
	// - no energy sharing
	// - energy sharing to ONE contiguos (left OR right) strip
	vector<double> multi_hit = bsts.FindStrip(layer-1, sector-1, isensor, Lxyz);

	// n_multi_hits: number of contiguous strips
	// it is either single strip (n_multi_hits = 1) or both strips are hit (n_multi_hits=2)
	int n_multi_hits = multi_hit.size()/2;

	// closest strip: this identifier is always set
	// assigning the closest strip number to identifier id
	yid[4].id = (int) multi_hit[0];

	// assigning variable id_sharing to all the identifiers
	// for the first strip
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
	yid[4].id_sharing = multi_hit[1];


	// additional strip
	if(n_multi_hits == 2) {
		// the first four identifiers are
		// 0. superlayer
		// 1. region
		// 2. sector
		// 3. sensor
		// if n_multi_hits has size two:
		// creating ONE additional identifier vector.
		// the first 4 identifier are the same as the original strip, so copying them
		// the sharing percentage for that strip is multi_hit[3]
		for(int j=0; j<4; j++) {
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

		// last identifier for the additional hit is is the strip id
		// copying infos from original strip
		identifier this_id;
		this_id.name       = yid[4].name;
		this_id.rule       = yid[4].rule;
		this_id.time       = yid[4].time;
		this_id.TimeWindow = yid[4].TimeWindow;
		this_id.TrackId    = yid[4].TrackId;

		// now assigning strip ID and sharing percentage
		this_id.id         = (int) multi_hit[2];
		this_id.id_sharing = multi_hit[3];
		yid.push_back(this_id);

	}

	return yid;
}

void bst_HitProcess::initWithRunNumber(int runno)
{
	if(bstc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		bstc = initializeBSTConstants(runno);
		bstc.runNo = runno;
	}
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


	vector<identifier> fullID;

	identifier superLayerID;
	superLayerID.id = 1;
	superLayerID.name = "svtNoise";

	identifier regionID;
	regionID.id = 1;
	regionID.name = "svtNoise";

	identifier sectorID;
	sectorID.id = 1;
	sectorID.name = "svtNoise";

	identifier sensorID;
	sensorID.id = 1;
	sensorID.name = "svtNoise";

	identifier stripID;
	stripID.id = 1;
	stripID.name = "svtNoise";

	fullID.push_back(superLayerID);
	fullID.push_back(regionID);
	fullID.push_back(sectorID);
	fullID.push_back(sensorID);
	fullID.push_back(stripID);

	double energy = 0.5;
	double time = 5.5;
	int pid = 123;

	MHit* thisNoiseHit = new MHit(energy, time, fullID, pid);
	noiseHits.push_back(thisNoiseHit);



	return noiseHits;
}


map< string, vector <int> >  bst_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;

	return MH;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > bst_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double bst_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}


// this static function will be loaded first thing by the executable
bstConstants bst_HitProcess::bstc = initializeBSTConstants(-1);










