// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// gemc headers
#include "rich_hitprocess.h"

static richConstants initializeRICHConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	richConstants richc;
	if(runno == -1) return richc;
	
	
	// no quantities read yet.
	
	//	// database
	//	richc.runNo = runno;
	//	richc.date       = "2016-03-15";
	//	if(getenv ("CCDB_CONNECTION") != NULL)
	//		richc.connection = (string) getenv("CCDB_CONNECTION");
	//	else
	//		richc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	//
	//	richc.variation  = "main";
	//	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(richc.connection));
	
	return richc;
}



map<string, double> rich_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;
	
	dgtz["hitn"]   = hitn;
	
	return dgtz;
}


// this routine needs to be modified
// no point drawing should be made here, but in MHit
// finding the PMT should be in a different function,
// with parameters coming from DB

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"

vector<identifier> rich_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

map< string, vector <int> >  rich_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


void rich_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;

	if(richc.runNo != runno) {
//		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		richc = initializeRICHConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		richc.runNo = runno;
	}
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> rich_HitProcess :: electronicNoise()
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



// - charge: returns charge/time digitized information / step
map< int, vector <double> > rich_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The input are charge value, time
double rich_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}



// this static function will be loaded first thing by the executable
richConstants rich_HitProcess::richc = initializeRICHConstants(-1);





