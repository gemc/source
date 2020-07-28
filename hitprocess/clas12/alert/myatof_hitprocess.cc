// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "myatof_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

static myatofConstants initializeMYATOFConstants(int runno, string digiVariation = "default") {
	myatofConstants atc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if (runno == -1) return atc;
	
	atc.runNo = runno;
	atc.date = "2020-04-20";
	if (getenv("CCDB_CONNECTION") != NULL)
		atc.connection = (string) getenv("CCDB_CONNECTION");
	else
		atc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	return atc;
}


// this methos is for implementation of digitized outputs and the first one that needs to be implemented.
map<string, double> myatof_HitProcess::integrateDgt(MHit* aHit, int hitn) {
	
	// digitized output
	map<string, double> dgtz;
	// hit ids
	vector<identifier> identity = aHit->GetId();

//	int sector;
//	int zregion;
//	int paddle;
	double adc;
	double time;

	if(aHit->isBackgroundHit == 1) {

		vector<double>        stepTime    = aHit->GetTime();
		cout << " This is a background hit with time " << stepTime[0] << endl;
		dgtz["sector"]     = 0;
		dgtz["superlayer"]      = 0;
		dgtz["layer"]      = 0;
		dgtz["paddle"]       = 0;
		dgtz["time"]        = stepTime[0];
		dgtz["hitn"]       = hitn;

		if(filterDummyBanks == false) {
			dgtz["adc"]       = 0;
		}
		return dgtz;
	}

	adc = 999.999;
	time = 888.888;

	dgtz["sector"] = identity[0].id;
	dgtz["superlayer"] = identity[1].id;
	dgtz["layer"] = identity[2].id;
	dgtz["paddle"] = identity[3].id;
	dgtz["adc"] = adc;
	dgtz["time"] = time;
	dgtz["hitn"] = hitn;		//(2202,99)

	cout << " start of the ATOF hit " << endl;
	cout << " value in identity[0].id = sector var: " << identity[0].id << endl;
	cout << " value in identity[1].id = superlayer var: " << identity[1].id << endl;
	cout << " value in identity[2].id = layer var: " << identity[2].id << endl;
	cout << " value in identity[3].id = paddle var: " << identity[3].id << endl;
	//cout << " value in identity[3].id var: " << identity[3].id << endl;
	cout << " value in adc var: " << adc << endl;
	cout << " value in time var: " << time << endl;
	cout << " value in hitn var: " << hitn << endl;
	cout << " end of the ATOF hit " << endl;



	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}

// this method is to locate the hit event, it returns a hitted wire or paddle; this is also the one that needs to be implemented at first.
vector<identifier> myatof_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector) {
	
	id[id.size()-1].id_sharing = 1;
	return id;


}

// - electronicNoise: returns a vector of hits generated / by electronics.

vector<MHit*> myatof_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;
	
	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);
	
	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)
	
	return noiseHits;
}

map< string, vector <int> > myatof_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map< string, vector <int> > MH;
	
	return MH;
}

// - charge: returns charge/time digitized information / step

map< int, vector <double> > myatof_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)

double myatof_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

void myatof_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;

	if (atc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		atc = initializeMYATOFConstants(runno, digiVariation);
		atc.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
myatofConstants myatof_HitProcess::atc = initializeMYATOFConstants(-1);




