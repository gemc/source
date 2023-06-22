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

        string timestamp = "";
        if(digiSnapshotTime != "no") {
                timestamp = ":"+digiSnapshotTime;
        }
	
	// database
	richc.runNo = runno;
	
	//	richc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != nullptr)
	  richc.connection = (string) getenv("CCDB_CONNECTION");
	else
	  richc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
        richc.variation  = "main";
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(richc.connection));

	vector<vector<double>> data;
	
	// Eventually: not limited to module1

	// read timewalk correction
	snprintf(richc.database, sizeof(richc.database), "/calibration/rich/module1/time_walk:%d:%s%s", richc.runNo, digiVariation.c_str(), timestamp.c_str());
	calib->GetCalib(data,richc.database);

	for(int row = 0; row<data.size(); row++){	  
	  int ipmt = data[row][1];
	  richc.timewalkCorr_D0[ipmt-1] = data[row][3];
	  richc.timewalkCorr_m1[ipmt-1]	= data[row][4];
	  richc.timewalkCorr_m2[ipmt-1]	= data[row][5];
	  richc.timewalkCorr_T0[ipmt-1]	= data[row][6];	  
	}	

        // read time offset
        snprintf(richc.database, sizeof(richc.database), "/calibration/rich/module1/time_offset:%d:%s%s", richc.runNo, digiVariation.c_str(), timestamp.c_str());
        calib->GetCalib(data,richc.database);

        for(int row = 0; row<data.size(); row++){
          int ipmt = data[row][1];
          richc.timeOffsetCorr[ipmt-1] = data[row][3];
        }

	// TODO: set dead pixels? hot pixels?
	return richc;
}


// digitized info integrated over hit
// should this match what was in rich_sector4/rich__bank.txt?
// or should it roughly match RICH::tdc from data.json?

map<string, double> rich_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	//cout << "hit number " << hitn << endl;
        vector<identifier> identity = aHit->GetId();
        int idsector = identity[0].id;
        int idpmt = identity[1].id;
        int idpixel = 0;//identity[2].id;
	
	rejectHitConditions = false;
	writeHit = true;

	// setting all values generically just to see how they get printed out
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = idsector;
	dgtz["layer"] = idpmt;
	dgtz["component"] = idpixel; // TODO: function for local position -> pixel #

	dgtz["TDC_TDC"]=0;
	dgtz["TDC_order"]=2;
	
	writeHit = true;	
	
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

// not really used in any others?
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
	return 1.0;
}



// this static function will be loaded first thing by the executable
richConstants rich_HitProcess::richc = initializeRICHConstants(-1);





