// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <math.h>
#include <random>

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "ahdc_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// V. Sergeyeva, started on 29 May 2020


// this method is for connection to calibration database and extraction of calibration parameters
static ahdcConstants initializeAHDCConstants(int runno, string digiVariation = "default") {
	ahdcConstants atc;
	
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
map<string, double> ahdc_HitProcess::integrateDgt(MHit* aHit, int hitn) {
	
	// digitized output
	map<string, double> dgtz;
	// hit ids
	vector<identifier> identity = aHit->GetId();

	// From here the implementation of what we consider as a hit
	// And dgtz variables calculation algorithms

	//int SECI  = identity[0].id - 1; // sector id
	//int SLI   = identity[1].id - 1; // superlayer id
	//int LAY   = identity[2].id - 1; // layer id
	//int nwire = identity[3].id;     // wire id

//	int superlayer;
//	int layer;
//	int wireNum;
//	double doca;
	double adc;
	double time;
	
	

	if(aHit->isBackgroundHit == 1) {

		vector<double>        stepTime    = aHit->GetTime();
		cout << " This is a background hit with time " << stepTime[0] << endl;
		dgtz["superlayer"]     = 0;
		dgtz["layer"]      = 0;
		dgtz["wireNum"]       = 0;
		dgtz["time"]        = stepTime[0];
		dgtz["hitn"]       = hitn;

		if(filterDummyBanks == false) {
			dgtz["adc"]       = 0;
		}
		return dgtz;
	}

	
	// true information
	
	trueInfos tInfos(aHit);

	vector<int>           stepTrackId = aHit->GetTIds();
	vector<double>        stepTime    = aHit->GetTime();
	vector<double>        mgnf        = aHit->GetMgnf();
	// energy at each step
	// for example tInfos.eTot is total energy deposited
	// tInfos.eTot is the sum of all steps s of Edep[s]
	vector<G4double>      Edep        = aHit->GetEdep();
	vector<G4ThreeVector> pos         = aHit->GetPos();
	// local variable for each step
	vector<G4ThreeVector> Lpos        = aHit->GetLPos();
	// take momentum for each step
	vector<G4ThreeVector> mom         = aHit->GetMoms();
	vector<double>        E           = aHit->GetEs();

	unsigned nsteps = Edep.size();
	
	cout << " AHDC hitprocess: number of steps in a hit: " << nsteps << endl;
	

	//
	// Get the information x,y,z and Edep at each ionization point
	//

	// -------------------------- TIME SHIFT for non-primary tracks ---------------------------
	//if(aHit->GetTId() != timeShift_map.cend()->first){
	//if(aHit->GetTId()< 3) timeShift_map[aHit->GetTId()] = 0.0;
	//else timeShift_map[aHit->GetTId()] = G4RandFlat::shoot(-8000.,8000.);
	//}

	double LposX=0.;
	double LposY=0.;
	double LposZ=0.;
	double DiffEdep=0.;

//	double shift_t = 0.0;


	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		LposX = Lpos[s].x();
		LposY = Lpos[s].y();
		LposZ = Lpos[s].z();

		DiffEdep = Edep[s];

		//z_cm = LposZ/10.0;

		/*
		 // rtpc example has code for drift length and drift time calculation
		 // all in cm
		 a_t=a_t1*(pow(z_cm,4)) + a_t2*(pow(z_cm,3)) + a_t3*(pow(z_cm,2)) + a_t4*z_cm + a_t5;
		 b_t=b_t1*(pow(z_cm,4)) + b_t2*(pow(z_cm,3)) + b_t3*(pow(z_cm,2)) + b_t4*z_cm + b_t5;
		 ...
		 */
		// For now, we use simple variables, without realistic formulae

		// time shift
		//shift_t = timeShift_map.find(aHit->GetTId())->second;
		// NO time shift
		//shift_t = 0.0;

		//tdc=t_s2pad+t_gap+shift_t;
		time = stepTime[s]++;

		//adc=DiffEdep;
		adc = DiffEdep++;
	}

	// Here are the dgtz varibles that we want to calculate using MC true info of a hit
	// They are visible in the gemc simulation output: integrated digitized bank (2302,0)
	dgtz["superlayer"] = identity[0].id;	//identity[0].id; 	//(2302,1)
	dgtz["layer"] = identity[1].id;		//identity[2].id;		//(2302,2)
	dgtz["wireNum"] = identity[2].id;	//identity[3].id;	//(2302,3)
	dgtz["adc"]    = adc;		//(2302,4)
	dgtz["time"]   = time;		//(2302,5)
	dgtz["hitn"] = hitn;		//(2302,99)

	cout << " start of the AHDC hit " << endl;
	cout << " value in identity[0].id = superlayer var: " << identity[0].id << endl;
	cout << " value in identity[1].id = layer var: " << identity[1].id << endl;
	cout << " value in identity[2].id = wireNum var: " << identity[2].id << endl;
	cout << " value in identity[3].id var: " << identity[3].id << endl;
	cout << " value in adc var: " << adc << endl;
	cout << " value in time var: " << time << endl;
	cout << " value in hitn var: " << hitn << endl;
	cout << " end of the AHDC hit " << endl;

	//cout << " value in superlayer var: " << identity[0].id << endl;
	//cout << " value in layer var: " << identity[1].id << endl;
	//cout << " value in wireNum var: " << identity[2].id << endl;


	
	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}



// this method is to locate the hit event, it returns a hitted wire or paddle; this is also the one that needs to be implemented at first.
vector<identifier> ahdc_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector) {

	//id[id.size()-1].id_sharing = 1;
	//return id;

	
	vector<identifier> yid = id;
	
	int nwire = 13;

	/*
	 G4StepPoint   *prestep   = aStep->GetPreStepPoint();
	 G4StepPoint   *poststep  = aStep->GetPostStepPoint();
	 G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction
	 G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
	 ->GetTopTransform().TransformPoint(xyz);

	 double ylength = Detector.dimensions[3];  ///< G4Trap Semilength
	 double deltay  = 0.9;
	 double loc_y   = Lxyz.y() + ylength;    ///< Distance from bottom of G4Trap. ministaggger does not affect it since the field/guardwires are fixed.

	 int nwire = (int) floor(loc_y/deltay);

	 // resetting nwire for extreme cases
	 if(nwire <= 0 )  nwire = 1;
	 if(nwire >= 31) nwire = 30;
	 */

	// setting wire number
	yid[3].id = nwire;

	// checking that the next wire is not the one closer
	//if(fabs( (nwire+1)*deltay - loc_y ) < fabs( nwire*deltay - loc_y ) && nwire != 112 )
	//yid[3].id = nwire + 1;

	// all energy to this wire (no energy sharing)
	yid[3].id_sharing = 1;

	return yid;
	
}

// - electronicNoise: returns a vector of hits generated / by electronics.
// additional method, can be implemented later
vector<MHit*> ahdc_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;
	
	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);
	
	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)
	
	return noiseHits;
}

map< string, vector <int> > ahdc_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map< string, vector <int> > MH;
	
	return MH;
}

// - charge: returns charge/time digitized information / step
// this method is implemented in ftof, but information from this bank is not translated into the root format right now (29/05/2020)
// the output is only visible in .txt output of gemc simulation + <option name="SIGNALVT" value="ftof"/> into gcard
map< int, vector <double> > ahdc_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)

double ahdc_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

void ahdc_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;

	if (atc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		atc = initializeAHDCConstants(runno, digiVariation);
		atc.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
ahdcConstants ahdc_HitProcess::atc = initializeAHDCConstants(-1);




