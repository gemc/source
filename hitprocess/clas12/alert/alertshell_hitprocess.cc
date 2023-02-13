// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <math.h>
#include <random>

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>

#define PI 3.1415926535

using namespace ccdb;

// gemc headers
#include "alertshell_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// V. Sergeyeva, started on 13 Oct. 2020


// this method is for connection to calibration database and extraction of calibration parameters
static alertshellConstants initializeALERTSHELLConstants(int runno, string digiVariation = "default") {
	alertshellConstants atc;
	
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
map<string, double> alertshell_HitProcess::integrateDgt(MHit* aHit, int hitn) {
	
	// digitized output
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	// From here the implementation of what we consider as a hit
	// And dgtz variables calculation algorithms
	
	//int superlayer;
	//int layer;
	//int wire;
	//double doca = 100.0;
	//double adc;
	//double time;
	
	
	
	if(aHit->isBackgroundHit == 1) {
		
		vector<double>        stepTime    = aHit->GetTime();
		cout << " This is a background hit with time " << stepTime[0] << endl;
		//dgtz["superlayer"]     = 0;
		//dgtz["layer"]      = 0;
		//dgtz["wire"]       = 0;
		//dgtz["time"]        = stepTime[0];
		dgtz["hitn"]       = hitn;
		
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
	
	// unsigned nsteps = Edep.size();
	
	// 	double signal_t = 0.0;
	// double signal_tTimesEdep = 0.0;
	
	//cout << " alertshell hitprocess: number of steps in a hit: " << nsteps << endl;
	
	//	double LposX=0.;
	//	double LposY=0.;
	//	double LposZ=0.;
	
	// this is for energy deposit calculation
	//	double E_wire = 0.0;
	//	double E_tot_wire = 0.0;
	// double attenuation = 10.0; // mm!!!
	// double adc_energy = 0.0;
	// double EYld = 10.0;
	
	// double totEdepMC = 0.0;
	
	//	for(unsigned int s=0; s<tInfos.nsteps; s++) {
	//		LposX = Lpos[s].x();
	//		LposY = Lpos[s].y();
	//		LposZ = Lpos[s].z();
	//
	//		totEdepMC = totEdepMC+Edep[s];
	//
	//		// time calculation
	//		signal_t = stepTime[s];
	//		cout << "signal_t: " << signal_t << endl;
	//	}
	//	cout << "First loop on steps ends" << endl;
	
	/*
	 dgtz["superlayer"] = identity[0].id;	//(2302,1)
	 dgtz["layer"] = identity[1].id;		//(2302,2)
	 dgtz["wire"] = identity[2].id;	//(2302,3)
	 dgtz["doca"]    = doca;		//(2302,4)
	 dgtz["subcell"]    = subcell;	// subcell 1 is on the right of the signal wire, subcell 2 is on the left of the signal wire
	 dgtz["adc_energy"]    = adc_energy;		//(2302,5)
	 dgtz["wire_energy"]    = E_tot_wire;		//(2302,5)
	 dgtz["totEdep_MC"]    = totEdepMC;
	 dgtz["signal"]   = signal;
	 dgtz["time"]   = time;		//(2302,6)
	 */
	dgtz["hitn"] = hitn;		//(2302,99)
	
	//	cout << " start of the alertshell hit " << endl;
	//	cout << " value in MC totEdep: " << totEdepMC << endl;
	//	cout << " value in signal var: " << signal_t << endl;
	//	cout << " value in hitn var: " << hitn << endl;
	//	cout << " ************** Hit ended! **************** " << endl;
	
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}



// this method is to locate the hit event, it returns a hitted wire or paddle; this is also the one that needs to be implemented at first.
vector<identifier> alertshell_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector) {
	
	//id[id.size()-1].id_sharing = 1;
	//return id;
	
	
	vector<identifier> yid = id;
	
	//int nwire = 13;
	
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
	//yid[3].id = nwire;
	
	// checking that the next wire is not the one closer
	//if(fabs( (nwire+1)*deltay - loc_y ) < fabs( nwire*deltay - loc_y ) && nwire != 112 )
	//yid[3].id = nwire + 1;
	
	// all energy to this wire (no energy sharing)
	//yid[3].id_sharing = 1;
	
	return yid;
	
}

// - electronicNoise: returns a vector of hits generated / by electronics.
// additional method, can be implemented later
vector<MHit*> alertshell_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;
	
	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);
	
	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)
	
	return noiseHits;
}

map< string, vector <int> > alertshell_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map< string, vector <int> > MH;
	
	return MH;
}

// - charge: returns charge/time digitized information / step
// this method is implemented in ftof, but information from this bank is not translated into the root format right now (29/05/2020)
// the output is only visible in .txt output of gemc simulation + <option name="SIGNALVT" value="ftof"/> into gcard
map< int, vector <double> > alertshell_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)

double alertshell_HitProcess::voltage(double charge, double time, double forTime) {
	return 0.0;
}

void alertshell_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	
	if (atc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		atc = initializeALERTSHELLConstants(runno, digiVariation);
		atc.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
alertshellConstants alertshell_HitProcess::atc = initializeALERTSHELLConstants(-1);




