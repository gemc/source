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
	if (getenv("CCDB_CONNECTION") != nullptr)
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

	//int superlayer;
        //int layer;
	//int wire;
	double doca = 100.0;
	//double adc;
	double time;
	
	

	if(aHit->isBackgroundHit == 1) {

		vector<double>        stepTime    = aHit->GetTime();
			cout << " This is a background hit with time " << stepTime[0] << endl;
		dgtz["superlayer"]     = 0;
		dgtz["layer"]      = 0;
		dgtz["wire"]       = 0;
		dgtz["time"]        = stepTime[0];
		dgtz["hitn"]       = hitn;

		if(filterDummyBanks == false) {
			dgtz["doca"]       = 0;
			dgtz["energy"] = 0;
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
	
	double signal_t = 0.0;
	double signal_tTimesEdep = 0.0;
	
	cout << " AHDC hitprocess: number of steps in a hit: " << nsteps << endl;
    
	double LposX=0.;
	double LposY=0.;
	double LposZ=0.;

	double driftVelocity = 0.026;  // mm/ns // drift velocity is 26 um/ns, taken from DC
	
	//vector<double> CellVertex;
	double CellVertex = 0.0;
	
	double xV0 = 0.0;
	double yV0 = 0.0;
	double xV3 = 0.0;
	double yV3 = 0.0;
	double xV4 = 0.0;
	double yV4 = 0.0;	
	double xV7 = 0.0;
	double yV7 = 0.0;

	double Y_sigwire_top = 0.0; // z=-150 mm
	double X_sigwire_top = 0.0;
	double Z_sigwire_top = -150.0;// Global coordiantes = -22.3 mm // Local coordinates = -150.0 mm!
	double Y_sigwire_bot = 0.0; // z=+150 mm
	double X_sigwire_bot = 0.0;
	double Z_sigwire_bot = 150.0;// Global coordiantes = 277.7 mm // Local coordinates = 150.0 mm!
	
	double dim_id_1, dim_id_2, dim_id_7, dim_id_8;

	int subcell = 0; // subcell value = 1 or 2, for one same cell, 1 is to the right of the signal wire, 2 is to the left of the signal wire.
	// Vertices #0 and #3 are always the first/last ones to define the top face of G4 generic trapezoide
	// Vertices #6 and #9 are always the first/last ones to define the bottom face of G4 generic trapezoide

	double SigWireRadiusTop = 0.0; 
	double SigWireRadiusBottom = 0.0;
	int wirealigned = 0;
	double SigWireAlphaTop = 0.0; 
	double SigWireAlphaBottom = 0.0;
	double SigWireDeltaAlpha = 0.0;
	
	// this is to find the distance from the hit to the signal wire and that is perpendicular to the wire! mm!
	double L_ab, L_ah, L_bh, H_abh; // for a triangle abh with: a = top sig. wire coordinates, b = bottom sig. wire coordinates, h = hit coordinates (mm);
	
	// this is for additional calculations for check
	double HitRadius, AlphaHit;

	// this is for energy deposit calculation
	double E_wire = 0.0;
	double E_tot_wire = 0.0;
	double attenuation = 1.0; // mm!!!
	double adc_energy = 0.0;
	double EYld = 1.0;
	
	for(int p=0; p<17; p++)
	{
		CellVertex = aHit->GetDetector().dimensions[p]; // G4 Generic Trapezoide dimensions G4 Generic Trapezoide dimensions
		cout << "Hitted cell parameter" << p << " value = " << CellVertex << endl;
	}
	
	dim_id_1 = aHit->GetDetector().dimensions[1];
	dim_id_2 = aHit->GetDetector().dimensions[2];

	dim_id_7 = aHit->GetDetector().dimensions[7];
	dim_id_8 = aHit->GetDetector().dimensions[8];

		yV3 = aHit->GetDetector().dimensions[8];
		xV3 = aHit->GetDetector().dimensions[7];
		yV0 = aHit->GetDetector().dimensions[2];
		xV0 = aHit->GetDetector().dimensions[1];

		yV7 = aHit->GetDetector().dimensions[16];
		xV7 = aHit->GetDetector().dimensions[15];
		yV4 = aHit->GetDetector().dimensions[10];
		xV4 = aHit->GetDetector().dimensions[9];
		
		if ( abs(dim_id_2) > abs(dim_id_8))
		{
 			subcell = 1;

			X_sigwire_top = xV3 + (xV0 - xV3)/2;
			Y_sigwire_top = yV3 + (yV0 - yV3)/2; // z=-150 mm
			X_sigwire_bot = xV7 + (xV4 - xV7)/2;
			Y_sigwire_bot = yV7 + (yV4 - yV7)/2; // z=+150 mm
						
			
		}
		else
		{
			subcell = 2;

			X_sigwire_top = xV0 + (xV3 - xV0)/2;
			Y_sigwire_top = yV0 + (yV3 - yV0)/2; // z=-150 mm
			X_sigwire_bot = xV4 + (xV7 - xV4)/2;
			Y_sigwire_bot = yV4 + (yV7 - yV4)/2; // z=+150 mm
		}
	
	cout << " shared side defined by points: (" << dim_id_1 << "," <<dim_id_2 << ") and (" << dim_id_7 << "," <<dim_id_8 << ") for the 1rst face; " << endl;
	cout << " and  by points: (" << xV4 << "," << yV4 << ") and (" << xV7 << "," << yV7 << ") for the 2nd face; all in mm!!!" << endl;
 	cout << "sub-cell number = " << subcell << "; hitted cell wire coord.: top face = (" << X_sigwire_top << ", " << Y_sigwire_top << "); bottom face = (" << X_sigwire_bot << ", " << Y_sigwire_bot << ");" << endl;
	
	SigWireRadiusTop = sqrt(pow(X_sigwire_top,2) + pow(Y_sigwire_top,2)); // mm!
	SigWireRadiusBottom = sqrt(pow(X_sigwire_bot,2) + pow(Y_sigwire_bot,2)); // mm!
	SigWireAlphaTop = acos(X_sigwire_top/SigWireRadiusTop); // * 180.0 / PI to pass into degrees, acos() gives the value in radians:
	SigWireAlphaBottom = acos(X_sigwire_bot/SigWireRadiusBottom); // * 180.0 / PI;
	SigWireDeltaAlpha = SigWireAlphaTop - SigWireAlphaBottom;
	
	// Plan ZX for a cell, signal wire orientation
	
	
	if ( SigWireRadiusTop == SigWireRadiusBottom )
	{
		wirealigned = 1;
	} 

	cout << " Signal wire radius (top face) = " << SigWireRadiusTop << "; Signal wire radius (bottom face) = " << SigWireRadiusBottom << "; sig. wire aligned = " << wirealigned << endl;
	cout << " Signal wire alpha (top face) degrees = " << (SigWireAlphaTop * 180.0 / PI)  << "; Signal wire alpha (bottom face) degrees = " << (SigWireAlphaBottom * 180.0 / PI) << "; sig. wire Delta Alpha (top - bot.) = " << (SigWireDeltaAlpha * 180.0 / PI) << endl;

	cout << " ************** Hit started! **************** " << endl;
	cout << "First loop on steps begins" << endl;

	L_ab = sqrt( pow((X_sigwire_bot - X_sigwire_top),2) + pow((Y_sigwire_bot - Y_sigwire_top),2) + pow((Z_sigwire_bot - Z_sigwire_top),2) );
	int local = 0;
	    	for(unsigned int s=0; s<tInfos.nsteps; s++)
		{
			LposX = Lpos[s].x();
			LposY = Lpos[s].y();
			LposZ = Lpos[s].z();

			if (LposZ>150.0) cout << " ######### Global coordinate is seen in Z!!! ########## " << endl;
			if (LposZ<-22.3) 
			{
				cout << " ######### Local coordinate is seen in Z!!! ########## " << endl;
				local = local + 1;
			}

			HitRadius = sqrt(pow(LposX,2) + pow(LposY,2)); // mm!
			AlphaHit = acos(LposX/HitRadius); // radians! and for a precize Z_hit = LposZ (mm)

			// Calculation of distance from hit to signal wire and perpendicular to the wire!
			L_ah = sqrt( pow((X_sigwire_top - LposX),2) + pow((Y_sigwire_top - LposY),2) + pow((Z_sigwire_top - LposZ),2) );
			L_bh = sqrt( pow((X_sigwire_bot - LposX),2) + pow((Y_sigwire_bot - LposY),2) + pow((Z_sigwire_bot - LposZ),2) );
			H_abh = L_ah * sqrt( 1 - pow(((L_ah*L_ah + L_ab*L_ab - L_bh*L_bh)/(2*L_ab*L_ah)),2) );
			
			// variables check for doca calculation
			cout << "Hit (X,Y,Z) location (mm) (" << LposX << ", " << LposY << ", " <<  LposZ << ")" << endl;
			cout << "Signal wire length (mm) = " << L_ab << endl;
			cout << "distance hit->signal wire &perpendicular to wire (mm) = " << H_abh << endl;
			cout << "Hit radius (mm) and alpha (deg.) = (" << HitRadius << ", " << (AlphaHit) * 180.0 / PI << ")" << endl;
			
			if ( H_abh <= doca )
			{
				doca = H_abh; // mm!!!
			}

			// energy deposit calculation
			E_wire = Edep[s] *exp(-H_abh/attenuation);
			E_tot_wire = E_tot_wire + E_wire;	

			// time calculation
			signal_t = stepTime[s] + (H_abh/driftVelocity);
			cout << "signal_t: " << signal_t << ", stepTime: " << stepTime[s] << endl;	
			signal_tTimesEdep = signal_tTimesEdep + (stepTime[s] + H_abh/driftVelocity) * E_wire;
		
			//time = stepTime[s]++;
            		//adc = Edep[s]++;	
			
		}
	cout << "First loop on steps ends" << endl;

			// energy adc value
			adc_energy = E_tot_wire * EYld;
			
			// Just to test, time is linear with doca
			double a = 3.5;
			double b = 4.0;
			time = a*doca+b + signal_tTimesEdep/E_tot_wire;
			
		
			// Here are the dgtz varibles that we want to calculate using MC true info of a hit
			// They are visible in the gemc simulation output: integrated digitized bank (2302,0)
			dgtz["superlayer"] = identity[0].id;	//(2302,1)
       			dgtz["layer"] = identity[1].id;		//(2302,2)
        		dgtz["wire"] = identity[2].id;	//(2302,3)
        		dgtz["doca"]    = doca;		//(2302,4)
			dgtz["energy"]    = adc_energy;		//(2302,5)
			dgtz["time"]   = time;		//(2302,6)	
			dgtz["hitn"] = hitn;		//(2302,99)
		
		cout << " start of the AHDC hit " << endl;
		cout << " value in identity[0].id = superlayer var: " << identity[0].id << endl;	
		cout << " value in identity[1].id = layer var: " << identity[1].id << endl;
		cout << " value in identity[2].id = wire var: " << identity[2].id << endl;
		cout << " value in identity[3].id var: " << identity[3].id << endl;
		cout << " doca value = dist. hit->sig.wire & perpendicular to signal wire) (mm!!!): " << doca << endl;
		cout << " value in doca output: " << doca << endl;
		cout << " value in energy output: " << adc_energy << endl;
		cout << " value in time var: " << time << endl;
		cout << " value in hitn var: " << hitn << endl;
		cout << " if local <> 0 than hit in local reference; local = " << local << endl;
		cout << " ************** Hit ended! **************** " << endl;
		
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




