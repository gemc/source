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
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	int sector    = 0;
	int layer     = 10 * identity[0].id + identity[1].id ; // 10*superlayer + layer
	int component = identity[2].id;

	if(aHit->isBackgroundHit == 1) {

		double totEdep  = aHit->GetEdep()[0];
		double stepTime = aHit->GetTime()[0];
		double tdc      = stepTime;

		dgtz["hitn"]      = hitn;
		dgtz["sector"]    = sector;
		dgtz["layer"]     = layer;
		dgtz["component"] = component;
		dgtz["ADC_order"] = 0;
		dgtz["ADC_ADC"]   = (int) totEdep;
		dgtz["ADC_time"]  = tdc;
		dgtz["ADC_ped"]   = 0;
		
		dgtz["TDC_order"] = 0;
		dgtz["TDC_TDC"]   = tdc;

		return dgtz;

	}
	
	double doca = 100.0;

	// true information
	
	trueInfos tInfos(aHit);
	
	vector<int>           stepTrackId = aHit->GetTIds();
	vector<double>        stepTime    = aHit->GetTime();
	vector<double>        mgnf        = aHit->GetMgnf();
	vector<G4double>      Edep        = aHit->GetEdep();
	vector<G4ThreeVector> pos         = aHit->GetPos(); // local position variable for each step
	vector<G4ThreeVector> Lpos        = aHit->GetLPos();
	vector<G4ThreeVector> mom         = aHit->GetMoms();
	vector<double>        E           = aHit->GetEs();
	
	
//	double signal_t = 0.0;
	double signal_tTimesEdep = 0.0;
	
	
	double LposX=0.;
	double LposY=0.;
	double LposZ=0.;
	
	double driftVelocity = 0.026;  // mm/ns // drift velocity is 26 um/ns, taken from DC
	
	//vector<double> CellVertex;
	// double CellVertex = 0.0;
	
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
	
	double dim_id_2, dim_id_8;
	// int subcell = 0; // subcell value = 1 or 2, for one same cell, 1 is to the right of the signal wire, 2 is to the left of the signal wire.

	// Vertices #0 and #3 are always the first/last ones to define the top face of G4 generic trapezoide
	// Vertices #6 and #9 are always the first/last ones to define the bottom face of G4 generic trapezoide
	
	//	double SigWireRadiusTop = 0.0;
	//	double SigWireRadiusBottom = 0.0;
	// int wirealigned = 0;
	// double SigWireAlphaTop = 0.0;
	// double SigWireAlphaBottom = 0.0;
	// double SigWireDeltaAlpha = 0.0;
	
	// this is to find the distance from the hit to the signal wire and that is perpendicular to the wire! mm!
	double L_ab, L_ah, L_bh, H_abh; // for a triangle abh with: a = top sig. wire coordinates, b = bottom sig. wire coordinates, h = hit coordinates (mm);
	
	// this is for additional calculations for check
	// double HitRadius;
	//, AlphaHit;
	
	// this is for energy deposit calculation
	double E_wire = 0.0;
	double E_tot_wire = 0.0;
	double attenuation = 10.0; // mm!!!
	double adc_energy = 0.0;
	double EYld = 10.0;
	
	double totEdepMC = 0.0;
	
//	for(int p=0; p<17; p++)
//	{
//		CellVertex = aHit->GetDetector().dimensions[p]; // G4 Generic Trapezoide dimensions G4 Generic Trapezoide dimensions
//		cout << "Hitted cell parameter" << p << " value = " << CellVertex << endl;
//	}
	
	//	dim_id_1 = aHit->GetDetector().dimensions[1];
	
	 dim_id_2 = aHit->GetDetector().dimensions[2];
	 dim_id_8 = aHit->GetDetector().dimensions[8];
	
	yV3 = aHit->GetDetector().dimensions[8];
	xV3 = aHit->GetDetector().dimensions[7];
	yV0 = aHit->GetDetector().dimensions[2];
	xV0 = aHit->GetDetector().dimensions[1];
	
	yV7 = aHit->GetDetector().dimensions[16];
	xV7 = aHit->GetDetector().dimensions[15];
	yV4 = aHit->GetDetector().dimensions[10];
	xV4 = aHit->GetDetector().dimensions[9];


	if ( abs(dim_id_2) > abs(dim_id_8)) {
		// subcell = 1;

		X_sigwire_top = xV3 + (xV0 - xV3)/2;
		Y_sigwire_top = yV3 + (yV0 - yV3)/2; // z=-150 mm
		X_sigwire_bot = xV7 + (xV4 - xV7)/2;
		Y_sigwire_bot = yV7 + (yV4 - yV7)/2; // z=+150 mm

	} else {
		// subcell = 2;

		X_sigwire_top = xV0 + (xV3 - xV0)/2;
		Y_sigwire_top = yV0 + (yV3 - yV0)/2; // z=-150 mm
		X_sigwire_bot = xV4 + (xV7 - xV4)/2;
		Y_sigwire_bot = yV4 + (yV7 - yV4)/2; // z=+150 mm
	}


	//cout << " shared side defined by points: (" << dim_id_1 << "," <<dim_id_2 << ") and ("  << "," <<dim_id_8 << ") for the 1rst face; " << endl;
	//cout << " and  by points: (" << xV4 << "," << yV4 << ") and (" << xV7 << "," << yV7 << ") for the 2nd face; all in mm!!!" << endl;
	//cout << "sub-cell number = "  << "; hitted cell wire coord.: top face = (" << X_sigwire_top << ", " << Y_sigwire_top << "); bottom face = (" << X_sigwire_bot << ", " << Y_sigwire_bot << ");" << endl;
	
	//	SigWireRadiusTop = sqrt(pow(X_sigwire_top,2) + pow(Y_sigwire_top,2)); // mm!
	//	SigWireRadiusBottom = sqrt(pow(X_sigwire_bot,2) + pow(Y_sigwire_bot,2)); // mm!
	// SigWireAlphaTop = acos(X_sigwire_top/SigWireRadiusTop); // * 180.0 / PI to pass into degrees, acos() gives the value in radians:
	// SigWireAlphaBottom = acos(X_sigwire_bot/SigWireRadiusBottom); // * 180.0 / PI;
	// SigWireDeltaAlpha = SigWireAlphaTop - SigWireAlphaBottom;
	
	// Plan ZX for a cell, signal wire orientation
	
	
	//	if ( SigWireRadiusTop == SigWireRadiusBottom )
	//	{
	//		wirealigned = 1;
	//	}
	
	//cout << " Signal wire radius (top face) = " << SigWireRadiusTop << "; Signal wire radius (bottom face) = " << SigWireRadiusBottom << "; sig. wire aligned = " << wirealigned << endl;
	//cout << " Signal wire alpha (top face) degrees = " << (SigWireAlphaTop * 180.0 / PI)  << "; Signal wire alpha (bottom face) degrees = " << (SigWireAlphaBottom * 180.0 / PI) << "; sig. wire Delta Alpha (top - bot.) = " << (SigWireDeltaAlpha * 180.0 / PI) << endl;
	
	// cout << " ************** Hit started! **************** " << endl;
	//cout << "First loop on steps begins" << endl;
	
	L_ab = sqrt( pow((X_sigwire_bot - X_sigwire_top),2) + pow((Y_sigwire_bot - Y_sigwire_top),2) + pow((Z_sigwire_bot - Z_sigwire_top),2) );
	// int local = 0;
	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		LposX = Lpos[s].x();
		LposY = Lpos[s].y();
		LposZ = Lpos[s].z();
		
//		if (LposZ>150.0) cout << " ######### Global coordinate is seen in Z!!! ########## " << endl;
//		if (LposZ<-22.3)
//		{
//			cout << " ######### Local coordinate is seen in Z!!! ########## " << endl;
//			local = local + 1;
//		}
		
		// HitRadius = sqrt(pow(LposX,2) + pow(LposY,2)); // mm!
		// AlphaHit = acos(LposX/HitRadius); // radians! and for a precize Z_hit = LposZ (mm)
		
		// Calculation of distance from hit to signal wire and perpendicular to the wire!
		L_ah = sqrt( pow((X_sigwire_top - LposX),2) + pow((Y_sigwire_top - LposY),2) + pow((Z_sigwire_top - LposZ), 2) );
		L_bh = sqrt( pow((X_sigwire_bot - LposX),2) + pow((Y_sigwire_bot - LposY),2) + pow((Z_sigwire_bot - LposZ), 2) );
		H_abh = L_ah * sqrt( 1 - pow(((L_ah*L_ah + L_ab*L_ab - L_bh*L_bh)/(2*L_ab*L_ah)), 2) );
		
		// variables check for doca calculation
//		cout << "Hit (X,Y,Z) location (mm) (" << LposX << ", " << LposY << ", " <<  LposZ << ")" << endl;
//		cout << "Signal wire length (mm) = " << L_ab << endl;
//		cout << "distance hit->signal wire &perpendicular to wire (mm) = " << H_abh << endl;
//		cout << "Hit radius (mm) and alpha (deg.) = (" << HitRadius << ", " << (AlphaHit) * 180.0 / PI << ")" << endl;
		
		if ( H_abh <= doca ) {
			doca = H_abh; // mm!!!
		}
		
		// energy deposit calculation
		E_wire = Edep[s] *exp(-H_abh/attenuation);
		E_tot_wire = E_tot_wire + E_wire;
		
		totEdepMC = totEdepMC+Edep[s];
		
		// time calculation
		// signal_t = stepTime[s] + (H_abh/driftVelocity);
		// cout << "signal_t: " << signal_t << ", stepTime: " << stepTime[s] << endl;
            
                // docasig is a fit to sigma vs distance plot. A second order pol used for the fit (p0+p1*x+p2*x*x).
                // drift velocity as a function of distance. pol2 fitted to t vs x plot and drift velocity is derived from the fit (1/(dt/dx)).
                // both sig vs dist and t vs dist plots are taken from  Lucien Causse's PhD thesis ("Development of a stereo drift chamber for the Jefferson Laboratory ALERT Experiment."). 
                // plots were digitized and then fitted to a pol2. 

                double driftP1=-16.17;
                double driftP2=24.81;
                double docasig = 337.3-210.3*doca+34.7*pow(doca,2);
                std::default_random_engine dseed(time(0)); //seed
                std::normal_distribution<double> ddist(doca, docasig); //a resuolution affect is added to doca.
                double doca_r =ddist(dseed);
	        driftVelocity = 1./(driftP1+2.*driftP2*doca_r);  // mm/ns // drift velocity as a function of distance. pol2 fitted to t vs x plot and drift velocity is then extracted from dx/dt, d/dt(p0+p1*x+p2*x^2)=p1+2*p2*x.
		signal_tTimesEdep = signal_tTimesEdep + (stepTime[s] + H_abh/driftVelocity) * E_wire;
		// cout << "signal_tTimesEdep: " << signal_tTimesEdep << endl;
		
		//time = stepTime[s]++;
		//adc = Edep[s]++;
		
	}
	//cout << "First loop on steps ends" << endl;
	
	// energy adc value
	adc_energy = E_tot_wire * EYld;
	
	// Just to test, time is linear with doca
	double a = 5.0;
	double b = 5.0;
	//double time = a*doca+b + signal_tTimesEdep/E_tot_wire;
	double time = signal_tTimesEdep/E_tot_wire;
//	double signal = 0.0;
//	signal = signal_tTimesEdep/E_tot_wire;
	
	
	// Here are the dgtz varibles that we want to calculate using MC true info of a hit
	// They are visible in the gemc simulation output: integrated digitized bank (2302,0)
//	dgtz["superlayer"]  = identity[0].id;	//(2302,1)
//	dgtz["layer"]       = identity[1].id;		//(2302,2)
//	dgtz["wire"]        = identity[2].id;	//(2302,3)
//	dgtz["doca"]        = doca;		//(2302,4)
//	dgtz["adc_energy"]  = adc_energy;		//(2302,5)
//	dgtz["wire_energy"] = E_tot_wire;		//(2302,5)
//	dgtz["totEdep_MC"]  = totEdepMC;
//	dgtz["signal"]      = signal;
//	dgtz["time"]        = time;		//(2302,6)
//	dgtz["hitn"]        = hitn;		//(2302,99)
	
	//cout << " start of the AHDC hit " << endl;
	//cout << " value in identity[0].id = superlayer var: " << identity[0].id << endl;
	//cout << " value in identity[1].id = layer var: " << identity[1].id << endl;
	//cout << " value in identity[2].id = wire var: " << identity[2].id << endl;
	//cout << " value in identity[3].id var: " << identity[3].id << endl;
	//cout << " doca value = dist. hit->sig.wire & perpendicular to signal wire) (mm!!!): " << doca << endl;
	//cout << " value in doca output: " << doca << endl;
	//cout << " value in MC totEdep: " << totEdepMC << endl;
	//cout << " value in wire energy deposit: " << E_tot_wire << endl;
	//cout << " value in adc energy deposit: " << adc_energy << endl;
	//cout << " value in signal var: " << signal << endl;
	//cout << " value in time var: " << time << endl;
	//cout << " value in hitn var: " << hitn << endl;
	//cout << " if local <> 0 than hit in local reference; local = " << local << endl;
	// cout << " ************** Hit ended! **************** " << endl;
	
	//cout << " value in superlayer var: " << identity[0].id << endl;
	//cout << " value in layer var: " << identity[1].id << endl;
	//cout << " value in wireNum var: " << identity[2].id << endl;
	
	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = sector;
	dgtz["layer"]     = layer;
	dgtz["component"] = component;
	dgtz["ADC_order"] = 0;
	dgtz["ADC_ADC"]   = (int) 100000*adc_energy;
	dgtz["ADC_time"]  = time;
	dgtz["ADC_ped"]   = doca*1000;
	
	dgtz["TDC_order"] = 0;
	dgtz["TDC_TDC"]   = time;

	
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}



vector<identifier> ahdc_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector) {

	id[id.size()-1].id_sharing = 1;
	return id;
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




