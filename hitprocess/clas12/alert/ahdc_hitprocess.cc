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

	int sector    = 1;
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
		dgtz["ADC_order"] = 1;
		dgtz["ADC_ADC"]   = (int) totEdep;
		dgtz["ADC_time"]  = tdc;
		dgtz["ADC_ped"]   = 0;
		
		dgtz["TDC_order"] = 0;
		dgtz["TDC_TDC"]   = tdc;

		return dgtz;

	}

	ahdcSignal *Signal = new ahdcSignal(aHit,hitn);
	Signal->SetTmin(0);
	Signal->SetTmax(6000);
	//Signal->SetDelay(1000); // ns
	//Signal->SetSamplingTime(44); // ns
	Signal->SetElectronYield(100000);
	Signal->Digitize();
	std::map<std::string,double> output = Signal->Decode();

	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = sector;
	dgtz["layer"]     = layer;
	dgtz["component"] = component;
	dgtz["ADC_order"] = 1;
	dgtz["ADC_ADC"]   = (int) output["max_value"]; // adc
	dgtz["ADC_time"]  = output["t_ovr"]; // ns
	dgtz["ADC_ped"]   = (int) output["noise_level"]; // adc
	dgtz["ADC_integral"] = (int) output["integral"]; // adc per 44 ns
	dgtz["ADC_timestamp"] = output["t_start"]; // ns
	dgtz["ADC_t_cfd"] = output["t_cfd"]; // ns
	dgtz["ADC_mctime"] = Signal->GetMCTime(); // ns
	dgtz["ADC_nsteps"] = Signal->Get_nsteps();
	dgtz["ADC_mcEtot"] = Signal->GetMCEtot(); // keV

	//dgtz["TDC_order"] = 0;
	//dgtz["TDC_TDC"]   = output["t_start"];
	std::vector<double> SDgtz = Signal->GetDgtz();
	for (int itr=1;itr<=136;itr++){
		std::ostringstream sEntry;
		sEntry << "wf136_s" << itr;
		dgtz[sEntry.str()] = (int) SDgtz.at(itr-1);
	}
	delete Signal;
	
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


// -------------
// ahdcSignal
// -------------

void ahdcSignal::ComputeDocaAndTime(MHit * aHit){
	vector<G4ThreeVector> Lpos        = aHit->GetLPos();
	int nsteps = Lpos.size();
	double LposX, LposY, LposZ;
	
	// ALERT geometry
	double X_sigwire_top = 0; // [mm]
	double Y_sigwire_top = 0;
	double Z_sigwire_top = -150; 
	double X_sigwire_bot = 0; // [mm]
	double Y_sigwire_bot = 0;
	double Z_sigwire_bot = 150;
	
	// Compute Y_sigwire_top, Z_sigwire_top, Y_sigwire_bot, Z_sigwire_bot
	double xV0 = 0.0;
	double yV0 = 0.0;
	double xV3 = 0.0;
	double yV3 = 0.0;
	double xV4 = 0.0;
	double yV4 = 0.0;
	double xV7 = 0.0;
	double yV7 = 0.0;
	double dim_id_2, dim_id_8;
	
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
	}
	else {
		// subcell = 2;
		X_sigwire_top = xV0 + (xV3 - xV0)/2;
		Y_sigwire_top = yV0 + (yV3 - yV0)/2; // z=-150 mm
		X_sigwire_bot = xV4 + (xV7 - xV4)/2;
		Y_sigwire_bot = yV4 + (yV7 - yV4)/2; // z=+150 mm
	}

	// Triangle abh
	// a (sigwire_top), b (sigwire_bot), h (hit position)
	// H_abh is the distance between hit and the wire and perpendicular to the wire
	double L_ab, L_ah, L_bh, H_abh;
	// Compute the distance between top and bottom of the wire
	L_ab = sqrt(pow(X_sigwire_top-X_sigwire_bot,2) + pow(Y_sigwire_top-Y_sigwire_bot,2) + pow(Z_sigwire_top-Z_sigwire_bot,2));
	for (int s=0;s<nsteps;s++) {
		// Load current hit positions
		LposX = Lpos[s].x();
		LposY = Lpos[s].y();
		LposZ = Lpos[s].z();
		// Compute distance
		L_ah = sqrt(pow(X_sigwire_top-LposX,2) + pow(Y_sigwire_top-LposY,2) + pow(Z_sigwire_top-LposZ,2));
		L_bh = sqrt(pow(X_sigwire_bot-LposX,2) + pow(Y_sigwire_bot-LposY,2) + pow(Z_sigwire_bot-LposZ,2));
		// Compute the height of a triangular (see documentation for the demonstration of the formula)
		H_abh = L_ah*sqrt(1 - pow((L_ah*L_ah + L_ab*L_ab - L_bh*L_bh)/(2*L_ah*L_ab),2)); // this is the d.o.c.a of a given hit (!= MHit)
		Doca.push_back(H_abh);
		// Add a resolution on doca
		double docasig = 337.3-210.3*H_abh+34.7*pow(H_abh,2); // um // fit sigma vs distance // Fig 4.14 (right), L. Causse's thesis
		docasig = docasig/1000; // mm
		std::default_random_engine dseed(time(0)); //seed
		std::normal_distribution<double> docadist(H_abh, docasig);
		// Compute time
		double driftTime = 7*H_abh + 7*pow(H_abh,2) + 4*pow(H_abh,3); // fit t vs distance //  Fig 4.12 (right), L. Causse's thesis
		DriftTime.push_back(driftTime);
	}
}

void ahdcSignal::GenerateNoise(double mean, double stdev){
	int Npts = (int) floor( (tmax-tmin)/samplingTime );
	std::random_device rd;      // Create a random device to seed the generator
	std::mt19937 gen(rd());     // Create a random number engine (e.g., Mersenne Twister)
	for (int i=0;i<Npts;i++){
		std::normal_distribution<double> draw(mean,stdev);
		double value = draw(gen);
		if (value < 0) value = 0;
		Noise.push_back(value);
	}
}

void ahdcSignal::Digitize(){
	this->GenerateNoise(300,30);
	int Npts = (int) floor( (tmax-tmin)/samplingTime );
	for (int i=0;i<Npts;i++) {
		double value = this->operator()(tmin + i*samplingTime); //in keV/ns
		value = (int) floor(electronYield*value + Noise.at(i)); //convert in ADC +  noise
		int adc = (value < adc_max) ? value : adc_max; // saturation effect 
		Dgtz.push_back(adc);
	}
}

std::map<std::string,double> ahdcSignal::Decode(){

	double t_start, t_ovr, t_max_value, max_value, integral;
	max_value = Dgtz.at(0);
	integral = 0;
	int Npts = Dgtz.size();
	int i_max = 0;
	
	// compute max_value
	for (int i=0;i<Npts;i++){
		if (max_value < Dgtz.at(i)) {
			max_value = Dgtz.at(i);
			i_max = i; // useful
		}
	}
	if (max_value == adc_max) { // there is a  plateau (saturation)
		int i_max2 = i_max;
		while (i_max2 < Npts-1){
			if (Dgtz.at(i_max2) == adc_max) {
				i_max2++;
			} 
			else {break;}
		}
		i_max = (int) (i_max+i_max2-1)/2;
	}
	else {  // normal case
		// averaging of max_value
		if ((i_max > 2) and (i_max < Npts-2)){
			max_value = 0;
			for (int i=-2;i<=2;i++){ max_value += Dgtz.at(i_max+i);}
			max_value = max_value/5; // done
		}
	}
	t_max_value = i_max*samplingTime; // done
	
	// define noise and threshold
	double noise = 0;
	for (int i=0;i<5;i++){ noise += Noise.at(i);} 
	noise = noise/5; 
	double threshold = (max_value+noise)/2.0;
	
	// compute t_start
	int i_start = 0;
	for (int i=0;i<i_max;i++){
		if (Dgtz.at(i) < threshold) {
			i_start = i; // last pass below threshold and before max_value
		}
	}	// at this stage : i_start < t_start/samplingTime < i_start+1
	int i1 = i_start; // 1 index below 
	int i2 = i_start+1; // 1 index above
	if (i1 < 0) {i1 = 0; } 
	if (i2 >= Npts) {i2 = Npts-1;}
	double slope = (Dgtz.at(i1) - Dgtz.at(i2))/(i1-i2); 
	t_start = tmin + samplingTime*(i1 + (threshold-Dgtz.at(i1))/slope); // done
	
	// compute t_ovr
	int i_ovr = i_max;
	while (i_ovr < Npts-1) {
		if (Dgtz.at(i_ovr) > threshold){
			i_ovr++; // first pass below threshold starting from max_value
		}
		else { break;}
	}      // at this stage : i_ovr-1 < t_ovr/samplingTime < i_ovr
	if (i_ovr < Npts-2) {
		i1 = i_ovr-1; 
		i2 = i_ovr;
		if (i1 < 1) {i1 = 0; }
		slope = (Dgtz.at(i1) - Dgtz.at(i2))/(i1-i2);
		t_ovr = tmin + samplingTime*(i1 + (threshold-Dgtz.at(i1))/slope) - t_start; // done // it's a time interval
	}
	else { t_ovr = samplingTime*i_ovr;}

	// compute integral
	double i_inf = t_start/samplingTime;
	double i_sup = (t_start+t_ovr)/samplingTime;
	integral = 0;
	for (int i=0;i<Npts;i++){
		if ((i >= i_inf) and (i <= i_sup)){
			integral += (Dgtz.at(i)-threshold);
		}
	}
	integral = integral/1; // done // adc per 44 ns
	// constant fraction discriminator 
	double t_cfd = this->Apply_CFD(0.3,5);
	// output
	std::map<std::string,double> output;
	output["t_start"] = t_start;
	output["t_ovr"] = t_ovr;
	output["integral"] = integral;
	output["max_value"] = max_value;
	output["t_max_value"] = t_max_value;
	output["threshold"] = threshold;
	output["noise_level"] = noise;
	output["t_cfd"] = t_cfd;
	
	return output;

}

double ahdcSignal::Apply_CFD(double CFD_fraction, int CFD_delay){
	int Npts = Dgtz.size();
	std::vector<double> Data = Dgtz;
	// Remove noise 
	double noise = 0;
	for (int i=0;i<5;i++){
		noise += Data.at(i);
	}
	noise = noise/5;
	double ymax = 0;
	for (int i=0;i<Npts;i++){
		Data[i] = Data.at(i) - noise;
		if (ymax < Data.at(i)) ymax = Data.at(i);
	}
	// Start CFD
	std::vector<double> signal(Npts,0.0);
	for (int i=0;i<Npts;i++){
		signal[i] += (1-CFD_fraction)*Data.at(i);
		if (i < Npts-CFD_delay){
			signal[i] += -1*CFD_fraction*Data.at(i+CFD_delay);
		}
	}
	int i_min=0, i_max=0;
	for (int i=0;i<Npts;i++){
		if (signal.at(i_max) < signal.at(i)) i_max = i;
	}
	for (int i=0;i<i_max;i++){ // add this loop to be sure that i_min < i_max
		if (signal.at(i_min) > signal.at(i)) i_min = i;
	}
	// Deternine t_cfd
	int i_ref = 0;
	for (int i=i_min;i<=i_max;i++){
		if (signal.at(i) < 0){
			i_ref = i;
		}
	} // last pass below zero
	int i1 = i_ref; // 1 index below
	int i2 = i_ref+1; // 1 index above
	if (i1 < 0) {i1 = 0; }
	if (i2 >= Npts) {i2 = Npts-1;}
	double slope = (signal.at(i1) - signal.at(i2))/(i1-i2);
	double i_cfd;
	i_cfd = i1 + (0-signal.at(i1))/slope; // DONE
	double t_cfd = tmin + i_cfd*samplingTime; // DONE
	
	return t_cfd;
}

double ahdcSignal::GetMCTime(){
	if (nsteps == 0){ return 0; }
	double mctime = 0;
	double Etot = 0;
	for (int s=0;s<nsteps;s++){
		mctime += DriftTime.at(s)*Edep.at(s);
		Etot += Edep.at(s);
	}
	mctime = mctime/Etot;
	return mctime;
}

double ahdcSignal::GetMCEtot(){
	if (nsteps == 0) { return 0;}
	double mcEtot = 0;
	for (int s=0;s<nsteps;s++){
		mcEtot += Edep.at(s);
	}
	return mcEtot;
}


