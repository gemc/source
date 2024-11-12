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

	ahdcSignal *Signal = new ahdcSignal(aHit,hitn,0,6000,1000,44,240);
	Signal->SetElectronYield(100000);
	Signal->Digitize();
	std::map<std::string,double> output = Signal->Extract();

	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = sector;
	dgtz["layer"]     = layer;
	dgtz["component"] = component;
	dgtz["ADC_order"] = 1;
	dgtz["ADC_ADC"]   = (int) output["adcMax"]; 
	dgtz["ADC_time"]  = output["timeMax"];
	dgtz["ADC_ped"]   = (int) output["adcOffset"]; 
	dgtz["ADC_integral"] = (int) output["integral"]; 
	dgtz["ADC_timestamp"] = 0;
	dgtz["ADC_timeRiseCFA"] = output["timeRiseCFA"]; 
	dgtz["ADC_timeCFD"] = output["timeCFD"]; 
	dgtz["ADC_timeOVR"] = output["timeOverThresholdCFA"];
	dgtz["ADC_mctime"] = Signal->GetMCTime(); 
	dgtz["ADC_nsteps"] = Signal->nsteps;
	dgtz["ADC_mcEtot"] = Signal->GetMCEtot(); 

	//dgtz["TDC_order"] = 0;
	//dgtz["TDC_TDC"]   = output["t_start"];
	
	dgtz["wf136_order"] = 1;
	dgtz["wf136_timestamp"] = 0;
	std::vector<short> SDgtz = Signal->GetDgtz();
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
		short adc = (value < ADC_LIMIT) ? value : ADC_LIMIT; // saturation effect 
		Dgtz.push_back(adc);
	}
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

#include <TString.h>
std::map<std::string,double> ahdcSignal::Extract(){
	ahdcExtractor T(samplingTime,0.5f,5,0.3f);
	T.adcOffset = (short) (Dgtz[0] + Dgtz[1] + Dgtz[2] + Dgtz[3] + Dgtz[4])/5;
	std::map<std::string,double> output = T.extract(Dgtz);
	if (nsteps > 10){
		T.Show(TString::Format("./output/SignalDecoded_%d_%d_%d_%d.pdf",hitn,sector,layer,component).Data());
		T.ShowCFD(TString::Format("./output/SignalCFD_%d_%d_%d_%d.pdf",hitn,sector,layer,component).Data());
	}
	return output;
}

std::map<std::string,double> ahdcExtractor::extract(const std::vector<short> samples){
	samplesCorr = samples;
	this->waveformCorrection();
	this->fitAverage();
	this->fitParabolic();
	this->computeTimeAtConstantFractionAmplitude();
	this->computeTimeUsingConstantFractionDiscriminator();
	//this->fineTimeStampCorrection();
	std::map<std::string,double> output;
	output["binMax"] = binMax;
	output["binOffset"] = binOffset;
	output["adcMax"] = adcMax;
	output["timeMax"] = timeMax;
	output["integral"] = integral;
	output["timeRiseCFA"] = timeRiseCFA;
	output["timeFallCFA"] = timeFallCFA;
	output["timeOverThresholdCFA"] = timeOverThresholdCFA;
	output["timeCFD"] = timeCFD;
	output["adcOffset"] = adcOffset;
	return output;

}

void ahdcExtractor::waveformCorrection(){
	binNumber = samplesCorr.size();
	binMax = 0;
	adcMax = (short) (samplesCorr[0] - adcOffset);
	integral = 0;
	for (int bin = 0; bin < binNumber; bin++){
		samplesCorr[bin] = (short) (samplesCorr[bin] - adcOffset);
		if (adcMax < samplesCorr[bin]){
			adcMax = samplesCorr[bin];
			binMax = bin;
		}
		integral += samplesCorr[bin];
	}
	/*
	 * If adcMax + adcOffset == ADC_LIMIT, that means there is saturation
	 * In that case, binMax is the middle of the first plateau
	 * This convention can be changed
	 */
	if ((short) adcMax + adcOffset == ADC_LIMIT) {
		int binMax2 = binMax;
		for (int bin = binMax; bin < binNumber; bin++){
			if (samplesCorr[bin] + adcOffset == ADC_LIMIT) {
				binMax2 = bin;
			}
			else {
				break;
			}
		}
		binMax = (binMax + binMax2)/2;
	}
	binOffset = sparseSample*binMax;
	timeMax = (binMax + binOffset)*samplingTime;
}


void ahdcExtractor::fitAverage(){
	if ((binMax - 2 >= 0) && (binMax + 2 <= binNumber - 1)){
		adcMax = 0;
		for (int bin = binMax - 2; bin <= binMax + 2; bin++){
			adcMax += samplesCorr[bin];
		}
		adcMax = adcMax/5;
	}
}

void ahdcExtractor::fitParabolic(){}

void ahdcExtractor::fineTimeStampCorrection(){}

void ahdcExtractor::computeTimeAtConstantFractionAmplitude(){
	float threshold = amplitudeFractionCFA*adcMax;
	// timeRiseCFA
	int binRise = 0;
	for (int bin = 0; bin < binMax; bin++){
		if (samplesCorr[bin] < threshold)
			binRise = bin;  // last pass below threshold and before adcMax
	} // at this stage : binRise < timeRiseCFA/samplingTime <= binRise + 1 // timeRiseCFA is determined by assuming a linear fit between binRise and binRise + 1
	float slopeRise = 0;
	if (binRise + 1 <= binNumber-1)
		slopeRise = samplesCorr[binRise+1] - samplesCorr[binRise];
	float fittedBinRise = (slopeRise == 0) ? binRise : binRise + (threshold - samplesCorr[binRise])/slopeRise;
	timeRiseCFA = (fittedBinRise + binOffset)*samplingTime; // binOffset is determined in wavefromCorrection() // must be the same for all time ? // or must be defined using fittedBinRise*sparseSample

	// timeFallCFA
	int binFall = binMax;
	for (int bin = binMax; bin < binNumber; bin++){
		if (samplesCorr[bin] > threshold){
				binFall = bin;
		}
		else {
				binFall = bin;
				break; // first pass below the threshold
		}
	} // at this stage : binFall - 1 <= timeRiseCFA/samplingTime < binFall // timeFallCFA is determined by assuming a linear fit between binFall - 1 and binFall
	float slopeFall = 0;
	if (binFall - 1 >= 0)
		slopeFall = samplesCorr[binFall] - samplesCorr[binFall-1];
	float fittedBinFall = (slopeFall == 0) ? binFall : binFall-1 + (threshold - samplesCorr[binFall-1])/slopeFall;
	timeFallCFA = (fittedBinFall + binOffset)*samplingTime;
	
	// timeOverThreshold
	timeOverThresholdCFA = timeFallCFA - timeRiseCFA;
}

void ahdcExtractor::computeTimeUsingConstantFractionDiscriminator(){
	std::vector<float> signal(binNumber,0.0);
	// signal generation
	for (int bin = 0; bin < binNumber; bin++){
		signal[bin] = (1 - fractionCFD)*samplesCorr[bin]; // we fill it with a fraction of the original signal
		if (bin < binNumber - binDelayCFD)
			signal[bin] += -1*fractionCFD*samplesCorr[bin + binDelayCFD]; // we advance and invert a complementary fraction of the original signal and superimpose it to the previous signal
	}
	// determine the two humps
	int binHumpSup = 0;
	int binHumpInf = 0;
	for (int bin = 0; bin < binNumber; bin++){
		if (signal[bin] > signal[binHumpSup])
			binHumpSup = bin;
	}
	for (int bin = 0; bin < binHumpSup; bin++){ // this loop has been added to be sure : binHumpInf < binHumpSup
		if (signal[bin] < signal[binHumpInf])
			binHumpInf = bin;
	}
	// research for zero
	int binZero = 0;
	for (int bin = binHumpInf; bin <= binHumpSup; bin++){
		if (signal[bin] < 0)
			binZero = bin; // last pass below zero
	} // at this stage : binZero < timeCFD/samplingTime <= binZero + 1 // timeCFD is determined by assuming a linear fit between binZero and binZero + 1
	float slopeCFD = 0;
	if (binZero + 1 <= binNumber)
		slopeCFD = signal[binZero+1] - signal[binZero];
	float fittedBinZero = (slopeCFD == 0) ? binZero : binZero + (0 - signal[binZero])/slopeCFD;
	timeCFD = (fittedBinZero + binOffset)*samplingTime;
	//
	samplesCFD = signal;
}

#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1.h"
#include "TGraphPolar.h"
#include "TGaxis.h"
#include <time.h>
#include "TLine.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TArrow.h"
#include <fstream>

void ahdcExtractor::Show(const char * filename){
	double tmin = 0;
	double tmax = 6000;
        double threshold = adcOffset + amplitudeFractionCFA*adcMax;

        // Main graph
        double ymax = 0;
        TGraph* gr1 = new TGraph(binNumber);
        for (int i=0;i<binNumber;i++){
                int adc = samplesCorr.at(i) + adcOffset;
                if (ymax < adc) ymax = adc;
                gr1->SetPoint(i,tmin + i*samplingTime,adc);
        }
        // Graph for filling
	int binNumberOVR = 0; // binNumberOverthresholdCFA
	for (int bin=0; bin< binNumber;bin++){
		if ((bin*samplingTime > timeRiseCFA) && (bin*samplingTime < timeFallCFA))
			binNumberOVR++;
	}
        TGraph* gr2 = new TGraph(binNumberOVR+2);
        gr2->SetPoint(0,timeRiseCFA,threshold);
        gr2->SetPoint(binNumberOVR+1,timeFallCFA, threshold);
	int binRiseCFA = (int) timeRiseCFA/samplingTime;
        for (int i=1;i<=binNumberOVR;i++){
                gr2->SetPoint(i,tmin+samplingTime*(binRiseCFA+i),samplesCorr.at(binRiseCFA+i)+adcOffset);
        }

        // Plot graph
        TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
        gr1->SetTitle("");
        //gr1->SetTitle(TString::Format("%s : p = #Box MeV, #theta = #Box deg, #phi = #Box deg",pid2name[pid].data()));
        gr1->GetXaxis()->SetTitle("Time (ns)");
        gr1->GetXaxis()->SetTitleSize(0.05);
        gr1->GetYaxis()->SetTitle("Charge (adc)");
        gr1->GetYaxis()->SetTitleSize(0.05);
        gr1->GetYaxis()->SetRangeUser(0,ymax+0.05*ymax);
        gr1->SetMarkerColor(kBlack);
        gr1->SetMarkerSize(5);
        gr1->SetLineColor(kBlue);
        gr1->Draw("APL");
        gr2->SetFillColorAlpha(kGreen, 1.0);
        gr2->Draw("F");

        // View decoding
        TLine* line1 = new TLine(timeRiseCFA,0,timeRiseCFA,threshold); line1->SetLineWidth(1); line1->SetLineColor(kRed); line1->SetLineStyle(2); line1->Draw(); // timeRiseCFA
        TLine* line2 = new TLine(0,threshold,timeRiseCFA,threshold); line2->SetLineWidth(1); line2->SetLineColor(kRed); line2->SetLineStyle(2); line2->Draw(); // timeRiseCFA
        TLine* line3 = new TLine(timeFallCFA,0,timeFallCFA,threshold); line3->SetLineWidth(1); line3->SetLineColor(kRed); line3->SetLineStyle(2); line3->Draw(); // timeFallCFA
        TArrow* arrow1 = new TArrow(timeRiseCFA,threshold,timeFallCFA,threshold,0.02,"<>"); arrow1->SetLineWidth(1); arrow1->SetLineColor(kRed); arrow1->Draw(); // timeOverThresholdCFA
        TLine* line4 = new TLine(tmin,adcOffset,tmax,adcOffset); line4->SetLineWidth(1); line4->SetLineColor(kRed); line4->SetLineStyle(2); line4->Draw(); // adcOffset
        TLine* line5 = new TLine(0,adcMax+adcOffset,timeMax,adcMax+adcOffset); line5->SetLineWidth(1); line5->SetLineColor(kRed); line5->SetLineStyle(2); line5->Draw(); // adcMax+adcOffset
        TLine* line6 = new TLine(timeMax,0,timeMax,adcMax+adcOffset); line6->SetLineWidth(1); line6->SetLineColor(kRed); line6->SetLineStyle(2); line6->Draw(); // adcMax+adcOffset

        TLatex data;
        data.SetTextSize(0.03);
        data.SetTextAlign(13);
        data.DrawLatexNDC(0.5,0.8,TString::Format("#bf{#bf{timeRiseCFA} =  %.2lf ns}",timeRiseCFA).Data());
        data.DrawLatexNDC(0.5,0.8-0.05,TString::Format("#bf{#bf{timeOverThresholdCFA} =  %.2lf ns}",timeOverThresholdCFA).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*2,TString::Format("#bf{#bf{adcMax+adcOffset} =  %.0lf adc }",adcMax+adcOffset).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*3,TString::Format("#bf{#bf{integral} =  %.0lf adc per 44 ns}",integral).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*4,TString::Format("#bf{#bf{adcOffset} =  %d adc}",adcOffset).Data());
        data.DrawLatexNDC(0.5,0.8-0.05*5,"#bf{1 adc =  10^{-5} keV/ns }");
        data.DrawLatexNDC(0.5,0.8-0.05*6,TString::Format("#bf{#bf{timeCFD} = %.2lf ns}",timeCFD));
        data.SetTextAlign(11);
        data.DrawLatex(timeRiseCFA,0+(adcMax+adcOffset)*0.02,"timeRiseCFA");
        data.DrawLatex(timeMax,(adcMax+adcOffset)+(adcMax+adcOffset)*0.02,"adcMax+adcOffset");
        data.DrawLatex(timeMax,threshold,"timeOverThresholdCFA");
        data.DrawLatex(0+tmax*0.02,threshold,"threshold");
        data.DrawLatex(tmax, adcOffset+(adcMax+adcOffset)*0.02,"adcOffset");

        canvas1->Print(filename);
        delete line1; delete line2; delete line3; delete line5;
        delete gr1; delete gr2; delete arrow1;
        delete canvas1;
}

void ahdcExtractor::ShowCFD(const char * filename){
	double tmin = 0;
	double tmax = 6000;
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	TGraph* gr1 = new TGraph(binNumber);
	TGraph* gr2 = new TGraph(binNumber);
	for (int i=0;i<binNumber;i++){
		gr1->SetPoint(i,tmin + i*samplingTime,samplesCorr.at(i)+adcOffset);
		gr2->SetPoint(i,tmin + i*samplingTime,samplesCFD.at(i));
	}
	gr2->SetTitle(TString::Format("#bf{CFD},  fraction = %.1lf, delay = %d index units",fractionCFD,binDelayCFD));
        gr2->GetXaxis()->SetTitle("Time (ns)");
        gr2->GetXaxis()->SetTitleSize(0.05);
        gr2->GetYaxis()->SetTitle("Charge (adc)");
        gr2->GetYaxis()->SetTitleSize(0.05);
        //gr2->SetLineStyle(1);
        gr2->SetLineColor(kRed);
        gr2->SetMarkerColor(kRed);
        gr2->SetMarkerSize(5);
        gr2->GetYaxis()->SetRangeUser(-fractionCFD*adcMax-0.1*adcMax,adcOffset+adcMax+0.1*adcMax);
        gr2->Draw("ALP");

        gr1->SetMarkerColor(kBlue);
        gr1->SetMarkerSize(5);
        gr1->SetLineColor(kBlue);
        //gr1->SetLineStyle(2);
        gr1->Draw("PL");

        TGaxis* axis1 = new TGaxis(tmin,0,tmax,0,tmin,tmax,510,"");
        axis1->SetLineColor(kGreen);
        axis1->SetLabelColor(kGreen);
        //axis1->SetTitle("index units");
        axis1->Draw();

        TLegend* legend = new TLegend(0.7,0.8,0.9,0.9);
        legend->AddEntry(gr1,"Digitized signal","l");
        legend->AddEntry(gr2,"CFD signal","l");
        legend->Draw();

        TLatex data;
        data.SetTextSize(0.04);
        data.SetTextAlign(13);
        data.DrawLatexNDC(0.7,0.6,TString::Format("#bf{#bf{timeCFD} =  %.2lf ns}",timeCFD).Data());


        //canvas1->Print(TString::Format("./output/SignalCFD_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	canvas1->Print(filename);
        delete gr1; delete gr2; delete canvas1;
        delete axis1; delete legend;
}
