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

	ahdcSignal *Signal = new ahdcSignal(aHit,hitn);
	// Set parameters for digitization
	Signal->SetTmin(0);
	Signal->SetTmax(6000);
	Signal->SetDelay(1000);
	Signal->SetSamplingTime(44); // ns
	Signal->SetElectronYield(100000);
	//Signal->SetAdcMax(10000);
	Signal->Digitize();
	std::map<std::string,double> output = Signal->Decode(Signal->Get_nsteps() > 10);
	//std::map<std::string,double> output = Signal->Decode(false);

	//if (Signal->Get_nsteps() >= 10) {
		//Signal->PrintBeforeProcessing();
		//Signal->PrintAllShapes();
		//Signal->PrintAfterProcessing();
		//Signal->PrintNoise();
		//Signal->PrintSignal();
		//Signal->ShowDecoding();
	//	output = Signal->Decode();
	//}
	delete Signal;

	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = sector;
	dgtz["layer"]     = layer;
	dgtz["component"] = component;
	dgtz["ADC_order"] = 0;
	dgtz["ADC_ADC"]   = (int) output["max_value"]; // adc
	dgtz["ADC_time"]  = output["t_ovr"]; // ns
	dgtz["ADC_ped"]   = (int) output["noise_level"]; // adc
	dgtz["ADC_integral"] = (int) output["integral"]; // adc per 44 ns
	dgtz["ADC_timestamp"] = output["t_start"]; // ns
	dgtz["ADC_t_cfd"] = output["t_cfd"];

	dgtz["TDC_order"] = 0;
	dgtz["TDC_TDC"]   = output["t_start"];
		
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
	
	// ******** cout alert geometry ***********
	// 
	// ****************************************
	/*if (nsteps > 10) {
		std::cout << "     ===> alert geometry, geant4 type :   " << aHit->GetDetector().type << std::endl;
		std::cout << "     ===> alert geometry, geant4 type :   " << aHit->GetDetector().dimensions.size() << std::endl;
		for (int i=0;i< (int) aHit->GetDetector().dimensions.size() ;i++) {
		std::cout << "     dim_%d" << i << " : "  << aHit->GetDetector().dimensions[i] << std::endl;
		}
	}*/
	// ******** end cout *********************

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
	// std::cout << "=======> Inside ComputeDocaAndTime" << std::endl;
	// std::cout << "   X_sigwire_top : " << X_sigwire_top << std::endl;
	// std::cout << "   X_sigwire_bot : " << X_sigwire_bot << std::endl;
	// std::cout << "   Y_sigwire_top : " << Y_sigwire_top << std::endl;
	// std::cout << "   Y_sigwire_bot : " << Y_sigwire_bot << std::endl;

	// Triangle abh
	// a (sigwire_top), b (sigwire_bot), h (hit position)
	// H_abh is the distance between hit and the wire and perpendicular to the wire
	double L_ab, L_ah, L_bh, H_abh;
	// Compute the distance between top and bottom of the wire
	L_ab = sqrt(pow(X_sigwire_top-X_sigwire_bot,2) + pow(Y_sigwire_top-Y_sigwire_bot,2) + pow(Z_sigwire_top-Z_sigwire_bot,2));
	// std::cout << "   L_ab : " << L_ab << std::endl;
	// double doca = DBL_MAX;
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
		//if (doca > H_abh) doca = H_abh; 
		//Height.at(s) = H_abh;
		Doca.push_back(H_abh);
		// Add a resolution on doca
		double docasig = 337.3-210.3*H_abh+34.7*pow(H_abh,2); // um // fit sigma vs distance // Fig 4.14 (right), L. Causse's thesis
		docasig = docasig/1000; // mm
		std::default_random_engine dseed(time(0)); //seed
		std::normal_distribution<double> docadist(H_abh, docasig);
		// Compute time
		double driftTime = 7*H_abh + 7*pow(H_abh,2) + 4*pow(H_abh,3); // fit t vs distance //  Fig 4.12 (right), L. Causse's thesis
		//Time.at(s) = driftTime; // ns
		DriftTime.push_back(driftTime);

		 // if ((s==0) and (nsteps > 10)) {
		 //	 std::cout << "=======> Inside ComputeDocaAndTime" << std::endl;		 	
		 //      LposX = Lpos[s].x();
		 //      LposY = Lpos[s].y();
		 //      LposZ = Lpos[s].z();
		 //      std::cout << "      L_ah  : " << L_ah << std::endl;
		 //      std::cout << "      L_bh  : " << L_bh << std::endl;
		 //      std::cout << "      H_abh : " << H_abh << std::endl;
		 //	 std::cout << "      docasig : " << docasig << std::endl;
		 // }

	}
}

namespace futils {
	bool cart2polar3D(double x, double y, double z, double & rho, double & theta, double & phi){
		rho = sqrt(x*x+y*y+z*z);
		if (rho <= __DBL_EPSILON__) {return false;} // if rho == 0
		theta = acos(z/rho);
		if (y >= 0){
			phi = acos(x/(rho*sin(theta)));
		} 
		else {
			phi = 2*PI - acos(x/(rho*sin(theta)));			
		}
		return true;
	}

}


void ahdcSignal::PrintBeforeProcessing(){
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	// Draw stems 
	TGraph* gr1 = new TGraph(nsteps);
	for (int s=0;s<nsteps;s++){
		// Draw points
		gr1->SetPoint(s,DriftTime.at(s),Edep.at(s));
	}
	// Draw graph
	//gr1->SetTitle("Deposited energy in each steps");
	gr1->SetTitle("");
	gr1->GetXaxis()->SetTitle("G4Time (ns)");
	gr1->GetYaxis()->SetTitle("Edep (keV)");
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(kRed);
	gr1->SetMarkerSize(2);
	gr1->Draw("AP");
	for (int s=0;s<nsteps;s++){
		// Draw lines
		TLine* line = new TLine(DriftTime.at(s),0,DriftTime.at(s),Edep.at(s));
		line->SetLineWidth(1);
		line->SetLineColor(kBlack);
		line->Draw();
	}
	
	canvas1->Print(TString::Format("./output/SignalBeforeProcessing_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	delete gr1; 
	delete canvas1;
}

void ahdcSignal::PrintAllShapes(){
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	TLegend* legend = new TLegend();
	// Draw all shapes
	int Npts = 1000;
	double ymax=0;
	int s_ref=0;
	// Define ymax and s_ref
	for (int s=0;s<nsteps;s++){ 
		double xRange[Npts], yRange[Npts];
		for (int i=0;i<Npts;i++){
			xRange[i] = tmin + i*(tmax-tmin)/Npts;
			yRange[i] =  Edep.at(s)*ROOT::Math::landau_pdf(xRange[i]-delay,600/2.5,DriftTime.at(s))*1000; // normalisation constant for better
		    	if (ymax < yRange[i]) {ymax = yRange[i]; s_ref = s;}
		}
	}
	for (int s=0;s<nsteps;s++){
		if (ymax < Edep.at(s)) 
			ymax = Edep.at(s);
	}
	// plot each distribation
	{ 	// In s_ref
		TGraph* gr2 = new TGraph(Npts);
		double xRange[Npts], yRange[Npts];
		for (int i=0;i<Npts;i++){
			xRange[i] = tmin + i*(tmax-tmin)/Npts;
			yRange[i] =  Edep.at(s_ref)*ROOT::Math::landau_pdf(xRange[i]-delay,600/2.5,DriftTime.at(s_ref))*1000; // normalisation constant for a better view
			gr2->SetPoint(i,xRange[i],yRange[i]);
		}
		gr2->SetLineColor(s_ref+1);
		gr2->SetFillColorAlpha(2+s_ref%38,1.0);
		gr2->SetFillStyle(3001);
		//gr2->SetTitle("Spreed of each Edep over the time using a Landau distribution");
		gr2->SetTitle("");
		gr2->GetXaxis()->SetTitle("Time (ns)");
		gr2->GetYaxis()->SetTitle("#frac{d Edep}{dt} (keV/ns)");
		gr2->GetYaxis()->SetRangeUser(0,ymax+0.05*ymax);
		gr2->Draw("AL"); // we use the axis of s_ref as reference

	}
	for (int s=0;s<nsteps;s++){ 
		if (s == s_ref) continue;
		else {
			TGraph* gr2 = new TGraph(Npts);
			double xRange[Npts], yRange[Npts];
			for (int i=0;i<Npts;i++){
				xRange[i] = tmin + i*(tmax-tmin)/Npts;
				yRange[i] =  Edep.at(s)*ROOT::Math::landau_pdf(xRange[i]-delay,600/2.5,DriftTime.at(s))*1000; // normalisation constant for a better view
				gr2->SetPoint(i,xRange[i],yRange[i]);
			}
			gr2->SetLineColor(s+1);
			gr2->SetFillColorAlpha(2+s%38,1.0);
			gr2->SetFillStyle(3001);
			gr2->Draw("L");
			legend->AddEntry(gr2,TString::Format("Shape %d",s),"l");
	      	}
	}
	// Draw stems
	TGraph* gr1 = new TGraph(nsteps);
	for (int s=0;s<nsteps;s++){
		// Draw points
		gr1->SetPoint(s,DriftTime.at(s)+delay,Edep.at(s));
	}
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(kRed);
	gr1->SetMarkerSize(2);
	gr1->Draw("P");
	for (int s=0;s<nsteps;s++){
		// Draw lines
		TLine* line = new TLine(DriftTime.at(s)+delay,0,DriftTime.at(s)+delay,Edep.at(s));
		line->SetLineWidth(1);
		line->SetLineColor(kBlack);
		line->Draw();
	}
       	// Draw text
	TLatex latex2;
	latex2.SetTextSize(0.025);
	latex2.SetTextAlign(13);
	//latex2.DrawLatex(tmax/2, 2*ymax/3,"#bf{#splitline{All shapes are nomalised x 2000}{for a better view}}");
	latex2.DrawLatex(tmax/2, 2*ymax/3,TString::Format("#bf{#splitline{All shapes are nomalised x 1000}{A delay of #bf{%.1lf ns} as been added}}",delay).Data());
	// Draw legend
	legend->SetX1(0.82);
	legend->SetY1(0.3); 
	legend->SetX2(0.95);
	legend->SetY2(0.95);
	legend->AddEntry(gr1,"Edep in each steps","p");
	legend->Draw();
	// Print file
	canvas1->Print(TString::Format("./output/SignalAllShapes_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	delete gr1; delete legend;
	delete canvas1;
}

void ahdcSignal::PrintAfterProcessing(){
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768); 
	// Draw graph
	double ymax = 0;
	int Npts = 1000;
    	TGraph* gr1 = new TGraph(Npts);
    	for (int i=0;i<Npts;i++){
		double x_ = tmin + i*(tmax-tmin)/Npts;
		double y_ = this->operator()(x_);
		if (ymax < y_) ymax = y_;
		gr1->SetPoint(i,x_,y_);
    	}
	gr1->SetLineColor(kBlue);
	gr1->SetFillColorAlpha(kBlue,1.0);
	gr1->SetFillStyle(3001);
	//gr1->SetTitle("");
	gr1->SetTitle(TString::Format("%s : p = □  MeV, #theta = □  deg, #phi = □  deg",pid2name[pid].data()));
	gr1->GetXaxis()->SetTitle("Time (ns)");
	gr1->GetYaxis()->SetTitle("#frac{d Edep}{dt} (keV/ns)");
	//gr1->Draw("ALF"); 
	gr1->Draw("AL");
	// Draw text
	TLatex latex2;
	latex2.SetTextSize(0.03);
	latex2.SetTextAlign(13);
	//latex2.DrawLatex(tmax/2, 2*ymax/3,TString::Format("#bf{A delay of #bf{%.1lf ns} as been added}",delay).Data());
	// Print file
	canvas1->Print(TString::Format("./output/SignalAfterProcessing_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	delete gr1; 
	delete canvas1;
}

void ahdcSignal::GenerateNoise(double mean, double stdev){
	int Npts = (int) floor( (tmax-tmin)/samplingTime );
	// define de 1st value
	//std::default_random_engine dseed(time(0)); //seed
	std::random_device rd;      // Create a random device to seed the generator
	std::mt19937 gen(rd());     // Create a random number engine (e.g., Mersenne Twister)
	std::normal_distribution<double> draw1(mean,stdev);
	// double value = draw1(dseed);
	double value = draw1(gen);
	if (value < 0) value = 0;
	Noise.push_back(value);
	for (int i=1;i<Npts;i++){
		//std::normal_distribution<double> draw(Noise.at(i-1),stdev);
		std::normal_distribution<double> draw(mean,stdev);
		value = draw(gen);
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

void ahdcSignal::PrintNoise(){
	//int Npts = Noise.size();
	int Npts = 100;
	// Define canvas
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	// Draw graph
	double dt = (tmax-tmin)/Npts;
	TGraph* gr1 = new TGraph(Npts);
	for (int i=0;i<Npts;i++){
		gr1->SetPoint(i,tmin + i*dt,Noise.at(i));
	}
	gr1->SetLineColor(kBlue);
	gr1->SetMarkerColor(kRed);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerSize(1);
	gr1->SetTitle("");
	gr1->GetXaxis()->SetTitle("Time (ns)");
	gr1->GetYaxis()->SetTitle("Noise (adc)");
	gr1->Draw("APL");
	// Print file
	canvas1->Print(TString::Format("./output/SignalNoise_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	delete gr1;
	delete canvas1;
}

void ahdcSignal::PrintSignal(){
	int Npts = Dgtz.size(); 
	// Histogram
	TH1D * hist = new TH1D("hist_adc","hist_adc",Npts,tmin,tmax);
	TGraph* gr1 = new TGraph(Npts);
	double ymax = 0;
	for (int i=0;i<Npts;i++){
		int adc = Dgtz.at(i); // in ADC
		for (int j=0;j<adc;j++)
			hist->Fill(tmin + i*samplingTime);
		if (ymax < adc) ymax = adc;
		gr1->SetPoint(i,tmin + i*samplingTime,adc);
	}
	
	// Plot graph 
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	/*gStyle->SetOptStat("nemruo");
	hist->GetXaxis()->SetTitle("Time (ns)");
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetTitle("Charge (adc)");
	hist->GetYaxis()->SetTitleSize(0.05);
	hist->Draw();*/
	
	//gr1->SetTitle("");
	gr1->SetTitle(TString::Format("%s : p = #Box MeV, #theta = #Box deg, #phi = #Box deg",pid2name[pid].data()));
	gr1->GetXaxis()->SetTitle("Time (ns)");
	gr1->GetYaxis()->SetTitle("Charge (adc)");
	//gr1->GetYaxis()->SetRangeUser(0,ymax+0.05*ymax);
	gr1->SetLineColor(kBlue);
	gr1->SetMarkerStyle(1);
	gr1->SetMarkerSize(5);
	gr1->SetMarkerColor(kRed);
	//gr1->SetFillColorAlpha(kBlue,1.0);
	//gr1->SetFillStyle(3001);
	gr1->Draw("APL");

	// Draw text
	TLatex latex1;
	latex1.SetTextSize(0.04);
	latex1.SetTextAlign(23);
	//latex1.DrawLatex(tmax/2, 2*ymax/3,TString::Format("#bf{A delay of #bf{%.1lf ns} as been added}",delay).Data());
	// Print file
	canvas1->Print(TString::Format("./output/SignalDigitized_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	delete hist;
	delete gr1;
	delete canvas1;

}


std::map<std::string,double> ahdcSignal::Decode(bool printFigure){

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
	int Npts2=0;
	for (int i=0;i<Npts;i++){
		if ((i >= i_inf) and (i <= i_sup)){
			integral += (Dgtz.at(i)-noise);
			Npts2++;
		}
	}
	integral = integral/1; // done // adc per 44 ns
	// constant fraction discriminator 
	double t_cfd = this->Apply_CFD(0.3,60,printFigure);
	// output
	std::map<std::string,double> output;
	output["t_start"] = t_start;
	output["t_ovr"] = t_ovr;
	output["integral"] = integral;
	output["max_value"] = max_value;
	output["t_max_value"] = t_max_value;
	output["threshold"] = threshold;
	output["noise_level"] = noise;
	output["Npts2"] = Npts2;
	output["t_cfd"] = t_cfd;
	
	if (printFigure) this->ShowDecoding(output);
	return output;

}

void ahdcSignal::ShowDecoding(std::map<std::string,double> output){
	int Npts = Dgtz.size(); 
	int Npts2 = output["Npts2"];
	double t_start = output["t_start"];
	double t_ovr = output["t_ovr"];
	double integral = output["integral"];
	double max_value = output["max_value"];
	double t_max_value = output["t_max_value"];
	double threshold = output["threshold"];
	double noise = output["noise_level"];
	int i_inf = (int) t_start/samplingTime;
	double t_cfd = output["t_cfd"];

	// Main graph
	double ymax = 0;
	TGraph* gr1 = new TGraph(Npts);
	for (int i=0;i<Npts;i++){
		int adc = Dgtz.at(i); 
		if (ymax < adc) ymax = adc;
		gr1->SetPoint(i,tmin + i*samplingTime,adc);
	}
	// Graph for filling
	TGraph* gr2 = new TGraph(Npts2+2);
	gr2->SetPoint(0,t_start,threshold);
	gr2->SetPoint(Npts2+1,t_start+t_ovr, threshold);
	for (int i=1;i<=Npts2;i++){
		gr2->SetPoint(i,tmin+samplingTime*(i_inf+i),Dgtz.at(i_inf+i));
	}

	// Plot graph 
	TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
	gr1->SetTitle(TString::Format("%s : p = #Box MeV, #theta = #Box deg, #phi = #Box deg",pid2name[pid].data()));
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
	TLine* line1 = new TLine(t_start,0,t_start,threshold); line1->SetLineWidth(1); line1->SetLineColor(kRed); line1->SetLineStyle(2); line1->Draw(); // t_start
	TLine* line2 = new TLine(0,threshold,t_start,threshold); line2->SetLineWidth(1); line2->SetLineColor(kRed); line2->SetLineStyle(2); line2->Draw(); // t_start
	TLine* line3 = new TLine(t_start+t_ovr,0,t_start+t_ovr,threshold); line3->SetLineWidth(1); line3->SetLineColor(kRed); line3->SetLineStyle(2); line3->Draw(); // t_ovr
	TArrow* arrow1 = new TArrow(t_start,threshold,t_start+t_ovr,threshold,0.02,"<>"); arrow1->SetLineWidth(1); arrow1->SetLineColor(kRed); arrow1->Draw(); // t_ovr
	TLine* line4 = new TLine(tmin,noise,tmax,noise); line4->SetLineWidth(1); line4->SetLineColor(kRed); line4->SetLineStyle(2); line4->Draw(); // noise level
	TLine* line5 = new TLine(0,max_value,t_max_value,max_value); line5->SetLineWidth(1); line5->SetLineColor(kRed); line5->SetLineStyle(2); line5->Draw(); // max_value
	TLine* line6 = new TLine(t_max_value,0,t_max_value,max_value); line6->SetLineWidth(1); line6->SetLineColor(kRed); line6->SetLineStyle(2); line6->Draw(); // max_value

	TLatex data;
	data.SetTextSize(0.03);
	data.SetTextAlign(13);
	data.DrawLatexNDC(0.5,0.8,TString::Format("#bf{#bf{t_start} =  %.2lf ns}",t_start).Data());
	data.DrawLatexNDC(0.5,0.8-0.05,TString::Format("#bf{#bf{t_ovr} =  %.2lf ns}",t_ovr).Data());
	data.DrawLatexNDC(0.5,0.8-0.05*2,TString::Format("#bf{#bf{max_value} =  %.0lf adc }",max_value).Data());
	data.DrawLatexNDC(0.5,0.8-0.05*3,TString::Format("#bf{#bf{integral} =  %.0lf adc per 44 ns}",integral).Data());
	data.DrawLatexNDC(0.5,0.8-0.05*4,TString::Format("#bf{#bf{noise level} =  %.0lf adc}",noise).Data());
	data.DrawLatexNDC(0.5,0.8-0.05*5,"#bf{1 adc =  10^{-5} keV/ns }");
	data.DrawLatexNDC(0.5,0.8-0.05*6,TString::Format("#bf{#bf{t_cfd} = %.2lf ns}",t_cfd));
	data.SetTextAlign(11);	
	data.DrawLatex(t_start,0+max_value*0.02,"t_start");
	data.DrawLatex(t_max_value,max_value+max_value*0.02,"max_value");
	data.DrawLatex(t_max_value,threshold,"t_ovr");
	data.DrawLatex(0+tmax*0.02,threshold,"threshold");
	data.DrawLatex(tmax, noise+max_value*0.02,"noise level");

	canvas1->Print(TString::Format("./output/SignalDecoded_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
	delete line1; delete line2; delete line3; delete line5;
	delete gr1; delete gr2; delete arrow1;
	delete canvas1;
}



double ahdcSignal::Apply_CFD(double CFD_fraction, int CFD_delay, bool printFigure){
	int Npts = Dgtz.size();
	std::vector<double> Data = Dgtz;
	// Remove noise 
	double noise = 0;
	for (int i=0;i<5;i++){
		noise += Data.at(i);
	}
	noise = noise/5;
	for (int i=0;i<Npts;i++){
		Data[i] = Data.at(i) - noise;
	}
	// Start CFD
	std::vector<double> signal(Npts,0.0);
	for (int i=0;i<Npts;i++){
		signal[i] += -1*CFD_fraction*Data.at(i);
		if (i >= CFD_delay) {
			signal[i] += Data.at(i-CFD_delay);
		}
	}
	// Deternine t_cfd
	int i_ref = 0;
	for (int i=0;i<Npts-1;i++){
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
	
	if (printFigure) {
		std::cout << "--------------------------------------" << std::endl;
		std::cout << "|       t_CFD  :  "  << t_cfd << std::endl;
		std::cout << "--------------------------------------" << std::endl;
		
		TCanvas* canvas1 = new TCanvas("c1","c1 title",1366,768);
		TGraph* gr1 = new TGraph(Npts);
		TGraph* gr2 = new TGraph(Npts);
		for (int i=0;i<Npts;i++){
			gr1->SetPoint(i,tmin + i*samplingTime,Dgtz.at(i));
			gr2->SetPoint(i,tmin + i*samplingTime,signal.at(i));
		}
		gr2->SetTitle(TString::Format("#bf{CFD},  fraction = %.0lf, delay = %d index units",CFD_fraction,CFD_delay));
		gr2->GetXaxis()->SetTitle("Time (ns)");
		gr2->GetXaxis()->SetTitleSize(0.05);
		gr2->GetYaxis()->SetTitle("Charge (adc)");
		gr2->GetYaxis()->SetTitleSize(0.05);
		//gr2->SetLineStyle(1);
		gr2->SetLineColor(kRed);
		gr2->SetMarkerColor(kRed);
		gr2->SetMarkerSize(5);
		gr2->Draw("ALP");

		gr1->SetMarkerColor(kBlue);
		gr1->SetMarkerSize(5);
		gr1->SetLineColor(kBlue);
		//gr1->SetLineStyle(2);
		gr1->Draw("PL");

		TGaxis* axis1 = new TGaxis(0,0,tmax,0,0,tmax,510,"");
		axis1->SetLineColor(kGreen);
		axis1->SetLabelColor(kGreen);
		axis1->Draw();

		TLatex data;
		data.SetTextSize(0.04);
		data.SetTextAlign(13);
		data.DrawLatexNDC(0.7,0.8,TString::Format("#bf{#bf{t_cfd} =  %.2lf ns}",t_cfd).Data());
		
		canvas1->Print(TString::Format("./output/SignalCFD_%d_%d_%d_%d.pdf",hitn,sector,layer,component));
		delete gr1; delete gr2; delete canvas1;
		delete axis1;
		
		//std::ofstream fichier(TString::Format("./output/Data_CFD_%d_%d_%d_%d.txt",hitn,sector,layer,component));
		//if (fichier.is_open()) {
		//	fichier << Npts << std::endl;
		//	for (int i=0;i<Npts;i++){
		//		fichier << Dgtz.at(i) << " ";
		//	}
		//	fichier.close();
		//}
	}
	
	return t_cfd;
}

//double ahdcSignal::Apply_CFD(double CFD_fraction, double CFD_delay, bool printFigure){ // CFD_delay in ns
//	int i_CFD_delay = (int) CFD_delay/samplingTime;
//	return this->Apply_CFD(CFD_fraction,i_CFD_delay,printFigure);
//}

//void ahdcSignal::Show_CFD(){
//
//
//}
