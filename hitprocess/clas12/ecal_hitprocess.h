#ifndef ECAL_HITPROCESS_H
#define ECAL_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// ec constants
// these are loaded with initWithRunNumber
class ecConstants
{
public:
	// runNo is mandatory variable to keep track of run number changes
	int    runNo;
	string date;
	string connection;
	char   database[80];
        int    outputRAW=0;  //bypass calibration, resolution and efficiency 
	
	static const int nsect  = 6;  // Number of sectors
	static const int nlayer = 9;  // layer=1-3 (PCAL) 4-6 (ECinner) 7-9 (ECouter)
	static const int nview  = 3;  // Number of views, U,V and W
	
	// For strip dependent constants read from CCDB
	// Array [6][9][3] -> sector,layer,view sector=1-6 layer=1-3 (PCAL) 4-6 (ECinner) 7-9 (ECouter) view=1-3 (U,V,W)
	
	//attlen: attenuation length
	vector<double> attlen[nsect][nlayer][10];
	
	//gain: pmt gain
	vector<double> gain[nsect][nlayer];
	
	//timing: TDC calibration constants
	vector<double>  dtime[nsect][nlayer][9];
	vector<double>  ftime[nsect][nlayer][7];
  
	vector<double>  fadc_offset[nsect][nlayer];
	vector<double>   tmf_offset[nsect][nlayer];
        vector<double>  global_time_walk[nsect][nlayer];
  
  	double  tdc_global_offset;
	double fadc_global_offset;
	
	//veff: effective velocity (cm/ns)
	vector<double> dveff[nsect][nlayer];
	vector<double> fveff[nsect][nlayer];

        // DSC/TSC efficiency
        vector<double> deff[nsect][nlayer][3];

        // FADC threshold
        vector<double> fthr[nsect][nlayer];
  
        // Timiing resolution
        vector<double> dtres[nsect][nlayer][4];
        vector<double> ftres[nsect][nlayer][4];
	
	// status:
	//	0 - fully functioning
	//	1 - noADC
	//	2 - noTDC
	//	3 - noADC, noTDC (PMT is dead)
	// 5 - any other reconstruction problem
	vector<int> status[nsect][nlayer];	
	
	// ======== FADC Pedestals and sigmas ===========
	double pedestal[nsect][nlayer][nview] = {};
	double pedestal_sigm[nsect][nlayer][nview] = {};
	
  //	double TDC_time_to_evio;    // Conversion from time (ns) to TDC format
  	double ADC_GeV_to_evio;     // Conversion from energy (GeV) to FADC250 format
	double pmtQE;               // Quantum efficiency of PMT
	double pmtDynodeGain;       // PMT dynode gain
	double pmtDynodeK;          // PMT dynode secondary emission statistics factor: K=0 (Poisson) K=1 (exponential)
	double pmtFactor;           // Contribution to FWHM from PMT statistical fluctuations.
	
	// voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
	double vpar[4];
	
	// translation table
	TranslationTable TT;
};


// Class definition
class ecal_HitProcess : public HitProcess
{
public:
	
	~ecal_HitProcess(){;}
	
	// constants initialized with initWithRunNumber
	static ecConstants ecc;
	
	void initWithRunNumber(int runno);
	
	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);
	
	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);
	
	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);
	
	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
	// creates the HitProcess
	static HitProcess *createHitClass() {return new ecal_HitProcess;}
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
private:
	
	double fadc_precision = 0.0625;  // 62 picoseconds resolution
	double convert_to_precision(double time) {
		return (int( time / fadc_precision ) * fadc_precision);
	}
  
        double getTRES(double x, double p0, double p1, double p2, double p3) {
	  return (p0*exp(pow(x,p3)/120000) + p1/x + p2/pow(x,0.5));
        }
	
};

#endif
