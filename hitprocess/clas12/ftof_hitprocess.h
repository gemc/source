#ifndef FTOF_HITPROCESS_H
#define FTOF_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"


class ftofConstants
{
public:

	// Database parameters
	int    runNo;
	string variation;
	string date;
	string connection;
	char   database[80];
        
	// For paddle dependent constants read from CCDB
	// Array [6][3][2] -> sector,panel,LR

	// status:
	//	0 - fully functioning
	//	1 - noADC
	//	2 - noTDC
	//	3 - noADC, noTDC (PMT is dead)
	//      5 - any other reconstruction problem
	vector<int> status[6][3][2];

    // tdc_conc: tdc conversion factors
    vector<double> tdcconv[6][3][2];
    
	// veff: effective velocity
	vector<double> veff[6][3][2];

	// attlen: attenuation length
	vector<double> attlen[6][3][2];

	// countsForMIP: Desired ADC channel for MIP peak calibration
	vector<double> countsForMIP[6][3][2];

	// twlk: Time walk correction, 3 constants each for L and R
	vector<double> twlk[6][3][6];

	// toff_LR and tof_P2P: time offsets for Left-Right and Paddle-to-Paddle
        vector<double> toff_LR[6][3];
        vector<double> toff_RFpad[6][3];
        vector<double> toff_P2P[6][3];

    // ======== FADC Pedestals and sigmas ===========
    double pedestal[6][3][62][2] = {};
    double pedestal_sigm[6][3][62][2] = {};
    
    // tres: Gaussian sigma for smearing time resolution
    vector<double> tres[3];
    
	int    npaddles[3];  // Number of paddles for Panel 1A, 1B and 2.
	int    thick[3];     // Thickness of paddles (cm) for Panel 1A, 1B and 2.
	double dEdxMIP;      // Nominal MIP specific energy loss (MeV/gm/cm2)
	double dEMIP[3];     // Nominal MIP energy loss (MeV) for Panel 1A, 1B and 2.

	double pmtPEYld;      // Photoelectron yield (p.e./MeV)
	double pmtQE;         // Quantum efficiency of PMT
	double pmtDynodeGain; // PMT dynode gain
	double pmtDynodeK;    // PMT dynode secondary emission statistics factor: K=0 (Poisson) K=1 (exponential)
	double pmtFactor;     // Contribution to FWHM from PMT statistical fluctuations.
//	double tdcLSB;        // Conversion from ns to TDC channel.

	//	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
	double vpar[4];

 	// translation table
	TranslationTable TT;        
       
        
};


// Class definition
/// \class ftof_HitProcess
/// <b> Forward Time of Flight Hit Process Routine</b>\n\n

class ftof_HitProcess : public HitProcess
{
public:

	~ftof_HitProcess(){;}

	// constants initialized with initWithRunNumber
	static ftofConstants ftc;

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
	static HitProcess *createHitClass() {return new ftof_HitProcess;}

	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

};



#endif




